import numpy as np
from scipy.integrate import cumulative_trapezoid
import PyIRI
import PyIRI.main_library as ml
import pymsis

from config import (yr, mn, dy, F107, F107A, AP, lat0, lon0, kB, g0, RE, omega, kh, Pr, z_ref)


def iri_profiles(ut_h, lats, lons, z_km):
    F2, F1, E, Es, sun, mag, EDP = ml.IRI_density_1day(yr, mn, dy, np.atleast_1d(ut_h).astype(float),
        np.atleast_1d(lons).astype(float), np.atleast_1d(lats).astype(float),z_km, F107, PyIRI.coeff_dir, 0)
    N0 = np.transpose(EDP, (0, 2, 1))
    return np.maximum(N0, 1e6), F2


def msis_profiles(ut_h, z_km):
    
    ut_h = np.atleast_1d(ut_h)
    dates = np.array([np.datetime64(f"{yr:04d}-{mn:02d}-{dy:02d}T00:00") + np.timedelta64(int(round(h * 3600)), "s") for h in ut_h])
    nt = len(dates)
    out = pymsis.calculate(dates, lon0, lat0, z_km,f107s=np.full(nt, F107), f107as=np.full(nt, F107A),aps=np.full((nt, 7), AP))
    out = np.nan_to_num(out.reshape(nt, len(z_km), -1))

    V = pymsis.Variable
    ms = {k: out[..., getattr(V, k)] for k in
         ["MASS_DENSITY", "N2", "O2", "O", "HE", "H", "AR", "N", "ANOMALOUS_O", "NO", "TEMPERATURE"]}
    ms["rho"], ms["Tn"] = ms.pop("MASS_DENSITY"), ms.pop("TEMPERATURE")
    return ms


def gw_struct(ms, z_km):

    z = z_km * 1e3
    Tn, rho = ms["Tn"], ms["rho"]
    mono = ms["O"] + ms["HE"] + ms["H"] + ms["AR"] + ms["N"] + ms["ANOMALOUS_O"]
    di = ms["N2"] + ms["O2"] + ms["NO"]
    n = mono + di
    mbar = rho / n
    g = g0 * (RE / (RE + z)) ** 2
    cp = kB * (2.5 * mono + 3.5 * di) / n / mbar
    cv = kB * (1.5 * mono + 2.5 * di) / n / mbar
    gamma = cp / cv

    H = kB * Tn / (mbar * g)
    cs2 = gamma * g * H
    N2 = g / Tn * (np.gradient(Tn, z, axis=-1) + g / cp)
    wa2 = cs2 / (4 * H ** 2)
    m2 = kh ** 2 * (N2 / omega ** 2 - 1) + (omega ** 2 - wa2) / cs2   # inviscid (Hines 1960)
    m_inv = -np.sqrt(m2)                            # m < 0: upward energy
    cgz = -m_inv * omega ** 3 / (kh ** 2 * N2 - omega ** 4 / cs2)

    # Dissipative dispersion relation (Vadas & Fritts 2005), q = kh^2 + m^2:
    # (omega + i nu q)(omega + i ka q)(q + 1/4H^2) = N^2 kh^2
    nu = 3.34e-7 * Tn ** 0.71 / rho
    ka = nu / Pr
    s = 1 / (4 * H ** 2)
    a3 = -nu * ka + 0j
    a2 = 1j * omega * (nu + ka) - nu * ka * s
    a1 = omega ** 2 + 1j * omega * (nu + ka) * s
    a0 = omega ** 2 * s - N2 * kh ** 2 + 0j
    comp = np.zeros(Tn.shape + (3, 3), complex)
    comp[..., 0, 0], comp[..., 0, 1], comp[..., 0, 2] = -a2 / a3, -a1 / a3, -a0 / a3
    comp[..., 1, 0] = comp[..., 2, 1] = 1.0
    roots = np.linalg.eigvals(comp)                      

    # Follow the gravity-wave branch upward by continuity
    q = np.empty(Tn.shape, complex)
    prev = kh ** 2 + m_inv[..., 0] ** 2
    for k in range(Tn.shape[-1]):
        r = roots[..., k, :]
        j = np.argmin(np.abs(r - prev[..., None]), axis=-1)
        q[..., k] = np.take_along_axis(r, j[..., None], axis=-1)[..., 0]
        prev = q[..., k]

    m = np.sqrt(q - kh ** 2)
    m = np.where(m.imag < 0, -m, m)                       # Im m >= 0: decay upward
    iref = np.argmin(np.abs(z_km - z_ref))
    ln_a_nodamp = -0.5 * np.log(rho)
    lna = ln_a_nodamp - cumulative_trapezoid(m.imag, z, axis=-1, initial=0)
    a = np.exp(lna - lna[..., [iref]])
    a_nodamp = np.exp(ln_a_nodamp - ln_a_nodamp[..., [iref]])
    phi = cumulative_trapezoid(m.real, z, axis=-1, initial=0)
    phi = phi - phi[..., [iref]]
    return dict(z_km=z_km, H=H, N2=N2, T_BV=2 * np.pi / np.sqrt(np.maximum(N2, 1e-12)),
                cs2=cs2, gamma=gamma, m2=m2, m_inv=m_inv, cgz=cgz, m=m,
                lamZ=2 * np.pi / np.abs(m.real), lamZ_inv=2 * np.pi / np.abs(m_inv),
                nu=nu, diss_ratio=nu * (kh ** 2 + m_inv ** 2) / omega,
                a=a, a_nodamp=a_nodamp, phi=phi, mbar=mbar, Tn=Tn, rho=rho)
