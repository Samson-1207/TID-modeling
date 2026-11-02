"""
Hooke (1968) TID model: ion response to a gravity-wave neutral wind in the
geomagnetic field, with chemical loss and field-aligned ambipolar diffusion.
All parameters come from config.py.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
import ppigrf

from config import (lat0, lon0, t0, altF, kB, qe, m_Op, omega, kh, az, kx, ky,
                    TeOverTn, TiOverTn, smooth_km, recomb_loss, diff)

# Geomagnetic field at 300 km (IGRF). ppigrf returns east, north, up [nT].
Be, Bn, Bu = ppigrf.igrf(lon0, lat0, 300, t0.replace(tzinfo=None))
Bvec = np.array([Be, Bn, Bu]).ravel() * 1e-9   # [T]
B0 = np.linalg.norm(Bvec)
bhat = Bvec / B0
bx, by, bz = bhat
I_deg = np.rad2deg(np.arcsin(-bz))             # inclination (positive downward)
D_deg = np.rad2deg(np.arctan2(bx, by))         # declination (positive east)
kb_h = kx * bx + ky * by                       # horizontal part of k.b


def ion_rates(ms, N0):
    """Perturbation loss rate L [1/s], ambipolar diffusion coefficient
    Da [m^2/s] and ion-neutral collision frequency nu_in [1/s]."""
    Tn = ms["Tn"]
    if N0.ndim == 3 and Tn.ndim == 2:
        ms = {k: v[:, None, :] for k, v in ms.items()}
        Tn = ms["Tn"]

    Te, Ti = TeOverTn * Tn, TiOverTn * Tn
    r = np.clip(Tn, 300, 1700) / 300
    k1 = 1.533e-12 - 5.92e-13 * r + 8.6e-14 * r ** 2                       # O+ + N2 [cm^3/s]
    k2 = 2.82e-11 - 7.74e-12 * r + 1.073e-12 * r ** 2 - 5.17e-14 * r ** 3 + 9.65e-16 * r ** 4  # O+ + O2
    beta = (k1 * ms["N2"] + k2 * ms["O2"]) * 1e-6       # [1/s]
    alpha = 3e-7 * (300 / Te) ** 0.7 * 1e-6             # [m^3/s]
    aN = alpha * N0
    L = beta * aN * (2 * beta + aN) / (beta + aN) ** 2  # linearised beta/alpha loss
    Tr = 0.5 * (Ti + Tn)
    nu_in = (6.82e-10 * ms["N2"] + 6.64e-10 * ms["O2"]) * 1e-6 \
        + 3.67e-11 * ms["O"] * 1e-6 * np.sqrt(Tr) * (1 - 0.064 * np.log10(Tr)) ** 2
    Da = kB * (Te + Ti) / (m_Op * nu_in)
    return L, Da, nu_in


def neutral_wind(gw):
    """Complex neutral wind (uE, uN, uZ) for |u_h| = 1 m/s at z_ref."""
    uh = gw["a"] * np.exp(1j * gw["phi"])
    w = -kh * uh / (gw["m"] + 1j / (2 * gw["H"]))       # anelastic continuity
    Uvec = np.stack([uh * np.sin(az), uh * np.cos(az), w])
    return uh, Uvec


def ion_velocity(Uvec, kappa):
    """Steady ion drift driven by the neutral wind; kappa = Omega_i / nu_in."""
    b = bhat.reshape((3,) + (1,) * (Uvec.ndim - 1))
    Upar = np.sum(Uvec * b, axis=0)
    Uperp = Uvec - Upar * b
    Uxb = np.stack([Uperp[1] * b[2] - Uperp[2] * b[1],
                    Uperp[2] * b[0] - Uperp[0] * b[2],
                    Uperp[0] * b[1] - Uperp[1] * b[0]])
    return Upar * b + (Uperp + kappa * Uxb) / (1 + kappa ** 2)


def hooke_delta(N0, gw, ms):
    """Complex relative perturbation dNe/Ne for a 1 m/s wind at z_ref,
    and the complex horizontal neutral wind."""
    uh, Uvec = neutral_wind(gw)
    m = gw["m"]
    if N0.ndim == 3:
        Uvec, m = Uvec[:, :, None, :], m[:, None, :]
    L, Da, nu_in = ion_rates(ms, N0)
    kappa = (qe * B0 / m_Op) / nu_in
    vi = ion_velocity(Uvec, kappa)
    dz = (altF[1] - altF[0]) * 1e3
    dvz = np.gradient(vi[2], dz, axis=-1)
    lnN0_s = gaussian_filter1d(np.log(N0), smooth_km / (altF[1] - altF[0]), axis=-1, mode="nearest")
    dlnN0 = np.gradient(lnN0_s, dz, axis=-1)
    kpar = kb_h + m.real * bz
    denom = -1j * omega + (L if recomb_loss else 0) + (Da * kpar ** 2 if diff else 0)
    delta_hat = -((1j * (kx * vi[0] + ky * vi[1]) + dvz) + vi[2] * dlnN0) / denom
    return delta_hat, uh
