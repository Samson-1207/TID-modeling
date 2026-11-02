"""
Single source of all model parameters.

Every module (iri_msis_gw.py, numerical_model.py) and the notebook import
their settings from here, so a change made in this file is seen everywhere.
"""
import numpy as np
from datetime import datetime, timedelta, timezone

# ---------------------------------------------------------------------------
# Site, date and space-weather drivers
# ---------------------------------------------------------------------------
lat0, lon0 = 69.12, 20.0          # Tromsø / EISCAT [deg N, deg E]
yr, mn, dy = 2020, 12, 8          # date
F107, F107A, AP = 82.4, 85.4, 4   # daily F10.7, 81-day mean F10.7, Ap

# ---------------------------------------------------------------------------
# Grids
# ---------------------------------------------------------------------------
alt = np.arange(80, 500.01, 5)              # display altitude grid [km]
altF = np.arange(80, 500.01, 1)             # computation altitude grid [km]
ialtF_out = np.searchsorted(altF, alt)      # indices of alt inside altF
xh = np.arange(-500, 500.01, 10)            # horizontal east grid [km]
yh = np.arange(-500, 500.01, 10)            # horizontal north grid [km]
dt_min = 1                                  # time step [min]
ut = np.arange(0, 24, dt_min / 60)          # UT [h]
t0 = datetime(yr, mn, dy, tzinfo=timezone.utc)
times = np.array([t0 + timedelta(hours=float(h)) for h in ut])

# ---------------------------------------------------------------------------
# TID / gravity-wave parameters
# ---------------------------------------------------------------------------
period = 40            # wave period [min]
lambda_h = 200         # horizontal wavelength [km]
azimuth = 170          # propagation azimuth, clockwise from north [deg]
amp_Fpeak = 0.1        # |dNe/Ne| imposed at hmF2
z_ref = 250            # reference altitude for wind amplitude / phase [km]
sigma_efolding = 60    # reserved: packet e-folding time [min] (not used yet)
lambda_z = None        # reserved: prescribed vertical wavelength (not used yet)

# ---------------------------------------------------------------------------
# Hooke-model switches and plasma parameters
# ---------------------------------------------------------------------------
recomb_loss = True     # include chemical loss L in the denominator
diff = True            # include field-aligned ambipolar diffusion Da*k_par^2
TeOverTn, TiOverTn = 1, 1
Pr = 0.7               # Prandtl number (thermal diffusivity = nu / Pr)
smooth_km = 5          # Gaussian smoothing of ln N0 before differentiation [km]

fradar = 224e6         # radar frequency [Hz] (EISCAT VHF)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
kB = 1.380649e-23
amu = 1.66053907e-27
qe = 1.602176634e-19
g0 = 9.80665
RE = 6371e3
c0 = 299792458
m_Op = 16 * amu

# ---------------------------------------------------------------------------
# Derived wave quantities
# ---------------------------------------------------------------------------
omega = 2 * np.pi / (period * 60)          # [rad/s]
kh = 2 * np.pi / (lambda_h * 1e3)          # [rad/m]
az = np.deg2rad(azimuth)
kx, ky = kh * np.sin(az), kh * np.cos(az)  # x = east, y = north
c_phase = omega / kh                        # horizontal phase speed [m/s]
