"""
Single source of all model parameters.
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
alt = np.arange(80, 500.01, 5)              
altF = np.arange(80, 500.01, 1)            
ialtF_out = np.searchsorted(altF, alt)      
xh = np.arange(-500, 500.01, 10)            
yh = np.arange(-500, 500.01, 10)           
dt_min = 1                                  
ut = np.arange(0, 24, dt_min / 60)         
t0 = datetime(yr, mn, dy, tzinfo=timezone.utc)
times = np.array([t0 + timedelta(hours=float(h)) for h in ut])

# TID / gravity-wave parameter

period = 40            # mstid period
lambda_h = 200         # horizontal wavelength
azimuth = 170          # propagation azimuth
amp_Fpeak = 0.1        # 10%  imposed at hmF2
z_ref = 250            # reference altitude for wind amplitude and phase explanation

# AGW-plasma paraeters

recomb_loss = True     
diff = True           
TeOverTn = 1
TiOverTn = 1
Pr = 0.7               
smooth_km = 5         

fradar = 224e6         # radar frequency [Hz] (EISCAT VHF)

# Constants

kB = 1.380649e-23
amu = 1.66053907e-27
qe = 1.602176634e-19
g0 = 9.80665
RE = 6371e3
c0 = 299792458
m_Op = 16 * amu

# Derived wave quantities

omega = 2 * np.pi / (period * 60)          # [rad/s]
kh = 2 * np.pi / (lambda_h * 1e3)          # [rad/m]
az = np.deg2rad(azimuth)
kx, ky = kh * np.sin(az), kh * np.cos(az)  # x = east, y = north
c_phase = omega / kh                        # horizontal phase speed [m/s]
