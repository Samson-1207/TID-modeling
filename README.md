# TID-modeling

Modelling medium-scale travelling ionospheric disturbances (MSTIDs) with
**Hooke's (1968) TID model** on top of **IRI** electron-density profiles,
driven by a dissipating gravity wave in an **NRLMSIS** thermosphere.
The default case is Tromsø / EISCAT (69.12° N, 20° E) on 8 December 2020.

## Model chain

```
MSIS (neutral air)                  IRI (background Ne, hmF2, NmF2)
      │                                          │
      ▼                                          │
gw_struct: gravity-wave vertical structure       │
(Hines 1960 inviscid + Vadas & Fritts 2005       │
 viscous/conductive dispersion relation)         │
      │                                          │
      ▼                                          ▼
neutral_wind → ion_velocity (IGRF B field) → hooke_delta → δNe/Ne (z, t) and (x, y, t)
                                  (chemical loss L, ambipolar diffusion Da)
```

The ion continuity equation solved for the relative perturbation δ = n′/N₀ is

    δ̂ = −[ i k_h·v + ∂v_z/∂z + v_z ∂lnN₀/∂z ] / ( −iω + L + D_a k∥² )

where **v** is the ion drift driven by the gravity-wave wind in the geomagnetic field.

## Repository contents

| File | Purpose |
|---|---|
| `config.py` | **All parameters** (site, date, grids, wave, model switches, constants). Edit here only. |
| `iri_msis_gw.py` | IRI background, MSIS neutral atmosphere, gravity-wave vertical structure |
| `numerical_model.py` | IGRF field, ion chemistry/collision rates, ion drift, Hooke continuity solution |
| `mstid_model.ipynb` | Runs the model and produces all figures (already executed, so figures show on GitHub) |
| `tests/test_model.py` | Smoke tests (`pytest`) |

## Installation

Python ≥ 3.10.

```bash
git clone https://github.com/<your-username>/TID-modeling.git
cd TID-modeling
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

```bash
pytest -q tests                    # quick check (a few seconds)
jupyter lab mstid_model.ipynb      # run the full model (about 20 s)
```

Change the scenario (period, wavelength, azimuth, amplitude, date, location, ...)
in `config.py`, then restart the kernel and run all cells.

## Main parameters (defaults)

| Parameter | Value | Meaning |
|---|---|---|
| `period` | 40 min | wave period |
| `lambda_h` | 200 km | horizontal wavelength (phase speed 83 m/s) |
| `azimuth` | 170° | propagation direction, clockwise from north |
| `amp_Fpeak` | 0.10 | \|δNe/Ne\| imposed at hmF2 |
| `z_ref` | 250 km | reference height for wind amplitude and phase |
| `Pr` | 0.7 | Prandtl number |
| `recomb_loss`, `diff` | True | include chemical loss / field-aligned diffusion |

## Assumptions and limitations

- Linear, monochromatic plane wave; WKB vertical structure (the |m|^½ refraction factor is omitted).
- Neutral atmosphere taken from MSIS at the site only (horizontally uniform).
- No electric fields (polarisation or high-latitude convection), no particle precipitation,
  no production term, no wave-induced neutral composition changes.
- The wind amplitude is rescaled at every time step so that |δNe/Ne| = `amp_Fpeak` at hmF2.

## References

- Hooke, W. H. (1968), Ionospheric irregularities produced by internal atmospheric gravity waves, *J. Atmos. Terr. Phys.*, 30, 795–823.
- Hines, C. O. (1960), Internal atmospheric gravity waves at ionospheric heights, *Can. J. Phys.*, 38, 1441–1481.
- Vadas, S. L., & Fritts, D. C. (2005), Thermospheric responses to gravity waves: Influences of increasing viscosity and thermal diffusivity, *J. Geophys. Res.*, 110, D15103.
- Schunk, R. W., & Nagy, A. F. (2009), *Ionospheres: Physics, Plasma Physics, and Chemistry*, Cambridge University Press.
- PyIRI: https://github.com/victoriyaforsythe/PyIRI · pymsis: https://github.com/SWxTREC/pymsis · ppigrf: https://github.com/IAGA-VMOD/ppigrf
