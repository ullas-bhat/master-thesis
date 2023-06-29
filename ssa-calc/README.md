# SSA-Calc: Satellite Sun Angles Calculator

ssa-calc is a Python package to calculate the position of a satellite in Earth orbit and its angles with respect to the Sun. This can be used in analyzing satellite properties that are influenced by solar illumination.

The package is based on work by Ziqi Zhang.

## Usage

The package contains a `Satellite` class that can be used to calculate the position and angles of a satellite with respect to the Sun. The position of the satellite is calculated using [SatMAD](https://satmad.readthedocs.io/en/stable/about.html) interpolators stored as either pickle or compressed pickle files. The interpolators can be initialized using an orbital propogator of your choice. The angles are calculated using the [JPL Ephemerides](https://ssd.jpl.nasa.gov/orbits.html) of the Sun and the calculated satellite position.

An example of using the package to calculate the satellite coordinates and its illumination status is shown below:

```python
import numpy as np
import astropy.units as u
from astropy.time import Time

import ssa_calc

POS_INTERP_PATH = r"./data/pos_interp.pkl   # Path to the pickle file containing the position interpolator
VEL_INTERP_PATH = r"./data/vel_interp.pkl   # Path to the pickle file containing the velocity interpolator
satellite = ssa_calc.Satellite(POS_INTERP_PATH, VEL_INTERP_PATH)

# Calculate the illumination status of the satellite at a given time
time = Time("2021-01-01 00:00:00", format="iso")    # time of interest
sat_coord = satellite.get_sat_coord(time)           # satellite coordinates at the time of interest
beta, fe, theta, theta_T = satellite.calc_beta_fe_theta(time)   # satellite position with respect to the sun

if theta_T < 2 * np.pi * (1 - fe):
    print("The satellite is not in eclipse")
else:
    print("The satellite is in eclipse")
```
