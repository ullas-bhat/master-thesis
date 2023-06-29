"""
SSA-Calc: Satellite Sun Angles Calculator
===========================================

Copyright (c) 2022 Ullas Bhat

License: MIT
For full license details, see LICENSE.rst

Module used to create a Satellite class that can be used to calculate the
position and angle of a satellite in Earth orbit with respect to the Sun.
"""

import pickle
import bz2
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from satmad.core.celestial_bodies_lib import SUN, EARTH


class Satellite:
    """
    Class with methods to calculate satellite position and velocity at a given time
    :attribute pos_interp: SatMAD Cartesian3D interpolator for satellite position
    :attribute vel_interp: SatMAD Cartesian3D interpolator for satellite velocity

    :method coord_from_interp: function to generate SkyCoord object of satellite
    position at input time
    :method __angle_between_vectors: internal function to calculate the angle
    between two input vectors
    :method calc_beta_fe_theta: function to calculate the beta angle, fraction
    of eclipse, theta and transformed theta angles
    """

    def __init__(self, pos_interp_path, vel_interp_path):
        """
        :param pos_interp_path: path to the SatMAD Cartesian3D interpolator for
        satellite position
        :param vel_interp_path: path to the SatMAD Cartesian3D interpolator for
        satellite velocity
        """

        if pos_interp_path[-3:] == "pkl" and vel_interp_path[-3:] == "pkl":
            print("Loading interpolators from pickle files...")
            with open(pos_interp_path, "rb") as f:
                self.pos_interp = pickle.load(f)
            with open(vel_interp_path, "rb") as f:
                self.vel_interp = pickle.load(f)
            print("Interpolators loaded.")
        elif pos_interp_path[-4:] == "pbz2" and vel_interp_path[-4:] == "pbz2":
            print("Loading interpolators from compressed pickle files...")
            with bz2.BZ2File(pos_interp_path, "rb") as f:
                self.pos_interp = pickle.load(f)
            with bz2.BZ2File(vel_interp_path, "rb") as f:
                self.vel_interp = pickle.load(f)
            print("Interpolators loaded.")
        else:
            raise ValueError("Invalid file type for interpolator. Must be .pkl or .pbz2")

    def coord_from_interp(self, time):
        """
        Function to generate SkyCoord object of satellite position at input time
        :param time: astropy Time object

        :return: SkyCoord object of satellite position at input time
        """

        # Get position and velocity from interpolators:
        pos = self.pos_interp(time.jd)
        vel = self.vel_interp(time.jd)

        # Convert to astropy.units.Quantity objects:
        pos = u.Quantity(pos, u.km)
        vel = u.Quantity(vel, u.km / u.s)

        # Generate SkyCoord object from position and velocity:
        coord = SkyCoord(
            representation_type="cartesian",
            x=pos[0],
            y=pos[1],
            z=pos[2],
            v_x=vel[0],
            v_y=vel[1],
            v_z=vel[2],
            obstime=time.iso,
            frame="gcrs",
        )
        return coord

    def __angle_between_vectors(self, vec_1, vec_2):
        """Function to calculate the angle between two input vectors
        :param vec_1: vector 1
        :param vec_2: vector 2

        :return: angle between vectors in radians
        """

        return np.arccos(
            np.dot(vec_1, vec_2) / (np.linalg.norm(vec_1) * np.linalg.norm(vec_2))
        )

    def calc_beta_fe_theta(self, time):
        """
        Function to calculate the beta, theta, theta_T and fe values for a given time
        :param time: astropy Time object

        :return: beta, theta, theta_T and fe values for a given time
        """

        # Get posision of satellite in GCRS frame at input time:
        sat_coord = self.coord_from_interp(time)

        # Get position of SUN in GCRS frame at input time:
        sun_coord = SUN.get_coord_list(time, ephemeris="jpl").gcrs

        # Decompose coordinates into position and velocity vectors:
        sat_pos = sat_coord.cartesian.xyz  # satellite position vector
        sat_vel = sat_coord.velocity.d_xyz  # satellite velocity vector
        sun_pos = sun_coord.cartesian.xyz  # sun position vector

        sat_h = np.cross(sat_pos, sat_vel)  # satellite angular momentum vector
        sat_h_unit = sat_h / np.linalg.norm(sat_h)  # unit normal vector to the orbital plane

        # Calculating beta angle:
        beta = np.pi / 2 * u.rad - self.__angle_between_vectors(sun_pos, sat_h_unit)
        if beta < 0 * u.rad:
            beta = -beta  # beta angle is always positive

        # Calculating fe value:
        Re = EARTH.ellipsoid.re  # Earth radius
        h = np.linalg.norm(sat_pos) - Re  # satellite altitude
        beta_str = np.arcsin(Re / (Re + h))  # critical beta angle

        if beta < beta_str:
            fe = np.arccos(np.sqrt(h**2 + 2 * Re * h) / ((Re + h) * np.cos(beta))) / (
                np.pi * u.rad
            )
        else:
            fe = 0  # no eclipse

        # Calculating theta angle:
        N = (
            sat_h_unit * np.linalg.norm(sun_pos) * np.sin(beta)
        )  # normal vector from orbital plane to sun
        if self.__angle_between_vectors(sun_pos, sat_h_unit) > np.pi / 2 * u.rad:
            N = -N  # N vector always points to the sun
            print("transform N")
        sun_proj = sun_pos - N  # sun vector projection on the orbital plane
        theta = self.__angle_between_vectors(sun_proj, sat_pos)
        if (
            np.dot(sat_h_unit, np.cross(sun_proj, sat_pos))
            / (np.linalg.norm(np.cross(sun_proj, sat_pos)))
            < 0
        ):
            theta = 2 * np.pi * u.rad - theta

        # Calculating theta_T angle:
        theta_T = theta + np.pi * u.rad * (1 - fe)
        if theta > np.pi * u.rad * (1 + fe):
            theta_T = theta - np.pi * u.rad * (1 + fe)

        return beta, fe, theta, theta_T
