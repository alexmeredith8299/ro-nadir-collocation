"""
This module contains useful constants relating to Earth and conversions
between different time standards. It also contains utility functions 
related to these conversions and properties.
"""
import numpy as np

### Earth-related constants #############

Earth_equatorial_radius = 6378137.0/1000 #km
Earth_polar_radius = 6356752.3142/1000 #km
km_to_degree = 1/111 #approx. equal to latitude-dependent value of 360/(2*math.pi*Earth radius) [degree/km]
mu = 3.986004418e5 #standard gravitational parameter of Earth [km^3/s^2]
sec_to_sidereal_day = 1/86164.0905 #sidereal day/sec

##########################################

### Simple time constants ###############

hours_to_sec = 3600 #sec/hr
sec_to_hours = 1/hours_to_sec #hr/sec
min_to_sec = 60 #sec/min
microsec_to_sec = 1e-6 #sec/microsec
day_to_sec = 24*hours_to_sec #sec/day

###########################################

### Time conversions #####################

utc_to_tai = 37 #UTC lags TAI because it has leap seconds and TAI doesn't [seconds]
tai_to_utc = -utc_to_tai #seconds
j2000_epoch_to_unix = 946684800 #offset between start of Unix time (1/1/1970) and start of J2000 time (1/1/2000) [seconds]
gps_epoch_to_unix = 315964800 #offset between start of Unix time (1/1/1970) and start of GPS time (6/1/1980) [seconds]
tai_epoch_to_unix = -378691200 #offset between start of Unix time (1/1/1970) and start of TAI epoch (1/1/1958) [seconds]
y2021_epoch_to_unix = 1609459200 #offset between start of Unix time (1/1/1970) and start of 2021

###########################################

## Utility functions #######################

def calculate_radius_of_earth(lat_geodetic):
    """
    This function calculates the latitude-dependent radius of the Earth.

    Parameters 
    ------------
        lat_geodetic: float
            Geodetic latitude[radians]

    Returns 
    --------
        r_e: float
            Latitude-adjusted radius of Earth, in km
    """
    a = Earth_equatorial_radius
    b = Earth_polar_radius
    r_e = np.sqrt(((a**2*np.cos(lat_geodetic))**2 + (b**2*np.sin(lat_geodetic))**2)/((a*np.cos(lat_geodetic))**2+(b*np.sin(lat_geodetic))**2))
    return r_e
