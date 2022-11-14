# -*- coding: utf-8 -*-
import os
import numpy as np
from sgp4.ext import rv2coe
from sgp4.api import Satrec
import math
from constants_and_utils import tai_to_utc, km_to_degree, hours_to_sec, min_to_sec

def RTCD(mw_sat,ro_sat,time_tolerance,spatial_tolerance,N, times, lats, lons):
    """
    This method finds colocations between radio occultations and microwave soundings for two satellites over
    one day of data. It considers colocations to occur when an occultation is colocated with the center of a footprint
    for a microwave sounding, and uses spatial and time checks to ensure that colocations occur within reasonable
    bounds for both time and space.

    Arguments
    -----------
        mw_sat: MW_G
            Microwave satellite
        ro_sat: RO
            RO satellite
        time_tolerance: float
            Time tolerance [seconds]
        N: int
            Number of pseudo-occultations (use 2 for linearized)
        times: np array
            Occultation times
        lats: np array
            Occultation latitudes
        lons: np array
            Occultation longitudes
    Returns
    ---------
        colocs: int
            Number of colocations
        clats: list
            List of colocation latitudes, corresponding to occultation latitude in radians
        clons: list
            List of colocation longitudes, corresponding to occultation longitude in radians
        ctimes: list
            List of colocation times, corresponding to occultation time in Unix hours
    """
    if mw_sat.day != ro_sat.day:
        raise ValueError(f'ERROR: MW satellite day {mw_sat.day} != RO satellite day {ro_sat.day}')

    p_matrix, times, lats, lons,rises = mw_sat.calc_p_matrix_fast_RTCD(time_tolerance, N, times, lats, lons)
    N = p_matrix.shape[1]
    if N == 1:
        dt = 0
    else:
        dt = 2*time_tolerance/(N-1)
    num_occ = p_matrix.shape[0]

    clats, clons, ctimes = [], [], []
    colocs = 0

    #Get sat period (we assume it doesn't change much over the day)
    s, tline = mw_sat.get_current_tle(mw_sat.day,times[0]*hours_to_sec+tai_to_utc)
    sat = Satrec.twoline2rv(s, tline)
    sat_period = 2*math.pi/sat.nm

    for j in range(num_occ):
        tj = times[j]
        lat = lats[j]
        lon = lons[j]
        #Get max scan angle
        ds = mw_sat.calculate_ds(lat)
        max_scan_angle = ((abs(ds)) + np.deg2rad(spatial_tolerance*km_to_degree))

        coloc = 0

        for k in range(1,N):
            p_prev_jk = p_matrix[j,k-1, :, 0]
            t_prev_jk = tj*hours_to_sec + (k-1-(N-1)/2)*dt
            p_jk = p_matrix[j, k, :, 0]
            t_jk = tj*hours_to_sec + (k-(N-1)/2)*dt

            x_prev, y_prev, z_prev = p_prev_jk[0], p_prev_jk[1], p_prev_jk[2]
            x_j, y_j, z_j = p_jk[0], p_jk[1], p_jk[2]

            num_revs = (t_jk-t_prev_jk)/(sat_period*min_to_sec)
            is_coloc = RTCD_core(x_prev, y_prev, z_prev, x_j, y_j, z_j, max_scan_angle, num_revs, spatial_tolerance)
            if is_coloc:
                colocs += 1
                ctimes.append(tj)
                clats.append(lat)
                clons.append(lon)
                break
    return colocs, clats,clons, ctimes 

def RTCD_core(x_prev, y_prev, z_prev, x_j, y_j, z_j, max_scan_angle, num_revs, spatial_tolerance):
    """
    This function checks if two points in the MW satellite frame are colocated.

    Arguments
    ----------
        x_prev : float
            x-coordinate of start point in MW satellite frame
        y_prev : float
            y-coordinate of start point in MW satellite frame
        z_prev : float
            z-coordinate of start point in MW satellite frame
        x_j : float
            x-coordinate of end point in MW satellite frame
        y_j : float
            y-coordinate of end point in MW satellite frame
        z_j : float
            z-coordinate of end point in MW satellite frame
        max_scan_angle : float
            maximum scan angle of the MW satellite
        num_revs: float
            fractional number of orbits completed between t_prev and t_j
        spatial_tolerance: float
            spatial tolerance for colocation in km

    Returns
    ----------
        is_coloc: bool
            True if the points are colocated, False otherwise
    """
    is_coloc = False

    #Find start and end arglat
    start_arglat = np.arctan2(y_prev,x_prev)
    end_arglat = np.arctan2(y_j,x_j)
    start_delta_arglat = start_arglat
    end_delta_arglat = end_arglat

    #Shift delta arglats such that end < start and |end-start| < 2pi 
    end_delta_arglat, start_delta_arglat = constrain_for_rtcd(end_delta_arglat, start_delta_arglat, num_full_revs=num_revs)
    min_delta_arglat = end_delta_arglat
    max_delta_arglat = start_delta_arglat
    start_left = False

    #Find start and end scan angle
    start_scan_angle = np.arcsin(z_prev)
    end_scan_angle = np.arcsin(z_j)

    #Shift scan angles to be within -pi and pi
    start_scan_angle = constrain_to_pi_range(start_scan_angle)
    end_scan_angle = constrain_to_pi_range(end_scan_angle)

    #Find scan angle in terms of arglat
    slope = (end_scan_angle-start_scan_angle)/(end_delta_arglat-start_delta_arglat)#rise over run

    #Need to check for 2pi crossing if num_full_revs is big enough, etc.
    i = 0 
    while i*2*math.pi < start_delta_arglat:
        y_int = end_scan_angle - slope*end_delta_arglat
        zero_crossing = y_int + i*slope*2*math.pi

        #Check for actually good zero crossing within scan range
        if min_delta_arglat < 0 and max_delta_arglat > 0 and np.abs(zero_crossing) < max_scan_angle:
            is_coloc = True
            break
        i +=1

    #Check for near-zero delta arglat crossing test
    if np.cos(start_delta_arglat) > np.cos(np.deg2rad(spatial_tolerance*km_to_degree)) and np.abs(start_scan_angle) < max_scan_angle:
        is_coloc = True
    if np.cos(end_delta_arglat) > np.cos(np.deg2rad(spatial_tolerance*km_to_degree)) and np.abs(end_scan_angle) < max_scan_angle:
        is_coloc = True

    return is_coloc

def constrain_for_rtcd(ang1, ang2, num_full_revs=0):
    """
    Constrains angles for the rotation method such that 
    -2pi < ang1 < 0 and ang1 < ang2 and such that 
    |ang2-ang1| < 2pi*(num_full_revs+1).

    Arguments
    -----------
        ang1: float
            angle in radians
        ang2: float
            angle in radians 
        num_full_revs: int (optional)
            number of full revolutions between ang1 and ang2 
            (e.g. if 0.75 revolutions between angles, num_full_revs=0)
    Returns
    --------
        ang1: float
            angle in radians constrained to be between -2pi and 0.
        ang2: float
            angle in radians constrained to be between ang1 and ang1+2pi*(num_full_revs+1).
    """
    while ang1 < -2*np.pi:
        ang1 = ang1+2*np.pi
    while ang1 > 0:
        ang1 = ang1-2*np.pi
    while ang2 < ang1+2*np.pi*num_full_revs-np.pi:
        ang2 = ang2+2*np.pi
    while ang2 > ang1+2*np.pi*(num_full_revs+1)-np.pi:
        ang2 = ang2-2*np.pi
    return ang1, ang2

def constrain_to_pi_range(ang):
    """
    Constrains an angle (in radians) to between -pi and +pi.

    Parameters
    -------------
        ang: float
            angle in radians
    Returns
    ---------
        ang: float
            angle in radians constrained to be between -pi and pi.
    """
    while ang < -np.pi:
        ang = ang+2*np.pi
    while ang > np.pi:
        ang = ang-2*np.pi
    return ang


