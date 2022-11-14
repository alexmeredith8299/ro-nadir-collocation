import os
import numpy as np
import bisect
from constants_and_utils import km_to_degree, sec_to_hours

def brute_force(mw_sat,ro_sat,time_tolerance, spatial_tolerance, mw_data, mw_times, ro_times, ro_lats, ro_lons, use_sorting=True):
    """
    This function is a wrapper that performs brute-force colocation finding. It loops over all occultation-microwave sounding pairs
    and checks each pair for cotemporality and spatial overlap, in order to find colocations that overlap in both
    time and space.

    Arguments
    ------------
        mw_sat: MWR_G
            Generalized microwave sounder object representing a microwave satellite with known sounding times and locations
        ro_sat: RO
            Radio occultation satellite object representing an RO sat with known occultation times and locations
        time_tolerance: float
            Time tolerance for colocations in seconds
        spatial_tolerance: float
            Spatial tolerance for colocations in km

    Returns
    -----------
        colocs: float
            Number of colocations
        clats: list
            Latitudes of radio occultations at each colocation (in radians)
        clons: list
            Longitudes of radio occultations at each colocation (in radians)
        ctimes: list
            Times of radio occultations at each colocation (in Unix hours)
    """
    if mw_sat.day != ro_sat.day:
        raise ValueError('ERROR: mw_sat.day != ro_sat.day')

    colocs = 0
    clats, clons, ctimes = [], [], []

    #for each occultation (each file in the folder represents a single occultation)
    for i in range(len(ro_times)):
        coloc = False
        occult_time, occult_lat, occult_lon = ro_times[i], ro_lats[i], ro_lons[i]
        if use_sorting:
            coloc, mw_lat, mw_lon, mw_time = sorted_brute_force_single_occultation(mw_times, mw_data, occult_time, occult_lat, occult_lon, time_tolerance, spatial_tolerance)
        else:
            coloc, mw_lat, mw_lon, mw_time = brute_force_single_occultation(mw_data, occult_time, occult_lat, occult_lon, time_tolerance, spatial_tolerance)

        if coloc:
            colocs += 1
            clats.append(occult_lat)
            clons.append(occult_lon)
            ctimes.append(occult_time)

    return colocs, clats,clons,ctimes

def sorted_brute_force_single_occultation(mw_times_sorted, mw_data, ro_time,ro_lat,ro_lon,time_tolerance, spatial_tolerance):
    """
    This function uses a brute-force method to find microwave soundings colocated with a single radio occultation.

    Arguments
    ------------
        mw_times_sorted: list
            Sorted list containing MW sounding times.
        mw_filelist: list
            Sorted list containing one entry for each file of MW data. Each entry contains a list of longitudes, latitudes,
            mid-sounding times, spacecraft positions, and spacecraft velocities.
        ro_time: float
            Time of radio occultation in Unix seconds
        ro_lat: float
            Latitude of radio occultation in radians
        ro_lon: float
            Longitude of radio occultation in radians
        time_tolerance: float
            Time tolerance for colocation in seconds
        spatial_tolerance: float
            Spatial tolerance for colocation in km
    Returns
    -----------
        coloc: boolean
            True iff colocation exists for this occultation
        lat: float
            MW beam footprint center latitude at colocation (if colocation exists)
        lon: float
            MW beam footprint center longitude at colocation (if colocation exists)
        time: float
            Time of MW sounding in Unix hours (if colocation exists)
    """
    start_ind = bisect.bisect_left(mw_times_sorted,ro_time-time_tolerance*sec_to_hours)
    end_ind = bisect.bisect_right(mw_times_sorted,ro_time+time_tolerance*sec_to_hours)
    for i in range(start_ind, end_ind): #iterate through each file
        lats, mid_time = mw_data[i][1], mw_data[i][2]
        for j in range(len(lats)):
            if abs(lats[j]-ro_lat) < np.deg2rad(spatial_tolerance*km_to_degree): #slight efficiency improvement to check if latitude is close before doing more calcs
                if np.sqrt((((mw_data[i][0][j] - ro_lon)**2)*np.cos(ro_lat)**2) + ((lats[j] - ro_lat)**2)) < np.deg2rad(spatial_tolerance*km_to_degree):
                    return True, lats[j], mw_data[i][0][j], mid_time
    return False, 0,0,0

def brute_force_single_occultation(mw_filelist, ro_time,ro_lat,ro_lon,time_tolerance, spatial_tolerance):
    """
    This function uses a brute-force method to find microwave soundings colocated with a single radio occultation.

    Arguments
    ------------
        mw_filelist: list
            List containing one entry for each file of MW data. Each entry contains a list of longitudes, latitudes,
            mid-sounding times, spacecraft positions, and spacecraft velocities.
        ro_time: float
            Time of radio occultation in Unix seconds
        ro_lat: float
            Latitude of radio occultation in radians
        ro_lon: float
            Longitude of radio occultation in radians
        time_tolerance: float
            Time tolerance for colocation in seconds
        spatial_tolerance: float
            Spatial tolerance for colocation in km
    Returns
    -----------
        coloc: boolean
            True iff colocation exists for this occultation
        lat: float
            MW beam footprint center latitude at colocation (if colocation exists)
        lon: float
            MW beam footprint center longitude at colocation (if colocation exists)
        time: float
            Time of MW sounding in Unix hours
    """
    for l in range(0, len(mw_filelist)): #iterate through each file
        longitudes, latitudes, mid_times= mw_filelist[l][0], mw_filelist[l][1], mw_filelist[l][2]
        for u in range(len(latitudes)):
            if abs(ro_time - mid_times[u]) < (time_tolerance*sec_to_hours): #slight efficiency improvement to move this out of innermost loop
                lat_u = latitudes[u]
                for c in range(len(lat_u)):
                    if abs(lat_u[c]-ro_lat) < np.deg2rad(spatial_tolerance*km_to_degree): #slight efficiency improvement to check if latitude is close before doing more calcs
                        if np.sqrt((((longitudes[u][c] - ro_lon)**2)*np.cos(ro_lat)**2) + ((lat_u[c] - ro_lat)**2)) < np.deg2rad(spatial_tolerance*km_to_degree):
                            return True, lat_u[c], longitudes[u][c], mid_times[u]
    return False, 0,0,0
