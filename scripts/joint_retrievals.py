from constants_and_utils import *
from satclass import *
from brute_force_clean import *
from rotation_RO import *
import csv
import sys
import argparse
import textwrap
import logging

def parse_args():
    """
    Parse command line arguments specifiyng the data to download.

    Returns
    --------
        args: argparse.Namespace
            parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Find collocations between an RO and MW satellite', formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
            Example usage:

            python3 joint_retrievals.py --mw_sat MetOpB --ro_sat MetOpB --day 21001 --time_tol 600 --spatial_tol 150 --n_suboccultations 21 --save_dir /Users/alexmeredith/ro-nadir-collocation/results --data_dir /Users/alexmeredith/ro-nadir-collocation/data --tle_dir /Users/alexmeredith/ro-nadir-collocation/tle

            The example above will find collocations between the MetOpB-GRAS instrument and the MetOpB-AMSU instrument on January 1, 2021 with a time tolerance of 10 minutes/600 seconds and a spatial tolerance of 150 km, using four different methods (two brute-force methods, linearized rotation method, and rotation method with sub-occultations with 21 sub-occultations -- equivalently, a sub-occultation every 30 seconds).

            It will then save the results in /Users/alexmeredith/ro-nadir-collocation/results.
            '''))
    parser.add_argument('--mw_sat', help='MW satellite name, must be `MetOpB`, `MetOpC`, `NOAA-20`, or `SNPP`', required=True)
    parser.add_argument('--ro_sat', help='RO satellite name, must be `MetOpB`, `MetOpC`, or `COSMIC-2`', required=True)
    parser.add_argument('--day', help='Date, should fall between 21001 and 21031', required=True)
    parser.add_argument('--time_tol', type=int, help='Time tolerance for collocation in seconds', required=True)
    parser.add_argument('--spatial_tol', type=int, help='Spatial tolerance for collocation in km', required=True)
    parser.add_argument('--n_suboccultations', type=int, help='Number of sub-occultations for rotation method with sub-occultations', required=True)
    parser.add_argument('--save_dir', help="Directory to save results in", required=True)
    parser.add_argument('--data_dir', help="Directory that MW and RO data folders are in", required=True)
    parser.add_argument('--tle_dir', help="Directory that TLEs are saved in", required=True)
    args = parser.parse_args()

    return args


def get_rot_false_pos_neg_timefails_spatialfails(rot_arr, bf_arr):
    """
    This function figures out which colocations failed the spatial check and time check in the rotation method,
    and which only failed the time check. It also correctly ID's false negatives and false positives from the
    rotation method.

    Parameters
    ------------
        rot_arr: list
            list of length-3 lists [lat, lon, time] corresponding to colocations found by rotation method
        bf_arr: list
            list of length-3 lists [lat, lon, time] corresponding to colocations found by brute-force method

    Returns
    ---------
        rot_false_pos: list
            list of length-3 lists [lat, lon, time] corresponding to colocations found by rotation method but
            NOT by brute-force method
        rot_false_neg: list
            list of length-3 lists [lat, lon, time] corresponding to colocations found by brute-force method but
            NOT by rotation method
    """
    rot_false_pos = []
    rot_false_neg = []

    for line in rot_arr:
        if line not in bf_arr:
            rot_false_pos.append(line)

    for line in bf_arr:
        if line not in rot_arr:
            rot_false_neg.append(line)

    rot_false_pos, rot_false_neg = eliminate_rounding_errors(rot_false_pos, rot_false_neg)
    return rot_false_pos, rot_false_neg

def eliminate_rounding_errors(rot_false_pos, rot_false_neg, eps=1e-7):
    """
    This function eliminates "rounding errors" from lists of rotation method false positives and false negatives. Basically, if a false positive and false negative both exist for the same exact occultation time, but the lat/lon of both are within floating point precision, it deletes both and considers those occultations to be the same. This fixes an issue where floating point precision issues with conversions between radians and degrees can lead to the rotation method and brute force method saving a slightly different lat/lon for the same collocation, which looks like disagreement between methods to the computer.

    Parameters
    -------------
        rot_false_pos: list
            list of length-3 lists [lat, lon, time] corresponding to colocations that the rotation method ID'd as
            false positives
        rot_false_neg: list
            list of length-3 lists [lat, lon, time] corresponding to colocations that the rotation method ID'd as
            false negatives

    Returns
    ---------
        fp_actual: list
            list of length-3 lists [lat, lon, time] corresponding to actual false positive colocations
        fn_actual: list
            list of length-3 lists [lat, lon, time] corresponding to actual false negative colocations
    """
    #No overlap if at least one of the lists has zero elements
    if len(rot_false_neg) == 0 or len(rot_false_pos) == 0:
        return rot_false_pos, rot_false_neg

    #Otherwise look for duplicates
    fp_actual, fn_actual = [], []
    for fp in rot_false_pos:
        occ_time, lat, lon = fp[2], fp[0], fp[1]
        match_time = False
        for fn in rot_false_neg:
            occ_time2, lat2, lon2 = fn[2], fn[0], fn[1]
            if occ_time2==occ_time and np.abs(lat2-lat) < eps and np.abs(lon2-lon) < eps:
                match_time = True
        if not match_time:
            fp_actual.append(fp)
            fn_actual.append(fn)
    return fp_actual, fn_actual

def get_coloc_list(clats, clons, ctimes):
    """
    Combines collocation latitudes, longitudes, and times into one list.

    Arguments 
    ----------
        clats: list
            List of floats with collocation latitudes
        clons: list
            List of floats with collocation longitudes
        ctimes: list
            List of collocation times

    Results 
    ---------
       coloc_list: list 
            List of length-3 lists containing lat, lon, time for 
            each collocation
    """
    coloc_list = []
    for i in range(len(clats)):
        coloc = [clats[i], clons[i], ctimes[i]]
        coloc_list.append(coloc)
    return coloc_list

def save_single_method_results(colocs, fname):
    """
    Saves results for a single method in a csv.

    Arguments 
    ----------
        colocs: list
            List of length-3 lists containing lat, lon, time for
            each collocation

        fname: str
            Name of file to save collocations in
    """
    with open(fname, 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(['lat', 'long', 'time'])
        csvwriter.writerows(colocs)

def retrieve_daily_colocations(args):
    day = args.day
    mw_sat_name = args.mw_sat
    ro_sat_name = args.ro_sat
    time_tolerance = args.time_tol
    spatial_tolerance = args.spatial_tol
    rtcd_N_1 = args.n_suboccultations
    save_dir = args.save_dir
    data_dir = args.data_dir
    tle_dir = args.tle_dir
    rtcd_N_2 = 2 

    summary_fname = f"{save_dir}/{day}_{mw_sat_name}_{ro_sat_name}_{time_tolerance}_summary_rtcd_only.txt"
    noaa_20 = MWR_G('NOAA20',7205,52.725,'idps',day, data_dir, tle_dir)
    MetOpC_AMSUA = MWR_G('MetOpC-AMSUA', 7198.5, 48.33, '1.1.2000',day, data_dir, tle_dir)
    MetOpB_AMSUA = MWR_G('MetOpB-AMSUA', 7198.5, 48.33, '1.1.2000',day, data_dir, tle_dir)
    snpp = MWR_G('SNPP', 7198.5, 52.725, 'idps',day, data_dir, tle_dir)

    MetOpC_GRAS = RO('MetOpC-GRAS', day, data_dir)
    MetOpB_GRAS = RO('MetOpB-GRAS', day, data_dir)
    cosmic2 = RO('COSMIC2',day, data_dir)

    ro_dict = {"COSMIC-2":cosmic2, "MetOpC":MetOpC_GRAS, "MetOpB":MetOpB_GRAS}
    mw_dict = {"MetOpC":MetOpC_AMSUA, "MetOpB":MetOpB_AMSUA, "NOAA-20":noaa_20, "SNPP":snpp}

    mw_satellite = mw_dict[mw_sat_name]
    ro_satellite = ro_dict[ro_sat_name]

    #Pre-load RO data
    start_ro_load = time.time()
    ro_times, ro_lats, ro_lons = ro_satellite.load_data()
    end_ro_load = time.time()
    ro_load_time = end_ro_load - start_ro_load
    ro_load_s = f"Loading RO data for {ro_sat_name} on day {day} took {ro_load_time} sec \n"
    logging.info(ro_load_s.strip())

    #Pre-load MW data
    start_mw_load = time.time()
    mw_data = mw_satellite.load_data()
    end_mw_load = time.time()
    start_mw_sorted_load = time.time()
    mw_data_sorted, mw_times = mw_satellite.load_sorted_data()
    end_mw_sorted_load = time.time()
    mw_load_time = end_mw_load - start_mw_load
    mw_sorted_load_time = end_mw_sorted_load - start_mw_sorted_load
    mw_load_s = f"Loading MW data for {mw_sat_name} on day {day} took {mw_load_time} sec \n"
    mw_sorted_load_s = f"Loading sorted MW data for {mw_sat_name} on day {day} took {mw_sorted_load_time} sec \n"
    logging.info(mw_load_s.strip())
    logging.info(mw_sorted_load_s.strip())

    #Compute old brute-force
    start_old_brute_force = time.time()
    colocs_old_bf, clats_old_bf,clons_old_bf,ctimes_old_bf = brute_force(mw_satellite, ro_satellite,time_tolerance,spatial_tolerance, mw_data, mw_times, ro_times, ro_lats, ro_lons, use_sorting=False)
    end_old_brute_force = time.time()
    old_bf_time = end_old_brute_force-start_old_brute_force
    old_bf_s = f"Brute force method #1 for colocations between {mw_sat_name} (MW) and {ro_sat_name} (RO) on day {day}, with time_tol={time_tolerance} completed in {old_bf_time} sec and found {colocs_old_bf} colocations \n"
    old_bf = get_coloc_list(clats_old_bf, clons_old_bf, ctimes_old_bf)
    old_bf_fname = f"{save_dir}/{day}_{mw_sat_name}_{ro_sat_name}_{time_tolerance}_brute_force.csv"
    save_single_method_results(old_bf, old_bf_fname)
    logging.info(old_bf_s.strip())

    start_new_brute_force = time.time()
    colocs_new_bf, clats_new_bf,clons_new_bf,ctimes_new_bf = brute_force(mw_satellite, ro_satellite,time_tolerance, spatial_tolerance, mw_data_sorted, mw_times, ro_times, ro_lats, ro_lons, use_sorting=True)
    end_new_brute_force = time.time()
    new_bf_time = end_new_brute_force-start_new_brute_force
    new_bf_s = f"Brute force method #2 for colocations between {mw_sat_name} (MW) and {ro_sat_name} (RO) on day {day}, with time_tol={time_tolerance} completed in {new_bf_time} sec and found {colocs_new_bf} colocations \n"
    new_bf = get_coloc_list(clats_new_bf, clons_new_bf, ctimes_new_bf)
    new_bf_fname = f"{save_dir}/{day}_{mw_sat_name}_{ro_sat_name}_{time_tolerance}_brute_force_sorted.csv"
    save_single_method_results(new_bf, new_bf_fname)
    logging.info(new_bf_s.strip())
    bf_colocs = get_coloc_list(clats_old_bf, clons_old_bf, ctimes_old_bf)

    N=rtcd_N_1
    first_rtcd_N = N
    start_old_rtcd = time.time()
    colocs, clats,clons,ctimes,_,_ = RTCD(mw_satellite, ro_satellite,time_tolerance, spatial_tolerance, N, ro_times, ro_lats, ro_lons)
    end_old_rtcd = time.time()
    old_rtcd_time = end_old_rtcd-start_old_rtcd
    old_rot = get_coloc_list(clats, clons, ctimes)
    old_rot_fname = f"{save_dir}/{day}_{mw_sat_name}_{ro_sat_name}_{time_tolerance}_rotation_suboccs_N_{rtcd_N_1}.csv"
    save_single_method_results(old_rot, old_rot_fname)
    fp_old_rot, fn_old_rot = get_rot_false_pos_neg_timefails_spatialfails(old_rot, bf_colocs)
    num_bf_colocs = len(bf_colocs)
    if num_bf_colocs > 0:
        old_acc = (num_bf_colocs-len(fp_old_rot)-len(fn_old_rot))*100.0/num_bf_colocs
    elif num_bf_colocs== 0 and len(fp_old_rot) == 0 and len(fn_old_rot) == 0:
        old_acc = 100.0
    else:
        old_acc = 0.0
    old_rtcd_colocs = colocs
    old_rtcd_s = f"Rotation method with sub-occultations ran with {first_rtcd_N} sub-occultations in {old_rtcd_time} seconds and found {colocs} colocations, with {len(fp_old_rot)} false positives and {len(fn_old_rot)} false negatives for an accuracy of {old_acc}%"
    logging.info(old_rtcd_s.strip())

    N=rtcd_N_2
    second_rtcd_N = N
    start_new_rtcd = time.time()
    colocs, clats,clons,ctimes,_,_ = RTCD(mw_satellite, ro_satellite,time_tolerance, spatial_tolerance, N, ro_times, ro_lats, ro_lons)
    end_new_rtcd = time.time()
    new_rtcd_time = end_new_rtcd-start_new_rtcd
    new_rot = get_coloc_list(clats, clons, ctimes)
    new_rot_fname = f"{save_dir}/{day}_{mw_sat_name}_{ro_sat_name}_{time_tolerance}_rotation_linearized.csv"
    save_single_method_results(new_rot, new_rot_fname)
    fp_new_rot, fn_new_rot = get_rot_false_pos_neg_timefails_spatialfails(new_rot, bf_colocs)
    num_bf_colocs = len(bf_colocs)
    if num_bf_colocs > 0:
        new_acc = (num_bf_colocs-len(fp_new_rot)-len(fn_new_rot))*100.0/num_bf_colocs
    elif num_bf_colocs== 0 and len(fp_new_rot) == 0 and len(fn_new_rot) == 0:
        new_acc = 100.0
    else:
        new_acc = 0.0
    new_rtcd_colocs = colocs
    new_rtcd_s = f"Linearized rotation method ran in {new_rtcd_time} seconds and found {colocs} colocations, with {len(fp_new_rot)} false positives and {len(fn_new_rot)} false negatives for an accuracy of {new_acc}%"
    logging.info(new_rtcd_s.strip())

    with open(summary_fname, "w") as f:
        f.write(ro_load_s)
        f.write(mw_load_s)
        f.write(mw_sorted_load_s)
        f.write(old_bf_s)
        f.write(new_bf_s)
        f.write(old_rtcd_s)
        f.write(new_rtcd_s)

def main():
    args = parse_args()
    retrieve_daily_colocations(args)

if __name__ == "__main__":
    main()
