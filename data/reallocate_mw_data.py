import os
import logging
import argparse
import textwrap
import datetime

def parse_args():
    """
    Parse command line arguments specifying the data to move to folders.

    Returns
    --------
        args: argparse.Namespace
            parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Find collocations between an RO and MW satellite', formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
            Example usage:

            python3 reallocate_mw_data.py --data_dir /Users/alexmeredith/ro-nadir-collocation/SNPP-Data --start_day 21001 --end_day 21031 --time_tol 600 --sat_type jpss 

            The example above creates folders called `21001`...`21031` and allocates SNPP data for each day into the appropriate folder (+ data from 600 seconds before and after the day in question).
            '''))
    parser.add_argument('--data_dir', help='Path to data to move', required=True)
    parser.add_argument('--start_day', help='Day in TLE format (ex. 21001)', required=True)
    parser.add_argument('--end_day', help='Day in TLE format (ex. 21031)', required=True)
    parser.add_argument('--time_tol', type=int, help='Time tolerance', required=True)
    parser.add_argument('--sat_type', help='Satellite type. Supported types are (`jpss` and `metop`)', required=True)
    return parser.parse_args()

def yyddd_to_datetime(yy_jd):
    """
    Convert day in format {yy}{ddd} (e.g. 21001 = January 1, 2021)
    to Python datetime object.
    Arguments 
    -----------
        yy_jd: int or str 
            julian date with 2-digit year
    Returns 
    --------
        dt: datetime
            datetime.datetime corresponding to 0:00:00 on the represented date
    """
    year,jd = int(str(yy_jd)[0:2]), int(str(yy_jd)[2:])
    dt = datetime.datetime(year+2000, 1, 1)+datetime.timedelta(days=jd-1)
    return dt

def year_and_day_of_year_to_jd(year, day_of_year):
    """
    Converts a year and day of year to a yyddd formatted Julian date.
    Example usage:
        year = 2021 and day_of_year = 1  -> 21001
        year = 1999 and day_of_year = 2 -> year out of range error
        year = 2000 and day_of_year = 366 -> 00366
    Arguments 
    ----------
        year: int
            year (will throw error if not between 2000 and 2099)
        day_of_year: int
            day of year
    """
    if year < 2000 or year > 2099:
        raise ValueError('Year must be between 2000 and 2099 (inclusive).')
    year_str = str(year)[2:4]
    day_of_year_str = str(day_of_year)
    while len(day_of_year_str) < 3:
        day_of_year_str = '0' + day_of_year_str
    return year_str + day_of_year_str

def datetime_to_yyddd(dt, buffer):
    """
    Given a datetime.datetime object, return all Julian dates within buffer 
    seconds of the datetime in yyddd format.
    Example usage:
        datetime.datetime(2021, 1, 1, 23, 50) with a buffer of 1200 seconds (20 min)
        matches Jan 1 2021 and Jan 2 2021, and would return ['21001', '21002']
    Arguments 
    ----------
        dt: datetime.datetime
            time object to classify
        buffer: int or float
            buffer in seconds
    Returns 
    -------
        jd_list: list
            list of strings representing matching Julian dates
    """
    jd_set = set() 

    #Ensure that we sample at least one datetime per day between -buffer/86400 to buffer/86400
    buffer_days = int(buffer/86400)
    dt_list = [dt - datetime.timedelta(days=x) for x in range(-buffer_days, buffer_days)]

    #Get endpoints of time window from dt-buffer to dt + buffer
    dt_list.append(dt + datetime.timedelta(seconds=buffer))
    dt_list.append(dt - datetime.timedelta(seconds=buffer))

    for dt_sample in dt_list:
        day_of_year = dt_sample.timetuple().tm_yday
        year = dt_sample.year

        jd = year_and_day_of_year_to_jd(year, day_of_year)
        jd_set.add(jd)

    return sorted(list(jd_set))

def get_end_of_prev_data(root_path, day, prev_day, prev_date, extra_min, verbose=False):
    """
    Get all files that contain data within extra_min of the
    start of the current day.
    Arguments
    ----------
        root_path: str
            root directory
        day: str
            current day
        prev_day: str
            previous day
        extra_min: int
            extra minutes to add to the start of the current day
    Returns
    -------
        list of files that contain data within extra_min of the
        start of the current day
    """
    files = []
    for file in os.listdir(os.path.join(root_path, prev_day)):
        file_tokens = file.split('_')
        date = file_tokens[2]
        hh = file_tokens[4][1:3]
        mm = file_tokens[4][3:5]
        if date == prev_date and 60 - int(mm) + 60 * (23-int(hh)) < extra_min:
            files.append(os.path.join(root_path, prev_day, file))
    if verbose:
        for file in files:
            logging.info(f"Moving {file} to {day}")
    return files

def get_start_of_next_data(root_path, day, next_day, next_date, extra_min, verbose=False):
    """
    Get all files that contain data within extra_min of the
    end of the current day.
    Arguments
    ----------
        root_path: str
            root directory
        day: str
            current day
        next_day: str
            previous day
        extra_min: int
            extra minutes to add to the start of the current day
    Returns
    -------
        list of files that contain data within extra_min of the
        start of the current day
    """
    files = []
    for file in os.listdir(os.path.join(root_path, next_day)):
        file_tokens = file.split('_')
        hh = file_tokens[3][1:3]
        mm = file_tokens[3][3:5]
        date = file_tokens[2]
        if date == next_date and int(mm) + 60 * (int(hh)) < extra_min:
            files.append(os.path.join(root_path, next_day, file))
    if verbose:
        for file in files:
            logging.info(f"Moving {file} to {day}")
    return files

def metop_fname_to_start_time(fname):
    """
    Transforms a filename for a MetOp sounding file into a datetime 
    representing the time of the first sounding in the file.
    MetOp sounding files have the following format:
    `W_XX-EUMETSAT-Darmstadt, SOUNDING+SATELLITE, METOP{A, B, C}+
        AMSUA_C_EUMP_{yyyy}{mm}{dd}{hh}{mm}{ss}_{ffff}_eps_o_l1.nc`
    where fffff is the five-digit file number.
    Example usage:
     `W_XX-EUMETSAT-Darmstadt, SOUNDING+SATELLITE, METOPC+
        AMSUA_C_EUMP_20210119223423_43278_eps_o_l1.nc` ->
        datetime.datetime(2021, 1, 19, 22, 34, 23)
    Arguments 
    ----------
        fname: str
            name of a MetOp sounding file in the specified format
    Returns 
    --------
        dt: datetime.datetime
            datetime object corresponding to the first sounding time 
                in the file
    """
    try:
        fname_tokens = fname.split('_')
        date_info = fname_tokens[4]
        year = int(date_info[0:4])
        month = int(date_info[4:6])
        day = int(date_info[6:8])
        hour = int(date_info[8:10])
        minute = int(date_info[10:12])
        second = int(date_info[12:14])
    except:
        raise ValueError(f'File {fname} was not in expected MetOp data format and could not be parsed.')

    return datetime.datetime(year, month, day, hour, minute, second)

def jpss_fname_to_start_time(fname):
    """
    Transforms a filename for a GATMO file containing soundings for
    a satellite in the JPSS program into a datetime representing the
    time of the first sounding in the file.
    GATMO files have the following format:
    `GATMO_{nnn}_d{yyyy}{mm}{dd}_t{hh}{mm}{sss}_e{hh}{mm}{sss}
        _b{ffff}_c{gggggggggggggggggggg}_nobc_ops.h5`
    where ffff and ggg.. are irrelevant numbers for our purposes. nnn 
    is a code identifying the satellite and is j01 for NOAA-20 and npp
    for Suomi-NPP. The second hh-mm-ss time corresponds to the last 
    sounding time and isn't used in this function.
    Example usage:
        `GATMO_npp_d20210101_t2356346_e0004343_b47579
            _c20210102040435288988_nobc_ops.h5` ->
        datetime.datetime(2021, 1, 1, 23, 56, 35)
    Arguments 
    ----------
        fname: str
            name of a GATMO file in the specified format
    Returns 
    ---------
        dt: datetime.datetime
            datetime object corresponding to the first sounding time 
                in the file
    """
    try:
        fname_tokens = fname.split('_')
        date_info = fname_tokens[2]
        year = int(date_info[1:5])
        month = int(date_info[5:7])
        day = int(date_info[7:9])
        time_info = fname_tokens[3]
        hour = int(time_info[1:3])
        minute = int(time_info[3:5])
        second = int(round(int(time_info[5:7]) + int(time_info[7])/10))
        if second > 59:
            second = 59
    except:
        raise ValueError(f'File {fname} was not in expected JPSS data format and could not be parsed.')

    return datetime.datetime(year, month, day, hour, minute, second)

def group_sounding_files_by_jd(sounding_files, min_jd, max_jd, buffer, fname_to_time_func):
    """
    Given a list of sounding files, a min and max Julian date, and a 
    time buffer in seconds, find a file-to-Julian-date matchup such that 
    each day from min_jd to max_jd (inclusive) corresponds to the minimal
    list of sounding files fully covering the time window (day start - buffer 
    sec, day end + buffer sec).
    It is assumed that there are no gaps between files, so the end time of a file 
    is assumed to be the same as the start time of the next file.
    Example:
        ['W_XX-..._20210101235000', 'W_XX-...20210102010000'], 21001, 21002, 600 
        yields '21001' matched to the first file and '21002' matched to both files
    Arguments 
    ----------
        sounding_files: list of str
            list of sounding files to consider
        min_jd: str 
            minimum Julian date to match files to (files before this date won't be 
                allocated)
        max_jd: str
            maximum Julian date to match files to (files after this date won't be
                allocated)
        buffer: int
            time buffer in seconds
        fname_to_time_func: function
            function mapping str -> datetime which converts the start time of a file 
                to a datetime given the file name
    Returns 
    --------
        jd_to_files: dict
            dict with str reps of Julian dates as keys and lists of file names as values
    """
    jds_to_files = {}
    min_jd_int = int(min_jd)
    max_jd_int = int(max_jd)
    sounding_files = sorted(sounding_files)
    for file in sounding_files:
        jds = datetime_to_yyddd(fname_to_time_func(file), buffer)
        for jd in jds:
            if int(jd) >= min_jd_int and int(jd) <= max_jd_int:
                if jd in jds_to_files:
                    jds_to_files[jd].append(file)
                else:
                    jds_to_files[jd] = [file]

    #Need to check if file before first file in day ends after (day-buffer)
    #(same as first file starts after day-buffer)
    for jd in jds_to_files:
        file_list = sorted(jds_to_files[jd])
        first_file_ind = sounding_files.index(file_list[0])
        if first_file_ind > 0:
            first_file_time = fname_to_time_func(sounding_files[first_file_ind])
            if (yyddd_to_datetime(jd) - first_file_time).seconds <= buffer:
                file_list.append(sounding_files[first_file_ind-1])
                jds_to_files[jd] = sorted(file_list)
    return jds_to_files

def allocate_files(root_path, min_jd, max_jd, buffer, fname_to_time_func):
    """
    Wrapper around group_sounding_files_by_jd. Given a folder containing 
    sounding files, find a jd -> file mapping using group_sounding_files_by_jd,
    and for each jd create a folder and copy all sounding files corresponding 
    to that jd into the jd folder.
    Example usage:
        allocate_files('/Users/alexmeredith/gpsro-data/NOAA20-Data/raw', '21001', 
            '21031', 10800, jpss_fname_to_start_time)
    This will create folders 21001...21031 in `raw` and fill them appropriately.
    Arguments
    ----------
        root_path: str
            path to folder containing sounding files
        min_jd: str
            minimum Julian date to match files to
        max_jd: str
            maximum Julian date to match files to
        buffer: int
            time buffer in seconds
        fname_to_time_func: function
            function mapping str -> datetime which converts the start time of a file
                to a datetime given the file name
    """
    sounding_files = os.listdir(root_path)
    jds_to_files = group_sounding_files_by_jd(sounding_files, min_jd, max_jd, buffer, fname_to_time_func)
    for jd in jds_to_files:
        os.makedirs(f'{root_path}/{jd}', exist_ok=True)
        for file in jds_to_files[jd]:
            os.system(f'cp {root_path}/{file} {root_path}/{jd}/{file}')

def main():
    args = parse_args()
    name_to_time = jpss_fname_to_start_time if args.sat_type == 'jpss' else metop_fname_to_start_time if args.sat_type == 'metop' else None
    if name_to_time == None:
        raise ValueError(f'Satellite type {args.sat_type} not supported. Only `jpss` and `metop` are supported.')
    allocate_files(args.data_dir, args.start_day, args.end_day, args.time_tol, name_to_time)

if __name__ == "__main__":
    main()

