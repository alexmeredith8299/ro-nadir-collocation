import ssl
import os
import sys
import argparse
import textwrap

def parse_args():
    """
    Parse command line arguments specifiyng the data to download.

    Returns
    --------
        args: argparse.Namespace
            parsed command line arguments
    """
    parser = argparse.ArgumentParser(description='Download RO satellite data', formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
            Example usage:
                python3 get_ro_data.py --sat metopb --start 21001 --end 21031 --postProc 1 --level 2 --format wetPrf --filepath /Users/alexmeredith/ro-data/MetOpB-GRAS-Data

            The example above will download occultation data for Metop-B for the entire month of 2021 from UCAR. This will only work on Unix systems where `tar`, `mv`, `wget`, and `cp` are valid system calls.
            '''))
    parser.add_argument('--sat', help='RO satellite name, should be cosmic2, metopb, or metopc', required=True)
    parser.add_argument('--start', type=int, help='Start date, should fall between 21001 and 21031', required=True)
    parser.add_argument('--end', type=int, help='End date, should fall between 21001 and 21031', required=True)
    parser.add_argument('--postProc', type=int, help='0 for near real-time date (nrt), 1 for post-processed data (postProc)', required=True)
    parser.add_argument('--level', help='Level of data, should be 1b or 2', required=True)
    parser.add_argument('--format', help='UCAR data format (ex. wetPrf)', required=True)
    parser.add_argument('--filepath', help='Path to folder where satellite data is or should be stored', required=True)
    parser.add_argument('--download', dest='download', action='store_const', const=download, default=do_nothing)
    parser.add_argument('--unpack', dest='unpack', action='store_const', const=unpack, default=do_nothing)
    args = parser.parse_args()

    return args

def check_args(sat, start, end, postProc, level, format, filepath):
    """
    Validate args and raise an exception if anything is invalid.

    Arguments 
    ----------
        sat: str
            Satellite name
        start: int
            Start date
        end: int
            End date
        postProc: int
            Whether or not to use post-processed data or near-real-time retrievals
        level: str
            Level of data (1a, 1b, 2)
        format: str
            UCAR data format (ex. wetPrf or atmPrf)
        filepath: str
            Path to save data. Should be absolute path.
    """
    if start < 21001 or start > 21031:
        raise ValueError("Start date given by START must be between 21001 and 21031")
    if end < 21001 or end > 21031:
        raise ValueError("End date given by END must be between 21001 and 21031")
    if end < start:
        raise ValueError("End date given by END is before start date given by START")
    try:
        os.system(f"cd {filepath}")
    except:
        raise OSError("Cannot open RO data directory given by FILEPATH")
    if sat != 'cosmic2' and sat != 'metopb' and sat != 'metopc':
        raise ValueError("RO sat given by SAT must be equal to 'metopb', 'metopc', or 'cosmic2'")
    if level != '1a' and level !='1b' and level!='2':
        raise ValueError("Data level given by LEVEL must be equal to '1a', '1b', or '2'")
    if postProc != 0 and postProc != 1:
        raise ValueError("POSTPROC must be 1 (if postproc) or 0 (if nrt)")

def download(sat, start, end, postProc, level, format, output_path):
    """
    Downloads data.

    Arguments 
    ----------
        sat: str
            Satellite name
        start: int
            Start date
        end: int
            End date
        postProc: int
            Whether or not to use post-processed data or near-real-time retrievals
        level: str
            Level of data (1a, 1b, 2)
        format: str
            UCAR data format (ex. wetPrf or atmPrf)
        output_path: str
            Path to save data. Should be absolute path.
    """
    nrt_or_postProc = 'nrt'
    if postProc:
        nrt_or_postProc = 'postProc'
    print("Downloading")
    for date in (range(start, end + 1)):
        date_stem = str(date)[2:5]
        url = f"https://data.cosmic.ucar.edu/gnss-ro/{sat}/{nrt_or_postProc}/level{level}/2021/{date_stem}/{format}_{nrt_or_postProc}_2021_{date_stem}.tar.gz"
        filename = f"{output_path}/{format}_{nrt_or_postProc}_2021_{date_stem}.tar"
        try:
            print(f"Downloading data for day {date} from {url}...")
            ssl._create_default_https_context = ssl._create_stdlib_context
            os.system(f"wget -O {filename} --show-progress {url}")
        except:
            raise RuntimeError(f"Could not download data. Please navigate to {url} in a browser to make sure the file {filename} exists")

def unpack(sat, start, end, postProc, level, format, filepath):
    """
    Unpack .tar archive containing data.

    Arguments 
    ----------
        sat: str
            Satellite name
        start: int
            Start date
        end: int
            End date
        postProc: int
            Whether or not to use post-processed data or near-real-time retrievals
        level: str
            Level of data (1a, 1b, 2)
        format: str
            UCAR data format (ex. wetPrf or atmPrf)
        filepath: str
            Path to save data. Should be absolute path.
    """
    print("Unpacking")
    nrt_or_postProc = 'nrt'
    if postProc:
        nrt_or_postProc = 'postProc'
    os.system(f"mkdir {filepath}/Tarballs")
    for date in (range(start, end + 1)):
        if not os.path.exists(f"{filepath}/{str(date)}"):
            os.mkdir(f"{filepath}/{str(date)}")
        print(f"Unpacking data for day {start}...")
        date_stem = str(date)[2:5]
        filename = f"{format}_{nrt_or_postProc}_2021_{date_stem}.tar"
        os.system(f"mv {filepath}/{filename} {filepath}/{str(date)}")
        os.system(f"cd {filepath}/{date} && tar -xvzf {filepath}/{date}/{filename}")
        os.system(f"cd {filepath}/{date} && mv {filepath}/{date}/*.tar {filepath}/Tarballs")

def do_nothing(sat, start, end, postProc, level, format, filepath):
    pass

def main():
    """
    Parse args, download data, then unpack data.
    """
    args = parse_args()
    check_args(args.sat, args.start, args.end, args.postProc, args.level, args.format, args.filepath)
    download(args.sat, args.start, args.end, args.postProc, args.level, args.format, args.filepath)
    unpack(args.sat, args.start, args.end, args.postProc, args.level, args.format, args.filepath)

if __name__ == "__main__":
    main()
