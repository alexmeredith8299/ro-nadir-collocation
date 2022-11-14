# Efficient Collocation of GNSS Radio Occultation Soundings with Passive Nadir Microwave Soundings
This repository is associated with the paper "Efficient Collocation of GNSS Radio Occultation Soundings with Passive Nadir Microwave Soundings". It contains code and [instructions on how to download open-source data](#data) associated with the paper.

## Software

All code relevant to the paper is in the [scripts](https://github.com/alexmeredith8299/ro-nadir-collocation/tree/main/scripts) folder. The brute-force algorithm implementation is in `brute_force_clean.py`. The rotation-collocation algorithm implementation is in `rotation_RO.py` and depends on the satellite classes defined in `satclass.py`. All constants and one utility function related to calculating the latitude-dependent radius of Earth are in `constants_and_utils.py`. The code to retrieve a day's worth of collocations is in `joint_retrievals.py`.

### Setup and Dependencies

The code has some basic dependencies, namely the `sgp4` and `astropy` Python libraries for orbit propagation and coordinate conversions, `numpy` for various mathematical operations, `h5py` for reading HDF5 files, and `netCDF4` for reading netCDF files. To install the dependencies:

1. Install Python with a version >= 3.9 if not already installed
2. Install [pip](https://pypi.org/project/pip/) if not already installed, for package management
3. Run `python3 -m pip install -r requirements.txt` to install all dependencies

### Example Retrievals

Assuming that Metop-B GRAS and Metop-B AMSU data for January 1, 2021 are in /path/to/data/MetOpB-GRAS-Data/21001 and /path/to/data/MetOpB-AMSUA-Data/21001 respectively, and that Metop-B TLEs are available in /path/to/tle/MetOpB-AMSUA.txt (see [this section](#organizing-data) for a guide to organizing data):

`python3 joint_retrievals.py --mw_sat MetOpB --ro_sat MetOpB --day 21001 --time_tol 600 --spatial_tol 150 --n_suboccultations 21 --save_dir /path/to/results --data_dir /path/to/data --tle_dir /path/to/tle`

will find all collocations between Metop-B GRAS and Metop-B AMSU for January 1, 2021 with both brute-force methods and both rotation-collocation methods with a time tolerance of 10 minutes (600 seconds) and a spatial tolerance of 150 km. 

For the rotation-collocation method with sub-occultations, 21 sub-occultations will be used. All results will be saved to /path/to/results -- four .csv results files will be saved with the latitude, longitude, and time of each collocated MetopB-GRAS occultation for each method, and one summary .txt file will be saved summarizing the accuracy and time taken by each method.

## Data

All data files used for this paper are available from [NOAA](https://www.avl.class.noaa.gov/saa/products/welcome), [EUMETSAT](https://www.eumetsat.int/eumetsat-data-centre), or [UCAR](https://data.cosmic.ucar.edu/gnss-ro/).

### RO Data

The [data](https://github.com/alexmeredith8299/ro-nadir-collocation/tree/main/data) folder contains two helper scripts for downloading and sorting data. The first script, `get_ro_data.py`, can be used to download all the COSMIC-2, Metop-B, and Metop-C RO data from UCAR. For example, the following line will download Metop-B data for the month of January 2021:

`python3 get_ro_data.py --sat metopb --start 21001 --end 21031 --postProc 1 --level 2 --format wetPrf --filepath /path/to/data/MetOpB-GRAS-Data`

### Microwave Radiance Data

To download microwave sounding data for NOAA-20 and SNPP, use the [NOAA CLASS data ordering system](https://www.avl.class.noaa.gov/saa/products/welcome) to order "JPSS ATMS Sensor Data Record Operational (ATMS_SDR)" data. In the "Advanced Search" section, select the "ATMS SDR Ellipsoid Geolocation (GATMO) (public 12/10/2011)" checkbox. Highlight both NOAA-20 and SNPP in the "Satellite" box and choose "Either" for "Node". 

To download microwave sounding data for Metop-B and Metop-C, use the [EUMETSAT Data Centre](https://www.eumetsat.int/eumetsat-data-centre) to order "AMSU-A Level 1B - Metop - Global" Data. Select "METOP-B" and "METOP-C" for the satellites of interest, and select "netCDF" for the data format.

After downloading microwave sounding data, use `reallocate_mw_data.py` to subdivide the microwave sounding data into folders for each day. For example, the following line will reallocate SNPP data for the month of January 2021 using a time tolerance of 10 minutes/600 seconds:

`python3 reallocate_mw_data.py --data_dir /path/to/data/SNPP-Data --start_day 21001 --end_day 21031 --time_tol 600 --sat_type jpss`

It's important that the time tolerance matches the time tolerance used for retrieving collocations -- otherwise, the brute-force method will either have extra data to consider and will run more slowly, or will not have all necessary data and will miss some collocations.

### TLEs

All TLEs used for this paper are available from [Celestrak](https://celestrak.org/). Simply order TLEs from Celestrak for the time period of interest and download them. 

### Organizing Data

This section is extremely important. The code in `satclass.py` expects that there will be some folder /path/to/data containing folders called "NOAA20-Data", "MetOpB-AMSUA-Data", "MetOpC-AMSUA-Data", and "SNPP-Data" with the microwave radiance data, and "COSMIC2-Data", "MetOpB-GRAS-Data", and "MetOpC-GRAS-Data" with the radio occultation data. Each of these satellite folders will contain subfolders for each day of data labeled "21001", "21002", etc.

It also expects a folder /path/to/tle containing TLEs for all the microwave satellites in files called "MetOpB-AMSUA.txt", "MetOpC-AMSUA.txt", "NOAA-20.txt", and "SNPP.txt".

Using a different folder structure will not work out-of-the-box and will require some light modifications to `satclass.py`.
