import numpy as np
import os
from netCDF4 import Dataset
import logging
import time
import h5py
from astropy.coordinates import TEME, ITRS
from astropy import coordinates as coord
from constants_and_utils import *
from sgp4.ext import rv2coe
from sgp4.api import Satrec
import erfa
from astropy.coordinates.builtin_frames.utils import get_polar_motion
import math

class RO():
    """
    GNSS-RO satellite object, which represents an RO satellite on a specific day.

    Parameters
    -------------
        name: string
            Name of the satellite (ex. NOAA-20)
        day: string
            Day of year [in TLE form, ie, '20289' for the 289th day of 2020]

    Attributes
    -----------
        name: string
            Name of the satellite (ex. NOAA-20)
        day: string
            Day of year [in TLE form, ie, '20289' for the 289th day of 2020]
        data: string
            Filepath to occultation data
    """
    def __init__(self,name,day, data_dir):
        """
        Initializes RO satellite object (see class definition)
        """
        self.name = name
        self.day = day
        self.data = f'{data_dir}/{name}-Data/{day}/'

    def parse_occultation_file(self,filename):
        """
        This function parses data for a single occultation.

        Parameters
        -----------
            filename: string
                Path to file containing information about one occultation

        Returns
        --------
            occult_time: float
                Time of occultation point [Unix time, hours]
            occult_lat: float
                Latitude of perigee point at occultation point [Radians]
            occult_lon: float
                Longitude of perigee point at occultation point [Radians]
        """
        ncin = Dataset(self.data + filename, 'r', format='NETCDF3')
        occult_start_time = (ncin.start_time + gps_epoch_to_unix)*sec_to_hours
        occult_time = occult_start_time
        occult_lat = np.deg2rad(ncin.lat)
        occult_lon = np.deg2rad(ncin.lon)
        return occult_time, occult_lat, occult_lon

    def load_data(self):
        """
        This function loads times, latitudes, and longitudes for RO soundings
        from an RO satellite object.

        Returns
        --------
            times: list
                list of sounding times in Unix seconds
            lats: list
                list of sounding latitudes in radians
            lons: list
                list of sounding longitudes in radians
        """
        times, lats, lons = [], [], []
        for file in os.listdir(self.data):
            occult_time, occult_lat, occult_lon = self.parse_occultation_file(file)
            times.append(occult_time)
            lats.append(occult_lat)
            lons.append(occult_lon)
        return times, lats, lons

class MWR_G():
    """
    Generalized MWR satellite object. This represents a MW satellite with orbital properties and
    microwave soundings over a specific day, with a generalized orbit that isn't necessarily sun-synchronous.

    Parameters
    ------------
        name: string
            Name of satellite (ex. NOAA-20)
        inclination: float
            Orbital inclination [degrees]
        semi_major_axis: float
            Semi-major axis [km]
        max_scan_angle: float
            MWR instrument maximum scan angle [degrees]
        time_unit: string
            Time units of MWR data ['iet' or '1.1.2000']
        day: string
            Day of year [in TLE form, ie, '20289' for the 289th day of 2020]

    Attributes
    ------------
        name: string
            Name of satellite (ex. NOAA-20)
        a: float
            Semi-major axis [km]
        xi: float
            MWR instrument maximum scan angle [radians]
        time_unit: string
            Time units of MWR data ['iet' or '1.1.2000']
        day: string
            Day of year [in TLE form, ie, '20289' for the 289th day of 2020]
        tle: string
            Filepath to file containing TLEs
        data: string
            Filepath to sounding data
    """
    def __init__(self,name,semimajor_axis, max_scan_angle,time_unit,day, data_dir, tle_dir):
        """
        This function initializes the MWR_G object. See the class definition for further information.
        """
        self.name = name
        self.a = semimajor_axis
        self.xi = np.deg2rad(max_scan_angle)
        self.day = day

        self.tle = f'{tle_dir}/{name}.txt'
        self.data = f'{data_dir}/{name}-Data/{day}/'
        self.time_unit = time_unit

    def parse_data(self, filename):
        """
        Function to parse full MWR datasets for NOAA-20 and MetOp-C

        Parameters
        ------------
            filename: string
                MWR data filename
        Returns
        ---------
            longitudes: list
                List of sounding latitudes [radians]
            latitudes: list
                List of sounding longitudes [radians]
            record_start_times: list
                List of times (one time per cross track scan) [Unix Time seconds]
            sc_positions: list
                List of spacecraft positions [km] in ITRS coordinate system
            sc_velocities: list
                List of spacecraft velocities [km] in ITRS coordinate system
        """
        if self.name == 'NOAA20' or self.name=="SNPP":
            file = h5py.File(filename,'r')
            longitudes = np.deg2rad(list(file['All_Data/ATMS-SDR-GEO_All/Longitude']))
            latitudes = np.deg2rad(list(file['All_Data/ATMS-SDR-GEO_All/Latitude']))
            mid_times = self.time_conversion(np.array(file['All_Data/ATMS-SDR-GEO_All/MidTime']))

            #Avoid bad data (time = 0)
            good_times, good_lon, good_lat= [], [], []

            for i in range(len(mid_times)):
                if mid_times[i] == 0.0:
                    logging.warning(f"Skipping some data with no sounding time available for {self.name} on {self.day}")
                else:
                    good_times.append(mid_times[i])
                    good_lon.append(longitudes[i])
                    good_lat.append(latitudes[i])
            return good_lon, good_lat, good_times
        elif self.name == 'MetOpC-AMSUA' or self.name == 'MetOpB-AMSUA':
            test = Dataset(filename, 'r', format='NETCDF3')
            latitudes = np.deg2rad(test.variables['lat'][:])
            longitudes = np.deg2rad(test.variables['lon'][:])
            record_start_times = self.time_conversion(test.variables['record_start_time'][:])
            return longitudes, latitudes, record_start_times
        else:
            raise ValueError(f'ERROR: parsing data is only supported for NOAA-20 (`NOAA20`), Suomi-NPP (`SNPP`), Metop-B (`MetOpB-AMSUA`), and Metop-C (`MetOpC-AMSUA`). Satellite with name {self.name} is not supported.')

    def load_sorted_data(self):
        """
        This function loads a day's worth of soundings for a MW satellite and sorts them by time.

        Returns
        ---------
            file_list: list
                List of soundings (sorted by time, increasing). Each elt is a list that contains lon, lat, and sounding time
            mw_times: list
                List of MW times sorted by time, increasing
        """
        data_mw = os.listdir(self.data) #list of mw data files (each file represents multiple soundings)
        file_list = []
        mw_times = []
        for l in range(0,len(data_mw)): #iterate through each file
            f = self.data + str(os.listdir(self.data)[l])
            if self.name!='MetOpC-AMSUA' or f.split('/')[-1][0] != '.':
                longitudes, latitudes, mid_times = self.parse_data(f)
                for i in range(len(longitudes)):
                    file_list.append([longitudes[i], latitudes[i], mid_times[i]])
                    mw_times.append(mid_times[i])
        mw_times = sorted(mw_times)
        file_list = sorted(file_list, key=lambda x:x[2])
        return file_list, mw_times

    def load_data(self):
        """
        This function loads a day's worth of sounding data for a microwave satellite. Reading from files is slow!!!
        This speeds up brute force colocation search roughly 20-40x (informally tested on Macbook Pro).

        Returns
        ------------
            file_list: list
                List containing one entry for each file of MW data. Each entry contains a list of longitudes, latitudes,
                mid-sounding times, spacecraft positions, and spacecraft velocities.
        """
        data_mw = os.listdir(self.data) #list of mw data files (each file represents multiple soundings)
        file_list = []
        for l in range(0,len(data_mw)): #iterate through each file
            f = self.data + str(os.listdir(self.data)[l])
            if self.name!='MetOpC-AMSUA' or f.split('/')[-1][0] != '.':
                longitudes, latitudes, mid_times = self.parse_data(f)
                file_list.append((longitudes, latitudes, mid_times))
        return file_list

    def time_conversion(self,times):
        """
        Function to convert MWR sounding times to Unix Epoch Time

        Parameters
        ------------
            times: list
                List of MWR sounding times in IDPS Epoch Time (IET) [Microseconds since
                epoch of 1 Jan 1958, excludes leap seconds] or 1.1.2000 epoch time [seconds since
                1 Jan 2000, includes leap seconds]
        Returns
        ----------
            converted_times: list
                List of MWR sounding times in Unix time converted from seconds to hours,
                *excluding leap seconds*. This is equivalent to TAI time, but measured from
                the start of the Unix epoch, in hours.
                [hours since 1 Jan 1970, excludes leap seconds]
        """
        if self.time_unit == 'idps' or self.time_unit == 'iet':
            for i in range(len(times)):
                if times[i] > 0:
                    times[i] = times[i]*microsec_to_sec + tai_epoch_to_unix
                else:
                    times[i] = 0
                    logging.warning(f"Sounding time was <= 0")
        elif self.time_unit == '1.1.2000':
            for i in range(len(times)):
                times[i] = times[i] + j2000_epoch_to_unix + utc_to_tai 
        else:
            raise ValueError('Error: time units not defined. Only TAI time in microseconds (`iet` or `idps`) and J2000 epoch time in seconds (`1.1.2000`) are supported.')

        return times*sec_to_hours

    def get_current_tle(self, day, occultation_time):
        """
        Function to parse MWR TLE files into the relevant TLE for the current epoch

        Parameters
        ------------
            day: string
                Day of interest [in TLE form, ie, '20289' for the 289th day of 2020]
            occultation_time: float
                Time of RO sounding [Unix time, seconds]

        Returns
        ---------
            tle_1: string
                first line of TLE
            tle_2: string
                second line of TLE
        """
        if day[0:2] != '21':
            raise ValueError(f'TLE parsing is only valid for the year of 2021. Day {day} not supported.')
        file = open(self.tle,'r').read()
        tle_lines = file.strip().splitlines()
        lines = []
        for i in range(len(tle_lines)):
            if tle_lines[i].split()[0] == '1' and tle_lines[i].split()[3].split('.')[0] == day:
                lines.append(i)

        lines.append(lines[-1] + 2)
        for a in range(len(lines)):
            current_tle = tle_lines[lines[a]]
            tle_1 = tle_lines[lines[a]-2] #previous "1" TLE line
            tle_2 = tle_lines[lines[a]-1] #previous "2" TLE line
            current_tle_time = y2021_epoch_to_unix + (float(current_tle.split()[3][2:14])-1)*day_to_sec
            if occultation_time <= current_tle_time:
                return tle_1, tle_2 

        raise ValueError(f'ERROR: could not find TLE for day {day} and occultation time {occultation_time}')

    def calculate_ds(self, lat_geodetic, a=None):
        """
        This function calculates the latitude-dependent scan angle for use in the rotation method.

        Parameters
        -----------
            lat_geodetic: float
                Geodetic latitude [radians]
            a: float
                Altitude of satellite [km]. If None, constant semimajor axis is used.

        Returns
        ---------
            ds: float
                Latitude-dependent scan angle [radians]
        """
        return self.calculate_ds_from_xi(lat_geodetic, self.xi, a=a)

    def calculate_ds_from_xi(self, lat_geodetic, xi, a=None):
        """
        This function calculates the latitude-dependent scan angle for use in the rotation method.

        Parameters
        -----------
            lat_geodetic: float
                Geodetic latitude [radians]
            xi: float
                Scan distance [radians]
            a: float
                Altitude of satellite [km]. If None, constant semimajor axis is used.

        Returns
        ---------
            ds: float
                Latitude-dependent scan angle [radians]
        """
        if a == None:
            a = self.a
        r = calculate_radius_of_earth(lat_geodetic)
        ds = np.arcsin((a/r)*np.sin(xi)) - xi
        return ds

    def r_vector_from_jd_and_polarmat(self, jd, polarmat, lon, lat):
        """
        Function to find the inertial latitude and longitude of an occultation.

        Parameters
        -----------
            jd: float
                Julian date representing occultation time
            polarmat: np array
                Matrix representing polar motion
            lat: float
                Latitude of occultation [radians]
            lon: float
                Longitude of occultation [radians]
        Returns
        --------
            lon_inert: float
                Inertial longitude of occultation [radians]
            lat_inert: float
                Inertial latitude of occultation [radians]
        """

        #r vector in ECEF coordinates
        r_j = np.array([[np.cos(lon)*np.cos(lat)],[np.sin(lon)*np.cos(lat)],[np.sin(lat)]])

        #convert to TEME
        gst = erfa.gmst82(jd,0)
        teme_to_itrs_mat = erfa.c2tcio(np.eye(3), gst, polarmat)
        occ_TEME = teme_to_itrs_mat.T@r_j

        #Find inertial latitude and longitude from TEME
        lon_inert = np.arctan2(occ_TEME[1][0],occ_TEME[0][0])
        lat_inert = np.arcsin(occ_TEME[2][0])
        return lon_inert, lat_inert

    def calc_p_matrix_fast_RTCD(self, time_tolerance, N, ro_times, ro_lats, ro_lons):
        """
        Function to perform the RTCD-G rotational transformation and compute the p matrix
            for the given GNSS-RO satellite, by relaxing the assumption of contemporaneous
            RO and MW soundings.

        Parameters
        -----------
            time_tolerance: float
                Time matchup tolerance [seconds]
            N: int
                Number of pseudo-RO occultations to create (ex. 201 -- this will result in pseudo-soundings being created every 6 seconds if time_tolerance = 10 minutes)

        Returns
        -----------
            p_matrix: numpy array (m x N x 3 x 1)
                Position rotated into the MW frame for each occultation and pseudo-occultation
            times: list
                List of occultation times for GNSS-RO satellite object [Unix time hours]
            lats: list
                List of occultation latitudes for GNSS-RO satellite object [Radians]
            lons: list
                List of occultation longitudes for GNSS-RO satellite object [Radians]
            rises: list
                List of whether occultations are rising or setting (1 for setting, -1 for rising)
        """
        times, lats, lons, rises = ro_times, ro_lats, ro_lons, []
        if N == 1:
            dt = 0
        else:
            dt = 2*time_tolerance/(N-1)

        count = len(ro_times) #Number of RO soundings
        p_matrix_big = np.zeros((count, N, 3, 1))

        s, tline = "", ""
        polarmat = np.zeros((3,3))
        for i in range(len(times)):
            occult_time, occult_lat, occult_lon = times[i], lats[i], lons[i]

            #Convert the time (ensuring that you account for conversion to UTC)
            p = time.gmtime(occult_time*hours_to_sec+tai_to_utc)
            jd1_start, jd2_start = erfa.dtf2d("UTC", p.tm_year, p.tm_mon, p.tm_mday,
                                          p.tm_hour, p.tm_min, p.tm_sec)
            jd = jd1_start+jd2_start

            #Get polar motion (only once per day)
            if i == 0:
                xp, yp = get_polar_motion(jd)
                polarmat = erfa.pom00(xp, yp, 0)

            #Find inertial lat/long for occultation
            lon_inert, lat_inert = self.r_vector_from_jd_and_polarmat(jd, polarmat, occult_lon, occult_lat)
            sidereal_spin_per_sec = 2*np.pi*sec_to_sidereal_day#86164.0905

            #Find TLE (only really need one per day)
            if i == 0:
                s, tline = self.get_current_tle(self.day,occult_time*hours_to_sec+tai_to_utc)
                sat = Satrec.twoline2rv(s, tline)

            #Use sgp4 to get orbit position and then get new orbital elements
            e, r, v = sat.sgp4(jd, 0)
            (p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper) = rv2coe(r, v, mu)

            #Loop over k=0...N-1 to create pseudo-occultations
            for k in range(0, N):
                #find tjk
                tjk = occult_time*hours_to_sec + (k-((N-1)/2))*dt
                p = time.gmtime(tjk+tai_to_utc)
                jd1_start, jd2_start = erfa.dtf2d("UTC", p.tm_year, p.tm_mon, p.tm_mday,
                                          p.tm_hour, p.tm_min, p.tm_sec)
                jd = jd1_start + jd2_start
                #Use sgp4 to get orbit position and then get new orbital elements
                e, r, v = sat.sgp4(jd, 0)
                (p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper) = rv2coe(r, v, mu)
                #find lonjk
                lonjk = lon_inert+sidereal_spin_per_sec*(tjk-occult_time*hours_to_sec)
                #Try setting I and omega to new value from SGP4
                I = incl
                omega = node
                #Construct p matrix
                incl_matrix = np.array([[1, 0, 0],
                                    [0, np.cos(I), np.sin(I)],
                                    [0, -np.sin(I), np.cos(I)]])
                raan_matrix = np.array([[np.cos(lonjk - omega)*np.cos(lat_inert)],
                                [np.sin(lonjk - omega)*np.cos(lat_inert)],
                                [np.sin(lat_inert)]])
                arglat_matrix = np.array([[np.cos(argp+nu), np.sin(argp+nu), 0],
                                          [-np.sin(argp+nu), np.cos(argp+nu), 0],
                                          [0, 0, 1]])
                p_matrix = arglat_matrix @ incl_matrix @ raan_matrix 

                p_matrix_big[i, k, :, :] = p_matrix
        return p_matrix_big, times, lats, lons,rises
