/*
Class for loading, storing and interpolating solar wind data.

Copyright 2010, 2011 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOLAR_WIND_HPP
#define SOLAR_WIND_HPP

#include "algorithm"
#include "boost/array.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "cmath"
#include "cstring"
#include "fstream"
#include "iostream"
#include "string"
#include "utility"
#include "vector"

namespace solar_wind {

/*!
Data type used by the Solar_wind class.
*/
typedef boost::array<double, 8> solar_wind_t;

/*!
Data is stored at these indices in solar_wind_t
*/
enum {
	density = 0,
	n = 0,
	N = 0,

	temperature = 1,
	T = 1,

	vx = 2,
	VX = 2,

	vy = 3,
	VY = 3,

	vz = 4,
	VZ = 4,

	Bx = 5,
	BX = 5,

	By = 6,
	BY = 6,

	Bz = 7,
	BZ = 7
};


/*!
Same as boost::posix_time::from_iso_string() but for an extended iso string.
*/
boost::posix_time::ptime from_iso_extended_string(const char* iso_extended_string) {
	// remove separators and use boost::from_iso_string
	std::string iso_string;
	int length = strlen(iso_extended_string);
	iso_string.reserve(length);
	int i = 0;
	while (i < length) {
		if (iso_extended_string[i] == '.') {
			break;
		}

		if (iso_extended_string[i] != '-' and iso_extended_string[i] != ':') {
			iso_string.push_back(iso_extended_string[i]);
		}

		i++;
	}

	boost::posix_time::ptime result(boost::posix_time::from_iso_string(iso_string));
	return result;
}

/*!
Same as the one for char* but for a string.
*/
boost::posix_time::ptime from_iso_extended_string(const std::string iso_extended_string) {
	return from_iso_extended_string(iso_extended_string.c_str());
}


/*!
Solar wind data is given in SI units and GSE coordinate system.
*/
class Solar_wind
{
public:

	/*!
	Initializes this class instance with an average solar wind near the Earth.
	*/
	Solar_wind() {
		boost::posix_time::ptime default_time(from_iso_extended_string("1970-01-01T00:00:00Z"));

		// assume spiral angle at the Earth is 45 degrees
		solar_wind_t default_data = { 6.5e6, 1e5, -4e5, 29, 0, -6e-9 * M_SQRT1_2, 6e-9 * M_SQRT1_2 };

		this->solar_wind_data.push_back(std::make_pair(default_time, default_data));
	}

	/*!
	Initializes the class instance with a solar wind from the given file, see the load function for details.
	*/
	Solar_wind(const char* filename) {
		this->load(filename);
	}

	/*!
	Same as the one for char* but for a string.
	*/
	Solar_wind(const std::string filename) {
		this->load(filename);
	}

	/*!
	Loads solar wind data from given file.

	One solar wind value consists of a date and time in extended format of ISO 8601 and density, temperature, velocity components and magnetic field components in SI units and GSE coordinates.
	Fractional seconds are ignored at the moment.
	For example:
	1970-01-01T00:00:00Z 6.5e6 1e5 -4e5 29 9 -4 4 0
	*/
	void load(const char* filename) {
		std::ifstream infile(filename, std::ifstream::in);
		if (!infile.good()) {
			std::cerr << "Error opening file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		this->solar_wind_data.clear();

		while (true) {

			std::string time_string;
			infile >> time_string;
			if (infile.eof()) {
				break;
			}

			if (!infile.good()) {
				std::cerr << "Error reading time from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			boost::posix_time::ptime time(from_iso_extended_string(time_string));

			// warn the user if this time is earlier than previous
			if (this->solar_wind_data.size() > 0 && this->solar_wind_data[this->solar_wind_data.size() - 1].first > time) {
				std::cout << "Given time " << boost::posix_time::to_iso_extended_string(time) << " is earlier than previous " << boost::posix_time::to_iso_extended_string(this->solar_wind_data[this->solar_wind_data.size() - 1].first) << std::endl;
			}

			solar_wind_t data;
			infile >> data[0];
			if (!infile.good()) {
				std::cerr << "Error reading desity from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[1];
			if (!infile.good()) {
				std::cerr << "Error reading temperature from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[2];
			if (!infile.good()) {
				std::cerr << "Error reading vx from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[3];
			if (!infile.good()) {
				std::cerr << "Error reading vy from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[4];
			if (!infile.good()) {
				std::cerr << "Error reading vz from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[5];
			if (!infile.good()) {
				std::cerr << "Error reading Bx from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[6];
			if (!infile.good()) {
				std::cerr << "Error reading By from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}
			infile >> data[7];
			if (!infile.good()) {
				std::cerr << "Error reading Bz from file " << filename << std::endl;
				exit(EXIT_FAILURE);
			}

			this->solar_wind_data.push_back(std::make_pair(time, data));
		}

		sort(this->solar_wind_data.begin(), this->solar_wind_data.end());
	}

	/*!
	Same as the one for char* but for a string.
	*/
	void load(const std::string filename) {
		this->load(filename.c_str());
	}


	/*!
	Returns the solar wind at given time.

	Solar wind data is interpolated linearly if between existing data or the first / last value is returned if given time is outside of the range of existing data.
	*/
	solar_wind_t get_solar_wind(const boost::posix_time::ptime given_time) {

		if (given_time < this->solar_wind_data[0].first) {
			return this->solar_wind_data[0].second;
		}

		if (given_time > this->solar_wind_data[this->solar_wind_data.size() - 1].first) {
			return this->solar_wind_data[this->solar_wind_data.size() - 1].second;
		}

		int i = 0;
		while (given_time > this->solar_wind_data[i].first) {
			i++;
		}

		boost::posix_time::time_duration before, after;
		before = given_time - this->solar_wind_data[i - 1].first;
		after = this->solar_wind_data[i].first - given_time;

		// TODO: implement interpolation
		if (after > before) {
			return this->solar_wind_data[i - 1].second;
		} else {
			return this->solar_wind_data[i].second;
		}
	}


private:

	std::vector<std::pair<boost::posix_time::ptime, solar_wind_t> > solar_wind_data;


};

} // namespace solar_wind

#endif

