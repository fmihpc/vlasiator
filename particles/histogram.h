#pragma once
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "vectorclass.h"
#include "vector3d.h"

#define ERROR(format, ...) fprintf (stderr, "E: " format, ##__VA_ARGS__)

/* Histograms of 1D data */
class Histogram1D {
	protected:
		size_t num_bins;
		double* bins;

	public:
		Histogram1D(size_t n) : num_bins(n) {
			bins = new double[num_bins];
			memset(bins, 0, sizeof(double)*num_bins);
		}
		~Histogram1D() {
			delete[] bins;
		}

		/* Using MPI_Reduce, sum up all CPUs. */
		void mpi_reduce() {
			double* targetbins = new double[num_bins];
			if(!targetbins) {
				ERROR("allocation failed while mpi-reducing histogram.\n");
				return;
			}

			MPI::COMM_WORLD.Allreduce(bins,targetbins,num_bins,MPI::DOUBLE,MPI_SUM);

			/* Throw away the old bins. */
			delete[] bins;
			bins = targetbins;
		}

		void save(const char* filename) const;
		void load(const char* filename);
		virtual void saveAscii(const char* filename) const;
		virtual void addValue(double value) = 0;
};

class LinearHistogram1D : public Histogram1D {

	private:
		/* Low and high bound of the histogram */
		double low, high;

	public:
		LinearHistogram1D(size_t n, double _low, double _high) :
			Histogram1D(n), low(_low), high(_high) {};

		virtual void addValue(double value) {
			value -= low;
			value /= high-low;

			int histogram_bin = value * num_bins;
			
			if(histogram_bin < 0) {
				histogram_bin = 0;
			} else if (histogram_bin+1 > (ssize_t)num_bins) {
				histogram_bin = num_bins-1;
			}
			bins[histogram_bin]++;
		}

		virtual void saveAscii(const char* filename) const;
};

class LogHistogram1D : public Histogram1D {
	private:
		/* low and high bound of the histogram */
		double low, high;
		
	public:
		LogHistogram1D(size_t n, double _low, double _high) :
			Histogram1D(n), low(_low), high(_high) {};

		virtual void addValue(double value) {
			value /= low;
			value = log(value);
			value /= log(high/low);

			int histogram_bin = value * num_bins;
			
			if(histogram_bin < 0) {
				histogram_bin = 0;
			} else if (histogram_bin+1 > (ssize_t)num_bins) {
				histogram_bin = num_bins-1;
			}
			bins[histogram_bin]++;
		}
};


/* Histograms of 2D data */
class Histogram2D {
	protected:
		size_t num_bins[2];
		double* bins;

	public:
		Histogram2D(size_t n[2]) {
			num_bins[0] = n[0];
			num_bins[1] = n[1];

			bins = new double[num_bins[0]*num_bins[1]];
			memset(bins, 0, sizeof(double)*num_bins[0]*num_bins[1]);
		}
		~Histogram2D() {
			delete[] bins;
		}

		/* Using MPI_Reduce, sum up all CPUs. */
		void mpi_reduce() {
			double* targetbins = new double[num_bins[0]*num_bins[1]];
			if(!targetbins) {
				ERROR("allocation failed while mpi-reducing histogram.\n");
				return;
			}

			MPI::COMM_WORLD.Allreduce(bins,targetbins,num_bins[0]*num_bins[1],MPI::DOUBLE,MPI_SUM);

			/* Throw away the old bins. */
			delete[] bins;
			bins = targetbins;
		}

		void save(const char* filename) const;
		void load(const char* filename);
		virtual void addValue(Vec2d value, double weight=1.) = 0;

		// Access bins
		double operator()(int x, int y) {
			return bins[x + num_bins[0]*y];
		}
};

class LinearHistogram2D : public Histogram2D {

	private:
		/* Low and high bound of the histogram */
		Vec2d low, high;

	public:
		LinearHistogram2D(size_t n[2], Vec2d _low, Vec2d _high) :
			Histogram2D(n),
			low(_low), high(_high) {};
		LinearHistogram2D(size_t nx, size_t ny, Vec2d _low, Vec2d _high) :
			Histogram2D(&nx),
			low(_low), high(_high) {};

		virtual void addValue(Vec2d value, double weight=1.) {
			value -= low;
			value /= high-low;

			int histogram_bin[2];
			histogram_bin[0] = value[0] * num_bins[0];
			histogram_bin[1] = value[1] * num_bins[1];

			for(int i=0; i<2;i++) {
				if(histogram_bin[i] < 0) {
					histogram_bin[i] = 0;
				} else if (histogram_bin[i]+1 > (ssize_t)num_bins[i]) {
					histogram_bin[i] = num_bins[i]-1;
				}
			}
			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1]]+= weight;
		}

	  void addValueLinearInterpolate(Vec2d value, double weight=1.) {

			value -= low;
			value /= high-low;

			double v[2];
			v[0]=value[0]*num_bins[0];
			v[1]=value[1]*num_bins[1];

			// Interpolation parameter
			double a[2];
			a[0] = v[0] - floor(v[0]);
			a[1] = v[1] - floor(v[1]);

			// Ensure we are within bounds.
			if(v[0] < 0) {
				v[0] = 0;
				a[0] = 0;
			}
			if(v[1] < 0) {
				v[1] = 0;
				a[1] = 0;
			}

			if(v[0] > num_bins[0]-1) {
				v[0] = num_bins[0]-2;
				a[0] = 1;
			}
			if(v[1] > num_bins[1]-1) {
				v[1] = num_bins[1]-2;
				a[1] = 1;
			}

			// Assign.
			int histogram_bin[2];
			histogram_bin[0] = floor(v[0]);
			histogram_bin[1] = floor(v[1]);

			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1]]+= weight * (1.-a[0]) * (1.-a[1]);
			bins[histogram_bin[0] + num_bins[0]*(histogram_bin[1]+1)]+= weight * (1.-a[0]) * a[1];
			bins[histogram_bin[0]+1 + num_bins[0]*histogram_bin[1]]+= weight * a[0] * (1.-a[1]);
			bins[histogram_bin[0]+1 + num_bins[0]*(histogram_bin[1]+1)]+= weight * a[0] * a[1];
		}

		/* Bin-wise arithmetic on histograms */
		void operator+=(LinearHistogram2D& other);
		void operator-=(LinearHistogram2D& other);

		// Write and read a ASCII metadata file containing BOV data.
		void writeBovAscii(const char* filename, int index, const char* datafilename);
};

// Histogram with one linear and one logarithmic axis
class LinLogHistogram2D : public Histogram2D {
	private:
		/* Low and high bound of the histogram */
		Vec2d low, high;

	public:
		LinLogHistogram2D(size_t n[2], Vec2d _low, Vec2d _high) :
			Histogram2D(n),
			low(_low), high(_high) {};

		virtual void addValue(Vec2d value, double weight=1.) {
			value -= low;
			value /= high-low;

			int histogram_bin[2];
			histogram_bin[0] = value[0] * num_bins[0];
			histogram_bin[1] = value[1] * num_bins[1];

			for(int i=0; i<2;i++) {
				if(histogram_bin[i] < 0) {
					histogram_bin[i] = 0;
				} else if (histogram_bin[i]+1 > (ssize_t)num_bins[i]) {
					histogram_bin[i] = num_bins[i]-1;
				}
			}
			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1]]+= weight;
		}

	  void addValueLinearInterpolate(Vec2d value, double weight=1.) {

			double v[2];
			v[0] = value[0]- low[0];
			v[0] /= high[0]-low[0];
			v[1] /= value[1] - low[1];
			v[1] = log(v[1])/log(high[1]/low[1]);

			v[0]*=num_bins[0];
			v[1]*=num_bins[1];

			// Interpolation parameter
			double a[2];
			a[0] = value[0] - floor(value[0]);
			a[1] = value[1] - floor(value[1]);

			// Ensure we are within bounds.
			if(v[0] < 0) {
				v[0] = 0;
				a[0] = 0;
			}
			if(v[1] < 0) {
				v[1] = 0;
				a[1] = 0;
			}

			if(v[0] > num_bins[0]-1) {
				v[0] = num_bins[0]-2;
				a[0] = 1;
			}
			if(v[1] > num_bins[1]-1) {
				v[1] = num_bins[1]-2;
				a[1] = 1;
			}

			// Assign.
			int histogram_bin[2];
			histogram_bin[0] = floor(v[0]);
			histogram_bin[1] = floor(v[1]);

			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1]]+= weight * (1.-a[0]) * (1.-a[1]);
			bins[histogram_bin[0] + num_bins[0]*(histogram_bin[1]+1)]+= weight * (1.-a[0]) * a[1];
			bins[histogram_bin[0]+1 + num_bins[0]*histogram_bin[1]]+= weight * a[0] * (1.-a[1]);
			bins[histogram_bin[0]+1 + num_bins[0]*(histogram_bin[1]+1)]+= weight * a[0] * a[1];
		}
};

class LogHistogram2D : public Histogram2D {

	private:
		/* Low and high bound of the histogram */
		Vec2d low, high;

		LogHistogram2D(size_t n[2], Vec2d _low, Vec2d _high) :
			Histogram2D(n),
			low(_low), high(_high) {};

		virtual void addValue(Vec2d value, double weight=1.) {
			value /= low;

			double v[2];
			v[0] = log(value[0])/log(high[0]/low[0]);
			v[1] = log(value[1])/log(high[1]/low[1]);

			int histogram_bin[2];
			histogram_bin[0] = v[0] * num_bins[0];
			histogram_bin[1] = v[1] * num_bins[1];

			for(int i=0; i<2;i++) {
				if(histogram_bin[i] < 0) {
					histogram_bin[i] = 0;
				} else if (histogram_bin[i]+1 > (ssize_t)num_bins[i]) {
					histogram_bin[i] = num_bins[i]-1;
				}
			}
			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1]]+=weight;
		}

		/* Bin-wise arithmetic on histograms */
		void operator+=(LogHistogram2D& other);
		void operator-=(LogHistogram2D& other);
};

/* Histograms of 3D data */
class Histogram3D {
	protected:
		size_t num_bins[3];
		float* bins;

	public:
		Histogram3D(size_t n[3]) {
			num_bins[0] = n[0];
			num_bins[1] = n[1];
			num_bins[2] = n[2];

			bins = new float[num_bins[0]*num_bins[1]*num_bins[2]];
			memset(bins, 0, sizeof(float)*num_bins[0]*num_bins[1]*num_bins[2]);
		}
		~Histogram3D() {
			delete[] bins;
		}

		/* Using MPI_Reduce, sum up all CPUs. */
		void mpi_reduce() {
			float* targetbins = new float[num_bins[0]*num_bins[1]*num_bins[2]];
			if(!targetbins) {
				ERROR("allocation failed while mpi-reducing histogram.\n");
				return;
			}

			MPI::COMM_WORLD.Allreduce(bins,targetbins,num_bins[0]*num_bins[1]*num_bins[2],MPI::FLOAT,MPI_SUM);

			/* Throw away the old bins. */
			delete[] bins;
			bins = targetbins;
		}

		void save(const char* filename) const;
		void load(const char* filename);
		virtual void addValue(Vec3d value) = 0;

		// Access bins
		double operator()(int x, int y, int z) {
			return bins[x + num_bins[0]*y + num_bins[0]*num_bins[1]*z];
		}
};

class LinearHistogram3D : public Histogram3D {

	private:
		/* Low and high bound of the histogram */
		Vec3d low, high;

	public:
		LinearHistogram3D(size_t n[3], Vec3d _low, Vec3d _high) :
			Histogram3D(n),
			low(_low), high(_high) {};

		virtual void addValue(Vec3d value) {
			value -= low;
			value /= high-low;

			int histogram_bin[3];
			histogram_bin[0] = value[0] * num_bins[0];
			histogram_bin[1] = value[1] * num_bins[1];
			histogram_bin[2] = value[2] * num_bins[2];

			for(int i=0; i<3;i++) {
				if(histogram_bin[i] < 0) {
					histogram_bin[i] = 0;
				} else if (histogram_bin[i]+1 > (ssize_t)num_bins[i]) {
					histogram_bin[i] = num_bins[i]-1;
				}
			}
			bins[histogram_bin[0] + num_bins[0]*histogram_bin[1] + num_bins[0]*num_bins[1]*histogram_bin[2]]++;
		}

		Vec3d coords_for_cell(Vec3d cell) {
			Vec3d nx(num_bins[0],num_bins[1],num_bins[2]);
			return low + cell/nx * (high-low);
		}

		/* Bin-wise arithmetic on histograms */
		void operator+=(LinearHistogram3D& other);
		void operator-=(LinearHistogram3D& other);

		// Write and read a ASCII metadata file containing BOV data.
		void writeBovAscii(const char* filename, int index, const char* datafilename);
		void readBov(const char* filename);
};
