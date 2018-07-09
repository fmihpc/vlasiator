/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "histogram.h"

void Histogram1D::save(const char* filename) const {

   double* tempbuf;
   tempbuf = new double[num_bins];
   // MPI Reduce the histograms
   MPI_Allreduce(bins, tempbuf, num_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   int fd = open(filename, O_CREAT|O_TRUNC|O_WRONLY,0644);
   if(!fd || fd == -1) {
      ERROR("unable to write histogram file %s: %s\n", filename, strerror(errno));
      return;
   }

   for(ssize_t remain=sizeof(double) * num_bins; remain > 0;) {
      remain -= write(fd, ((char*)tempbuf) + remain - sizeof(double) * num_bins, remain);
   }

   close(fd);
}

void Histogram1D::saveAscii(const char* filename) const {

   FILE* f = fopen(filename, "w");
   if(!f) {
      ERROR("unable to write histogram file %s: %s\n", filename, strerror(errno));
      return;
   }

   /* Generic histograms can only be written without bin bound specifiers.
    * (Because how would we know them?) */
   for(unsigned int i = 0; i < num_bins; i++) {
      fprintf(f, "%lf\n", bins[i]);
   }

   fclose(f);
}

void LinearHistogram1D::saveAscii(const char* filename) const {

   FILE* f = fopen(filename, "w");
   if(!f) {
      ERROR("unable to write histogram file %s: %s\n", filename, strerror(errno));
      return;
   }

   /* Linear histograms are written with the center bin coordinate as
    * first field */
   for(unsigned int i = 0; i < num_bins; i++) {
      fprintf(f, "%lf %lf\n", low + (i + .5) * (high - low) / num_bins, bins[i]);
   }

   fclose(f);
}

void Histogram1D::load(const char* filename) {
   int fd = open(filename, O_RDONLY);
   int ret;
   if(!fd || fd==-1) {
      ERROR("unable to open histogram file %s for reading: %s\n", filename, strerror(errno));
      return;
   }

   /* TODO: Check filesize */
   size_t size = sizeof(double) * num_bins;

   for(ssize_t remain=size; remain > 0;) {
      ret = read(fd, ((char*)bins) + remain - size, remain);
      if(ret == 0 || !ret) {
         /* TODO: Proper errorhandling here. */
         break;
      }
   }

   close(fd);
}

void Histogram2D::save(const char* filename) const {

   double* tempbuf;
   int num_bins_tot = num_bins[0] * num_bins[1];
   tempbuf = new double[num_bins_tot];
   // MPI Reduce the histograms
   MPI_Allreduce(bins, tempbuf, num_bins_tot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   int fd = open(filename, O_CREAT|O_TRUNC|O_WRONLY,0644);
   if(!fd || fd == -1) {
      ERROR("unable to write histogram file %s: %s\n", filename, strerror(errno));
      return;
   }

   size_t size = sizeof(double) * num_bins_tot;

   for(ssize_t remain=size; remain > 0;) {
      remain -= write(fd, ((char*)tempbuf) + remain - size, remain);
   }

   close(fd);
}

void Histogram2D::load(const char* filename) {
   int fd = open(filename, O_RDONLY);
   int ret;
   if(!fd || fd == -1) {
      ERROR("unable to open histogram file %s for reading: %s\n", filename, strerror(errno));
      return;
   }

   /* TODO: Check filesize */
   size_t size = sizeof(double) * num_bins[0] * num_bins[1];

   for(ssize_t remain=size; remain > 0;) {
      ret = read(fd, ((char*)bins) + remain - size, remain);
      if(ret == 0 || !ret) {
         /* TODO: Proper errorhandling here. */
         ERROR("Read error: %s.\n", strerror(errno));
         break;
      }
      remain -= ret;
   }

   close(fd);
}

void Histogram3D::save(const char* filename) const {

   double* tempbuf;
   int num_bins_tot = num_bins[0] * num_bins[1] * num_bins[2];
   tempbuf = new double[num_bins_tot];
   // MPI Reduce the histograms
   MPI_Allreduce(bins, tempbuf, num_bins_tot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   int fd = open(filename, O_CREAT|O_TRUNC|O_WRONLY, 0644);
   if(!fd || fd == -1) {
      ERROR("unable to write histogram file %s: %s\n", filename, strerror(errno));
      return;
   }

   size_t size = sizeof(float) * num_bins_tot;

   for(ssize_t remain=size; remain > 0;) {
      remain -= write(fd, ((char*)tempbuf) + remain - size, remain);
   }

   close(fd);
}

void Histogram3D::load(const char* filename) {
   int fd = open(filename, O_RDONLY);
   int ret;
   if(!fd || fd == -1) {
      ERROR("unable to open histogram file %s for reading: %s\n", filename, strerror(errno));
      return;
   }

   /* TODO: Check filesize */
   size_t size = sizeof(float) * num_bins[0] * num_bins[1] * num_bins[2];

   for(ssize_t remain=size; remain > 0;) {
      ret = read(fd, ((char*)bins) + remain - size, remain);
      if(ret == 0 || !ret) {
         /* TODO: Proper errorhandling here. */
         ERROR("Read error: %s.\n", strerror(errno));
         break;
      }
      remain -= ret;
   }

   close(fd);
}

void LinearHistogram2D::writeBovAscii(const char* filename, int index, const char* datafilename) {

   FILE* f = fopen(filename, "w");
   if(!f) {
      ERROR("unable to write BOV ascii file %s: %s\n", filename, strerror(errno));
      return;
   }

   fprintf(f, "TIME: %i\n", index);
   fprintf(f, "DATA_FILE: %s\n", datafilename);
   fprintf(f, "DATA_SIZE: %lu %lu 1\n", num_bins[0], num_bins[1]);
   fprintf(f, "DATA_FORMAT: DOUBLE\nVARIABLE: f\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n");
   fprintf(f, "BRICK_ORIGIN: %lf %lf 0\n", low[0], low[1]);
   fprintf(f, "BRICK_SIZE: %lf %lf 1\n", high[0] - low[0], high[1] - low[1]);
   fprintf(f, "DATA_COMPONENTS: 1\n");

   fclose(f);
}

void LinearHistogram3D::writeBovAscii(const char* filename, int index, const char* datafilename) {

   FILE* f = fopen(filename, "w");
   if(!f) {
      ERROR("unable to write BOV ascii file %s: %s\n", filename, strerror(errno));
      return;
   }

   fprintf(f, "TIME: %i\n", index);
   fprintf(f, "DATA_FILE: %s\n", datafilename);
   fprintf(f, "DATA_SIZE: %lu %lu %lu\n", num_bins[0], num_bins[1], num_bins[2]);
   fprintf(f, "DATA_FORMAT: FLOAT\nVARIABLE: f\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n");
   fprintf(f, "BRICK_ORIGIN: %lf %lf %lf\n", low[0], low[1], low[2]);
   fprintf(f, "BRICK_SIZE: %lf %lf %lf\n", high[0] - low[0], high[1] - low[1], high[2] - low[2]);
   fprintf(f, "DATA_COMPONENTS: 1\n");

   fclose(f);
}

void LinearHistogram3D::readBov(const char* filename) {
   char buffer[256];
   char datafilename[256];
   Vec3d size;
   bool filenameread=false;
   FILE* f = fopen(filename, "r");
   if(!f) {
      ERROR("unable to read BOV ascii file %s: %s\n", filename, strerror(errno));
      return;
   }

   while(!feof(f)) {
      if(!fgets(buffer, 256, f)) {
         break;
      }

      // Parse line data
      if(!strncmp("DATA_FILE:", buffer, 10)) {
         sscanf(buffer, "DATA_FILE: %s\n", datafilename);
         filenameread = true;
         continue;
      }
      if(!strncmp("BRICK_ORIGIN:", buffer, 13)) {
         double l[3];
         sscanf(buffer, "BRICK_ORIGIN: %lf %lf %lf", &l[0], &l[1], &l[2]);
         low.load(l);
         continue;
      }
      if(!strncmp("BRICK_SIZE:", buffer, 11)) {
         double s[3];
         sscanf(buffer, "BRICK_SIZE: %lf %lf %lf\n", &s[0], &s[1], &s[2]);
         size.load(s);
         high=low + size;
         continue;
      }
   }

   if(!filenameread) {
      ERROR("BOV file %s did not contain a DATA_FILE statement!\n", filename);
      fclose(f);
      return;
   }

   // Load the actual histogram data from the data file
   load(datafilename);
   fclose(f);
}

/* --- Arithmetic operators for adding and/or substracting Histograms --- */

void LinearHistogram2D::operator+=(LinearHistogram2D& other) {
   /* Check precondition: Both Histograms have to be same size and energy range */
   if(num_bins[0] != other.num_bins[0] ||
         num_bins[1] != other.num_bins[1] ||
         low[0] != other.low[0] || low[1] != other.low[1] ||
         high[0] != other.high[0] || high[1] != other.high[1]) {
      ERROR("Attempted arithmetic on LinearHistogram2D with different bounds.\n");
      return;
   }

   for(size_t j=0; j<num_bins[1]; j++) {
      for(size_t i=0; i<num_bins[0]; i++) {
         bins[j * num_bins[0] + i] += other.bins[j * num_bins[0] +i];
      }
   }
}

void LinearHistogram2D::operator-=(LinearHistogram2D& other) {
   /* Check precondition: Both Histograms have to be same size and energy range */
   if(num_bins[0] != other.num_bins[0] ||
         num_bins[1] != other.num_bins[1] ||
         low[0] != other.low[0] || low[1] != other.low[1] ||
         high[0] != other.high[0] || high[1] != other.high[1]) {
      ERROR("Attempted arithmetic on LinearHistogram2D with different bounds.\n");
      return;
   }

   for(size_t j=0; j<num_bins[1]; j++) {
      for(size_t i=0; i<num_bins[0]; i++) {
         bins[j * num_bins[0] + i] -= other.bins[j * num_bins[0] +i];
      }
   }
}

void LogHistogram2D::operator+=(LogHistogram2D& other) {
   /* Check precondition: Both Histograms have to be same size and energy range */
   if(num_bins[0] != other.num_bins[0] ||
         num_bins[1] != other.num_bins[1] ||
         low[0] != other.low[0] || low[1] != other.low[1] ||
         high[0] != other.high[0] || high[1] != other.high[1]) {
      ERROR("Attempted arithmetic on LogHistogram2D with different bounds.\n");
      return;
   }

   for(size_t j=0; j<num_bins[1]; j++) {
      for(size_t i=0; i<num_bins[0]; i++) {
         bins[j * num_bins[0] + i] += other.bins[j * num_bins[0] +i];
      }
   }
}

void LogHistogram2D::operator-=(LogHistogram2D& other) {
   /* Check precondition: Both Histograms have to be same size and energy range */
   if(num_bins[0] != other.num_bins[0] ||
         num_bins[1] != other.num_bins[1] ||
         low[0] != other.low[0] || low[1] != other.low[1] ||
         high[0] != other.high[0] || high[1] != other.high[1]) {
      ERROR("Attempted arithmetic on LogHistogram2D with different bounds.\n");
      return;
   }

   for(size_t j=0; j<num_bins[1]; j++) {
      for(size_t i=0; i<num_bins[0]; i++) {
         bins[j * num_bins[0] + i] -= other.bins[j * num_bins[0] +i];
      }
   }
}
