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
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

#include "../fieldsolver/fs_limiters.h"
#include "vlsvreaderinterface.h"
#include "vlsv_util.h"

using namespace std;

struct Values {
   double t;
   double E[3];
   double phi;
};

void processFile(const size_t& index,const string& fname,vector<Values>& values0,vector<Values>& values1) {
   vlsvinterface::Reader vlsvReader;
   if (vlsvReader.open(fname) == false) {
      stringstream ss;
      ss << "ERROR, could not open file '" << fname << "' in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << ss.str();
      return;
   }
   
   if (vlsvReader.setCellIds() == false) {
      stringstream ss;
      ss << "ERROR, could not read Cell IDs from file '" << fname << "' in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << ss.str();
      return;
   }

   vector<uint64_t> cellIDs;
   cellIDs.push_back(63400);
   cellIDs.push_back(63401);
   cellIDs.push_back(63799);
   cellIDs.push_back(63800);
   cellIDs.push_back(63801);
   cellIDs.push_back(63802);
   cellIDs.push_back(64199);
   cellIDs.push_back(64200);
   cellIDs.push_back(64201);
   cellIDs.push_back(64202);
   cellIDs.push_back(64600);
   cellIDs.push_back(64601);

   map<uint64_t,array<double,3> > E_values;
   map<uint64_t,double> phi_values;
   for (size_t c=0; c<cellIDs.size(); ++c) {
   //for (map<uint64_t,array<double,3> >::iterator it=E_values.begin(); it!=E_values.end(); ++it) {
      array<double,3> arr;
      if (vlsvReader.getVariable("E_vol",cellIDs[c],arr) == false) {
      //if (vlsvReader.getVariable("E_vol",it->first,it->second) == false) {
         stringstream ss;
         ss << "ERROR, could not read E value from file '" << fname << "' in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
         return;
      } /*else {
         cout << "Read " << it->first << '\t' << (it->second)[0] << endl;
      }*/
      E_values[cellIDs[c]] = arr;
      
      array<double,1> phi_arr;
      if (vlsvReader.getVariable("poisson/potential",cellIDs[c],phi_arr) == false) {
         stringstream ss;
         ss << "ERROR, could not read phi value from file '" << fname << "' in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
      } else {
         phi_values[cellIDs[c]] = phi_arr[0];
      }
   }

   Values val;
   if (vlsvReader.readParameter("time",val.t) == false) {
      stringstream ss;
      ss << "ERROR, could not read 'time' value from file '" << fname << "' in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << ss.str();
      return;
   }
   
   // 0th order electric field
   val.E[0] = 0.25*(E_values[63800][0]+E_values[63801][0]+E_values[64200][0]+E_values[64201][0]);
   val.E[1] = 0.25*(E_values[63800][1]+E_values[63801][1]+E_values[64200][1]+E_values[64201][1]);
   val.E[2] = 0.25*(E_values[63800][2]+E_values[63801][2]+E_values[64200][2]+E_values[64201][2]);
   val.phi  = 0.25*(phi_values[63800]+phi_values[63801]+phi_values[64200]+phi_values[64201]);
   values0[index] = val;

   // 1st order electric field
   const double dx = 5.0;
   double d_E1x = (E_values[63800][0]-E_values[63799][0])/5.0;
   double d_E2x = (E_values[64200][0]-E_values[64199][0])/5.0;
   double d_E3x = (E_values[63802][0]-E_values[63801][0])/5.0;
   double d_E4x = (E_values[64202][0]-E_values[64201][0])/5.0;
   
   val.E[0]  = (E_values[63800][0] + 0.5*d_E1x*dx);
   val.E[0] += (E_values[64200][0] + 0.5*d_E2x*dx);
   val.E[0] += (E_values[63801][0] - 0.5*d_E3x*dx);
   val.E[0] += (E_values[64201][0] - 0.5*d_E4x*dx);
   val.E[0] *= 0.25;
   val.E[1] = 0;
   val.E[2] = 0;
   values1[index] = val;
}

int main(int argn,char* args[]) {
   
   if (argn != 2) {
      cerr << "USAGE: ./esail_intpol <file mask>" << endl;
      return 1;
   }
   
   //Get the file name
   const string mask = args[1];
   vector<string> fileList = toolutil::getFiles(mask);
   sort(fileList.begin(),fileList.end());
   //cout << "Found " << fileList.size() << " files" << endl;

   vector<Values> values0(fileList.size());
   vector<Values> values1(fileList.size());
   for (size_t f=0; f<fileList.size(); ++f) {
      processFile(f,fileList[f],values0,values1);
   }

   double charge = 100e9 * 1.602e-19;
   for (size_t f=0; f<values0.size(); ++f) {
      cout << values0[f].t << '\t' << charge*values0[f].E[0]/1e-9 << '\t' << charge*values1[f].E[0]/1e-9 << '\t' << values0[f].phi << endl;
   }

   return 0;
}
