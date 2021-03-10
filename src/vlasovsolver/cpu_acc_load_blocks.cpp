/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute, 
 * 2017 CSC - IT center for Science
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


#include "vec.h"
#include "cpu_acc_load_blocks.hpp"

void loadColumnBlockData(
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
   vmesh::GlobalID* blocks,
   vmesh::LocalID n_blocks,
   const int dimension,
   Vec* __restrict__ values) {
   // first set the 0 values for the two empty blocks 
   // we store above and below the existing blocks

   for (uint k=0; k<WID; ++k) {
      for (uint j = 0; j < WID; j += VECL/WID){ 
         values[i_pcolumnv(j, k, -1, n_blocks)] = Vec(0);
         values[i_pcolumnv(j, k, n_blocks, n_blocks)] = Vec(0);
      }
   }

   /*[[[cog
import cog

WID = 4
WID2 = WID * WID

for dimension in range(0, 2):
    if dimension == 0:    
      cell_indices_to_id = [ WID2, WID, 1]
    if dimension == 1:    
      cell_indices_to_id = [ 1, WID2, WID]
    if dimension == 2:    
      cell_indices_to_id = [ 1, WID, WID2]

    cellid_transpose=[]
    for k in range(0,WID):
        for j in range(0,WID):
            for i in range(0,WID):
                cellid_transpose.append(i * cell_indices_to_id[0] +  j * cell_indices_to_id[1] + k * cell_indices_to_id[2])
    
    cog.outl("if(dimension == %s ) {" % dimension)
    cog.outl("   for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {")
    cog.outl("      Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));")    
    for vecl in [4, 8, 16]:
        for accuracy in ["f", "d"]:
            cog.outl("#ifdef VEC%d%s_AGNER" % (vecl, accuracy.upper()))
            cell = 0
            for k in range(0, WID):
                for planeVector in range(0, WID2/vecl):
                    cog.out("      values[i_pcolumnv_b(%d, %d, block_k, n_blocks)] = gather%d%s<" % (planeVector, k, vecl, accuracy) )
                    for vi in range(0,vecl):
                        cog.out("%d" %cellid_transpose[cell])
                        if vi < vecl-1:
                            cog.out(" ,")
                        cell = cell + 1

                    cog.outl(">(data);")
            cog.outl("#endif //VEC%d%s_AGNER" % (vecl, accuracy.upper()))                    
    cog.outl("      //zero old output data")
    cog.outl("      for (uint i=0; i<WID3; ++i) {")
    cog.outl("         data[i]=0;")
    cog.outl("      }")
    cog.outl("   }")                
    cog.outl("}")                

    ]]]*/
   if(dimension == 0 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #ifdef VEC4F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4f<0 ,16 ,32 ,48>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4f<4 ,20 ,36 ,52>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4f<8 ,24 ,40 ,56>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4f<12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4f<1 ,17 ,33 ,49>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4f<5 ,21 ,37 ,53>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4f<9 ,25 ,41 ,57>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4f<13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4f<2 ,18 ,34 ,50>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4f<6 ,22 ,38 ,54>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4f<10 ,26 ,42 ,58>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4f<14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4f<3 ,19 ,35 ,51>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4f<7 ,23 ,39 ,55>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4f<11 ,27 ,43 ,59>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4f<15 ,31 ,47 ,63>(data);
   #endif //VEC4F_AGNER
   #ifdef VEC4D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4d<0 ,16 ,32 ,48>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4d<4 ,20 ,36 ,52>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4d<8 ,24 ,40 ,56>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4d<12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4d<1 ,17 ,33 ,49>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4d<5 ,21 ,37 ,53>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4d<9 ,25 ,41 ,57>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4d<13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4d<2 ,18 ,34 ,50>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4d<6 ,22 ,38 ,54>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4d<10 ,26 ,42 ,58>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4d<14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4d<3 ,19 ,35 ,51>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4d<7 ,23 ,39 ,55>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4d<11 ,27 ,43 ,59>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4d<15 ,31 ,47 ,63>(data);
   #endif //VEC4D_AGNER

#if defined(VEC4F_FALLBACK) || defined(VEC4D_FALLBACK)
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec(data[0], data[16], data[32], data[48]);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec(data[4], data[20], data[36], data[52]);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec(data[8], data[24], data[40], data[56]);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec(data[12], data[28], data[44], data[60]);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec(data[1], data[17], data[33], data[49]);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec(data[5], data[21], data[37], data[53]);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec(data[9], data[25], data[41], data[57]);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec(data[13], data[29], data[45], data[61]);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec(data[2], data[18], data[34], data[50]);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec(data[6], data[22], data[38], data[54]);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec(data[10], data[26], data[42], data[58]);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec(data[14], data[30], data[46], data[62]);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec(data[3], data[19], data[35], data[51]);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec(data[7], data[23], data[39], data[55]);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec(data[11], data[27], data[43], data[59]);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec(data[15], data[31], data[47], data[63]);
#endif //VEC4F_FALLBACK || VEC4D_FALLBACK

   #ifdef VEC8F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0 ,16 ,32 ,48 ,4 ,20 ,36 ,52>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<8 ,24 ,40 ,56 ,12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<1 ,17 ,33 ,49 ,5 ,21 ,37 ,53>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<9 ,25 ,41 ,57 ,13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<2 ,18 ,34 ,50 ,6 ,22 ,38 ,54>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<10 ,26 ,42 ,58 ,14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<3 ,19 ,35 ,51 ,7 ,23 ,39 ,55>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<11 ,27 ,43 ,59 ,15 ,31 ,47 ,63>(data);
   #endif //VEC8F_AGNER
   #ifdef VEC8D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0 ,16 ,32 ,48 ,4 ,20 ,36 ,52>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<8 ,24 ,40 ,56 ,12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<1 ,17 ,33 ,49 ,5 ,21 ,37 ,53>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<9 ,25 ,41 ,57 ,13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<2 ,18 ,34 ,50 ,6 ,22 ,38 ,54>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<10 ,26 ,42 ,58 ,14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<3 ,19 ,35 ,51 ,7 ,23 ,39 ,55>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<11 ,27 ,43 ,59 ,15 ,31 ,47 ,63>(data);
   #endif //VEC8D_AGNER

#if defined(VEC8F_FALLBACK) || defined(VEC8D_FALLBACK)
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec(data[0], data[16], data[32], data[48], data[4], data[20], data[36], data[52]);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec(data[8], data[24], data[40], data[56], data[12], data[28], data[44], data[60]);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec(data[1], data[17], data[33], data[49], data[5], data[21], data[37], data[53]);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec(data[9], data[25], data[41], data[57], data[13], data[29], data[45], data[61]);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec(data[2], data[18], data[34], data[50], data[6], data[22], data[38], data[54]);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec(data[10], data[26], data[42], data[58], data[14], data[30], data[46], data[62]);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec(data[3], data[19], data[35], data[51], data[7], data[23], data[39], data[55]);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec(data[11], data[27], data[43], data[59], data[15], data[31], data[47], data[63]);
#endif //VEC8F_FALLBACK || VEC8D_FALLBACK

   #ifdef VEC16F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0 ,16 ,32 ,48 ,4 ,20 ,36 ,52 ,8 ,24 ,40 ,56 ,12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<1 ,17 ,33 ,49 ,5 ,21 ,37 ,53 ,9 ,25 ,41 ,57 ,13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<2 ,18 ,34 ,50 ,6 ,22 ,38 ,54 ,10 ,26 ,42 ,58 ,14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<3 ,19 ,35 ,51 ,7 ,23 ,39 ,55 ,11 ,27 ,43 ,59 ,15 ,31 ,47 ,63>(data);
   #endif //VEC16F_AGNER
   #ifdef VEC16D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16d<0 ,16 ,32 ,48 ,4 ,20 ,36 ,52 ,8 ,24 ,40 ,56 ,12 ,28 ,44 ,60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16d<1 ,17 ,33 ,49 ,5 ,21 ,37 ,53 ,9 ,25 ,41 ,57 ,13 ,29 ,45 ,61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16d<2 ,18 ,34 ,50 ,6 ,22 ,38 ,54 ,10 ,26 ,42 ,58 ,14 ,30 ,46 ,62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16d<3 ,19 ,35 ,51 ,7 ,23 ,39 ,55 ,11 ,27 ,43 ,59 ,15 ,31 ,47 ,63>(data);
   #endif //VEC16D_AGNER
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
   }
   if(dimension == 1 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #ifdef VEC4F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4f<0 ,1 ,2 ,3>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4f<16 ,17 ,18 ,19>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4f<32 ,33 ,34 ,35>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4f<48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4f<4 ,5 ,6 ,7>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4f<20 ,21 ,22 ,23>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4f<36 ,37 ,38 ,39>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4f<52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4f<8 ,9 ,10 ,11>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4f<24 ,25 ,26 ,27>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4f<40 ,41 ,42 ,43>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4f<56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4f<12 ,13 ,14 ,15>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4f<28 ,29 ,30 ,31>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4f<44 ,45 ,46 ,47>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4f<60 ,61 ,62 ,63>(data);
   #endif //VEC4F_AGNER

#if defined(VEC4F_FALLBACK) || defined(VEC4D_FALLBACK)
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec(data[0], data[1], data[2], data[3]);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec(data[16], data[17], data[18], data[19]);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec(data[32], data[33], data[34], data[35]);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec(data[48], data[49], data[50], data[51]);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec(data[4], data[5], data[6], data[7]);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec(data[20], data[21], data[22], data[23]);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec(data[36], data[37], data[38], data[39]);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec(data[52], data[53], data[54], data[55]);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec(data[8], data[9], data[10], data[11]);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec(data[24], data[25], data[26], data[27]);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec(data[40], data[41], data[42], data[43]);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec(data[56], data[57], data[58], data[59]);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec(data[12], data[13], data[14], data[15]);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec(data[28], data[29], data[30], data[31]);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec(data[44], data[45], data[46], data[47]);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec(data[60], data[61], data[62], data[63]);
#endif //VEC4F_FALLBACK || VEC4D_FALLBACK

   #ifdef VEC4D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4d<0 ,1 ,2 ,3>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4d<16 ,17 ,18 ,19>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4d<32 ,33 ,34 ,35>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4d<48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4d<4 ,5 ,6 ,7>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4d<20 ,21 ,22 ,23>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4d<36 ,37 ,38 ,39>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4d<52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4d<8 ,9 ,10 ,11>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4d<24 ,25 ,26 ,27>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4d<40 ,41 ,42 ,43>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4d<56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4d<12 ,13 ,14 ,15>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4d<28 ,29 ,30 ,31>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4d<44 ,45 ,46 ,47>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4d<60 ,61 ,62 ,63>(data);
   #endif //VEC4D_AGNER
   #ifdef VEC8F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0 ,1 ,2 ,3 ,16 ,17 ,18 ,19>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<32 ,33 ,34 ,35 ,48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<4 ,5 ,6 ,7 ,20 ,21 ,22 ,23>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<36 ,37 ,38 ,39 ,52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<8 ,9 ,10 ,11 ,24 ,25 ,26 ,27>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<40 ,41 ,42 ,43 ,56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<12 ,13 ,14 ,15 ,28 ,29 ,30 ,31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<44 ,45 ,46 ,47 ,60 ,61 ,62 ,63>(data);
   #endif //VEC8F_AGNER
   #ifdef VEC8D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0 ,1 ,2 ,3 ,16 ,17 ,18 ,19>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<32 ,33 ,34 ,35 ,48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<4 ,5 ,6 ,7 ,20 ,21 ,22 ,23>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<36 ,37 ,38 ,39 ,52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<8 ,9 ,10 ,11 ,24 ,25 ,26 ,27>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<40 ,41 ,42 ,43 ,56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<12 ,13 ,14 ,15 ,28 ,29 ,30 ,31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<44 ,45 ,46 ,47 ,60 ,61 ,62 ,63>(data);
   #endif //VEC8D_AGNER

#if defined (VEC8F_FALLBACK) || defined (VEC8D_FALLBACK)
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec(data[0], data[1], data[2], data[3], data[16], data[17], data[18], data[19]);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec(data[32], data[33], data[34], data[35], data[48], data[49], data[50], data[51]);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec(data[4], data[5], data[6], data[7], data[20], data[21], data[22], data[23]);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec(data[36], data[37], data[38], data[39], data[52], data[53], data[54], data[55]);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec(data[8], data[9], data[10], data[11], data[24], data[25], data[26], data[27]);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec(data[40], data[41], data[42], data[43], data[56], data[57], data[58], data[59]);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec(data[12], data[13], data[14], data[15], data[28], data[29], data[30], data[31]);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec(data[44], data[45], data[46], data[47], data[60], data[61], data[62], data[63]);
#endif //VEC8F_FALLBACK || VEC8D_FALLBACK

   #ifdef VEC16F_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0 ,1 ,2 ,3 ,16 ,17 ,18 ,19 ,32 ,33 ,34 ,35 ,48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<4 ,5 ,6 ,7 ,20 ,21 ,22 ,23 ,36 ,37 ,38 ,39 ,52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<8 ,9 ,10 ,11 ,24 ,25 ,26 ,27 ,40 ,41 ,42 ,43 ,56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<12 ,13 ,14 ,15 ,28 ,29 ,30 ,31 ,44 ,45 ,46 ,47 ,60 ,61 ,62 ,63>(data);
   #endif //VEC16F_AGNER
   #ifdef VEC16D_AGNER
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16d<0 ,1 ,2 ,3 ,16 ,17 ,18 ,19 ,32 ,33 ,34 ,35 ,48 ,49 ,50 ,51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16d<4 ,5 ,6 ,7 ,20 ,21 ,22 ,23 ,36 ,37 ,38 ,39 ,52 ,53 ,54 ,55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16d<8 ,9 ,10 ,11 ,24 ,25 ,26 ,27 ,40 ,41 ,42 ,43 ,56 ,57 ,58 ,59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16d<12 ,13 ,14 ,15 ,28 ,29 ,30 ,31 ,44 ,45 ,46 ,47 ,60 ,61 ,62 ,63>(data);
   #endif //VEC16D_AGNER
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
   }
//[[[end]]]

   if (dimension == 2) {
      // copy block data for all blocks. Dimension 2 is easy, here
      // data is in the right order
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
         uint offset = 0;
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               values[i_pcolumnv_b(planeVector, k, block_k, n_blocks)].load(data + offset);
               offset += VECL;
            }
         }
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
   }
}
