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
WIDs = [4,8]
for wi, WID in enumerate(WIDs):
      #WID = 8
   WID2 = WID * WID
   if wi == 0:
      cog.outl("#if WID == %d " % WID)
   else:
      cog.outl("#elif WID == %d " % WID)      

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

      cog.outl(" if(dimension == %s ) {" % dimension)
      cog.outl("   for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {")
      cog.outl("      Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));")    
      init = True
      for vecl in [4, 8, 16, 32, 64]:
            # guard against WID2//vecl == 0 and more exotic truncations, and vecl < WID not supported
            if vecl*(WID2//vecl) != WID2 or vecl < WID:
                  continue
            # Agner supports up to 16elem f, 8elem d, so not going further than that!
            for accuracy in ["f", "d"]:
                  if accuracy == "f" and vecl > 16:
                        continue
                  if accuracy == "d" and vecl > 8:
                        continue
                  if init:
                     cog.outl("#if defined(VEC%d%s_AGNER) && VECL == %d" % (vecl, accuracy.upper(), vecl))
                     init = False
                  else:
                     cog.outl("#elif defined(VEC%d%s_AGNER) && VECL == %d" % (vecl, accuracy.upper(), vecl))
                  cell = 0
                  cog.outl("// WID %d, vecl %d" % (WID, vecl))
                  for k in range(0, WID):
                    for planeVector in range(0, WID2//vecl):
                        cog.out("      values[i_pcolumnv_b(%d, %d, block_k, n_blocks)] = gather%d%s<" % (planeVector, k, vecl, accuracy) )
                        for vi in range(0,vecl):
                              cog.out("%d" %cellid_transpose[cell])
                              if vi < vecl-1:
                                cog.out(", ")
                              cell = cell + 1

                        cog.outl(">(data);")
            cog.outl("#elif defined(VEC_FALLBACK_GENERIC) && VECL == %d " % vecl)
            cog.outl("// WID %d, vecl %d" % (WID, vecl))
            cell = 0
            for k in range(0, WID):
                  for planeVector in range(0, WID2//vecl):
                        cog.out("      values[i_pcolumnv_b(%d, %d, block_k, n_blocks)] = Vec({" % (planeVector, k) )
                        for vi in range(0,vecl):
                              cog.out("data[%d]" %cellid_transpose[cell])
                              if vi < vecl-1:
                                    cog.out(", ")
                              cell = cell + 1
                        cog.outl("});")
      cog.outl("#else // Fell through, never fall into this particular pit again")
      cog.outl("\tstd::cerr << \"Undefined VECTORCLASS flag or implementation missing in loadColumnBlockData() before \" << __FILE__ << \":\" << __LINE__ << std::endl;")
      cog.outl("\tabort();")
      cog.outl("#endif")  



      cog.outl("      //zero old output data")
      cog.outl("      for (uint i=0; i<WID3; ++i) {")
      cog.outl("         data[i]=0;")
      cog.outl("      }")
      cog.outl("   }")                
      cog.outl(" }")

cog.outl("#else")
cog.outl("  std::cerr << \"Undefined WID encountered in \" << __FILE__ << \":\" << __LINE__ << std::endl;")
cog.outl("  abort();")      
cog.outl("#endif")          

    ]]]*/
   #if WID == 4 
    if(dimension == 0 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #if defined(VEC4F_AGNER) && VECL == 4
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4f<0, 16, 32, 48>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4f<4, 20, 36, 52>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4f<8, 24, 40, 56>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4f<12, 28, 44, 60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4f<1, 17, 33, 49>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4f<5, 21, 37, 53>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4f<9, 25, 41, 57>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4f<13, 29, 45, 61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4f<2, 18, 34, 50>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4f<6, 22, 38, 54>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4f<10, 26, 42, 58>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4f<14, 30, 46, 62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4f<3, 19, 35, 51>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4f<7, 23, 39, 55>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4f<11, 27, 43, 59>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4f<15, 31, 47, 63>(data);
   #elif defined(VEC4D_AGNER) && VECL == 4
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4d<0, 16, 32, 48>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4d<4, 20, 36, 52>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4d<8, 24, 40, 56>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4d<12, 28, 44, 60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4d<1, 17, 33, 49>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4d<5, 21, 37, 53>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4d<9, 25, 41, 57>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4d<13, 29, 45, 61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4d<2, 18, 34, 50>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4d<6, 22, 38, 54>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4d<10, 26, 42, 58>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4d<14, 30, 46, 62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4d<3, 19, 35, 51>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4d<7, 23, 39, 55>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4d<11, 27, 43, 59>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4d<15, 31, 47, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 4 
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[16], data[32], data[48]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[4], data[20], data[36], data[52]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[8], data[24], data[40], data[56]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[12], data[28], data[44], data[60]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[17], data[33], data[49]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[5], data[21], data[37], data[53]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[9], data[25], data[41], data[57]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[13], data[29], data[45], data[61]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[18], data[34], data[50]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[6], data[22], data[38], data[54]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[10], data[26], data[42], data[58]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[14], data[30], data[46], data[62]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[19], data[35], data[51]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[7], data[23], data[39], data[55]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[11], data[27], data[43], data[59]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[15], data[31], data[47], data[63]});
   #elif defined(VEC8F_AGNER) && VECL == 8
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0, 16, 32, 48, 4, 20, 36, 52>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<8, 24, 40, 56, 12, 28, 44, 60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<1, 17, 33, 49, 5, 21, 37, 53>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<9, 25, 41, 57, 13, 29, 45, 61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<2, 18, 34, 50, 6, 22, 38, 54>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<10, 26, 42, 58, 14, 30, 46, 62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<3, 19, 35, 51, 7, 23, 39, 55>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<11, 27, 43, 59, 15, 31, 47, 63>(data);
   #elif defined(VEC8D_AGNER) && VECL == 8
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0, 16, 32, 48, 4, 20, 36, 52>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<8, 24, 40, 56, 12, 28, 44, 60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<1, 17, 33, 49, 5, 21, 37, 53>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<9, 25, 41, 57, 13, 29, 45, 61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<2, 18, 34, 50, 6, 22, 38, 54>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<10, 26, 42, 58, 14, 30, 46, 62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<3, 19, 35, 51, 7, 23, 39, 55>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<11, 27, 43, 59, 15, 31, 47, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 8 
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[16], data[32], data[48], data[4], data[20], data[36], data[52]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[8], data[24], data[40], data[56], data[12], data[28], data[44], data[60]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[17], data[33], data[49], data[5], data[21], data[37], data[53]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[9], data[25], data[41], data[57], data[13], data[29], data[45], data[61]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[18], data[34], data[50], data[6], data[22], data[38], data[54]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[10], data[26], data[42], data[58], data[14], data[30], data[46], data[62]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[19], data[35], data[51], data[7], data[23], data[39], data[55]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[11], data[27], data[43], data[59], data[15], data[31], data[47], data[63]});
   #elif defined(VEC16F_AGNER) && VECL == 16
   // WID 4, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0, 16, 32, 48, 4, 20, 36, 52, 8, 24, 40, 56, 12, 28, 44, 60>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<1, 17, 33, 49, 5, 21, 37, 53, 9, 25, 41, 57, 13, 29, 45, 61>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<2, 18, 34, 50, 6, 22, 38, 54, 10, 26, 42, 58, 14, 30, 46, 62>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<3, 19, 35, 51, 7, 23, 39, 55, 11, 27, 43, 59, 15, 31, 47, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 16 
   // WID 4, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[16], data[32], data[48], data[4], data[20], data[36], data[52], data[8], data[24], data[40], data[56], data[12], data[28], data[44], data[60]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[17], data[33], data[49], data[5], data[21], data[37], data[53], data[9], data[25], data[41], data[57], data[13], data[29], data[45], data[61]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[18], data[34], data[50], data[6], data[22], data[38], data[54], data[10], data[26], data[42], data[58], data[14], data[30], data[46], data[62]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[19], data[35], data[51], data[7], data[23], data[39], data[55], data[11], data[27], data[43], data[59], data[15], data[31], data[47], data[63]});
   #else // Fell through, never fall into this particular pit again
   	std::cerr << "Undefined VECTORCLASS flag or implementation missing in loadColumnBlockData() before " << __FILE__ << ":" << __LINE__ << std::endl;
   	abort();
   #endif
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
    }
    if(dimension == 1 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #if defined(VEC4F_AGNER) && VECL == 4
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4f<0, 1, 2, 3>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4f<16, 17, 18, 19>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4f<32, 33, 34, 35>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4f<48, 49, 50, 51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4f<4, 5, 6, 7>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4f<20, 21, 22, 23>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4f<36, 37, 38, 39>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4f<52, 53, 54, 55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4f<8, 9, 10, 11>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4f<24, 25, 26, 27>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4f<40, 41, 42, 43>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4f<56, 57, 58, 59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4f<12, 13, 14, 15>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4f<28, 29, 30, 31>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4f<44, 45, 46, 47>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4f<60, 61, 62, 63>(data);
   #elif defined(VEC4D_AGNER) && VECL == 4
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather4d<0, 1, 2, 3>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather4d<16, 17, 18, 19>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather4d<32, 33, 34, 35>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather4d<48, 49, 50, 51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather4d<4, 5, 6, 7>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather4d<20, 21, 22, 23>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather4d<36, 37, 38, 39>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather4d<52, 53, 54, 55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather4d<8, 9, 10, 11>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather4d<24, 25, 26, 27>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather4d<40, 41, 42, 43>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather4d<56, 57, 58, 59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather4d<12, 13, 14, 15>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather4d<28, 29, 30, 31>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather4d<44, 45, 46, 47>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather4d<60, 61, 62, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 4 
   // WID 4, vecl 4
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[16], data[17], data[18], data[19]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[48], data[49], data[50], data[51]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[4], data[5], data[6], data[7]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[20], data[21], data[22], data[23]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[36], data[37], data[38], data[39]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[52], data[53], data[54], data[55]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[24], data[25], data[26], data[27]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[56], data[57], data[58], data[59]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[12], data[13], data[14], data[15]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[28], data[29], data[30], data[31]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[44], data[45], data[46], data[47]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[60], data[61], data[62], data[63]});
   #elif defined(VEC8F_AGNER) && VECL == 8
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0, 1, 2, 3, 16, 17, 18, 19>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<32, 33, 34, 35, 48, 49, 50, 51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<4, 5, 6, 7, 20, 21, 22, 23>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<36, 37, 38, 39, 52, 53, 54, 55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<8, 9, 10, 11, 24, 25, 26, 27>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<40, 41, 42, 43, 56, 57, 58, 59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<12, 13, 14, 15, 28, 29, 30, 31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<44, 45, 46, 47, 60, 61, 62, 63>(data);
   #elif defined(VEC8D_AGNER) && VECL == 8
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0, 1, 2, 3, 16, 17, 18, 19>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<32, 33, 34, 35, 48, 49, 50, 51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<4, 5, 6, 7, 20, 21, 22, 23>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<36, 37, 38, 39, 52, 53, 54, 55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<8, 9, 10, 11, 24, 25, 26, 27>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<40, 41, 42, 43, 56, 57, 58, 59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<12, 13, 14, 15, 28, 29, 30, 31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<44, 45, 46, 47, 60, 61, 62, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 8 
   // WID 4, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[16], data[17], data[18], data[19]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35], data[48], data[49], data[50], data[51]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[4], data[5], data[6], data[7], data[20], data[21], data[22], data[23]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[36], data[37], data[38], data[39], data[52], data[53], data[54], data[55]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[24], data[25], data[26], data[27]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43], data[56], data[57], data[58], data[59]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[12], data[13], data[14], data[15], data[28], data[29], data[30], data[31]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[44], data[45], data[46], data[47], data[60], data[61], data[62], data[63]});
   #elif defined(VEC16F_AGNER) && VECL == 16
   // WID 4, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<4, 5, 6, 7, 20, 21, 22, 23, 36, 37, 38, 39, 52, 53, 54, 55>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<8, 9, 10, 11, 24, 25, 26, 27, 40, 41, 42, 43, 56, 57, 58, 59>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 16 
   // WID 4, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[16], data[17], data[18], data[19], data[32], data[33], data[34], data[35], data[48], data[49], data[50], data[51]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[4], data[5], data[6], data[7], data[20], data[21], data[22], data[23], data[36], data[37], data[38], data[39], data[52], data[53], data[54], data[55]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[24], data[25], data[26], data[27], data[40], data[41], data[42], data[43], data[56], data[57], data[58], data[59]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[12], data[13], data[14], data[15], data[28], data[29], data[30], data[31], data[44], data[45], data[46], data[47], data[60], data[61], data[62], data[63]});
   #else // Fell through, never fall into this particular pit again
   	std::cerr << "Undefined VECTORCLASS flag or implementation missing in loadColumnBlockData() before " << __FILE__ << ":" << __LINE__ << std::endl;
   	abort();
   #endif
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
    }
   #elif WID == 8 
    if(dimension == 0 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #if defined(VEC8F_AGNER) && VECL == 8
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0, 64, 128, 192, 256, 320, 384, 448>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<8, 72, 136, 200, 264, 328, 392, 456>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather8f<16, 80, 144, 208, 272, 336, 400, 464>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather8f<24, 88, 152, 216, 280, 344, 408, 472>(data);
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = gather8f<32, 96, 160, 224, 288, 352, 416, 480>(data);
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = gather8f<40, 104, 168, 232, 296, 360, 424, 488>(data);
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = gather8f<48, 112, 176, 240, 304, 368, 432, 496>(data);
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = gather8f<56, 120, 184, 248, 312, 376, 440, 504>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<1, 65, 129, 193, 257, 321, 385, 449>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<9, 73, 137, 201, 265, 329, 393, 457>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather8f<17, 81, 145, 209, 273, 337, 401, 465>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather8f<25, 89, 153, 217, 281, 345, 409, 473>(data);
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = gather8f<33, 97, 161, 225, 289, 353, 417, 481>(data);
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = gather8f<41, 105, 169, 233, 297, 361, 425, 489>(data);
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = gather8f<49, 113, 177, 241, 305, 369, 433, 497>(data);
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = gather8f<57, 121, 185, 249, 313, 377, 441, 505>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<2, 66, 130, 194, 258, 322, 386, 450>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<10, 74, 138, 202, 266, 330, 394, 458>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather8f<18, 82, 146, 210, 274, 338, 402, 466>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather8f<26, 90, 154, 218, 282, 346, 410, 474>(data);
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = gather8f<34, 98, 162, 226, 290, 354, 418, 482>(data);
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = gather8f<42, 106, 170, 234, 298, 362, 426, 490>(data);
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = gather8f<50, 114, 178, 242, 306, 370, 434, 498>(data);
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = gather8f<58, 122, 186, 250, 314, 378, 442, 506>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<3, 67, 131, 195, 259, 323, 387, 451>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<11, 75, 139, 203, 267, 331, 395, 459>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather8f<19, 83, 147, 211, 275, 339, 403, 467>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather8f<27, 91, 155, 219, 283, 347, 411, 475>(data);
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = gather8f<35, 99, 163, 227, 291, 355, 419, 483>(data);
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = gather8f<43, 107, 171, 235, 299, 363, 427, 491>(data);
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = gather8f<51, 115, 179, 243, 307, 371, 435, 499>(data);
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = gather8f<59, 123, 187, 251, 315, 379, 443, 507>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather8f<4, 68, 132, 196, 260, 324, 388, 452>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather8f<12, 76, 140, 204, 268, 332, 396, 460>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather8f<20, 84, 148, 212, 276, 340, 404, 468>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather8f<28, 92, 156, 220, 284, 348, 412, 476>(data);
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = gather8f<36, 100, 164, 228, 292, 356, 420, 484>(data);
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = gather8f<44, 108, 172, 236, 300, 364, 428, 492>(data);
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = gather8f<52, 116, 180, 244, 308, 372, 436, 500>(data);
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = gather8f<60, 124, 188, 252, 316, 380, 444, 508>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather8f<5, 69, 133, 197, 261, 325, 389, 453>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather8f<13, 77, 141, 205, 269, 333, 397, 461>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather8f<21, 85, 149, 213, 277, 341, 405, 469>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather8f<29, 93, 157, 221, 285, 349, 413, 477>(data);
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = gather8f<37, 101, 165, 229, 293, 357, 421, 485>(data);
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = gather8f<45, 109, 173, 237, 301, 365, 429, 493>(data);
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = gather8f<53, 117, 181, 245, 309, 373, 437, 501>(data);
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = gather8f<61, 125, 189, 253, 317, 381, 445, 509>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather8f<6, 70, 134, 198, 262, 326, 390, 454>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather8f<14, 78, 142, 206, 270, 334, 398, 462>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather8f<22, 86, 150, 214, 278, 342, 406, 470>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather8f<30, 94, 158, 222, 286, 350, 414, 478>(data);
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = gather8f<38, 102, 166, 230, 294, 358, 422, 486>(data);
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = gather8f<46, 110, 174, 238, 302, 366, 430, 494>(data);
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = gather8f<54, 118, 182, 246, 310, 374, 438, 502>(data);
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = gather8f<62, 126, 190, 254, 318, 382, 446, 510>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather8f<7, 71, 135, 199, 263, 327, 391, 455>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather8f<15, 79, 143, 207, 271, 335, 399, 463>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather8f<23, 87, 151, 215, 279, 343, 407, 471>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather8f<31, 95, 159, 223, 287, 351, 415, 479>(data);
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = gather8f<39, 103, 167, 231, 295, 359, 423, 487>(data);
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = gather8f<47, 111, 175, 239, 303, 367, 431, 495>(data);
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = gather8f<55, 119, 183, 247, 311, 375, 439, 503>(data);
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = gather8f<63, 127, 191, 255, 319, 383, 447, 511>(data);
   #elif defined(VEC8D_AGNER) && VECL == 8
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0, 64, 128, 192, 256, 320, 384, 448>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<8, 72, 136, 200, 264, 328, 392, 456>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather8d<16, 80, 144, 208, 272, 336, 400, 464>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather8d<24, 88, 152, 216, 280, 344, 408, 472>(data);
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = gather8d<32, 96, 160, 224, 288, 352, 416, 480>(data);
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = gather8d<40, 104, 168, 232, 296, 360, 424, 488>(data);
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = gather8d<48, 112, 176, 240, 304, 368, 432, 496>(data);
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = gather8d<56, 120, 184, 248, 312, 376, 440, 504>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<1, 65, 129, 193, 257, 321, 385, 449>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<9, 73, 137, 201, 265, 329, 393, 457>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather8d<17, 81, 145, 209, 273, 337, 401, 465>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather8d<25, 89, 153, 217, 281, 345, 409, 473>(data);
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = gather8d<33, 97, 161, 225, 289, 353, 417, 481>(data);
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = gather8d<41, 105, 169, 233, 297, 361, 425, 489>(data);
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = gather8d<49, 113, 177, 241, 305, 369, 433, 497>(data);
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = gather8d<57, 121, 185, 249, 313, 377, 441, 505>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<2, 66, 130, 194, 258, 322, 386, 450>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<10, 74, 138, 202, 266, 330, 394, 458>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather8d<18, 82, 146, 210, 274, 338, 402, 466>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather8d<26, 90, 154, 218, 282, 346, 410, 474>(data);
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = gather8d<34, 98, 162, 226, 290, 354, 418, 482>(data);
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = gather8d<42, 106, 170, 234, 298, 362, 426, 490>(data);
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = gather8d<50, 114, 178, 242, 306, 370, 434, 498>(data);
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = gather8d<58, 122, 186, 250, 314, 378, 442, 506>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<3, 67, 131, 195, 259, 323, 387, 451>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<11, 75, 139, 203, 267, 331, 395, 459>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather8d<19, 83, 147, 211, 275, 339, 403, 467>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather8d<27, 91, 155, 219, 283, 347, 411, 475>(data);
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = gather8d<35, 99, 163, 227, 291, 355, 419, 483>(data);
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = gather8d<43, 107, 171, 235, 299, 363, 427, 491>(data);
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = gather8d<51, 115, 179, 243, 307, 371, 435, 499>(data);
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = gather8d<59, 123, 187, 251, 315, 379, 443, 507>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather8d<4, 68, 132, 196, 260, 324, 388, 452>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather8d<12, 76, 140, 204, 268, 332, 396, 460>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather8d<20, 84, 148, 212, 276, 340, 404, 468>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather8d<28, 92, 156, 220, 284, 348, 412, 476>(data);
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = gather8d<36, 100, 164, 228, 292, 356, 420, 484>(data);
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = gather8d<44, 108, 172, 236, 300, 364, 428, 492>(data);
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = gather8d<52, 116, 180, 244, 308, 372, 436, 500>(data);
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = gather8d<60, 124, 188, 252, 316, 380, 444, 508>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather8d<5, 69, 133, 197, 261, 325, 389, 453>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather8d<13, 77, 141, 205, 269, 333, 397, 461>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather8d<21, 85, 149, 213, 277, 341, 405, 469>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather8d<29, 93, 157, 221, 285, 349, 413, 477>(data);
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = gather8d<37, 101, 165, 229, 293, 357, 421, 485>(data);
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = gather8d<45, 109, 173, 237, 301, 365, 429, 493>(data);
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = gather8d<53, 117, 181, 245, 309, 373, 437, 501>(data);
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = gather8d<61, 125, 189, 253, 317, 381, 445, 509>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather8d<6, 70, 134, 198, 262, 326, 390, 454>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather8d<14, 78, 142, 206, 270, 334, 398, 462>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather8d<22, 86, 150, 214, 278, 342, 406, 470>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather8d<30, 94, 158, 222, 286, 350, 414, 478>(data);
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = gather8d<38, 102, 166, 230, 294, 358, 422, 486>(data);
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = gather8d<46, 110, 174, 238, 302, 366, 430, 494>(data);
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = gather8d<54, 118, 182, 246, 310, 374, 438, 502>(data);
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = gather8d<62, 126, 190, 254, 318, 382, 446, 510>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather8d<7, 71, 135, 199, 263, 327, 391, 455>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather8d<15, 79, 143, 207, 271, 335, 399, 463>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather8d<23, 87, 151, 215, 279, 343, 407, 471>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather8d<31, 95, 159, 223, 287, 351, 415, 479>(data);
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = gather8d<39, 103, 167, 231, 295, 359, 423, 487>(data);
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = gather8d<47, 111, 175, 239, 303, 367, 431, 495>(data);
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = gather8d<55, 119, 183, 247, 311, 375, 439, 503>(data);
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = gather8d<63, 127, 191, 255, 319, 383, 447, 511>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 8 
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[64], data[128], data[192], data[256], data[320], data[384], data[448]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[8], data[72], data[136], data[200], data[264], data[328], data[392], data[456]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[16], data[80], data[144], data[208], data[272], data[336], data[400], data[464]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[24], data[88], data[152], data[216], data[280], data[344], data[408], data[472]});
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = Vec({data[32], data[96], data[160], data[224], data[288], data[352], data[416], data[480]});
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = Vec({data[40], data[104], data[168], data[232], data[296], data[360], data[424], data[488]});
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = Vec({data[48], data[112], data[176], data[240], data[304], data[368], data[432], data[496]});
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = Vec({data[56], data[120], data[184], data[248], data[312], data[376], data[440], data[504]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[65], data[129], data[193], data[257], data[321], data[385], data[449]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[9], data[73], data[137], data[201], data[265], data[329], data[393], data[457]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[17], data[81], data[145], data[209], data[273], data[337], data[401], data[465]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[25], data[89], data[153], data[217], data[281], data[345], data[409], data[473]});
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = Vec({data[33], data[97], data[161], data[225], data[289], data[353], data[417], data[481]});
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = Vec({data[41], data[105], data[169], data[233], data[297], data[361], data[425], data[489]});
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = Vec({data[49], data[113], data[177], data[241], data[305], data[369], data[433], data[497]});
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = Vec({data[57], data[121], data[185], data[249], data[313], data[377], data[441], data[505]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[66], data[130], data[194], data[258], data[322], data[386], data[450]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[10], data[74], data[138], data[202], data[266], data[330], data[394], data[458]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[18], data[82], data[146], data[210], data[274], data[338], data[402], data[466]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[26], data[90], data[154], data[218], data[282], data[346], data[410], data[474]});
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = Vec({data[34], data[98], data[162], data[226], data[290], data[354], data[418], data[482]});
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = Vec({data[42], data[106], data[170], data[234], data[298], data[362], data[426], data[490]});
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = Vec({data[50], data[114], data[178], data[242], data[306], data[370], data[434], data[498]});
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = Vec({data[58], data[122], data[186], data[250], data[314], data[378], data[442], data[506]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[67], data[131], data[195], data[259], data[323], data[387], data[451]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[11], data[75], data[139], data[203], data[267], data[331], data[395], data[459]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[19], data[83], data[147], data[211], data[275], data[339], data[403], data[467]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[27], data[91], data[155], data[219], data[283], data[347], data[411], data[475]});
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = Vec({data[35], data[99], data[163], data[227], data[291], data[355], data[419], data[483]});
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = Vec({data[43], data[107], data[171], data[235], data[299], data[363], data[427], data[491]});
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = Vec({data[51], data[115], data[179], data[243], data[307], data[371], data[435], data[499]});
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = Vec({data[59], data[123], data[187], data[251], data[315], data[379], data[443], data[507]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[4], data[68], data[132], data[196], data[260], data[324], data[388], data[452]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[12], data[76], data[140], data[204], data[268], data[332], data[396], data[460]});
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = Vec({data[20], data[84], data[148], data[212], data[276], data[340], data[404], data[468]});
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = Vec({data[28], data[92], data[156], data[220], data[284], data[348], data[412], data[476]});
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = Vec({data[36], data[100], data[164], data[228], data[292], data[356], data[420], data[484]});
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = Vec({data[44], data[108], data[172], data[236], data[300], data[364], data[428], data[492]});
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = Vec({data[52], data[116], data[180], data[244], data[308], data[372], data[436], data[500]});
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = Vec({data[60], data[124], data[188], data[252], data[316], data[380], data[444], data[508]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[5], data[69], data[133], data[197], data[261], data[325], data[389], data[453]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[13], data[77], data[141], data[205], data[269], data[333], data[397], data[461]});
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = Vec({data[21], data[85], data[149], data[213], data[277], data[341], data[405], data[469]});
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = Vec({data[29], data[93], data[157], data[221], data[285], data[349], data[413], data[477]});
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = Vec({data[37], data[101], data[165], data[229], data[293], data[357], data[421], data[485]});
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = Vec({data[45], data[109], data[173], data[237], data[301], data[365], data[429], data[493]});
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = Vec({data[53], data[117], data[181], data[245], data[309], data[373], data[437], data[501]});
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = Vec({data[61], data[125], data[189], data[253], data[317], data[381], data[445], data[509]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[6], data[70], data[134], data[198], data[262], data[326], data[390], data[454]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[14], data[78], data[142], data[206], data[270], data[334], data[398], data[462]});
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = Vec({data[22], data[86], data[150], data[214], data[278], data[342], data[406], data[470]});
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = Vec({data[30], data[94], data[158], data[222], data[286], data[350], data[414], data[478]});
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = Vec({data[38], data[102], data[166], data[230], data[294], data[358], data[422], data[486]});
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = Vec({data[46], data[110], data[174], data[238], data[302], data[366], data[430], data[494]});
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = Vec({data[54], data[118], data[182], data[246], data[310], data[374], data[438], data[502]});
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = Vec({data[62], data[126], data[190], data[254], data[318], data[382], data[446], data[510]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[7], data[71], data[135], data[199], data[263], data[327], data[391], data[455]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[15], data[79], data[143], data[207], data[271], data[335], data[399], data[463]});
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = Vec({data[23], data[87], data[151], data[215], data[279], data[343], data[407], data[471]});
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = Vec({data[31], data[95], data[159], data[223], data[287], data[351], data[415], data[479]});
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = Vec({data[39], data[103], data[167], data[231], data[295], data[359], data[423], data[487]});
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = Vec({data[47], data[111], data[175], data[239], data[303], data[367], data[431], data[495]});
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = Vec({data[55], data[119], data[183], data[247], data[311], data[375], data[439], data[503]});
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = Vec({data[63], data[127], data[191], data[255], data[319], data[383], data[447], data[511]});
   #elif defined(VEC16F_AGNER) && VECL == 16
   // WID 8, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0, 64, 128, 192, 256, 320, 384, 448, 8, 72, 136, 200, 264, 328, 392, 456>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather16f<16, 80, 144, 208, 272, 336, 400, 464, 24, 88, 152, 216, 280, 344, 408, 472>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather16f<32, 96, 160, 224, 288, 352, 416, 480, 40, 104, 168, 232, 296, 360, 424, 488>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather16f<48, 112, 176, 240, 304, 368, 432, 496, 56, 120, 184, 248, 312, 376, 440, 504>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<1, 65, 129, 193, 257, 321, 385, 449, 9, 73, 137, 201, 265, 329, 393, 457>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather16f<17, 81, 145, 209, 273, 337, 401, 465, 25, 89, 153, 217, 281, 345, 409, 473>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather16f<33, 97, 161, 225, 289, 353, 417, 481, 41, 105, 169, 233, 297, 361, 425, 489>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather16f<49, 113, 177, 241, 305, 369, 433, 497, 57, 121, 185, 249, 313, 377, 441, 505>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<2, 66, 130, 194, 258, 322, 386, 450, 10, 74, 138, 202, 266, 330, 394, 458>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather16f<18, 82, 146, 210, 274, 338, 402, 466, 26, 90, 154, 218, 282, 346, 410, 474>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather16f<34, 98, 162, 226, 290, 354, 418, 482, 42, 106, 170, 234, 298, 362, 426, 490>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather16f<50, 114, 178, 242, 306, 370, 434, 498, 58, 122, 186, 250, 314, 378, 442, 506>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<3, 67, 131, 195, 259, 323, 387, 451, 11, 75, 139, 203, 267, 331, 395, 459>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather16f<19, 83, 147, 211, 275, 339, 403, 467, 27, 91, 155, 219, 283, 347, 411, 475>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather16f<35, 99, 163, 227, 291, 355, 419, 483, 43, 107, 171, 235, 299, 363, 427, 491>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather16f<51, 115, 179, 243, 307, 371, 435, 499, 59, 123, 187, 251, 315, 379, 443, 507>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather16f<4, 68, 132, 196, 260, 324, 388, 452, 12, 76, 140, 204, 268, 332, 396, 460>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather16f<20, 84, 148, 212, 276, 340, 404, 468, 28, 92, 156, 220, 284, 348, 412, 476>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather16f<36, 100, 164, 228, 292, 356, 420, 484, 44, 108, 172, 236, 300, 364, 428, 492>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather16f<52, 116, 180, 244, 308, 372, 436, 500, 60, 124, 188, 252, 316, 380, 444, 508>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather16f<5, 69, 133, 197, 261, 325, 389, 453, 13, 77, 141, 205, 269, 333, 397, 461>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather16f<21, 85, 149, 213, 277, 341, 405, 469, 29, 93, 157, 221, 285, 349, 413, 477>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather16f<37, 101, 165, 229, 293, 357, 421, 485, 45, 109, 173, 237, 301, 365, 429, 493>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather16f<53, 117, 181, 245, 309, 373, 437, 501, 61, 125, 189, 253, 317, 381, 445, 509>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather16f<6, 70, 134, 198, 262, 326, 390, 454, 14, 78, 142, 206, 270, 334, 398, 462>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather16f<22, 86, 150, 214, 278, 342, 406, 470, 30, 94, 158, 222, 286, 350, 414, 478>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather16f<38, 102, 166, 230, 294, 358, 422, 486, 46, 110, 174, 238, 302, 366, 430, 494>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather16f<54, 118, 182, 246, 310, 374, 438, 502, 62, 126, 190, 254, 318, 382, 446, 510>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather16f<7, 71, 135, 199, 263, 327, 391, 455, 15, 79, 143, 207, 271, 335, 399, 463>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather16f<23, 87, 151, 215, 279, 343, 407, 471, 31, 95, 159, 223, 287, 351, 415, 479>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather16f<39, 103, 167, 231, 295, 359, 423, 487, 47, 111, 175, 239, 303, 367, 431, 495>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather16f<55, 119, 183, 247, 311, 375, 439, 503, 63, 127, 191, 255, 319, 383, 447, 511>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 16 
   // WID 8, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[64], data[128], data[192], data[256], data[320], data[384], data[448], data[8], data[72], data[136], data[200], data[264], data[328], data[392], data[456]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[16], data[80], data[144], data[208], data[272], data[336], data[400], data[464], data[24], data[88], data[152], data[216], data[280], data[344], data[408], data[472]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[32], data[96], data[160], data[224], data[288], data[352], data[416], data[480], data[40], data[104], data[168], data[232], data[296], data[360], data[424], data[488]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[48], data[112], data[176], data[240], data[304], data[368], data[432], data[496], data[56], data[120], data[184], data[248], data[312], data[376], data[440], data[504]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[65], data[129], data[193], data[257], data[321], data[385], data[449], data[9], data[73], data[137], data[201], data[265], data[329], data[393], data[457]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[17], data[81], data[145], data[209], data[273], data[337], data[401], data[465], data[25], data[89], data[153], data[217], data[281], data[345], data[409], data[473]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[33], data[97], data[161], data[225], data[289], data[353], data[417], data[481], data[41], data[105], data[169], data[233], data[297], data[361], data[425], data[489]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[49], data[113], data[177], data[241], data[305], data[369], data[433], data[497], data[57], data[121], data[185], data[249], data[313], data[377], data[441], data[505]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[66], data[130], data[194], data[258], data[322], data[386], data[450], data[10], data[74], data[138], data[202], data[266], data[330], data[394], data[458]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[18], data[82], data[146], data[210], data[274], data[338], data[402], data[466], data[26], data[90], data[154], data[218], data[282], data[346], data[410], data[474]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[34], data[98], data[162], data[226], data[290], data[354], data[418], data[482], data[42], data[106], data[170], data[234], data[298], data[362], data[426], data[490]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[50], data[114], data[178], data[242], data[306], data[370], data[434], data[498], data[58], data[122], data[186], data[250], data[314], data[378], data[442], data[506]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[67], data[131], data[195], data[259], data[323], data[387], data[451], data[11], data[75], data[139], data[203], data[267], data[331], data[395], data[459]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[19], data[83], data[147], data[211], data[275], data[339], data[403], data[467], data[27], data[91], data[155], data[219], data[283], data[347], data[411], data[475]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[35], data[99], data[163], data[227], data[291], data[355], data[419], data[483], data[43], data[107], data[171], data[235], data[299], data[363], data[427], data[491]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[51], data[115], data[179], data[243], data[307], data[371], data[435], data[499], data[59], data[123], data[187], data[251], data[315], data[379], data[443], data[507]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[4], data[68], data[132], data[196], data[260], data[324], data[388], data[452], data[12], data[76], data[140], data[204], data[268], data[332], data[396], data[460]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[20], data[84], data[148], data[212], data[276], data[340], data[404], data[468], data[28], data[92], data[156], data[220], data[284], data[348], data[412], data[476]});
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = Vec({data[36], data[100], data[164], data[228], data[292], data[356], data[420], data[484], data[44], data[108], data[172], data[236], data[300], data[364], data[428], data[492]});
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = Vec({data[52], data[116], data[180], data[244], data[308], data[372], data[436], data[500], data[60], data[124], data[188], data[252], data[316], data[380], data[444], data[508]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[5], data[69], data[133], data[197], data[261], data[325], data[389], data[453], data[13], data[77], data[141], data[205], data[269], data[333], data[397], data[461]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[21], data[85], data[149], data[213], data[277], data[341], data[405], data[469], data[29], data[93], data[157], data[221], data[285], data[349], data[413], data[477]});
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = Vec({data[37], data[101], data[165], data[229], data[293], data[357], data[421], data[485], data[45], data[109], data[173], data[237], data[301], data[365], data[429], data[493]});
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = Vec({data[53], data[117], data[181], data[245], data[309], data[373], data[437], data[501], data[61], data[125], data[189], data[253], data[317], data[381], data[445], data[509]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[6], data[70], data[134], data[198], data[262], data[326], data[390], data[454], data[14], data[78], data[142], data[206], data[270], data[334], data[398], data[462]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[22], data[86], data[150], data[214], data[278], data[342], data[406], data[470], data[30], data[94], data[158], data[222], data[286], data[350], data[414], data[478]});
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = Vec({data[38], data[102], data[166], data[230], data[294], data[358], data[422], data[486], data[46], data[110], data[174], data[238], data[302], data[366], data[430], data[494]});
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = Vec({data[54], data[118], data[182], data[246], data[310], data[374], data[438], data[502], data[62], data[126], data[190], data[254], data[318], data[382], data[446], data[510]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[7], data[71], data[135], data[199], data[263], data[327], data[391], data[455], data[15], data[79], data[143], data[207], data[271], data[335], data[399], data[463]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[23], data[87], data[151], data[215], data[279], data[343], data[407], data[471], data[31], data[95], data[159], data[223], data[287], data[351], data[415], data[479]});
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = Vec({data[39], data[103], data[167], data[231], data[295], data[359], data[423], data[487], data[47], data[111], data[175], data[239], data[303], data[367], data[431], data[495]});
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = Vec({data[55], data[119], data[183], data[247], data[311], data[375], data[439], data[503], data[63], data[127], data[191], data[255], data[319], data[383], data[447], data[511]});
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 32 
   // WID 8, vecl 32
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[64], data[128], data[192], data[256], data[320], data[384], data[448], data[8], data[72], data[136], data[200], data[264], data[328], data[392], data[456], data[16], data[80], data[144], data[208], data[272], data[336], data[400], data[464], data[24], data[88], data[152], data[216], data[280], data[344], data[408], data[472]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[32], data[96], data[160], data[224], data[288], data[352], data[416], data[480], data[40], data[104], data[168], data[232], data[296], data[360], data[424], data[488], data[48], data[112], data[176], data[240], data[304], data[368], data[432], data[496], data[56], data[120], data[184], data[248], data[312], data[376], data[440], data[504]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[65], data[129], data[193], data[257], data[321], data[385], data[449], data[9], data[73], data[137], data[201], data[265], data[329], data[393], data[457], data[17], data[81], data[145], data[209], data[273], data[337], data[401], data[465], data[25], data[89], data[153], data[217], data[281], data[345], data[409], data[473]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[33], data[97], data[161], data[225], data[289], data[353], data[417], data[481], data[41], data[105], data[169], data[233], data[297], data[361], data[425], data[489], data[49], data[113], data[177], data[241], data[305], data[369], data[433], data[497], data[57], data[121], data[185], data[249], data[313], data[377], data[441], data[505]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[66], data[130], data[194], data[258], data[322], data[386], data[450], data[10], data[74], data[138], data[202], data[266], data[330], data[394], data[458], data[18], data[82], data[146], data[210], data[274], data[338], data[402], data[466], data[26], data[90], data[154], data[218], data[282], data[346], data[410], data[474]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[34], data[98], data[162], data[226], data[290], data[354], data[418], data[482], data[42], data[106], data[170], data[234], data[298], data[362], data[426], data[490], data[50], data[114], data[178], data[242], data[306], data[370], data[434], data[498], data[58], data[122], data[186], data[250], data[314], data[378], data[442], data[506]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[67], data[131], data[195], data[259], data[323], data[387], data[451], data[11], data[75], data[139], data[203], data[267], data[331], data[395], data[459], data[19], data[83], data[147], data[211], data[275], data[339], data[403], data[467], data[27], data[91], data[155], data[219], data[283], data[347], data[411], data[475]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[35], data[99], data[163], data[227], data[291], data[355], data[419], data[483], data[43], data[107], data[171], data[235], data[299], data[363], data[427], data[491], data[51], data[115], data[179], data[243], data[307], data[371], data[435], data[499], data[59], data[123], data[187], data[251], data[315], data[379], data[443], data[507]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[4], data[68], data[132], data[196], data[260], data[324], data[388], data[452], data[12], data[76], data[140], data[204], data[268], data[332], data[396], data[460], data[20], data[84], data[148], data[212], data[276], data[340], data[404], data[468], data[28], data[92], data[156], data[220], data[284], data[348], data[412], data[476]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[36], data[100], data[164], data[228], data[292], data[356], data[420], data[484], data[44], data[108], data[172], data[236], data[300], data[364], data[428], data[492], data[52], data[116], data[180], data[244], data[308], data[372], data[436], data[500], data[60], data[124], data[188], data[252], data[316], data[380], data[444], data[508]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[5], data[69], data[133], data[197], data[261], data[325], data[389], data[453], data[13], data[77], data[141], data[205], data[269], data[333], data[397], data[461], data[21], data[85], data[149], data[213], data[277], data[341], data[405], data[469], data[29], data[93], data[157], data[221], data[285], data[349], data[413], data[477]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[37], data[101], data[165], data[229], data[293], data[357], data[421], data[485], data[45], data[109], data[173], data[237], data[301], data[365], data[429], data[493], data[53], data[117], data[181], data[245], data[309], data[373], data[437], data[501], data[61], data[125], data[189], data[253], data[317], data[381], data[445], data[509]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[6], data[70], data[134], data[198], data[262], data[326], data[390], data[454], data[14], data[78], data[142], data[206], data[270], data[334], data[398], data[462], data[22], data[86], data[150], data[214], data[278], data[342], data[406], data[470], data[30], data[94], data[158], data[222], data[286], data[350], data[414], data[478]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[38], data[102], data[166], data[230], data[294], data[358], data[422], data[486], data[46], data[110], data[174], data[238], data[302], data[366], data[430], data[494], data[54], data[118], data[182], data[246], data[310], data[374], data[438], data[502], data[62], data[126], data[190], data[254], data[318], data[382], data[446], data[510]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[7], data[71], data[135], data[199], data[263], data[327], data[391], data[455], data[15], data[79], data[143], data[207], data[271], data[335], data[399], data[463], data[23], data[87], data[151], data[215], data[279], data[343], data[407], data[471], data[31], data[95], data[159], data[223], data[287], data[351], data[415], data[479]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[39], data[103], data[167], data[231], data[295], data[359], data[423], data[487], data[47], data[111], data[175], data[239], data[303], data[367], data[431], data[495], data[55], data[119], data[183], data[247], data[311], data[375], data[439], data[503], data[63], data[127], data[191], data[255], data[319], data[383], data[447], data[511]});
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 64 
   // WID 8, vecl 64
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[64], data[128], data[192], data[256], data[320], data[384], data[448], data[8], data[72], data[136], data[200], data[264], data[328], data[392], data[456], data[16], data[80], data[144], data[208], data[272], data[336], data[400], data[464], data[24], data[88], data[152], data[216], data[280], data[344], data[408], data[472], data[32], data[96], data[160], data[224], data[288], data[352], data[416], data[480], data[40], data[104], data[168], data[232], data[296], data[360], data[424], data[488], data[48], data[112], data[176], data[240], data[304], data[368], data[432], data[496], data[56], data[120], data[184], data[248], data[312], data[376], data[440], data[504]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[1], data[65], data[129], data[193], data[257], data[321], data[385], data[449], data[9], data[73], data[137], data[201], data[265], data[329], data[393], data[457], data[17], data[81], data[145], data[209], data[273], data[337], data[401], data[465], data[25], data[89], data[153], data[217], data[281], data[345], data[409], data[473], data[33], data[97], data[161], data[225], data[289], data[353], data[417], data[481], data[41], data[105], data[169], data[233], data[297], data[361], data[425], data[489], data[49], data[113], data[177], data[241], data[305], data[369], data[433], data[497], data[57], data[121], data[185], data[249], data[313], data[377], data[441], data[505]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[2], data[66], data[130], data[194], data[258], data[322], data[386], data[450], data[10], data[74], data[138], data[202], data[266], data[330], data[394], data[458], data[18], data[82], data[146], data[210], data[274], data[338], data[402], data[466], data[26], data[90], data[154], data[218], data[282], data[346], data[410], data[474], data[34], data[98], data[162], data[226], data[290], data[354], data[418], data[482], data[42], data[106], data[170], data[234], data[298], data[362], data[426], data[490], data[50], data[114], data[178], data[242], data[306], data[370], data[434], data[498], data[58], data[122], data[186], data[250], data[314], data[378], data[442], data[506]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[3], data[67], data[131], data[195], data[259], data[323], data[387], data[451], data[11], data[75], data[139], data[203], data[267], data[331], data[395], data[459], data[19], data[83], data[147], data[211], data[275], data[339], data[403], data[467], data[27], data[91], data[155], data[219], data[283], data[347], data[411], data[475], data[35], data[99], data[163], data[227], data[291], data[355], data[419], data[483], data[43], data[107], data[171], data[235], data[299], data[363], data[427], data[491], data[51], data[115], data[179], data[243], data[307], data[371], data[435], data[499], data[59], data[123], data[187], data[251], data[315], data[379], data[443], data[507]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[4], data[68], data[132], data[196], data[260], data[324], data[388], data[452], data[12], data[76], data[140], data[204], data[268], data[332], data[396], data[460], data[20], data[84], data[148], data[212], data[276], data[340], data[404], data[468], data[28], data[92], data[156], data[220], data[284], data[348], data[412], data[476], data[36], data[100], data[164], data[228], data[292], data[356], data[420], data[484], data[44], data[108], data[172], data[236], data[300], data[364], data[428], data[492], data[52], data[116], data[180], data[244], data[308], data[372], data[436], data[500], data[60], data[124], data[188], data[252], data[316], data[380], data[444], data[508]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[5], data[69], data[133], data[197], data[261], data[325], data[389], data[453], data[13], data[77], data[141], data[205], data[269], data[333], data[397], data[461], data[21], data[85], data[149], data[213], data[277], data[341], data[405], data[469], data[29], data[93], data[157], data[221], data[285], data[349], data[413], data[477], data[37], data[101], data[165], data[229], data[293], data[357], data[421], data[485], data[45], data[109], data[173], data[237], data[301], data[365], data[429], data[493], data[53], data[117], data[181], data[245], data[309], data[373], data[437], data[501], data[61], data[125], data[189], data[253], data[317], data[381], data[445], data[509]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[6], data[70], data[134], data[198], data[262], data[326], data[390], data[454], data[14], data[78], data[142], data[206], data[270], data[334], data[398], data[462], data[22], data[86], data[150], data[214], data[278], data[342], data[406], data[470], data[30], data[94], data[158], data[222], data[286], data[350], data[414], data[478], data[38], data[102], data[166], data[230], data[294], data[358], data[422], data[486], data[46], data[110], data[174], data[238], data[302], data[366], data[430], data[494], data[54], data[118], data[182], data[246], data[310], data[374], data[438], data[502], data[62], data[126], data[190], data[254], data[318], data[382], data[446], data[510]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[7], data[71], data[135], data[199], data[263], data[327], data[391], data[455], data[15], data[79], data[143], data[207], data[271], data[335], data[399], data[463], data[23], data[87], data[151], data[215], data[279], data[343], data[407], data[471], data[31], data[95], data[159], data[223], data[287], data[351], data[415], data[479], data[39], data[103], data[167], data[231], data[295], data[359], data[423], data[487], data[47], data[111], data[175], data[239], data[303], data[367], data[431], data[495], data[55], data[119], data[183], data[247], data[311], data[375], data[439], data[503], data[63], data[127], data[191], data[255], data[319], data[383], data[447], data[511]});
   #else // Fell through, never fall into this particular pit again
   	std::cerr << "Undefined VECTORCLASS flag or implementation missing in loadColumnBlockData() before " << __FILE__ << ":" << __LINE__ << std::endl;
   	abort();
   #endif
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
    }
    if(dimension == 1 ) {
      for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
         Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));
   #if defined(VEC8F_AGNER) && VECL == 8
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8f<0, 1, 2, 3, 4, 5, 6, 7>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8f<64, 65, 66, 67, 68, 69, 70, 71>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather8f<128, 129, 130, 131, 132, 133, 134, 135>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather8f<192, 193, 194, 195, 196, 197, 198, 199>(data);
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = gather8f<256, 257, 258, 259, 260, 261, 262, 263>(data);
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = gather8f<320, 321, 322, 323, 324, 325, 326, 327>(data);
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = gather8f<384, 385, 386, 387, 388, 389, 390, 391>(data);
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = gather8f<448, 449, 450, 451, 452, 453, 454, 455>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8f<8, 9, 10, 11, 12, 13, 14, 15>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8f<72, 73, 74, 75, 76, 77, 78, 79>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather8f<136, 137, 138, 139, 140, 141, 142, 143>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather8f<200, 201, 202, 203, 204, 205, 206, 207>(data);
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = gather8f<264, 265, 266, 267, 268, 269, 270, 271>(data);
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = gather8f<328, 329, 330, 331, 332, 333, 334, 335>(data);
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = gather8f<392, 393, 394, 395, 396, 397, 398, 399>(data);
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = gather8f<456, 457, 458, 459, 460, 461, 462, 463>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8f<16, 17, 18, 19, 20, 21, 22, 23>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8f<80, 81, 82, 83, 84, 85, 86, 87>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather8f<144, 145, 146, 147, 148, 149, 150, 151>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather8f<208, 209, 210, 211, 212, 213, 214, 215>(data);
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = gather8f<272, 273, 274, 275, 276, 277, 278, 279>(data);
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = gather8f<336, 337, 338, 339, 340, 341, 342, 343>(data);
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = gather8f<400, 401, 402, 403, 404, 405, 406, 407>(data);
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = gather8f<464, 465, 466, 467, 468, 469, 470, 471>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8f<24, 25, 26, 27, 28, 29, 30, 31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8f<88, 89, 90, 91, 92, 93, 94, 95>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather8f<152, 153, 154, 155, 156, 157, 158, 159>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather8f<216, 217, 218, 219, 220, 221, 222, 223>(data);
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = gather8f<280, 281, 282, 283, 284, 285, 286, 287>(data);
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = gather8f<344, 345, 346, 347, 348, 349, 350, 351>(data);
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = gather8f<408, 409, 410, 411, 412, 413, 414, 415>(data);
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = gather8f<472, 473, 474, 475, 476, 477, 478, 479>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather8f<32, 33, 34, 35, 36, 37, 38, 39>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather8f<96, 97, 98, 99, 100, 101, 102, 103>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather8f<160, 161, 162, 163, 164, 165, 166, 167>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather8f<224, 225, 226, 227, 228, 229, 230, 231>(data);
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = gather8f<288, 289, 290, 291, 292, 293, 294, 295>(data);
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = gather8f<352, 353, 354, 355, 356, 357, 358, 359>(data);
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = gather8f<416, 417, 418, 419, 420, 421, 422, 423>(data);
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = gather8f<480, 481, 482, 483, 484, 485, 486, 487>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather8f<40, 41, 42, 43, 44, 45, 46, 47>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather8f<104, 105, 106, 107, 108, 109, 110, 111>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather8f<168, 169, 170, 171, 172, 173, 174, 175>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather8f<232, 233, 234, 235, 236, 237, 238, 239>(data);
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = gather8f<296, 297, 298, 299, 300, 301, 302, 303>(data);
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = gather8f<360, 361, 362, 363, 364, 365, 366, 367>(data);
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = gather8f<424, 425, 426, 427, 428, 429, 430, 431>(data);
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = gather8f<488, 489, 490, 491, 492, 493, 494, 495>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather8f<48, 49, 50, 51, 52, 53, 54, 55>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather8f<112, 113, 114, 115, 116, 117, 118, 119>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather8f<176, 177, 178, 179, 180, 181, 182, 183>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather8f<240, 241, 242, 243, 244, 245, 246, 247>(data);
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = gather8f<304, 305, 306, 307, 308, 309, 310, 311>(data);
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = gather8f<368, 369, 370, 371, 372, 373, 374, 375>(data);
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = gather8f<432, 433, 434, 435, 436, 437, 438, 439>(data);
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = gather8f<496, 497, 498, 499, 500, 501, 502, 503>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather8f<56, 57, 58, 59, 60, 61, 62, 63>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather8f<120, 121, 122, 123, 124, 125, 126, 127>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather8f<184, 185, 186, 187, 188, 189, 190, 191>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather8f<248, 249, 250, 251, 252, 253, 254, 255>(data);
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = gather8f<312, 313, 314, 315, 316, 317, 318, 319>(data);
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = gather8f<376, 377, 378, 379, 380, 381, 382, 383>(data);
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = gather8f<440, 441, 442, 443, 444, 445, 446, 447>(data);
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = gather8f<504, 505, 506, 507, 508, 509, 510, 511>(data);
   #elif defined(VEC8D_AGNER) && VECL == 8
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather8d<0, 1, 2, 3, 4, 5, 6, 7>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather8d<64, 65, 66, 67, 68, 69, 70, 71>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather8d<128, 129, 130, 131, 132, 133, 134, 135>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather8d<192, 193, 194, 195, 196, 197, 198, 199>(data);
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = gather8d<256, 257, 258, 259, 260, 261, 262, 263>(data);
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = gather8d<320, 321, 322, 323, 324, 325, 326, 327>(data);
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = gather8d<384, 385, 386, 387, 388, 389, 390, 391>(data);
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = gather8d<448, 449, 450, 451, 452, 453, 454, 455>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather8d<8, 9, 10, 11, 12, 13, 14, 15>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather8d<72, 73, 74, 75, 76, 77, 78, 79>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather8d<136, 137, 138, 139, 140, 141, 142, 143>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather8d<200, 201, 202, 203, 204, 205, 206, 207>(data);
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = gather8d<264, 265, 266, 267, 268, 269, 270, 271>(data);
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = gather8d<328, 329, 330, 331, 332, 333, 334, 335>(data);
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = gather8d<392, 393, 394, 395, 396, 397, 398, 399>(data);
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = gather8d<456, 457, 458, 459, 460, 461, 462, 463>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather8d<16, 17, 18, 19, 20, 21, 22, 23>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather8d<80, 81, 82, 83, 84, 85, 86, 87>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather8d<144, 145, 146, 147, 148, 149, 150, 151>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather8d<208, 209, 210, 211, 212, 213, 214, 215>(data);
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = gather8d<272, 273, 274, 275, 276, 277, 278, 279>(data);
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = gather8d<336, 337, 338, 339, 340, 341, 342, 343>(data);
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = gather8d<400, 401, 402, 403, 404, 405, 406, 407>(data);
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = gather8d<464, 465, 466, 467, 468, 469, 470, 471>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather8d<24, 25, 26, 27, 28, 29, 30, 31>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather8d<88, 89, 90, 91, 92, 93, 94, 95>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather8d<152, 153, 154, 155, 156, 157, 158, 159>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather8d<216, 217, 218, 219, 220, 221, 222, 223>(data);
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = gather8d<280, 281, 282, 283, 284, 285, 286, 287>(data);
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = gather8d<344, 345, 346, 347, 348, 349, 350, 351>(data);
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = gather8d<408, 409, 410, 411, 412, 413, 414, 415>(data);
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = gather8d<472, 473, 474, 475, 476, 477, 478, 479>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather8d<32, 33, 34, 35, 36, 37, 38, 39>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather8d<96, 97, 98, 99, 100, 101, 102, 103>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather8d<160, 161, 162, 163, 164, 165, 166, 167>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather8d<224, 225, 226, 227, 228, 229, 230, 231>(data);
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = gather8d<288, 289, 290, 291, 292, 293, 294, 295>(data);
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = gather8d<352, 353, 354, 355, 356, 357, 358, 359>(data);
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = gather8d<416, 417, 418, 419, 420, 421, 422, 423>(data);
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = gather8d<480, 481, 482, 483, 484, 485, 486, 487>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather8d<40, 41, 42, 43, 44, 45, 46, 47>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather8d<104, 105, 106, 107, 108, 109, 110, 111>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather8d<168, 169, 170, 171, 172, 173, 174, 175>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather8d<232, 233, 234, 235, 236, 237, 238, 239>(data);
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = gather8d<296, 297, 298, 299, 300, 301, 302, 303>(data);
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = gather8d<360, 361, 362, 363, 364, 365, 366, 367>(data);
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = gather8d<424, 425, 426, 427, 428, 429, 430, 431>(data);
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = gather8d<488, 489, 490, 491, 492, 493, 494, 495>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather8d<48, 49, 50, 51, 52, 53, 54, 55>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather8d<112, 113, 114, 115, 116, 117, 118, 119>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather8d<176, 177, 178, 179, 180, 181, 182, 183>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather8d<240, 241, 242, 243, 244, 245, 246, 247>(data);
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = gather8d<304, 305, 306, 307, 308, 309, 310, 311>(data);
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = gather8d<368, 369, 370, 371, 372, 373, 374, 375>(data);
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = gather8d<432, 433, 434, 435, 436, 437, 438, 439>(data);
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = gather8d<496, 497, 498, 499, 500, 501, 502, 503>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather8d<56, 57, 58, 59, 60, 61, 62, 63>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather8d<120, 121, 122, 123, 124, 125, 126, 127>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather8d<184, 185, 186, 187, 188, 189, 190, 191>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather8d<248, 249, 250, 251, 252, 253, 254, 255>(data);
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = gather8d<312, 313, 314, 315, 316, 317, 318, 319>(data);
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = gather8d<376, 377, 378, 379, 380, 381, 382, 383>(data);
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = gather8d<440, 441, 442, 443, 444, 445, 446, 447>(data);
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = gather8d<504, 505, 506, 507, 508, 509, 510, 511>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 8 
   // WID 8, vecl 8
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[64], data[65], data[66], data[67], data[68], data[69], data[70], data[71]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[128], data[129], data[130], data[131], data[132], data[133], data[134], data[135]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[192], data[193], data[194], data[195], data[196], data[197], data[198], data[199]});
         values[i_pcolumnv_b(4, 0, block_k, n_blocks)] = Vec({data[256], data[257], data[258], data[259], data[260], data[261], data[262], data[263]});
         values[i_pcolumnv_b(5, 0, block_k, n_blocks)] = Vec({data[320], data[321], data[322], data[323], data[324], data[325], data[326], data[327]});
         values[i_pcolumnv_b(6, 0, block_k, n_blocks)] = Vec({data[384], data[385], data[386], data[387], data[388], data[389], data[390], data[391]});
         values[i_pcolumnv_b(7, 0, block_k, n_blocks)] = Vec({data[448], data[449], data[450], data[451], data[452], data[453], data[454], data[455]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[72], data[73], data[74], data[75], data[76], data[77], data[78], data[79]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[136], data[137], data[138], data[139], data[140], data[141], data[142], data[143]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[200], data[201], data[202], data[203], data[204], data[205], data[206], data[207]});
         values[i_pcolumnv_b(4, 1, block_k, n_blocks)] = Vec({data[264], data[265], data[266], data[267], data[268], data[269], data[270], data[271]});
         values[i_pcolumnv_b(5, 1, block_k, n_blocks)] = Vec({data[328], data[329], data[330], data[331], data[332], data[333], data[334], data[335]});
         values[i_pcolumnv_b(6, 1, block_k, n_blocks)] = Vec({data[392], data[393], data[394], data[395], data[396], data[397], data[398], data[399]});
         values[i_pcolumnv_b(7, 1, block_k, n_blocks)] = Vec({data[456], data[457], data[458], data[459], data[460], data[461], data[462], data[463]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[80], data[81], data[82], data[83], data[84], data[85], data[86], data[87]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[144], data[145], data[146], data[147], data[148], data[149], data[150], data[151]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[208], data[209], data[210], data[211], data[212], data[213], data[214], data[215]});
         values[i_pcolumnv_b(4, 2, block_k, n_blocks)] = Vec({data[272], data[273], data[274], data[275], data[276], data[277], data[278], data[279]});
         values[i_pcolumnv_b(5, 2, block_k, n_blocks)] = Vec({data[336], data[337], data[338], data[339], data[340], data[341], data[342], data[343]});
         values[i_pcolumnv_b(6, 2, block_k, n_blocks)] = Vec({data[400], data[401], data[402], data[403], data[404], data[405], data[406], data[407]});
         values[i_pcolumnv_b(7, 2, block_k, n_blocks)] = Vec({data[464], data[465], data[466], data[467], data[468], data[469], data[470], data[471]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[88], data[89], data[90], data[91], data[92], data[93], data[94], data[95]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[152], data[153], data[154], data[155], data[156], data[157], data[158], data[159]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[216], data[217], data[218], data[219], data[220], data[221], data[222], data[223]});
         values[i_pcolumnv_b(4, 3, block_k, n_blocks)] = Vec({data[280], data[281], data[282], data[283], data[284], data[285], data[286], data[287]});
         values[i_pcolumnv_b(5, 3, block_k, n_blocks)] = Vec({data[344], data[345], data[346], data[347], data[348], data[349], data[350], data[351]});
         values[i_pcolumnv_b(6, 3, block_k, n_blocks)] = Vec({data[408], data[409], data[410], data[411], data[412], data[413], data[414], data[415]});
         values[i_pcolumnv_b(7, 3, block_k, n_blocks)] = Vec({data[472], data[473], data[474], data[475], data[476], data[477], data[478], data[479]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35], data[36], data[37], data[38], data[39]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[96], data[97], data[98], data[99], data[100], data[101], data[102], data[103]});
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = Vec({data[160], data[161], data[162], data[163], data[164], data[165], data[166], data[167]});
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = Vec({data[224], data[225], data[226], data[227], data[228], data[229], data[230], data[231]});
         values[i_pcolumnv_b(4, 4, block_k, n_blocks)] = Vec({data[288], data[289], data[290], data[291], data[292], data[293], data[294], data[295]});
         values[i_pcolumnv_b(5, 4, block_k, n_blocks)] = Vec({data[352], data[353], data[354], data[355], data[356], data[357], data[358], data[359]});
         values[i_pcolumnv_b(6, 4, block_k, n_blocks)] = Vec({data[416], data[417], data[418], data[419], data[420], data[421], data[422], data[423]});
         values[i_pcolumnv_b(7, 4, block_k, n_blocks)] = Vec({data[480], data[481], data[482], data[483], data[484], data[485], data[486], data[487]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43], data[44], data[45], data[46], data[47]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[104], data[105], data[106], data[107], data[108], data[109], data[110], data[111]});
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = Vec({data[168], data[169], data[170], data[171], data[172], data[173], data[174], data[175]});
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = Vec({data[232], data[233], data[234], data[235], data[236], data[237], data[238], data[239]});
         values[i_pcolumnv_b(4, 5, block_k, n_blocks)] = Vec({data[296], data[297], data[298], data[299], data[300], data[301], data[302], data[303]});
         values[i_pcolumnv_b(5, 5, block_k, n_blocks)] = Vec({data[360], data[361], data[362], data[363], data[364], data[365], data[366], data[367]});
         values[i_pcolumnv_b(6, 5, block_k, n_blocks)] = Vec({data[424], data[425], data[426], data[427], data[428], data[429], data[430], data[431]});
         values[i_pcolumnv_b(7, 5, block_k, n_blocks)] = Vec({data[488], data[489], data[490], data[491], data[492], data[493], data[494], data[495]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[48], data[49], data[50], data[51], data[52], data[53], data[54], data[55]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[112], data[113], data[114], data[115], data[116], data[117], data[118], data[119]});
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = Vec({data[176], data[177], data[178], data[179], data[180], data[181], data[182], data[183]});
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = Vec({data[240], data[241], data[242], data[243], data[244], data[245], data[246], data[247]});
         values[i_pcolumnv_b(4, 6, block_k, n_blocks)] = Vec({data[304], data[305], data[306], data[307], data[308], data[309], data[310], data[311]});
         values[i_pcolumnv_b(5, 6, block_k, n_blocks)] = Vec({data[368], data[369], data[370], data[371], data[372], data[373], data[374], data[375]});
         values[i_pcolumnv_b(6, 6, block_k, n_blocks)] = Vec({data[432], data[433], data[434], data[435], data[436], data[437], data[438], data[439]});
         values[i_pcolumnv_b(7, 6, block_k, n_blocks)] = Vec({data[496], data[497], data[498], data[499], data[500], data[501], data[502], data[503]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[56], data[57], data[58], data[59], data[60], data[61], data[62], data[63]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[120], data[121], data[122], data[123], data[124], data[125], data[126], data[127]});
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = Vec({data[184], data[185], data[186], data[187], data[188], data[189], data[190], data[191]});
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = Vec({data[248], data[249], data[250], data[251], data[252], data[253], data[254], data[255]});
         values[i_pcolumnv_b(4, 7, block_k, n_blocks)] = Vec({data[312], data[313], data[314], data[315], data[316], data[317], data[318], data[319]});
         values[i_pcolumnv_b(5, 7, block_k, n_blocks)] = Vec({data[376], data[377], data[378], data[379], data[380], data[381], data[382], data[383]});
         values[i_pcolumnv_b(6, 7, block_k, n_blocks)] = Vec({data[440], data[441], data[442], data[443], data[444], data[445], data[446], data[447]});
         values[i_pcolumnv_b(7, 7, block_k, n_blocks)] = Vec({data[504], data[505], data[506], data[507], data[508], data[509], data[510], data[511]});
   #elif defined(VEC16F_AGNER) && VECL == 16
   // WID 8, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = gather16f<0, 1, 2, 3, 4, 5, 6, 7, 64, 65, 66, 67, 68, 69, 70, 71>(data);
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = gather16f<128, 129, 130, 131, 132, 133, 134, 135, 192, 193, 194, 195, 196, 197, 198, 199>(data);
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = gather16f<256, 257, 258, 259, 260, 261, 262, 263, 320, 321, 322, 323, 324, 325, 326, 327>(data);
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = gather16f<384, 385, 386, 387, 388, 389, 390, 391, 448, 449, 450, 451, 452, 453, 454, 455>(data);
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = gather16f<8, 9, 10, 11, 12, 13, 14, 15, 72, 73, 74, 75, 76, 77, 78, 79>(data);
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = gather16f<136, 137, 138, 139, 140, 141, 142, 143, 200, 201, 202, 203, 204, 205, 206, 207>(data);
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = gather16f<264, 265, 266, 267, 268, 269, 270, 271, 328, 329, 330, 331, 332, 333, 334, 335>(data);
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = gather16f<392, 393, 394, 395, 396, 397, 398, 399, 456, 457, 458, 459, 460, 461, 462, 463>(data);
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = gather16f<16, 17, 18, 19, 20, 21, 22, 23, 80, 81, 82, 83, 84, 85, 86, 87>(data);
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = gather16f<144, 145, 146, 147, 148, 149, 150, 151, 208, 209, 210, 211, 212, 213, 214, 215>(data);
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = gather16f<272, 273, 274, 275, 276, 277, 278, 279, 336, 337, 338, 339, 340, 341, 342, 343>(data);
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = gather16f<400, 401, 402, 403, 404, 405, 406, 407, 464, 465, 466, 467, 468, 469, 470, 471>(data);
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = gather16f<24, 25, 26, 27, 28, 29, 30, 31, 88, 89, 90, 91, 92, 93, 94, 95>(data);
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = gather16f<152, 153, 154, 155, 156, 157, 158, 159, 216, 217, 218, 219, 220, 221, 222, 223>(data);
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = gather16f<280, 281, 282, 283, 284, 285, 286, 287, 344, 345, 346, 347, 348, 349, 350, 351>(data);
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = gather16f<408, 409, 410, 411, 412, 413, 414, 415, 472, 473, 474, 475, 476, 477, 478, 479>(data);
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = gather16f<32, 33, 34, 35, 36, 37, 38, 39, 96, 97, 98, 99, 100, 101, 102, 103>(data);
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = gather16f<160, 161, 162, 163, 164, 165, 166, 167, 224, 225, 226, 227, 228, 229, 230, 231>(data);
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = gather16f<288, 289, 290, 291, 292, 293, 294, 295, 352, 353, 354, 355, 356, 357, 358, 359>(data);
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = gather16f<416, 417, 418, 419, 420, 421, 422, 423, 480, 481, 482, 483, 484, 485, 486, 487>(data);
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = gather16f<40, 41, 42, 43, 44, 45, 46, 47, 104, 105, 106, 107, 108, 109, 110, 111>(data);
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = gather16f<168, 169, 170, 171, 172, 173, 174, 175, 232, 233, 234, 235, 236, 237, 238, 239>(data);
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = gather16f<296, 297, 298, 299, 300, 301, 302, 303, 360, 361, 362, 363, 364, 365, 366, 367>(data);
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = gather16f<424, 425, 426, 427, 428, 429, 430, 431, 488, 489, 490, 491, 492, 493, 494, 495>(data);
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = gather16f<48, 49, 50, 51, 52, 53, 54, 55, 112, 113, 114, 115, 116, 117, 118, 119>(data);
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = gather16f<176, 177, 178, 179, 180, 181, 182, 183, 240, 241, 242, 243, 244, 245, 246, 247>(data);
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = gather16f<304, 305, 306, 307, 308, 309, 310, 311, 368, 369, 370, 371, 372, 373, 374, 375>(data);
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = gather16f<432, 433, 434, 435, 436, 437, 438, 439, 496, 497, 498, 499, 500, 501, 502, 503>(data);
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = gather16f<56, 57, 58, 59, 60, 61, 62, 63, 120, 121, 122, 123, 124, 125, 126, 127>(data);
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = gather16f<184, 185, 186, 187, 188, 189, 190, 191, 248, 249, 250, 251, 252, 253, 254, 255>(data);
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = gather16f<312, 313, 314, 315, 316, 317, 318, 319, 376, 377, 378, 379, 380, 381, 382, 383>(data);
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = gather16f<440, 441, 442, 443, 444, 445, 446, 447, 504, 505, 506, 507, 508, 509, 510, 511>(data);
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 16 
   // WID 8, vecl 16
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[64], data[65], data[66], data[67], data[68], data[69], data[70], data[71]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[128], data[129], data[130], data[131], data[132], data[133], data[134], data[135], data[192], data[193], data[194], data[195], data[196], data[197], data[198], data[199]});
         values[i_pcolumnv_b(2, 0, block_k, n_blocks)] = Vec({data[256], data[257], data[258], data[259], data[260], data[261], data[262], data[263], data[320], data[321], data[322], data[323], data[324], data[325], data[326], data[327]});
         values[i_pcolumnv_b(3, 0, block_k, n_blocks)] = Vec({data[384], data[385], data[386], data[387], data[388], data[389], data[390], data[391], data[448], data[449], data[450], data[451], data[452], data[453], data[454], data[455]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[72], data[73], data[74], data[75], data[76], data[77], data[78], data[79]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[136], data[137], data[138], data[139], data[140], data[141], data[142], data[143], data[200], data[201], data[202], data[203], data[204], data[205], data[206], data[207]});
         values[i_pcolumnv_b(2, 1, block_k, n_blocks)] = Vec({data[264], data[265], data[266], data[267], data[268], data[269], data[270], data[271], data[328], data[329], data[330], data[331], data[332], data[333], data[334], data[335]});
         values[i_pcolumnv_b(3, 1, block_k, n_blocks)] = Vec({data[392], data[393], data[394], data[395], data[396], data[397], data[398], data[399], data[456], data[457], data[458], data[459], data[460], data[461], data[462], data[463]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23], data[80], data[81], data[82], data[83], data[84], data[85], data[86], data[87]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[144], data[145], data[146], data[147], data[148], data[149], data[150], data[151], data[208], data[209], data[210], data[211], data[212], data[213], data[214], data[215]});
         values[i_pcolumnv_b(2, 2, block_k, n_blocks)] = Vec({data[272], data[273], data[274], data[275], data[276], data[277], data[278], data[279], data[336], data[337], data[338], data[339], data[340], data[341], data[342], data[343]});
         values[i_pcolumnv_b(3, 2, block_k, n_blocks)] = Vec({data[400], data[401], data[402], data[403], data[404], data[405], data[406], data[407], data[464], data[465], data[466], data[467], data[468], data[469], data[470], data[471]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31], data[88], data[89], data[90], data[91], data[92], data[93], data[94], data[95]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[152], data[153], data[154], data[155], data[156], data[157], data[158], data[159], data[216], data[217], data[218], data[219], data[220], data[221], data[222], data[223]});
         values[i_pcolumnv_b(2, 3, block_k, n_blocks)] = Vec({data[280], data[281], data[282], data[283], data[284], data[285], data[286], data[287], data[344], data[345], data[346], data[347], data[348], data[349], data[350], data[351]});
         values[i_pcolumnv_b(3, 3, block_k, n_blocks)] = Vec({data[408], data[409], data[410], data[411], data[412], data[413], data[414], data[415], data[472], data[473], data[474], data[475], data[476], data[477], data[478], data[479]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35], data[36], data[37], data[38], data[39], data[96], data[97], data[98], data[99], data[100], data[101], data[102], data[103]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[160], data[161], data[162], data[163], data[164], data[165], data[166], data[167], data[224], data[225], data[226], data[227], data[228], data[229], data[230], data[231]});
         values[i_pcolumnv_b(2, 4, block_k, n_blocks)] = Vec({data[288], data[289], data[290], data[291], data[292], data[293], data[294], data[295], data[352], data[353], data[354], data[355], data[356], data[357], data[358], data[359]});
         values[i_pcolumnv_b(3, 4, block_k, n_blocks)] = Vec({data[416], data[417], data[418], data[419], data[420], data[421], data[422], data[423], data[480], data[481], data[482], data[483], data[484], data[485], data[486], data[487]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43], data[44], data[45], data[46], data[47], data[104], data[105], data[106], data[107], data[108], data[109], data[110], data[111]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[168], data[169], data[170], data[171], data[172], data[173], data[174], data[175], data[232], data[233], data[234], data[235], data[236], data[237], data[238], data[239]});
         values[i_pcolumnv_b(2, 5, block_k, n_blocks)] = Vec({data[296], data[297], data[298], data[299], data[300], data[301], data[302], data[303], data[360], data[361], data[362], data[363], data[364], data[365], data[366], data[367]});
         values[i_pcolumnv_b(3, 5, block_k, n_blocks)] = Vec({data[424], data[425], data[426], data[427], data[428], data[429], data[430], data[431], data[488], data[489], data[490], data[491], data[492], data[493], data[494], data[495]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[48], data[49], data[50], data[51], data[52], data[53], data[54], data[55], data[112], data[113], data[114], data[115], data[116], data[117], data[118], data[119]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[176], data[177], data[178], data[179], data[180], data[181], data[182], data[183], data[240], data[241], data[242], data[243], data[244], data[245], data[246], data[247]});
         values[i_pcolumnv_b(2, 6, block_k, n_blocks)] = Vec({data[304], data[305], data[306], data[307], data[308], data[309], data[310], data[311], data[368], data[369], data[370], data[371], data[372], data[373], data[374], data[375]});
         values[i_pcolumnv_b(3, 6, block_k, n_blocks)] = Vec({data[432], data[433], data[434], data[435], data[436], data[437], data[438], data[439], data[496], data[497], data[498], data[499], data[500], data[501], data[502], data[503]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[56], data[57], data[58], data[59], data[60], data[61], data[62], data[63], data[120], data[121], data[122], data[123], data[124], data[125], data[126], data[127]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[184], data[185], data[186], data[187], data[188], data[189], data[190], data[191], data[248], data[249], data[250], data[251], data[252], data[253], data[254], data[255]});
         values[i_pcolumnv_b(2, 7, block_k, n_blocks)] = Vec({data[312], data[313], data[314], data[315], data[316], data[317], data[318], data[319], data[376], data[377], data[378], data[379], data[380], data[381], data[382], data[383]});
         values[i_pcolumnv_b(3, 7, block_k, n_blocks)] = Vec({data[440], data[441], data[442], data[443], data[444], data[445], data[446], data[447], data[504], data[505], data[506], data[507], data[508], data[509], data[510], data[511]});
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 32 
   // WID 8, vecl 32
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[64], data[65], data[66], data[67], data[68], data[69], data[70], data[71], data[128], data[129], data[130], data[131], data[132], data[133], data[134], data[135], data[192], data[193], data[194], data[195], data[196], data[197], data[198], data[199]});
         values[i_pcolumnv_b(1, 0, block_k, n_blocks)] = Vec({data[256], data[257], data[258], data[259], data[260], data[261], data[262], data[263], data[320], data[321], data[322], data[323], data[324], data[325], data[326], data[327], data[384], data[385], data[386], data[387], data[388], data[389], data[390], data[391], data[448], data[449], data[450], data[451], data[452], data[453], data[454], data[455]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[72], data[73], data[74], data[75], data[76], data[77], data[78], data[79], data[136], data[137], data[138], data[139], data[140], data[141], data[142], data[143], data[200], data[201], data[202], data[203], data[204], data[205], data[206], data[207]});
         values[i_pcolumnv_b(1, 1, block_k, n_blocks)] = Vec({data[264], data[265], data[266], data[267], data[268], data[269], data[270], data[271], data[328], data[329], data[330], data[331], data[332], data[333], data[334], data[335], data[392], data[393], data[394], data[395], data[396], data[397], data[398], data[399], data[456], data[457], data[458], data[459], data[460], data[461], data[462], data[463]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23], data[80], data[81], data[82], data[83], data[84], data[85], data[86], data[87], data[144], data[145], data[146], data[147], data[148], data[149], data[150], data[151], data[208], data[209], data[210], data[211], data[212], data[213], data[214], data[215]});
         values[i_pcolumnv_b(1, 2, block_k, n_blocks)] = Vec({data[272], data[273], data[274], data[275], data[276], data[277], data[278], data[279], data[336], data[337], data[338], data[339], data[340], data[341], data[342], data[343], data[400], data[401], data[402], data[403], data[404], data[405], data[406], data[407], data[464], data[465], data[466], data[467], data[468], data[469], data[470], data[471]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31], data[88], data[89], data[90], data[91], data[92], data[93], data[94], data[95], data[152], data[153], data[154], data[155], data[156], data[157], data[158], data[159], data[216], data[217], data[218], data[219], data[220], data[221], data[222], data[223]});
         values[i_pcolumnv_b(1, 3, block_k, n_blocks)] = Vec({data[280], data[281], data[282], data[283], data[284], data[285], data[286], data[287], data[344], data[345], data[346], data[347], data[348], data[349], data[350], data[351], data[408], data[409], data[410], data[411], data[412], data[413], data[414], data[415], data[472], data[473], data[474], data[475], data[476], data[477], data[478], data[479]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35], data[36], data[37], data[38], data[39], data[96], data[97], data[98], data[99], data[100], data[101], data[102], data[103], data[160], data[161], data[162], data[163], data[164], data[165], data[166], data[167], data[224], data[225], data[226], data[227], data[228], data[229], data[230], data[231]});
         values[i_pcolumnv_b(1, 4, block_k, n_blocks)] = Vec({data[288], data[289], data[290], data[291], data[292], data[293], data[294], data[295], data[352], data[353], data[354], data[355], data[356], data[357], data[358], data[359], data[416], data[417], data[418], data[419], data[420], data[421], data[422], data[423], data[480], data[481], data[482], data[483], data[484], data[485], data[486], data[487]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43], data[44], data[45], data[46], data[47], data[104], data[105], data[106], data[107], data[108], data[109], data[110], data[111], data[168], data[169], data[170], data[171], data[172], data[173], data[174], data[175], data[232], data[233], data[234], data[235], data[236], data[237], data[238], data[239]});
         values[i_pcolumnv_b(1, 5, block_k, n_blocks)] = Vec({data[296], data[297], data[298], data[299], data[300], data[301], data[302], data[303], data[360], data[361], data[362], data[363], data[364], data[365], data[366], data[367], data[424], data[425], data[426], data[427], data[428], data[429], data[430], data[431], data[488], data[489], data[490], data[491], data[492], data[493], data[494], data[495]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[48], data[49], data[50], data[51], data[52], data[53], data[54], data[55], data[112], data[113], data[114], data[115], data[116], data[117], data[118], data[119], data[176], data[177], data[178], data[179], data[180], data[181], data[182], data[183], data[240], data[241], data[242], data[243], data[244], data[245], data[246], data[247]});
         values[i_pcolumnv_b(1, 6, block_k, n_blocks)] = Vec({data[304], data[305], data[306], data[307], data[308], data[309], data[310], data[311], data[368], data[369], data[370], data[371], data[372], data[373], data[374], data[375], data[432], data[433], data[434], data[435], data[436], data[437], data[438], data[439], data[496], data[497], data[498], data[499], data[500], data[501], data[502], data[503]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[56], data[57], data[58], data[59], data[60], data[61], data[62], data[63], data[120], data[121], data[122], data[123], data[124], data[125], data[126], data[127], data[184], data[185], data[186], data[187], data[188], data[189], data[190], data[191], data[248], data[249], data[250], data[251], data[252], data[253], data[254], data[255]});
         values[i_pcolumnv_b(1, 7, block_k, n_blocks)] = Vec({data[312], data[313], data[314], data[315], data[316], data[317], data[318], data[319], data[376], data[377], data[378], data[379], data[380], data[381], data[382], data[383], data[440], data[441], data[442], data[443], data[444], data[445], data[446], data[447], data[504], data[505], data[506], data[507], data[508], data[509], data[510], data[511]});
   #elif defined(VEC_FALLBACK_GENERIC) && VECL == 64 
   // WID 8, vecl 64
         values[i_pcolumnv_b(0, 0, block_k, n_blocks)] = Vec({data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[64], data[65], data[66], data[67], data[68], data[69], data[70], data[71], data[128], data[129], data[130], data[131], data[132], data[133], data[134], data[135], data[192], data[193], data[194], data[195], data[196], data[197], data[198], data[199], data[256], data[257], data[258], data[259], data[260], data[261], data[262], data[263], data[320], data[321], data[322], data[323], data[324], data[325], data[326], data[327], data[384], data[385], data[386], data[387], data[388], data[389], data[390], data[391], data[448], data[449], data[450], data[451], data[452], data[453], data[454], data[455]});
         values[i_pcolumnv_b(0, 1, block_k, n_blocks)] = Vec({data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[72], data[73], data[74], data[75], data[76], data[77], data[78], data[79], data[136], data[137], data[138], data[139], data[140], data[141], data[142], data[143], data[200], data[201], data[202], data[203], data[204], data[205], data[206], data[207], data[264], data[265], data[266], data[267], data[268], data[269], data[270], data[271], data[328], data[329], data[330], data[331], data[332], data[333], data[334], data[335], data[392], data[393], data[394], data[395], data[396], data[397], data[398], data[399], data[456], data[457], data[458], data[459], data[460], data[461], data[462], data[463]});
         values[i_pcolumnv_b(0, 2, block_k, n_blocks)] = Vec({data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23], data[80], data[81], data[82], data[83], data[84], data[85], data[86], data[87], data[144], data[145], data[146], data[147], data[148], data[149], data[150], data[151], data[208], data[209], data[210], data[211], data[212], data[213], data[214], data[215], data[272], data[273], data[274], data[275], data[276], data[277], data[278], data[279], data[336], data[337], data[338], data[339], data[340], data[341], data[342], data[343], data[400], data[401], data[402], data[403], data[404], data[405], data[406], data[407], data[464], data[465], data[466], data[467], data[468], data[469], data[470], data[471]});
         values[i_pcolumnv_b(0, 3, block_k, n_blocks)] = Vec({data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31], data[88], data[89], data[90], data[91], data[92], data[93], data[94], data[95], data[152], data[153], data[154], data[155], data[156], data[157], data[158], data[159], data[216], data[217], data[218], data[219], data[220], data[221], data[222], data[223], data[280], data[281], data[282], data[283], data[284], data[285], data[286], data[287], data[344], data[345], data[346], data[347], data[348], data[349], data[350], data[351], data[408], data[409], data[410], data[411], data[412], data[413], data[414], data[415], data[472], data[473], data[474], data[475], data[476], data[477], data[478], data[479]});
         values[i_pcolumnv_b(0, 4, block_k, n_blocks)] = Vec({data[32], data[33], data[34], data[35], data[36], data[37], data[38], data[39], data[96], data[97], data[98], data[99], data[100], data[101], data[102], data[103], data[160], data[161], data[162], data[163], data[164], data[165], data[166], data[167], data[224], data[225], data[226], data[227], data[228], data[229], data[230], data[231], data[288], data[289], data[290], data[291], data[292], data[293], data[294], data[295], data[352], data[353], data[354], data[355], data[356], data[357], data[358], data[359], data[416], data[417], data[418], data[419], data[420], data[421], data[422], data[423], data[480], data[481], data[482], data[483], data[484], data[485], data[486], data[487]});
         values[i_pcolumnv_b(0, 5, block_k, n_blocks)] = Vec({data[40], data[41], data[42], data[43], data[44], data[45], data[46], data[47], data[104], data[105], data[106], data[107], data[108], data[109], data[110], data[111], data[168], data[169], data[170], data[171], data[172], data[173], data[174], data[175], data[232], data[233], data[234], data[235], data[236], data[237], data[238], data[239], data[296], data[297], data[298], data[299], data[300], data[301], data[302], data[303], data[360], data[361], data[362], data[363], data[364], data[365], data[366], data[367], data[424], data[425], data[426], data[427], data[428], data[429], data[430], data[431], data[488], data[489], data[490], data[491], data[492], data[493], data[494], data[495]});
         values[i_pcolumnv_b(0, 6, block_k, n_blocks)] = Vec({data[48], data[49], data[50], data[51], data[52], data[53], data[54], data[55], data[112], data[113], data[114], data[115], data[116], data[117], data[118], data[119], data[176], data[177], data[178], data[179], data[180], data[181], data[182], data[183], data[240], data[241], data[242], data[243], data[244], data[245], data[246], data[247], data[304], data[305], data[306], data[307], data[308], data[309], data[310], data[311], data[368], data[369], data[370], data[371], data[372], data[373], data[374], data[375], data[432], data[433], data[434], data[435], data[436], data[437], data[438], data[439], data[496], data[497], data[498], data[499], data[500], data[501], data[502], data[503]});
         values[i_pcolumnv_b(0, 7, block_k, n_blocks)] = Vec({data[56], data[57], data[58], data[59], data[60], data[61], data[62], data[63], data[120], data[121], data[122], data[123], data[124], data[125], data[126], data[127], data[184], data[185], data[186], data[187], data[188], data[189], data[190], data[191], data[248], data[249], data[250], data[251], data[252], data[253], data[254], data[255], data[312], data[313], data[314], data[315], data[316], data[317], data[318], data[319], data[376], data[377], data[378], data[379], data[380], data[381], data[382], data[383], data[440], data[441], data[442], data[443], data[444], data[445], data[446], data[447], data[504], data[505], data[506], data[507], data[508], data[509], data[510], data[511]});
   #else // Fell through, never fall into this particular pit again
   	std::cerr << "Undefined VECTORCLASS flag or implementation missing in loadColumnBlockData() before " << __FILE__ << ":" << __LINE__ << std::endl;
   	abort();
   #endif
         //zero old output data
         for (uint i=0; i<WID3; ++i) {
            data[i]=0;
         }
      }
    }
   #else
     std::cerr << "Undefined WID encountered in " << __FILE__ << ":" << __LINE__ << std::endl;
     abort();
   #endif
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
