

inline __device__ real velocityFluxX(real avg_neg,real avg_pos,real* cellParams,real* blockParams) {

   creal VY = blockParams[BlockParams::VYCRD] + (threadIdx.z+0.5f)*blockParams[BlockParams::DVY];
   creal VZ = blockParams[BlockParams::VZCRD] + (threadIdx.y+0.5f)*blockParams[BlockParams::DVZ];
   creal EX = cellParams[CellParams::EX];
   creal BY = cellParams[CellParams::BY];
   creal BZ = cellParams[CellParams::BZ];
   
   creal AX = blockParams[BlockParams::Q_PER_M]*(EX + VY*BZ - VZ*BY);
   return 0.5f*AX*(avg_neg + avg_pos) - 0.5f*fabsf(AX)*(avg_pos-avg_neg);
}

inline __device__ real velocityFluxY(real avg_neg,real avg_pos,real* cellParams,real* blockParams) {
   /*
    creal VX = blockParams[BlockParams::VXCRD] + (threadIdx.x+0.5f)*blockParams[BlockParams::DVX];
    creal VZ = blockParams[BlockParams::VZCRD] + (threadIdx.y+0.5f)*blockParams[BlockParams::DVZ];
    creal EY = cellParams[CellParams::EY];
    creal BX = cellParams[CellParams::BX];
    creal BZ = cellParams[CellParams::BZ];
    * 
    creal AY = blockParams[BlockParams::Q_PER_M]*(EY + VZ*BX - VX*BZ);
    return 0.5f*AY*(avg_neg + avg_pos) - 0.5f*fabsf(AY)*(avg_pos-avg_neg);
    */
   return 0.0f;
}

inline __device__ real velocityFluxZ(real avg_neg,real avg_pos,real* cellParams,real* blockParams) {
   /*
    creal VX = blockParams[BlockParams::VXCRD] + (threadIdx.x+0.5f)*blockParams[BlockParams::DVX];
    creal VY = blockParams[BlockParams::VYCRD] + (threadIdx.y+0.5f)*blockParams[BlockParams::DVY];
    creal EZ = cellParams[CellParams::EZ];
    creal BX = cellParams[CellParams::BX];
    creal BY = cellParams[CellParams::BY];
 
    creal AZ = blockParams[BlockParams::Q_PER_M]*(EZ + VX*BY - VY*BX);
    return 0.5f*AZ*(avg_neg + avg_pos) - 0.5f*fabsf(AZ)*(avg_pos-avg_neg);
    */
   return 0.0f;
}

inline __device__ real spatialFluxX(real avg_neg,real avg_pos,real* blockParams) {
   creal VX = blockParams[BlockParams::VXCRD] + (threadIdx.x+0.5f)*blockParams[BlockParams::DVX];
   return 0.5f*VX*(avg_neg + avg_pos) - 0.5f*fabsf(VX)*(avg_pos-avg_neg);
}

inline __device__ real spatialFluxY(real avg_neg,real avg_pos,real* blockParams) {
   creal VY = blockParams[BlockParams::VYCRD] + (threadIdx.y+0.5f)*blockParams[BlockParams::DVY];
   return 0.5f*VY*(avg_neg + avg_pos) - 0.5f*fabsf(VY)*(avg_pos-avg_neg);
}

inline __device__ real spatialFluxZ(real avg_neg,real avg_pos,real* blockParams,uint z_ind) {
   creal VZ = blockParams[BlockParams::VZCRD] + (z_ind+0.5f)*blockParams[BlockParams::DVZ];
   return 0.5f*VZ*(avg_neg + avg_pos) - 0.5f*fabsf(VZ)*(avg_pos-avg_neg);
}


