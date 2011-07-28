#ifndef DATAREDUCTIONOPERATOR_H
#define DATAREDUCTIONOPERATOR_H

#include <vector>
#include "definitions.h"
#include "cell_spatial.h"

namespace DRO {

   /** DRO::DataReductionOperator defines a base class for reducing simulation data
    * (six-dimensional distribution function) into more compact variables, e.g. 
    * scalar fields, which can be written into file(s) and visualized.
    * 
    * The intention is that each DRO::DataReductionOperator stores the reduced data 
    * into internal variables, whose values are written into a byte array when 
    * DRO::DataReductionOperator::appendReducedData is called.
    * 
    * If needed, a user can write his or her own DRO::DataReductionOperators, which 
    * are loaded when the simulation initializes.
    */
   class DataReductionOperator {
    public:
      DataReductionOperator();
      virtual ~DataReductionOperator();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      virtual bool setSpatialCell(const SpatialCell& cell);
      
    protected:
	
   };
   
   class MPIrank: public DataReductionOperator {
    public:
      MPIrank();
      ~MPIrank();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rank;
      int mpiRank;
   };

   class VariableB: public DataReductionOperator {
    public:
      VariableB();
      ~VariableB();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real Bx;
      Real By;
      Real Bz;
      Real* B;
   };

   class VariableVolB: public DataReductionOperator {
    public:
      VariableVolB();
      ~VariableVolB();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real* B;
   };
   
   class VariableE: public DataReductionOperator {
    public:
      VariableE();
      ~VariableE();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real Ex;
      Real Ey;
      Real Ez;
      Real* E;
   };

   class VariableVolE: public DataReductionOperator {
    public:
      VariableVolE();
      ~VariableVolE();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real* E;
   };
   
   class VariableRho: public DataReductionOperator {
    public:
      VariableRho();
      ~VariableRho();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rho;
   };

   class VariableRhoV: public DataReductionOperator {
    public:
      VariableRhoV();
      ~VariableRhoV();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rhovx;
      Real rhovy;
      Real rhovz;
      Real* rhov;
   };
} // namespace DRO

#endif

