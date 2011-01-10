#ifndef DATAREDUCTIONOPERATOR_H
#define DATAREDUCTIONOPERATOR_H

#include "definitions.h"
#include "cell_spatial.h"

namespace DRO {

   /** DataReductionOperator defines a base class for reducing simulation data
    * (six-dimensional distribution function) into more compact variables, e.g. 
    * scalar fields, which can be written into disk and visualized.
    * 
    * The intention is that each DataReductionOperator stores the reduced data 
    * into internal variables, whose values are appended into a byte array when 
    * appendReducedData is called.
    * 
    * If needed, the user can write his or her own DataReductionOperators, which 
    * are loaded when the simulation initializes.
    */
   class DataReductionOperator {
    public:
      DataReductionOperator();
      virtual ~DataReductionOperator();
      
      virtual bool appendReducedData(unsigned char* const byteArray);
      virtual unsigned char getElementByteSize() const;
      virtual std::string getName() const;
      virtual unsigned int getOutputByteSize() const;
      virtual unsigned char getVariableType() const;
      virtual bool reduceData(const Real* const avgs,const Real* const blockParams);
      virtual bool setSpatialCell(const SpatialCell& cell);
   
    protected:
      
   };
   
   class MPIrank: public DataReductionOperator {
    public:
      MPIrank();
      ~MPIrank();
      
      bool appendReducedData(unsigned char* const byteArray);
      unsigned char getElementByteSize() const;
      std::string getName() const;
      unsigned char getVariableType() const;
      bool reduceData(const Real* const avgs,const Real* const blockParams);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rank;
   };

   class VariableB: public DataReductionOperator {
    public:
      VariableB();
      ~VariableB();
      
      bool appendReducedData(unsigned char* const byteArray);
      unsigned char getElementByteSize() const;
      std::string getName() const;
      unsigned char getVariableType() const;
      bool reduceData(const Real* const avgs,const Real* const blockParams);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real Bx;
      Real By;
      Real Bz;
   };
   
   class VariableE: public DataReductionOperator {
    public:
      VariableE();
      ~VariableE();
      
      bool appendReducedData(unsigned char* const byteArray);
      unsigned char getElementByteSize() const;
      std::string getName() const;
      unsigned char getVariableType() const;
      bool reduceData(const Real* const avgs,const Real* const blockParams);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real Ex;
      Real Ey;
      Real Ez;
   };
   
   class VariableRho: public DataReductionOperator {
    public:
      VariableRho();
      ~VariableRho();
      
      bool appendReducedData(unsigned char* const byteArray);
      unsigned char getElementByteSize() const;
      std::string getName() const;
      unsigned char getVariableType() const;
      bool reduceData(const Real* const avgs,const Real* const blockParams);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rho;
   };

   class VariableRhoV: public DataReductionOperator {
    public:
      VariableRhoV();
      ~VariableRhoV();
      
      bool appendReducedData(unsigned char* const byteArray);
      unsigned char getElementByteSize() const;
      std::string getName() const;
      unsigned char getVariableType() const;
      bool reduceData(const Real* const avgs,const Real* const blockParams);
      bool setSpatialCell(const SpatialCell& cell);
      
    protected:
      Real rhovx;
      Real rhovy;
      Real rhovz;
   };
   
} // namespace DRO

#endif
   