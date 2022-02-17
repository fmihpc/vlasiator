#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cstring>
#include <cmath>
#include <errno.h>

using namespace std;

typedef double Real;
constexpr static int productionNumAccEnergies = 60;
constexpr static int productionNumTemperatures = 60;
constexpr static int productionNumParticleEnergies = 100;
constexpr static Real productionMinAccEnergy = 0.1; // keV
constexpr static Real productionMaxAccEnergy = 100.; // keV
constexpr static Real productionMinTemperature = 0.1; // keV
constexpr static Real productionMaxTemperature = 100.; // keV
constexpr static int numAtmosphereLevels = 20;

static const Real kB = 1.380649e-23;
static const Real CHARGE = 1.602176634e-19;
static const Real MASS_ELECTRON = 9.1093837015e-31;
static const Real MASS_PROTON = 1.67262158e-27;
static const Real recombAlpha = 2.4e-13; // m³/s

std::array< Real, productionNumParticleEnergies+1 > particle_energy;
std::array< Real, productionNumParticleEnergies > differentialFlux;
struct AtmosphericLayer {
   Real altitude;
   Real nui;
   Real nue;
   Real density;
   Real depth; // integrated density from the top of the atmosphere
   Real pedersencoeff;
   Real hallcoeff;
   Real parallelcoeff;
};
std::array<AtmosphericLayer, numAtmosphereLevels> atmosphere;

std::array< Real, numAtmosphereLevels > productionTable;


// Fractional energy dissipation rate for a isotropic beam, based on Rees (1963), figure 1
static Real ReesIsotropicLambda(Real x) {
   static const Real P[7] = { -11.639, 32.1133, -30.8543, 14.6063, -6.3375, 0.6138, 1.4946};
   Real lambda = (((((P[0] * x + P[1])*x +P[2])*x+ P[3])*x + P[4])*x +P[5])* x+P[6];
   if(x > 1. || lambda < 0) {
      return 0;
   }
   return lambda;
}


int main(int argc, char** argv) {

   if(argc < 3) {
      std::cerr << "Syntax: sigmaProfiles <Density (1/m³)> <Temperature (K)>" << std::endl;
      return 1;
   }
   Real rhon = atof(argv[1]);
   Real T = atof(argv[2]);
   std::string filename = "MSIS.dat";

   // -------------- Build atmosphere model --------------
   // These are the only height values (in km) we are actually interested in
   static const float alt[numAtmosphereLevels] = {
      66, 68, 71, 74, 78, 82, 87, 92, 98, 104, 111,
      118, 126, 134, 143, 152, 162, 172, 183, 194
   };

   // Open file, read in
   ifstream in(filename);
   if(!in) {
      cerr << "(ionosphere) WARNING: Atmospheric Model file " << filename << " could not be opened: " <<
         strerror(errno) << endl
         << "(ionosphere) All atmospheric values will be zero, and there will be no ionization!" << endl;
   }
   int altindex = 0;
   Real integratedDensity = 0;
   Real prevDensity = 0;
   Real prevAltitude = 0;
   while(in) {
      Real altitude, density, c1, c2, c3, extra;
      in >> altitude >>  c1 >> c2 >> c3 >> density >> extra;

      integratedDensity += (altitude - prevAltitude) *1000 * 0.5 * (density + prevDensity);
      prevAltitude = altitude;
      prevDensity = density;

      // When we encounter one of our reference layers, record its values
      if(altitude >= alt[altindex] && altindex < numAtmosphereLevels) {
         AtmosphericLayer newLayer;
         newLayer.altitude = altitude; // in km
         newLayer.density = density; // kg/m^3
         newLayer.depth = integratedDensity; // kg/m^2
         // Ion-neutral scattering frequencies (from Schunck and Nagy, 2009, Table 4.5)
         newLayer.nui = 1e-16*(2*c1 + 3.8*c2 + 5*c3);
         // Elctron-neutral scattering frequencies (Same source, Table 4.6)
         newLayer.nue = 1e-17*(2.33*c1 + 18.2*c2 + 8.9*c3);
         atmosphere[altindex++] = newLayer;
      }
   }

   // Now we have integrated density from the bottom of the atmosphere in the depth field.
   // Flip it around.
   for(int h=0; h<numAtmosphereLevels; h++) {
      atmosphere[h].depth = integratedDensity - atmosphere[h].depth;
   }

   // Calculate Hall and Pedersen conductivity coefficient based on charge carrier density
   const Real Bval = 5e-5; // TODO: Hardcoded B strength here?
   const Real NO_gyroFreq = CHARGE * Bval / (31*MASS_PROTON); // Ion (NO+) gyration frequency
   const Real e_gyroFreq = CHARGE * Bval / (MASS_ELECTRON); // Elctron gyration frequency
   for(int h=0; h<numAtmosphereLevels; h++) {
      // Vlasiator version
      Real sigma_i = CHARGE*CHARGE / ((31. * MASS_PROTON)  * atmosphere[h].nui);
      Real sigma_e = CHARGE*CHARGE / (MASS_ELECTRON  * atmosphere[h].nue);

      atmosphere[h].pedersencoeff = sigma_i * (atmosphere[h].nui * atmosphere[h].nui)/(atmosphere[h].nui*atmosphere[h].nui + NO_gyroFreq*NO_gyroFreq)
         + sigma_e *(atmosphere[h].nue * atmosphere[h].nue)/(atmosphere[h].nue*atmosphere[h].nue + e_gyroFreq*e_gyroFreq);
      atmosphere[h].hallcoeff = -sigma_i * (atmosphere[h].nui * NO_gyroFreq)/(atmosphere[h].nui*atmosphere[h].nui + NO_gyroFreq*NO_gyroFreq)
         + sigma_e *(atmosphere[h].nue * e_gyroFreq)/(atmosphere[h].nue*atmosphere[h].nue + e_gyroFreq*e_gyroFreq);

      atmosphere[h].parallelcoeff = sigma_e;


      // GUMICS version
      const double gyro = (1.6e-19*Bval/(31*1.66e-27));
      const double rho = atmosphere[h].nui/gyro;
      //atmosphere[h].pedersencoeff = (1.6e-19/Bval)*(rho/(1+rho*rho));
      //atmosphere[h].hallcoeff = rho * atmosphere[h].pedersencoeff;
   }


   // Energies of particles that sample the production array
   // are logspace-distributed from 10^-1 to 10^2.3 keV
   for(int e=0; e<productionNumParticleEnergies; e++) {
      particle_energy[e] = pow(10.0, -1.+e*(2.3+1.)/(productionNumParticleEnergies-1));
   }
   particle_energy[productionNumParticleEnergies] = 2*particle_energy[productionNumParticleEnergies-1] - particle_energy[productionNumParticleEnergies-2];


   // Precalculate scattering rates
   const Real eps_ion_keV = 0.035; // Energy required to create one ion
   std::array< std::array< Real, numAtmosphereLevels >, productionNumParticleEnergies > scatteringRate;
   for(int e=0;e<productionNumParticleEnergies; e++) {

      // From Rees, M. H. (1989), q 3.4.4
      const Real electronRange = 4.3e-6 + 5.36e-5 * pow(particle_energy[e], 1.67); // kg m^-2
      Real rho_R=0.;
      // Integrate downwards through the atmosphre to find density at depth=1
      for(int h=numAtmosphereLevels-1; h>=0; h--) {
         if(atmosphere[h].depth / electronRange > 1) {
            rho_R = atmosphere[h].density;
            break;
         }
      }
      if(rho_R == 0.) {
         rho_R = atmosphere[0].density;
      }

      for(int h=0; h<numAtmosphereLevels; h++) {
         const Real lambda = ReesIsotropicLambda(atmosphere[h].depth/electronRange);
         // Rees et al 1963, eq. 1
         //const Real rate = particle_energy[e] / (electronRange / rho_R) / eps_ion_keV *   lambda   *   atmosphere[h].density / integratedDensity; 
         // Rees 1989, eq. 3.3.7 / 3.3.8
         const Real rate = particle_energy[e] * lambda * atmosphere[h].density / electronRange / eps_ion_keV;
         scatteringRate[h][e] = max(0., rate); // m^-1
      }
   }



   // -------------- Build differential flux --------------
   Real tempenergy = kB * T / CHARGE / 1000;
   Real accenergy = productionMinAccEnergy;
   std::cerr << "# Temperature of " << T << " K == Thermal energy of " << tempenergy << " keV" << std::endl;
   Real integralFlux = 0;

   for(int p=0; p<productionNumParticleEnergies; p++) {
      // TODO: Kappa distribution here? Now only going for maxwellian
      Real energyparam = (particle_energy[p]-accenergy)/tempenergy;

      if(particle_energy[p] > accenergy) {
         Real deltaE = (particle_energy[p+1] - particle_energy[p])* 1e3*CHARGE;  // dE in J

         differentialFlux[p] = sqrt(1. / (2. * M_PI * MASS_ELECTRON))
            * particle_energy[p] / tempenergy / sqrt(tempenergy * 1e3 *CHARGE)
            * deltaE * exp(-energyparam);
         integralFlux += rhon * differentialFlux[p];
      } else {
         differentialFlux[p] = 0;
      }
   }

   // -------------- Fill production table --------------
   for(int h=0; h < numAtmosphereLevels; h++) {
      productionTable[h] = 0;
      for(int p=0; p<productionNumParticleEnergies; p++) {
         productionTable[h] += scatteringRate[h][p]*differentialFlux[p];
      }
   }

   // Calculate and output electron density and conductivities
   std::cout << "# Altitude (m)\tn_e (1/m³)\tsigmaP (mho)\tsigmaH (mho)\tsigmaParallel (mho)\tProduction rate (m^-3 s^-1)" << std::endl;
   std::array<Real, numAtmosphereLevels> electronDensity;
   Real SigmaH=0;
   Real SigmaP=0;
   Real SigmaParallel=0;

   for(int h=1; h < numAtmosphereLevels; h++) {
      Real qref = rhon*productionTable[h];

      // Get equilibrium electron density
      electronDensity[h] = sqrt(qref/recombAlpha);

      // Calculate conductivities
      Real halfdx = 1000 * 0.5 * (atmosphere[h].altitude -  atmosphere[h-1].altitude);
      Real halfCH = halfdx * 0.5 * (atmosphere[h-1].hallcoeff + atmosphere[h].hallcoeff);
      Real halfCP = halfdx * 0.5 * (atmosphere[h-1].pedersencoeff + atmosphere[h].pedersencoeff);
      Real halfCpara = halfdx * 0.5 * (atmosphere[h-1].parallelcoeff + atmosphere[h].parallelcoeff);

      Real sigmap = (electronDensity[h]+electronDensity[h-1]) * halfCP;
      Real sigmah = (electronDensity[h]+electronDensity[h-1]) * halfCH;
      Real sigmaParallel = (electronDensity[h]+electronDensity[h-1]) * halfCpara;

      SigmaP += sigmap;
      SigmaH += sigmah;
      SigmaParallel += sigmaParallel;

      std::cout << atmosphere[h].altitude << "\t" << 0.5*(electronDensity[h]+electronDensity[h-1]) << "\t" << sigmap/(2.*halfdx) << "\t" << sigmah/(2.*halfdx) << "\t" << sigmaParallel/(2.*halfdx) << "\t" << qref << std::endl;
   }
   std::cerr << std::endl;

   std::cerr << "Integral energy flux: " << integralFlux << " J/m²/s" << std::endl;
   std::cerr << "Height integrated conductivities:" << std::endl;
   std::cerr << "  SigmaH = " << SigmaH << std::endl;
   std::cerr << "  SigmaP = " << SigmaP << std::endl;
   std::cerr << "  Sigma∥ = " << SigmaParallel << std::endl;

   return 0;
}
