#include <iostream>
#include <array>
#include <cmath>

typedef double Real;
constexpr static int productionNumAccEnergies = 60;
constexpr static int productionNumTemperatures = 60;
constexpr static int productionNumParticleEnergies = 100;
constexpr static Real productionMinAccEnergy = 0.1; // keV
constexpr static Real productionMaxAccEnergy = 100.; // keV
constexpr static Real productionMinTemperature = 0.1; // keV
constexpr static Real productionMaxTemperature = 100.; // keV

static const Real kB = 1.380649e-23;
static const Real CHARGE = 1.602176634e-19;
static const Real MASS_ELECTRON = 9.1093837015e-31;

std::array< Real, productionNumParticleEnergies+1 > particle_energy;
std::array< Real, productionNumParticleEnergies > differentialFlux;


int main(int argc, char** argv) {

	if(argc < 3) {
		std::cerr << "Syntax: differentialFlux <Density (1/m³)> <Temperature (K)>" << std::endl;
		return 1;
	}
	Real rhon = atof(argv[1]);
	Real T = atof(argv[2]);

	// Energies of particles that sample the production array
	// are logspace-distributed from 10^-1 to 10^2.3 keV
	for(int e=0; e<productionNumParticleEnergies; e++) {
		particle_energy[e] = pow(10.0, -1.+e*(2.3+1.)/(productionNumParticleEnergies-1));
	}
	particle_energy[productionNumParticleEnergies] = 2*particle_energy[productionNumParticleEnergies-1] - particle_energy[productionNumParticleEnergies-2];

   Real tempenergy = kB * T / CHARGE / 1000;
	Real accenergy = productionMinAccEnergy;
	std::cerr << "# Temperature of " << T << " K == Thermal energy of " << tempenergy << " keV" << std::endl;

   for(int p=0; p<productionNumParticleEnergies; p++) {
      // TODO: Kappa distribution here? Now only going for maxwellian
      Real energyparam = (particle_energy[p]-accenergy)/tempenergy;

      if(particle_energy[p] > accenergy) {
         Real deltaE = (particle_energy[p+1] - particle_energy[p])* 1e3*CHARGE;  // dE in J

         differentialFlux[p] = sqrt(1. / (2. * M_PI * MASS_ELECTRON))
            * particle_energy[p] / tempenergy / sqrt(tempenergy * 1e3 *CHARGE)
            * deltaE * exp(-energyparam);
      } else {
         differentialFlux[p] = 0;
      }
   }

	std::cout << "#Energy (keV)\tFlux (1/m²/s)" << std::endl;
	for(int p=0; p<productionNumParticleEnergies; p++) {
		std::cout << particle_energy[p] << "\t" << rhon*differentialFlux[p] << std::endl;
	}

	return 0;
}
