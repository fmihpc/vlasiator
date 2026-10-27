#!/bin/bash
# GNU + OpenMPI + netCDF stack for GITM/TIEGCM
set -e   # stop on first error 

export EASYBUILD_PREFIX=/turso/group/spacephysics/applications/easybuild
module load PrgEnv-gnu/15.2.0
module unload MPICH gmpich 2>/dev/null || true
module use $EASYBUILD_PREFIX/modules/all
module use $EASYBUILD_PREFIX/modules/all/Compiler/GCCcore/15.2.0
module use $EASYBUILD_PREFIX/modules/all/Compiler/GCC/15.2.0
module use $EASYBUILD_PREFIX/modules/all/MPI/GCC/15.2.0/OpenMPI/5.0.10
module load OpenMPI/5.0.10 netCDF-Fortran

# Clone into a fixed name
git clone https://github.com/GITMCode/GITM.git GITM_testbuild
cd GITM_testbuild

./Config.pl -install -earth -compiler=gfortran
make

# Verification: the build is only "good" if these are true
test -f ext/Electrodynamics/src/Makefile || { echo "FAIL: Electrodynamics extension missing"; exit 1; }
find . -name "GITM.exe" | grep -q . || { echo "FAIL: GITM.exe not built"; exit 1; }
echo "BUILD OK"
(venv) [jpettit@turso03 ~]$ cat GITM_build_script_v2.sh 
#!/bin/bash
# GNU + OpenMPI + netCDF stack for GITM/TIEGCM
set -e   # stop on first error — "script finished" then means "build succeeded"

export EASYBUILD_PREFIX=/turso/group/spacephysics/applications/easybuild
module load PrgEnv-gnu/15.2.0
module unload MPICH gmpich 2>/dev/null || true
module use $EASYBUILD_PREFIX/modules/all
module use $EASYBUILD_PREFIX/modules/all/Compiler/GCCcore/15.2.0
module use $EASYBUILD_PREFIX/modules/all/Compiler/GCC/15.2.0
module use $EASYBUILD_PREFIX/modules/all/MPI/GCC/15.2.0/OpenMPI/5.0.10
module load OpenMPI/5.0.10 netCDF-Fortran

# Clone into a fixed name and USE THE SAME NAME afterwards
git clone https://github.com/GITMCode/GITM.git GITM_testbuild
cd GITM_testbuild

./Config.pl -install -earth -compiler=gfortran
make

# Verification: the build is only "good" if these are true
test -f ext/Electrodynamics/src/Makefile || { echo "FAIL: Electrodynamics extension missing"; exit 1; }
find . -name "GITM.exe" | grep -q . || { echo "FAIL: GITM.exe not built"; exit 1; }
echo "BUILD OK"
