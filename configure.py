#!/usr/bin/env python
# ------------------------------------------------------------------------------
# configure.py: Vlasiator configuration script in python.
#
# When configure.py is run, it uses the command line options and default 
# settings to create custom versions of Makefile from the template file 
# Makefile.in.
# Original version by CJW. Modified by Hongyang Zhou.
# ------------------------------------------------------------------------------

# Modules
import argparse
import glob
import re
import subprocess
import os

# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# --- Step 1. Prepare parser, add each of the arguments ------------------
vlasiator_description = (
    "Prepare custom Makefile for compiling Vlasiator"
)
vlasiator_epilog = (
    "Full documentation of options available at "
    "https://github.com/fmihpc/vlasiator/wiki/Installing-Vlasiator"
)
parser = argparse.ArgumentParser(description=vlasiator_description, epilog=vlasiator_epilog)

# -amr
parser.add_argument('-install',
                    action='store_true',
                    default=False,
                    help='install all the dependencies')

# --coord=[name]
parser.add_argument(
    '--coord',
    default='cartesian',
    choices=[
        'cartesian',
        'gr_user'],
    help='select coordinate system')

# --nx=[value]
parser.add_argument('--nx',
                    type=int,
                    default=4,
                    help='set number of cells in the block x dimension')

# --ny=[value]
parser.add_argument('--ny',
                    type=int,
                    default=4,
                    help='set number of cells in the block y dimension')

# --nz=[value]
parser.add_argument('--nz',
                    type=int,
                    default=4,
                    help='set number of cells in the block z dimension')

# --velocityorder=[value]
parser.add_argument('--velocityorder',
                    type=int,
                    default=5,
                    help='set order of velocity space acceleration')

# --spatialorder=[value]
parser.add_argument('--spatialorder',
                    type=int,
                    default=3,
                    help='set order of spatial translation')

# --field=[name]
parser.add_argument('--field',
                    default='londrillo_delzanna',
                    help='select field solver')

# --fieldorder=[value]
parser.add_argument('--fieldorder',
                    type=int,
                    default=2,
                    help='select field solver order')

# -amr
parser.add_argument('-amr',
                    action='store_true',
                    default=False,
                    help='enable AMR in velocity space')

# -debug
parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

# -debugsolver
parser.add_argument('-debugsolver',
                    action='store_true',
                    default=False,
                    help='enable debug flags for field solver')

# -debugionosphere
parser.add_argument('-debugionosphere',
                    action='store_true',
                    default=False,
                    help='enable debug flags for ionosphere module')

# -debugfloat
parser.add_argument('-debugfloat',
                    action='store_true',
                    default=False,
                    help='enable catching floating point exceptions')

# -float
parser.add_argument('-float',
                    action='store_true',
                    default=False,
                    help='enable single precision')

# -distfloat
parser.add_argument('-distfloat',
                    action='store_true',
                    default=False,
                    help='enable single precision for distribution function')

# -profile
parser.add_argument('-profile', 
                    dest='profile', 
                    action='store_true',
                    help='enable profiler')
parser.add_argument('-noprofile',
                    dest='profile',
                    action='store_false',
                    help='disable profiler')
parser.set_defaults(profile=True)

# -mpi
parser.add_argument('-mpi',
                    action='store_true',
                    default=True,
                    help='enable parallelization with MPI')

# -omp
parser.add_argument('-omp',
                    action='store_true',
                    default=False,
                    help='enable parallelization with OpenMP')

# --boost_path=[name]
parser.add_argument('--boost_path',
                    default='/usr/lib/x86_64-linux-gnu',
                    help='path to Boost libraries')

# --vlsv_path=[name]
parser.add_argument('--vlsv_path',
                    default='',
                    help='path to VLSV library')

# --vlsv_path=[name]
parser.add_argument('--dccrg_path',
                    default='',
                    help='path to DCCRG library')

# --phiprof_path=[name]
parser.add_argument('--fsgrid_path',
                    default='',
                    help='path to fsgrid library')

# --zoltan_path=[name]
parser.add_argument('--zoltan_path',
                    default='',
                    help='path to Zoltan library')

# --phiprof_path=[name]
parser.add_argument('--phiprof_path',
                    default='',
                    help='path to phiprof library')

# -jemalloc
parser.add_argument('-jemalloc',
                    action='store_true',
                    default=True,
                    help='enable Jemalloc memory allocator')

# --jemalloc_path=[name]
parser.add_argument('--jemalloc_path',
                    default='',
                    help='path to jemalloc library')

# -papi
parser.add_argument('-papi',
                    action='store_true',
                    default=False,
                    help='enable Papi memory profiler')

# The main choices for --cxx flag, using "ctype[-suffix]" formatting, where 
# "ctype" is the major family/suite/group of compilers and "suffix" may 
# represent variants of the compiler version and/or predefined sets of compiler
# options. The C++ compiler front ends are the main supported/documented options
# and are invoked on the command line, but the C front ends are also acceptable 
# selections and are mapped to the matching C++ front end:
# gcc -> g++, clang -> clang++, icc-> icpc
cxx_choices = [
    'g++',
    'g++-simd',
    'icpc',
    'icpc-debug',
    'icpc-phi',
    'cray',
    'clang++',
    'clang++-simd',
    'clang++-apple',
]


def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)

    if arg == 'clang':
        arg = 'clang++'
    else:
        arg = arg.replace('clang-', 'clang++-', 1)
    return arg


# --cxx=[name]
parser.add_argument(
    '--cxx',
    default='g++',
    type=c_to_cpp,
    choices=cxx_choices,
    help='select C++ compiler and default set of flags')

# --cflag=[string]
parser.add_argument('--cflag',
                    default=None,
                    help='additional string of flags to append to compiler/linker calls')

# --include=[name]
parser.add_argument(
    '--include',
    default=[],
    action='append',
    help=('extra path for included header files (-I<path>); can be specified multiple '
          'times'))

# --lib_path=[name]
parser.add_argument(
    '--lib_path',
    default=[],
    action='append',
    help=('extra path for linked library files (-L<path>); can be specified multiple '
          'times'))

# --lib=[name]
parser.add_argument(
    '--lib',
    default=[],
    action='append',
    help='name of library to link against (-l<lib>); can be specified multiple times')

# Parse command-line inputs
args = vars(parser.parse_args())

# --- Step 2. Test for incompatible arguments ----------------------------


# --- Step 3. Set definitions and Makefile options based on above argument

# Prepare dictionaries of substitutions to be made
definitions = {}
makefile_options = {}
makefile_options['LOADER_FLAGS'] = ''

# --cxx=[name] argument
if args['cxx'] == 'g++':
    # GCC is C++11 feature-complete since v4.8.1 (2013-05-31)
    definitions['COMPILER_CHOICE'] = 'g++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'g++-simd':
    # GCC version >= 4.9, for OpenMP 4.0; version >= 6.1 for OpenMP 4.5 support
    definitions['COMPILER_CHOICE'] = 'g++-simd'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -std=c++11 -fopenmp-simd -fwhole-program -flto -ffast-math '
        '-march=native -fprefetch-loop-arrays'
        # -march=skylake-avx512, skylake, core-avx2
        # -mprefer-vector-width=128  # available in gcc-8, but not gcc-7
        # -mtune=native, generic, broadwell
        # -mprefer-avx128
        # -m64 (default)
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc':
    # ICC is C++11 feature-complete since v15.0 (2014-08-26)
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -inline-forceinline -qopenmp-simd -qopt-prefetch=4 '
      '-qoverride-limits'  # -qopt-report-phase=ipo (does nothing without -ipo)
    )
    # -qopt-zmm-usage=high'  # typically harms multi-core performance on Skylake Xeon
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc-debug':
    # Disable IPO, forced inlining, and fast math. Enable vectorization reporting.
    # Useful for testing symmetry, SIMD-enabled functions and loops with OpenMP 4.5
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -xhost -qopenmp-simd -fp-model precise -qopt-prefetch=4 '
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits'
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc-phi':
    # Cross-compile for Intel Xeon Phi x200 KNL series (unique AVX-512ER and AVX-512FP)
    # -xMIC-AVX512: generate AVX-512F, AVX-512CD, AVX-512ER and AVX-512FP
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xMIC-AVX512 -inline-forceinline -qopenmp-simd '
      '-qopt-prefetch=4 -qoverride-limits'
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'cray':
    # Cray Compiling Environment 8.4 (2015-09-24) introduces C++11 feature completeness
    # (except "alignas"). v8.6 is C++14 feature-complete
    definitions['COMPILER_CHOICE'] = 'cray'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -h std=c++11 -h aggress -h vector3 -hfp3'
    makefile_options['LINKER_FLAGS'] = '-hwp -hpl=obj/lib'
    makefile_options['LIBRARY_FLAGS'] = '-lm'
if args['cxx'] == 'clang++':
    # Clang is C++11 feature-complete since v3.3 (2013-06-17)
    definitions['COMPILER_CHOICE'] = 'clang++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++-simd':
    # LLVM/Clang version >= 3.9 for most of OpenMP 4.0 and 4.5 (still incomplete; no
    # offloading, target/declare simd directives). OpenMP 3.1 fully supported in LLVM 3.7
    definitions['COMPILER_CHOICE'] = 'clang++-simd'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -fopenmp-simd'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++-apple':
    # Apple LLVM/Clang: forked version of the open-source LLVM project bundled in macOS
    definitions['COMPILER_CHOICE'] = 'clang++-apple'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''

# -float argument
if args['float']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DSP'
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDP'

# Distribution function precision
if args['distfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DSPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4F_FALLBACK' # vector backend type
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4D_FALLBACK' # vector backend type

if args['profile']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPROFILE'
 
if args['velocityorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PLM'       
elif args['velocityorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PPM'
elif args['velocityorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DACC_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for velocity space acceleration')

if args['spatialorder'] == 2:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PLM'       
elif args['spatialorder'] == 3:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PPM'
elif args['spatialorder'] == 5:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DTRANS_SEMILAG_PQM' 
else:
    raise SystemExit('### CONFIGURE ERROR: unknown semilag solver order for spatial translation')

# Make the field solver first-order in space and time
if args['fieldorder'] == 1:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_SPACE'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DFS_1ST_ORDER_TIME'

# Turn on AMR
if args['amr']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DAMR'

# Add -DNDEBUG to turn debugging off.
if args['debug']:
    pass
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DNDEBUG'

# Debug for field solver
if args['debugsolver']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_SOLVERS'

# Debug for ionosphere module
if args['debugionosphere']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDEBUG_IONOSPHERE'

# Catch floating point exceptions and stop execution
if args['debugfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DCATCH_FPE'

# -debug argument
if args['debug']:
    definitions['DEBUG_OPTION'] = 'DEBUG'
    # Completely replace the --cxx= sets of default compiler flags, disable optimization,
    # and emit debug symbols in the compiled binaries
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g'  # -Og
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] = '-O0 -h std=c++11'
    if args['cxx'] == 'bgxlc++':
        makefile_options['COMPILER_FLAGS'] = '-O0 -g -qlanglvl=extended0x'
    if args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g -xMIC-AVX512'
else:
    definitions['DEBUG_OPTION'] = 'NOT_DEBUG'

# -mpi argument
if args['mpi']:
    definitions['MPI_OPTION'] = 'MPI_PARALLEL'
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpc-phi' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -h mpi1'
else:
    definitions['MPI_OPTION'] = 'NOT_MPI_PARALLEL'

# -omp argument
if args['omp']:
    definitions['OPENMP_OPTION'] = 'OPENMP_PARALLEL'
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd' or args['cxx'] == 'clang++'
            or args['cxx'] == 'clang++-simd'):
        makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    if (args['cxx'] == 'clang++-apple'):
        # Apple Clang disables the front end OpenMP driver interface; enable it via the
        # preprocessor. Must install LLVM's OpenMP runtime library libomp beforehand
        makefile_options['COMPILER_FLAGS'] += ' -Xpreprocessor -fopenmp'
        makefile_options['LIBRARY_FLAGS'] += ' -lomp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] += ' -qopenmp'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -homp'
else:
    definitions['OPENMP_OPTION'] = 'NOT_OPENMP_PARALLEL'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -hnoomp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi':
        # suppressed messages:
        #   3180: pragma omp not recognized
        makefile_options['COMPILER_FLAGS'] += ' -diag-disable 3180'

# --cflag=[string] argument
if args['cflag'] is not None:
    makefile_options['COMPILER_FLAGS'] += ' '+args['cflag']

# --include=[name] arguments
for include_path in args['include']:
    makefile_options['COMPILER_FLAGS'] += ' -I'+include_path

# --lib_path=[name] arguments
for library_path in args['lib_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+library_path

# --lib=[name] arguments
for library_name in args['lib']:
    makefile_options['LIBRARY_FLAGS'] += ' -l'+library_name

# Assemble all flags of any sort given to compiler
definitions['COMPILER_FLAGS'] = ' '.join(
    [makefile_options[opt+'_FLAGS'] for opt in
     ['PREPROCESSOR', 'COMPILER', 'LINKER', 'LIBRARY']])

makefile_options['DCCRG_PATH'] = args['dccrg_path']
makefile_options['FSGRID_PATH'] = args['fsgrid_path']

# Add include paths to header-only libraries
if args['dccrg_path']:
    makefile_options['COMPILER_FLAGS'] += ' -I'+args['dccrg_path']
if args['fsgrid_path']:
    makefile_options['COMPILER_FLAGS'] += ' -I'+args['fsgrid_path']

makefile_options['PHIPROF_PATH'] = args['phiprof_path']
makefile_options['BOOST_PATH'] = args['boost_path']
makefile_options['ZOLTAN_PATH'] = args['zoltan_path']

# Add lib paths
if args['phiprof_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+args['phiprof_path']
if args['zoltan_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+args['zoltan_path']
if args['vlsv_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+args['vlsv_path']
if args['jemalloc_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+args['jemalloc_path']
if args['boost_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+args['boost_path']


# --- Step 4. Check dependencies -----------------------------------------

# Install dependencies
if args['install']:
    # Boost is skipped as it is too large to install here
    if not os.path.isdir("lib"):
        subprocess.check_call(["mkdir", "lib"]) 
    if not os.path.isdir("lib/vectorclass"):
        subprocess.check_call(["git", "clone", "https://github.com/vectorclass/version1.git"])
        subprocess.check_call(["git", "clone", "https://github.com/vectorclass/add-on.git"])
        subprocess.check_call(["cp", "add-on/vector3d/vector3d.h", "version1/"])
        subprocess.check_call(["mv", "version1", "lib/vectorclass"])
    if not os.path.isdir("lib/fsgrid"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/fsgrid.git"])
        subprocess.check_call(["mv", "fsgrid", "lib"])
    if not os.path.isdir("lib/dccrg"): 
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/dccrg.git"])
        subprocess.check_call(["mv", "dccrg", "lib"])
        os.chdir("lib/dccrg")
        subprocess.check_call(["git", "checkout", "01482cfba8"])
        os.chdir("../..")
    if not os.path.isdir("lib/phiprof"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/phiprof.git"])
        subprocess.check_call(["mv", "phiprof", "lib"])
        os.chdir("lib/phiprof/src")
        subprocess.check_call(["make"])
        os.chdir("../../..")
    if not os.path.isdir("lib/vlsv"):
        subprocess.check_call(["git", "clone", "https://github.com/fmihpc/vlsv.git"])
        subprocess.check_call(["mv", "vlsv", "lib"])
        os.chdir("lib/vlsv")
        subprocess.check_call(["make"])
        os.chdir("../..")
    if not os.path.isdir("lib/jemalloc"):
        subprocess.check_call(["wget", \
            "https://github.com/jemalloc/jemalloc/releases/download/4.0.4/jemalloc-4.0.4.tar.bz2"])
        subprocess.check_call(["tar", "-xf", "jemalloc-4.0.4.tar.bz2"])
        os.chdir("jemalloc-4.0.4")
        subprocess.check_call(["./configure","--prefix="+os.getcwd()+"/jemalloc", \
        "--with-jemalloc-prefix=je_"])
        subprocess.check_call(["make"])
        subprocess.check_call(["make", "install"])
        subprocess.check_call(["mv", "jemalloc", "../lib"])
        os.chdir("..")
    if not os.path.isdir("lib/Eigen"):
        subprocess.check_call(["wget", \
        "https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.bz2"])
        subprocess.check_call(["tar", "-xf", "eigen-3.2.8.tar.bz2"])
        subprocess.check_call(["cp", "-r", "eigen-3.2.8/Eigen", "lib"])
        subprocess.check_call(["wget", \
        "http://cs.sandia.gov/Zoltan/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz"])
        subprocess.check_call(["tar", "-xf", "zoltan_distrib_v3.83.tar.gz"])
        subprocess.check_call(["mkdir", "lib/zoltan"])
        os.chdir("lib/zoltan")
        subprocess.check_call(["../../Zoltan_v3.83/configure","--prefix="+os.getcwd(), \
        "--enable-mpi", "--with-mpi-compilers", "--with-gnumake", \
        "--with-id-type=ullong"])
        subprocess.check_call(["make", "-j", "4"])
        subprocess.check_call(["make", "install"])

        os.chdir("../..")

    for f in ["add-on", "eigen-3.2.8.tar.bz2","eigen-3.2.8",\
        "zoltan_distrib_v3.83.tar.gz", "Zoltan_v3.83", \
        "jemalloc-4.0.4", "jemalloc-4.0.4.tar.bz2"]:
        subprocess.check_call(["rm", "-rf", f])

    # Create path names for generating version.cpp in Makefile
    makefile_options['DCCRG_PATH'] = "lib/dccrg"
    makefile_options['FSGRID_PATH'] = "lib/fsgrid"
    makefile_options['PHIPROF_PATH'] = "lib/phiprof"
    makefile_options['ZOLTAN_PATH'] = "lib/zoltan"

    makefile_options['LINKER_FLAGS'] += " -Llib/zoltan/lib"
    makefile_options['LINKER_FLAGS'] += " -Llib/jemalloc/lib"
    makefile_options['LINKER_FLAGS'] += " -Llib/phiprof/lib"
    makefile_options['LINKER_FLAGS'] += " -Llib/vlsv"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/vectorclass"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/zoltan/include"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/fsgrid"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/dccrg"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib" # Eigen
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/phiprof/include"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/vlsv"
    makefile_options['COMPILER_FLAGS'] += ' -I'+"lib/jemalloc/include/jemalloc"

    for f in ["boost_program_options", "zoltan", "vlsv", "jemalloc", "phiprof"]:
        makefile_options['LIBRARY_FLAGS'] += ' -l'+f
else:
    # Check dependencies
    if any(os.path.isfile(os.path.join(p, "vectorf512.h")) for p in args['include']):
        if not any(os.path.isfile(os.path.join(p, "vector3d.h")) for p in args['include']):
            raise SystemExit('### CONFIGURE ERROR: vector3d.h not found!')
    else:
        raise SystemExit('### CONFIGURE ERROR: unknown vectorclass location!')

    if not os.path.isfile(os.path.join(args['fsgrid_path'], "fsgrid.hpp")):
        raise SystemExit('### CONFIGURE ERROR: unknown fsgrid location!')

    if not os.path.isfile(os.path.join(args['dccrg_path'], "dccrg.hpp")):
        raise SystemExit('### CONFIGURE ERROR: unknown dccrg location!')

    if not os.path.isfile(os.path.join(args['boost_path'], "libboost_program_options.a")):
        raise SystemExit('### CONFIGURE ERROR: unknown Boost location!')

    if not os.path.isfile(os.path.join(args['zoltan_path'], "libzoltan.a")):
        raise SystemExit('### CONFIGURE ERROR: unknown Zoltan location!')

    if not os.path.isfile(os.path.join(args['vlsv_path'], "conv_mtx_vlsv")):
        raise SystemExit('### CONFIGURE ERROR: vlsv not found or compiled!')

    if args['jemalloc']:
        if not os.path.isfile(os.path.join(args['jemalloc_path'], "libjemalloc.a")):
            raise SystemExit('### CONFIGURE ERROR: unknown jemalloc location!')

    if args['profile']:
        if not os.path.isfile(os.path.join(args['phiprof_path'], "libphiprof.a")):
            raise SystemExit('### CONFIGURE ERROR: unknown phiprof location!')


# --- Step 5. Create new files, finish up --------------------------------

# Read templates
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()

# Make substitutions
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

# Redirect field solver folder
if args['amr']:
    makefile_template = re.sub('vlasovsolver', 'vlasovsolver_amr', makefile_template)

# Write output files
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

# Finish with diagnostic output

print('Your Vlasiator distribution has now been configured with the following options:')
print('  Coordinate system:          ' + args['coord'])
print('  Linker flags:               ' + makefile_options['LINKER_FLAGS'] + ' '
      + makefile_options['LIBRARY_FLAGS'])
print('  Floating-point precision:   ' + ('single' if args['float'] else 'double'))
print('  Distribution precision:     ' + ('single' if args['distfloat'] else 'double'))
print('  Block size:                 ' + str(args['nx']) + ' ' \
                                       + str(args['ny']) + ' ' \
                                       + str(args['nz']))
print('  MPI parallelism:            ' + ('ON' if args['mpi'] else 'OFF'))
print('  OpenMP parallelism:         ' + ('ON' if args['omp'] else 'OFF'))
print('  Order of field solver:      ' + str(args['fieldorder']))
print('  Order of semilog velocity:  ' + str(args['velocityorder']))
print('  Order of semilog spatial:   ' + str(args['spatialorder']))
print('  AMR:                        ' + ('ON' if args['amr'] else 'OFF'))
print('  Profiler:                   ' + ('ON' if args['profile'] else 'OFF'))
print('  Memory tracker:             ' + ('ON' if args['papi'] else 'OFF'))
print('  Debug flags:                ' + ('ON' if args['debug'] else 'OFF'))
print('  Compiler:                   ' + args['cxx'])
print('  Compilation command:        ' + makefile_options['COMPILER_COMMAND'] + ' '
      + makefile_options['PREPROCESSOR_FLAGS'] + ' ' + makefile_options['COMPILER_FLAGS'])
