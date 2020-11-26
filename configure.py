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
import pkg_resources
import sys
import warnings
import fileinput


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

# -install
parser.add_argument('-install',
                    action='store_true',
                    default=False,
                    help='install all the dependencies')

# -save
parser.add_argument('-save',
                    action='store_true',
                    default=False,
                    help='save library paths into a customized Makefile')

# --save_machine=[name]
parser.add_argument('--save_machine',
                    default='new',
                    help='customized Makefile name')

# --coord=[name]
parser.add_argument(
    '--coord',
    default='cartesian',
    choices=[
        'cartesian',
        'user'],
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


# --machine=[name]
parser.add_argument('--machine',
                    default='',
                    help='choose machine name for specific setup')

# --boost_path=[string]
parser.add_argument('--boost_path',
                    default='/usr/lib/x86_64-linux-gnu',
                    help='path to Boost libraries')

# --vlsv_path=[string]
parser.add_argument('--vlsv_path',
                    default='',
                    help='path to VLSV library')

# --dccrg_path=[string]
parser.add_argument('--dccrg_path',
                    default='',
                    help='path to DCCRG library')

# --fsgrid_path=[string]
parser.add_argument('--fsgrid_path',
                    default='',
                    help='path to fsgrid library')

# --vectorclass_path=[string]
parser.add_argument('--vectorclass_path',
                    default='',
                    help='path to vectorclass library')

# --zoltan_path=[string]
parser.add_argument('--zoltan_path',
                    default='',
                    help='path to Zoltan library')

# --profile_path=[string]
parser.add_argument('--profile_path',
                    default='',
                    help='path to profile library')

# -jemalloc
parser.add_argument('-jemalloc',
                    action='store_true',
                    default=True,
                    help='enable Jemalloc memory allocator')

# --jemalloc_path=[string]
parser.add_argument('--jemalloc_path',
                    default='',
                    help='path to jemalloc library')

# -papi (WIP)
parser.add_argument('-papi',
                    action='store_true',
                    default=False,
                    help='enable Papi memory profiler')

# -silo
parser.add_argument('-silo',
                    action='store_true',
                    default=False,
                    help='enable silo format converter')

# --silo_path=[string]
parser.add_argument('--silo_path',
                    default='',
                    help='path to silo library')

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

if args['install'] and args['machine']:
    raise SystemExit('### CONFIGURE ERROR: does not support fresh installation with a preset machine makefile')

if args['save'] and args['save_machine'] == args['machine']:
    warnings.warn("Overwrite existings MAKE/Makefile."+args['machine']+"...")


# --- Step 3. Set Makefile options based on above argument

# Prepare dictionaries of substitutions to be made
makefile_options = {}
makefile_options['LOADER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''
makefile_options['LINKER_FLAGS'] = ''
makefile_options['LIBRARY_FLAGS'] = ''

# --cxx=[name] argument
if args['cxx'] == 'g++':
    # GCC is C++11 feature-complete since v4.8.1 (2013-05-31)
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
if args['cxx'] == 'g++-simd':
    # GCC version >= 4.9, for OpenMP 4.0; version >= 6.1 for OpenMP 4.5 support
    makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -std=c++11 -fopenmp-simd -fwhole-program -flto -ffast-math '
        '-march=native -fprefetch-loop-arrays'
        # -march=skylake-avx512, skylake, core-avx2
        # -mprefer-vector-width=128  # available in gcc-8, but not gcc-7
        # -mtune=native, generic, broadwell
        # -mprefer-avx128
        # -m64 (default)
    )
if args['cxx'] == 'icpc':
    # ICC is C++11 feature-complete since v15.0 (2014-08-26)
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -inline-forceinline -qopenmp-simd -qopt-prefetch=4 '
      '-qoverride-limits'  # -qopt-report-phase=ipo (does nothing without -ipo)
    )
    # -qopt-zmm-usage=high'  # typically harms multi-core performance on Skylake Xeon
if args['cxx'] == 'icpc-debug':
    # Disable IPO, forced inlining, and fast math. Enable vectorization reporting.
    # Useful for testing symmetry, SIMD-enabled functions and loops with OpenMP 4.5
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -xhost -qopenmp-simd -fp-model precise -qopt-prefetch=4 '
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits'
    )
if args['cxx'] == 'icpc-phi':
    # Cross-compile for Intel Xeon Phi x200 KNL series (unique AVX-512ER and AVX-512FP)
    # -xMIC-AVX512: generate AVX-512F, AVX-512CD, AVX-512ER and AVX-512FP
    makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xMIC-AVX512 -inline-forceinline -qopenmp-simd '
      '-qopt-prefetch=4 -qoverride-limits'
    )
if args['cxx'] == 'cray':
    # Cray Compiling Environment 8.4 (2015-09-24) introduces C++11 feature completeness
    # (except "alignas"). v8.6 is C++14 feature-complete
    makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['COMPILER_FLAGS'] = '-O3 -h std=c++11 -h aggress -h vector3 -hfp3'
    makefile_options['LINKER_FLAGS'] = '-hwp -hpl=obj/lib'
    makefile_options['LIBRARY_FLAGS'] = '-lm'
if args['cxx'] == 'clang++':
    # Clang is C++11 feature-complete since v3.3 (2013-06-17)
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
if args['cxx'] == 'clang++-simd':
    # LLVM/Clang version >= 3.9 for most of OpenMP 4.0 and 4.5 (still incomplete; no
    # offloading, target/declare simd directives). OpenMP 3.1 fully supported in LLVM 3.7
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -fopenmp-simd'
if args['cxx'] == 'clang++-apple':
    # Apple LLVM/Clang: forked version of the open-source LLVM project bundled in macOS
    makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'

# -float argument
if args['float']:
    precision = 'SP'
else:
    precision = 'DP'

makefile_options['PREPROCESSOR_FLAGS'] += ' -D'+precision

for key in ('vlsvdiff','vlsvextract','vlsv2silo'):
    makefile_options[key.upper()] = key+'_'+precision

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

# -mpi argument
if args['mpi']:
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpc-phi' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -h mpi1'
else:
    raise SystemExit('### CONFIGURE ERROR: -mpi is required for compilation!')

# -omp argument
if args['omp']:
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
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -hnoomp'
    if args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi':
        # suppressed messages:
        #   3180: pragma omp not recognized
        makefile_options['COMPILER_FLAGS'] += ' -diag-disable 3180'

# --cflag=[string]
if args['cflag'] is not None:
    makefile_options['COMPILER_FLAGS'] += ' '+args['cflag']

# Install dependencies
if args['install']:
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
    if not os.path.isdir("lib/zoltan"):
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
    makefile_options['BOOST_PATH'] = args['boost_path']
    makefile_options['INC_DCCRG'] = "lib/dccrg"
    makefile_options['INC_FSGRID'] = "lib/fsgrid"
    makefile_options['PROFILE_PATH'] = "lib/phiprof/lib"
    makefile_options['ZOLTAN_PATH'] = "lib/zoltan/lib"
    makefile_options['INC_VECTORCLASS'] = "lib/vectorclass"
    makefile_options['JEMALLOC_PATH'] = "lib/jemalloc/lib"
    makefile_options['VLSV_PATH'] = "lib/vlsv"
    makefile_options['SILO_PATH'] = ""

    # Boost is skipped as it is too large to install here
    if not os.path.isfile(os.path.join(makefile_options['BOOST_PATH'], "libboost_program_options.a")):
        errmsghead = """Boost not found: try to set the correct path, or """
        if 'Ubuntu' in os.uname()[3]:
            raise SystemExit(errmsghead+"install it manually as follows:"
            '\n sudo apt update\n sudo apt install libboost-all-dev')
        else:
            raise SystemExit(errmsghead+
            'search for how to install Boost!')

    makefile_options['INC_BOOST']       = ""
    makefile_options['INC_ZOLTAN']      = "lib/zoltan/include"
    makefile_options['INC_VLSV']        = "lib/vlsv"
    makefile_options['INC_SILO']        = ""
    makefile_options['INC_JEMALLOC']    = "lib/jemalloc/include/jemalloc"
    makefile_options['INC_PROFILE']     = "lib/phiprof/include"
    makefile_options['INC_EIGEN']       = "lib"
    makefile_options['INC_VECTORCLASS'] = "lib/vectorclass"

    for f in ["boost_program_options", "zoltan", "vlsv", "jemalloc", "phiprof"]:
        makefile_options['LIBRARY_FLAGS'] += ' -l'+f
else:
    # --include=[name]
    for include_path in args['include']:
        makefile_options['COMPILER_FLAGS'] += ' -I'+include_path

    # --lib_path=[name]
    for library_path in args['lib_path']:
        makefile_options['LINKER_FLAGS'] += ' -L'+library_path

    # --lib=[name]
    for library_name in args['lib']:
        makefile_options['LIBRARY_FLAGS'] += ' -l'+library_name

    makefile_options['PROFILE_PATH']    = args['profile_path']
    makefile_options['BOOST_PATH']      = args['boost_path']
    makefile_options['ZOLTAN_PATH']     = args['zoltan_path']
    makefile_options['JEMALLOC_PATH']   = args['jemalloc_path']
    makefile_options['VLSV_PATH']       = args['vlsv_path']
    makefile_options['SILO_PATH']       = args['silo_path']

    makefile_options['INC_DCCRG']       = args['dccrg_path']
    makefile_options['INC_FSGRID']      = args['fsgrid_path']
    makefile_options['INC_VECTORCLASS'] = args['vectorclass_path']
    makefile_options['INC_BOOST']       = ""
    makefile_options['INC_ZOLTAN']      = ""
    makefile_options['INC_VLSV']        = ""
    makefile_options['INC_SILO']        = ""
    makefile_options['INC_JEMALLOC']    = ""
    makefile_options['INC_PROFILE']     = ""
    makefile_options['INC_EIGEN']       = ""

if args['machine']:
    with open("MAKE/Makefile."+args['machine']) as f:
        for line in f:
            line = line.strip()
            if line.startswith("CXXFLAGS"):
                makefile_options['COMPILER_FLAGS'] = line.split(' = ')[1]
            elif line.startswith("LDFLAGS"):
                makefile_options['LINKER_FLAGS'] = line.split(' ')[1]
            elif line.startswith("LIBFLAGS"):
                makefile_options['LIBRARY_FLAGS'] = line.split(' ')[1]
            if line.startswith("INC_DCCRG"):
                makefile_options['INC_DCCRG'] = line.split(' = -I')[1]
            elif line.startswith("INC_FSGRID"):
                makefile_options['INC_FSGRID'] = line.split(' = -I')[1]
            elif line.startswith("LIB_PROFILE"):
                lib_path = line.split(' = -L')[1].split(' ')
                makefile_options['LIBRARY_FLAGS'] += ' '+lib_path[1]
                makefile_options['PROFILE_PATH'] = lib_path[0]
            elif line.startswith("LIB_ZOLTAN"):
                lib_path = line.split(' = -L')[1].split(' ')
                makefile_options['LIBRARY_FLAGS'] += ' '+lib_path[1]
                makefile_options['ZOLTAN_PATH'] = lib_path[0]
            elif line.startswith("LIB_VLSV"):
                lib_path = line.split(' = -L')[1].split(' ')
                makefile_options['LIBRARY_FLAGS'] += ' '+lib_path[1]
                makefile_options['VLSV_PATH'] = lib_path[0]
            elif line.startswith("LIB_JEMALLOC"):
                lib_path = line.split(' = -L')[1].split(' ')
                makefile_options['LIBRARY_FLAGS'] += ' '+lib_path[1]
                makefile_options['JEMALLOC_PATH'] = lib_path[0]
            elif line.startswith("INC_EIGEN"):
                makefile_options['INC_EIGEN'] = line.split(' = -I')[1]  
            elif line.startswith("INC_VECTORCLASS"):
                makefile_options['INC_VECTORCLASS'] = line.split(' = -I')[1]
            elif line.startswith("LIB_BOOST"):
                lib_path = line.split(' = -L')[1].split(' ')
                makefile_options['LIBRARY_FLAGS'] += ' '+lib_path[1]
                makefile_options['BOOST_PATH'] = lib_path[0]
            elif line.startswith("INC_BOOST"):
                makefile_options['INC_BOOST'] = line.split(' = -I')[1]
            elif line.startswith("INC_JEMALLOC"):
                makefile_options['INC_JEMALLOC'] = line.split(' = -I')[1]
            elif line.startswith("INC_VLSV"):
                makefile_options['INC_VLSV'] = line.split(' = -I')[1]
            elif line.startswith("INC_ZOLTAN"):
                makefile_options['INC_ZOLTAN'] = line.split(' = -I')[1]
            elif line.startswith("INC_PROFILE"):
                makefile_options['INC_PROFILE'] = line.split(' = -I')[1]
            elif line.startswith("LIB_SILO"):
                makefile_options['SILO_PATH'] = line.split(' = -L')[1]
            elif line.startswith("INC_SILO"):
                makefile_options['INC_SILO'] = line.split(' = -I')[1]


# --- Step 4. Check dependencies -----------------------------------------

# Check dependencies
if os.path.isfile(os.path.join(makefile_options['INC_VECTORCLASS'], "vectorf512.h")):
    if not os.path.isfile(os.path.join(makefile_options['INC_VECTORCLASS'], "vector3d.h")):
        raise SystemExit('### CONFIGURE ERROR: vector3d.h not found!')
else:
    raise SystemExit('### CONFIGURE ERROR: unknown vectorclass location!')

if not os.path.isfile(os.path.join(makefile_options['INC_FSGRID'], "fsgrid.hpp")):
    raise SystemExit('### CONFIGURE ERROR: unknown fsgrid location!')

if not os.path.isfile(os.path.join(makefile_options['INC_DCCRG'], "dccrg.hpp")):
    raise SystemExit('### CONFIGURE ERROR: unknown dccrg location!')

if not os.path.isfile(os.path.join(makefile_options['ZOLTAN_PATH'], "libzoltan.a")):
    raise SystemExit('### CONFIGURE ERROR: unknown Zoltan location!')

if not os.path.isfile(os.path.join(makefile_options['VLSV_PATH'], "conv_mtx_vlsv")):
    raise SystemExit('### CONFIGURE ERROR: vlsv not found or compiled!')

if not os.path.isfile(os.path.join(makefile_options['BOOST_PATH'], "libboost_program_options.a")):
    raise SystemExit('### CONFIGURE ERROR: unknown Boost location!')

if args['jemalloc']:
    if not os.path.isfile(os.path.join(makefile_options['JEMALLOC_PATH'], "libjemalloc.a")):
        raise SystemExit('### CONFIGURE ERROR: unknown jemalloc location!')

if args['profile']:
    if not os.path.isfile(os.path.join(makefile_options['PROFILE_PATH'], "libphiprof.a")):
        raise SystemExit('### CONFIGURE ERROR: unknown phiprof location!')


# --- Step 5. Create new files, finish up --------------------------------

if args['save']:
    makefile_name = 'MAKE/Makefile.'+args['save_machine']
    print('Saving the library paths and compiler flags in '+ makefile_name+'...')
    if os.path.isfile(makefile_name):
        for line in fileinput.input(files=makefile_name, inplace=True):
            isfind = False
            for key, val in makefile_options.items():
                if line.startswith(key):
                    isfind = True
                    if key.startswith('INC'):
                        print(key+' = -I'+val)
                    elif 'PATH' in key:
                        print(key+' = -L'+val)
            if line.startswith('CXXFLAGS'):
                print('CXXFLAGS = '+makefile_options['COMPILER_FLAGS'])
            elif line.startswith('LDFLAGS'):
                print('LDFLAGS = '+makefile_options['LINKER_FLAGS'])
            elif line.startswith('LDLIBS'):
                print('LDLIBS = '+makefile_options['LIBRARY_FLAGS'])
            elif not isfind:
                print(line, end='')

    else:
        with open(makefile_name, 'w') as f:
            f.write('CXXFLAGS = '+makefile_options['COMPILER_FLAGS']+'\n')
            f.write('LDFLAGS = '+makefile_options['LINKER_FLAGS']+'\n')
            f.write('LDLIBS = '+makefile_options['LIBRARY_FLAGS']+'\n\n')

            f.write('#======== Libraries ===========\n\n')
            f.write('INC_BOOST = -I'+makefile_options['INC_BOOST']+'\n')
            f.write('BOOST_PATH = -L'+makefile_options['BOOST_PATH']+'\n\n')

            f.write('INC_ZOLTAN = -I'+makefile_options['INC_ZOLTAN']+'\n')
            f.write('ZOLTAN_PATH = -L'+makefile_options['ZOLTAN_PATH']+'\n\n')

            f.write('INC_VLSV = -I'+makefile_options['INC_VLSV']+'\n')
            f.write('VLSV_PATH = -L'+makefile_options['VLSV_PATH']+'\n\n')

            f.write('INC_VECTORCLASS = -I'+makefile_options['INC_VECTORCLASS']+'\n\n')

            f.write('INC_DCCRG = -I'+makefile_options['INC_DCCRG']+'\n\n')

            f.write('INC_FSGRID = -I'+makefile_options['INC_FSGRID']+'\n\n')

            f.write('INC_JEMALLOC = -I'+makefile_options['INC_JEMALLOC']+'\n')
            f.write('JEMALLOC_PATH = -L'+makefile_options['JEMALLOC_PATH']+'\n\n')

            f.write('INC_PROFILE = -I'+makefile_options['INC_PROFILE']+'\n')
            f.write('PROFILE_PATH = -L'+makefile_options['PROFILE_PATH']+'\n\n')

            f.write('INC_EIGEN = -I'+makefile_options['INC_EIGEN']+'\n\n')

            f.write('INC_SILO = -I'+makefile_options['INC_SILO']+'\n')
            f.write('SILO_PATH = -L'+makefile_options['SILO_PATH']+'\n\n')
else:
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
print('  Machine:                    ' + (args['machine'] if args['machine'] else 'new'))
print('  Coordinate system:          ' + args['coord'])
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
      + makefile_options['COMPILER_FLAGS'])
