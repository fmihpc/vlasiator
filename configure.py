#!/usr/bin/env python
# ---------------------------------------------------------------------------------------
# configure.py: Athena++ configuration script in python. Original version by CJW.
#
# When configure.py is run, it uses the command line options and default settings to
# create custom versions of the files Makefile and src/defs.hpp from the template files
# Makefile.in and src/defs.hpp.in respectively.
#
# The following options are implememted:
#   -h  --help        help message
#   --prob=name       use src/pgen/name.cpp as the problem generator
#   --coord=xxx       use xxx as the coordinate system
#   --nghost=xxx      set NGHOST=xxx
#   --velocityorder=x set semilag solver order for velocity space acceleration (2,3,5)
#   --spatialorder=x  set semilag solver order for spatial acceleration (2,3,5)
#   -debug            enable debug flags (-g -O0); override other compiler options
#   -float            enable single precision (default double)
#   -distfloat        enable single precision for distribution function (default double)
#   -mpi              enable parallelization with MPI
#   -omp              enable parallelization with OpenMP
#   -profile          enable profiling (default on)
#   -noprofile        disable profiling
#   --vlsv_path=path  path to VLSV libraries (requires the VLSV library)
#   --cxx=xxx         use xxx as the C++ compiler
#   --cflag=string    append string whenever invoking compiler/linker
#   --include=path    use -Ipath when compiling
#   --lib_path=path   use -Lpath when linking
#   --lib=xxx         use -lxxx when linking
# ----------------------------------------------------------------------------------------

# Modules
import argparse
import glob
import re


# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'
defsfile_input = 'src/defs.hpp.in'
defsfile_output = 'src/defs.hpp'

# --- Step 1. Prepare parser, add each of the arguments ------------------
vlasiator_description = (
    "Prepare custom Makefile for compiling Vlasiator"
)
vlasiator_epilog = (
    "Full documentation of options available at "
    "https://github.com/PrincetonUniversity/athena-public-version/wiki/Configuring"
)
parser = argparse.ArgumentParser(description=vlasiator_description, epilog=vlasiator_epilog)

# --prob=[name]
pgen_directory = 'testpackage/tests'
# set pgen_choices to list of .cpp files in src/pgen/
pgen_choices = glob.glob(pgen_directory + '*.cpp')
# remove test folder prefix and '.cpp' extension from each filename
pgen_choices = [choice[len(pgen_directory):-4] for choice in pgen_choices]
parser.add_argument('--prob',
                    default='shock_tube',
                    choices=pgen_choices,
                    help='select problem generator')

# --coord=[name]
parser.add_argument(
    '--coord',
    default='cartesian',
    choices=[
        'cartesian',
        'cylindrical',
        'spherical_polar',
        'minkowski',
        'sinusoidal',
        'tilted',
        'schwarzschild',
        'kerr-schild',
        'gr_user'],
    help='select coordinate system')

# --nghost=[value]
parser.add_argument('--nghost',
                    default='2',
                    help='set number of ghost zones')

# --velocityorder=[value]
parser.add_argument('--velocityorder',
                    default=5,
                    help='set order of velocity space acceleration')

# --spatialorder=[value]
parser.add_argument('--spatialorder',
                    default=3,
                    help='set order of spatial translation')

# -debug
parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

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
                    default=False,
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
                    help='path to VLSV libraries')

# --zoltan_path=[name]
parser.add_argument('--zoltan_path',
                    default='',
                    help='path to Zoltan libraries')

# --jemalloc_path=[name]
parser.add_argument('--jemalloc_path',
                    default='',
                    help='path to jemalloc libraries')

# --phiprof_path=[name]
parser.add_argument('--phiprof_path',
                    default='',
                    help='path to phiprof libraries')

# The main choices for --cxx flag, using "ctype[-suffix]" formatting, where "ctype" is the
# major family/suite/group of compilers and "suffix" may represent variants of the
# compiler version and/or predefined sets of compiler options. The C++ compiler front ends
# are the main supported/documented options and are invoked on the command line, but the C
# front ends are also acceptable selections and are mapped to the matching C++ front end:
# gcc -> g++, clang -> clang++, icc-> icpc
cxx_choices = [
    'g++',
    'g++-simd',
    'icpc',
    'icpc-debug',
    'icpc-phi',
    'cray',
    'bgxlc++',
    'clang++',
    'clang++-simd',
    'clang++-apple',
]


def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)
    if arg == 'bgxl' or arg == 'bgxlc':
        arg = 'bgxlc++'

    if arg == 'clang':
        arg = 'clang++'
    else:
        arg = arg.replace('clang-', 'clang++-', 1)
    return arg


# --cxx=[name] argument
parser.add_argument(
    '--cxx',
    default='g++',
    type=c_to_cpp,
    choices=cxx_choices,
    help='select C++ compiler and default set of flags')

# --cflag=[string] argument
parser.add_argument('--cflag',
                    default=None,
                    help='additional string of flags to append to compiler/linker calls')

# --include=[name] arguments
parser.add_argument(
    '--include',
    default=[],
    action='append',
    help=('extra path for included header files (-I<path>); can be specified multiple '
          'times'))

# --lib_path=[name] arguments
parser.add_argument(
    '--lib_path',
    default=[],
    action='append',
    help=('extra path for linked library files (-L<path>); can be specified multiple '
          'times'))

# --lib=[name] arguments
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

# --prob=[name] argument
definitions['PROBLEM'] = makefile_options['PROBLEM_FILE'] = args['prob']

# --coord=[name] argument
definitions['COORDINATE_SYSTEM'] = makefile_options['COORDINATES_FILE'] = args['coord']

# --nghost=[value] argument
definitions['NUMBER_GHOST_CELLS'] = args['nghost']

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
if args['cxx'] == 'bgxlc++':
    # IBM XL C/C++ for BG/Q is NOT C++11 feature-complete as of v12.1.0.15 (2017-12-22)
    # suppressed messages:
    #   1500-036:  The NOSTRICT option has the potential to alter the program's semantics
    #   1540-1401: An unknown "pragma simd" is specified
    #   1586-083:  ld option ignored by IPA
    #   1586-233:  Duplicate definition of symbol ignored
    #   1586-267:  Inlining of specified subprogram failed due to the presence of a C++
    #                exception handler
    definitions['COMPILER_CHOICE'] = 'bgxlc++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'bgxlc++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -qhot=level=1:vector -qinline=level=5:auto -qipa=level=1:noobject'
      ' -qstrict=subnormals -qmaxmem=150000 -qlanglvl=extended0x -qsuppress=1500-036'
      ' -qsuppress=1540-1401 -qsuppress=1586-083 -qsuppress=1586-233'
      ' -qsuppress=1586-267'
    )
    makefile_options['LINKER_FLAGS'] = makefile_options['COMPILER_FLAGS']
    makefile_options['LIBRARY_FLAGS'] = ''
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

# distribution function precision
if args['distfloat']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DSPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4F_FALLBACK' # vector backend type
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DDPF'
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DVEC4D_FALLBACK' # vector backend type

if args['profile']:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DPROFILE'

#
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

#Add -DNDEBUG to turn debugging off.
if args['debug']:
    pass
else:
    makefile_options['PREPROCESSOR_FLAGS'] += ' -DNDEBUG'

# CXXFLAGS += -DDEBUG_SOLVERS
# CXXFLAGS += -DDEBUG_IONOSPHERE


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

# --- Step 4. Create new files, finish up --------------------------------

# Read templates
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()

# Make substitutions
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

# Write output files
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

# Finish with diagnostic output

print('Your Vlasiator distribution has now been configured with the following options:')
print('  Problem generator:          ' + args['prob'])
print('  Coordinate system:          ' + args['coord'])
print('  Debug flags:                ' + ('ON' if args['debug'] else 'OFF'))
print('  Linker flags:               ' + makefile_options['LINKER_FLAGS'] + ' '
      + makefile_options['LIBRARY_FLAGS'])
print('  Floating-point precision:   ' + ('single' if args['float'] else 'double'))
print('  Distribution precision:     ' + ('single' if args['distfloat'] else 'double'))
print('  Number of ghost cells:      ' + args['nghost'])
print('  MPI parallelism:            ' + ('ON' if args['mpi'] else 'OFF'))
print('  OpenMP parallelism:         ' + ('ON' if args['omp'] else 'OFF'))
print('  Semilog velocity order:     ' + str(args['velocityorder']))
print('  Semilog spatial order:      ' + str(args['spatialorder']))
print('  Profiler:                   ' + ('ON' if args['profile'] else 'OFF'))
print('  Compiler:                   ' + args['cxx'])
print('  Compilation command:        ' + makefile_options['COMPILER_COMMAND'] + ' '
      + makefile_options['PREPROCESSOR_FLAGS'] + ' ' + makefile_options['COMPILER_FLAGS'])