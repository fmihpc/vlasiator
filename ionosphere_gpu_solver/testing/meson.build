project(
    'ionosphere_gpu_testing',
    'cpp',
    default_options : ['warning_level=3', 'werror=true', 'optimization=3', 'cpp_std=c++20']
)

####################################################################################################
# Variables tracking sources and libraries
####################################################################################################

# Directories where compiler is able to include from
project_include_directories = []

dependencies = []

####################################################################################################
# Implementation 
####################################################################################################
project_include_directories += include_directories('include', '../include')

# Cuda modules can be added, eg: dependencies += dependency('cuda', version : '>=10', modules : ['cublas'])
dependencies += dependency('cuda', version : '>=10')

dependencies += meson.get_compiler('cpp').find_library(
    'ionosphere_gpu_solver',
    dirs : ['/wrk-vakka/users/pamiro/vlasiator/vlasiator/ionosphere_gpu_solver/']
    )

####################################################################################################
# Tests
####################################################################################################
subdir('tests')
