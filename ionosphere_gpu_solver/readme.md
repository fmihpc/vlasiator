# Compilation

Command to load required modules:

`
module purge; module load CUDA
`

Compiles into **libionosphere_gpu_solver.a** using:

`
make
`

Makefile also supports following commands:

Recompile:

`
make re
`

Clean all compilation artifacts:

`
make clean
`

Same as clean but also removes the library itself:

`
make cleanf
`


# Linking
Linking with **libionosphere_gpu_solver.a** one needs to also link with cuda runtime. Following compiler options should be enought:

`
-lcudart -L"path to inosphere_gpu_solver" -lionosphere_gpu_solver 
`

# Testing

Usually one do not need to interact with testing part. It is used as tool for development.

All the test related code is in **testing** folder. Testing uses Meson build system. One can obtain meson by loading it as module

`
module load Meson
`

Now one can create a build folder in **testing** folder with (meson is so called out of source build tool so all files that meson generates are inside build folder and it can always be deleted without worrying. It can always be regenerated)

`
meson setup build
`

Meson will create build folder and in it one can run the tests using 

`
meson test
`

## Writing tests

In **testsing/tests** folder tests are written in **prototypes** folder. In each test there is line

`
using test_type = TEST_TYPE_PROTOTYPE
`

when running the tests script will be run that will generate two source files where in one TEST_TYPE_PROTOTYPE is replaced with float and another where it is replace with double.  