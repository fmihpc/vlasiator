# Usage

1. Create a folder my_testcase_name in this directory
2. Add sw1.dat and my_testcase_name.cfg to the folder
3. Use a pre existing job script or create your own with the same format
4. If necessary, change which dimensions to scale and which partition to use in weakScaling.sh
5. Execute the weak scaling script using: 
```
./weakScaling.sh my_testcase_name job_script_name.sh MIN_EXPONENT MAX_EXPONENT
```
Where you need to replace MIN_EXPONENT and MAX_EXPONENT with appropriate values. The script will then be executed using 2^MIN_EXPONENT to 2^MAX_EXPONENT GPUs.