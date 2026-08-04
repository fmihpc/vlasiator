#!/bin/sh

echo "--------------------------------------------------------------------------------------------"
echo "   Testing write path                                                                       "
echo "--------------------------------------------------------------------------------------------"

filepath_test="./testing/foo/bar/bulk.0000000.vlsv" 
if [[ ! -f $filepath_test ]]; then 
   echo "File was not written to correct path at $filepath_test"
   exit 1
fi
mv ./testing/foo/bar/* ./
echo "Write path test successful!"
