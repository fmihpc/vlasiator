#!/bin/bash

#sed -i'' 's/a_z\^2/pC[a_z]*pC[a_z]/g' codifiedCoeffs8.txt
#sed -i'' 's/a_xz\^2/pC[a_xz]*pC[a_xz]/g' codifiedCoeffs8.txt

#sed -i'' 's/a_y\^2/pC[a_y]*pC[a_y]/g' codifiedCoeffs8.txt
#sed -i'' 's/a_xy\^2/pC[a_xy]*pC[a_xy]/g' codifiedCoeffs8.txt

sed -i'' 's/dx\^2/dx2/g' codifiedCoeffs8.txt
sed -i'' 's/dy\^2/dy2/g' codifiedCoeffs8.txt
sed -i'' 's/dz\^2/dz2/g' codifiedCoeffs8.txt
#sed -i'' 's/dx\^3/dx3/g' codifiedCoeffs8.txt
#sed -i'' 's/dx\^4/dx4/g' codifiedCoeffs8.txt
#sed -i'' 's/dx\^5/dx5/g' codifiedCoeffs8.txt
#sed -i'' 's/dx\^6/dx6/g' codifiedCoeffs8.txt


sed -i'' 's/a_/pC[a_/g' codifiedCoeffs8.txt
sed -i'' 's/b_/pC[b_/g' codifiedCoeffs8.txt
sed -i'' 's/c_/pC[c_/g' codifiedCoeffs8.txt

sed -i'' 's/x\*/x]*/g' codifiedCoeffs8.txt
sed -i'' 's/y\*/y]*/g' codifiedCoeffs8.txt
sed -i'' 's/z\*/z]*/g' codifiedCoeffs8.txt
sed -i'' 's/0\*/0]*/g' codifiedCoeffs8.txt

sed -i'' 's/x)/x])/g' codifiedCoeffs8.txt
sed -i'' 's/y)/y])/g' codifiedCoeffs8.txt
sed -i'' 's/z)/z])/g' codifiedCoeffs8.txt
sed -i'' 's/0)/0])/g' codifiedCoeffs8.txt

sed -i'' 's/x\-/x]\-/g' codifiedCoeffs8.txt
sed -i'' 's/y\-/y]\-/g' codifiedCoeffs8.txt
sed -i'' 's/z\-/z]\-/g' codifiedCoeffs8.txt

sed -i'' 's/BGBX/cp[CellParams::BGBX]/g' codifiedCoeffs8.txt
sed -i'' 's/BGBY/cp[CellParams::BGBY]/g' codifiedCoeffs8.txt
sed -i'' 's/BGBZ/cp[CellParams::BGBZ]/g' codifiedCoeffs8.txt

sed -i'' 's/pC\[pC\[/pC[/g' codifiedCoeffs8.txt

sed -i'' 's/dx\]/dx/g' codifiedCoeffs8.txt
sed -i'' 's/dy\]/dy/g' codifiedCoeffs8.txt
sed -i'' 's/dz\]/dz/g' codifiedCoeffs8.txt
