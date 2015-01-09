#!/bin/bash

### run on laptop
# /home/kempf/visit/bin/visit -v 2.7.3 -lb-random -nowin -cli -s background_dots.py
# /home/kempf/visit/bin/visit -debug 1 -lb-random -nowin -cli -s background.py


### run on voima
# pwd
# python parallel_mosaic.py 17 1833


### run anywhere, e.g. crystal
# for i in $( seq -f "%07g" 0 1 )
# do
#    echo $i
#    montage ${i}/ACB_3dv_0*.png -tile 13x -geometry 100x100+50+50 -texture ACB_background_rho_dotted_50_2664x2332_${i}.png tile_${i}.png
# done

# seq -f "%07g" 0 20  | parallel -j 20 --progress montage {}/ACB_3dv_0*.png -tile 13x -geometry 100x100+50+50 -texture ACB_background_rho_dotted_50_2664x2332_{}.png tile_{}.png
