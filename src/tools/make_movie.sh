#!/bin/bash


aba_rho_lin_with_inset_scale() 
{
    index=$1
    t=$(echo $index|gawk '{printf "%.1f", $1*0.5}')
    in_file=$(printf "stornext_vlasiator/visualizations/ABA/movies/rho_lin/rho_lin_%07d.png" $index)
    if [ ! -e $in_file ]
    then
        return
    fi
    if [ ! -e makemovie_scale.png ]
    then
        convert  stornext_vlasiator/visualizations/ABA/movies/rho_lin/rho_lin_0000200.png -crop 750x500+2120+560 -resize 50%  makemovie_scale.png
    fi
    convert $in_file -crop 600x600+600+1900 -border 2  makemovie_zoom_${index}.png
    convert $in_file -pointsize 100 -fill black  -gravity SouthWest -draw "text 20,20 \"t=$t s\"" -gravity South -crop 3060x3292+0+0 -resize 50% makemovie_global_${index}.png
    convert  makemovie_global_${index}.png makemovie_zoom_${index}.png -gravity NorthEast -geometry +10+10 -composite  \
        makemovie_scale.png  -gravity NorthEast -geometry +650+60 -composite \
        -gravity NorthEast -pointsize 50 -fill white -draw "text 750,120 \"Density\"" \
        -gravity NorthWest -pointsize 70 -fill white -draw "text 50,150 \"VLASIATOR\""  \
        makemovie_${index}.png
}



aba_rho_log_with_inset() 
{
    index=$1
    t=$(echo $index|gawk '{printf "%.1f", $1*0.5}')
    in_file=$(printf "stornext_vlasiator/visualizations/ABA/movies/rho_log/rho_log_%07d.png" $index)
    if [ ! -e $in_file ]
    then
        return
    fi

    convert $in_file -crop 600x600+600+1900 -border 2  makemovie_zoom_${index}.png 
    convert $in_file -pointsize 100 -fill black  -gravity SouthWest -draw "text 20,20 \"t=$t s\"" \
        -gravity South -crop 3060x3292+0+0 -resize 50% \
        makemovie_global_${index}.png 
    convert makemovie_global_${index}.png makemovie_zoom_${index}.png -gravity NorthEast -geometry +10+10 -composite  \
        -gravity NorthEast -pointsize 50 -fill black -draw "text 700,150 \"Density, log scale\"" \
        -gravity NorthWest -pointsize 70 -fill black -draw "text 50,150 \"VLASIATOR\"" \
        makemovie_${index}.png
}

aba_B_log_with_inset() 
{
    index=$1
    t=$(echo $index|gawk '{printf "%.1f", $1*0.5}')
    in_file=$(printf "stornext_vlasiator/visualizations/ABA/movies/B_mag_log/ABA_B_mag_log_%07d.png" $index)
    if [ ! -e $in_file ]
    then
        return
    fi

    convert $in_file -crop 600x600+600+1900 -border 2  makemovie_zoom_${index}.png 
    convert $in_file -pointsize 100 -fill white  -gravity SouthWest -draw "text 20,20 \"t=$t s\"" \
        -gravity South -crop 3060x3292+0+0 -resize 50% \
        makemovie_global_${index}.png 
    convert makemovie_global_${index}.png makemovie_zoom_${index}.png -gravity NorthEast -geometry +10+10 -composite  \
        -gravity NorthWest -pointsize 50 -fill white -draw "text 50,250 \"B magnitude, log scale\"" \
        -gravity NorthWest -pointsize 70 -fill white -draw "text 50,150 \"VLASIATOR\"" \
        makemovie_${index}.png
}


aba_rho_B_log() 
{
    index=$1
    t=$(echo $index|gawk '{printf "%.1f", $1*0.5}')
    rho=$(printf "stornext_vlasiator/visualizations/ABA/movies/rho_log/rho_log_%07d.png" $index)
    B=$(printf "stornext_vlasiator/visualizations/ABA/movies/B_mag_log/ABA_B_mag_log_%07d.png" $index)
    if [ ! -e $rho ]
    then
        return
    fi
    if [ ! -e $B ]
    then
        return
    fi

    convert $rho \
        -gravity SouthWest -pointsize 150 -fill black -draw "text 20,3100 \"VLASIATOR\"" \
        -gravity SouthWest -pointsize 100 -fill black -draw "text 1100,3100 \"Density\"" \
        -gravity SouthWest -pointsize 100 -fill black -draw "text 20,20 \"t=$t s\"" \
        -gravity SouthWest -crop 1532x3292+0+0 \
        -resize 50% \
        makemovie_rho_${index}.png 
    convert $B \
        -gravity SouthWest -pointsize 100 -fill black -draw "text 900,3100 \"B magnitude\"" \
        -gravity SouthWest -crop 1532x3292+0+0  \
        -resize 50% \
        makemovie_B_${index}.png 
    convert makemovie_rho_${index}.png makemovie_B_${index}.png +append \
        makemovie_${index}.png

}





export -f aba_rho_lin_with_inset_scale
export -f aba_rho_log_with_inset
export -f aba_B_log_with_inset
export -f aba_rho_B_log

#main script starts

if [[ $# -lt 3 ]]
then
    echo "Usage: $0 picture_generator_function first_frame last_frame [movie_name]"
    echo "  picture_generator_function  Bash function for generating the different frame images. Modify script to add a new one if needed. "
    echo "  first_frame                 Index of first frame"
    echo "  last_frame                  Index of last frame"
    echo "  movie_name                  Name for movie, frames converted to movie only of this is defined"
    echo "  "
    echo "List of picture_generator_functions:"
    declare -x -F 
    exit 
fi

frame_generator=$1
first_frame=$2
last_frame=$3
movie_name=$4

#All temporary files should be prefixed by makemovie_, these are all nuked when starting a new movie generation
rm makemovie_*.png

parallel --progress $frame_generator :::  $(seq $first_frame $last_frame) 

if [[ $# -eq 4 ]]
then

#    avconv -threads auto -f image2 -framerate 10  -i frame_%d.png  -pix_fmt yuv420p  -b 4000k  $movie_name
#MAC safe movie
    avconv -threads 8 -f  image2  -start_number $first_frame -i makemovie_%d.png  -b 4000k -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -f mp4 -vcodec libx264  -pix_fmt yuv420p $movie_name
#    rm frame*.png
fi
