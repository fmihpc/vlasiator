#!/bin/bash
#general code for processing images. Convert command in user supplied function
process_images() 
{
    in_file=$1
    first_frame=$2
    last_frame=$3
    convert_function=$4
    file_index=$(echo $in_file |  gawk -F[_.] '{printf "%d",$(NF-1)}')
    if [[ $file_index -ge $first_frame ]]
    then
        if [[ $file_index -lt $last_frame ]]
        then
            #FIXME, hardcoded time vs index
            t=$(echo $file_index|gawk '{printf "%.1f", $1*0.5}')
            out_file=frame_$((file_index - first_frame)).png            
            $4 $in_file $out_file $t
        fi
    fi
}

#Define here convert functions

convert_global_image() 
{
    in_file=$1
    out_file=$2
    t=$3
    convert $in_file -pointsize 100 -fill white  -gravity SouthWest -draw "text 20,20 \"With Hall t=$t s\"" -gravity South -crop 100%x94% -resize 50% $out_file 
}

convert_zoom_qp() 
{
    in_file=$1
    out_file=$2
    t=$3
    convert $in_file -pointsize 20 -fill black  -gravity SouthWest -draw "text 510,1760 \"With Hall t=$t s\"" -crop 800x800+500+1000 $out_file 
}

convert_with_inset() 
{
    in_file=$1
    out_file=$2
    t=$3
    convert $in_file -crop 600x600+600+1900 -border 2 zoom_$out_file 
    convert $in_file -pointsize 100 -fill white  -gravity SouthWest -draw "text 20,20 \"With Hall t=$t s\"" -gravity South -crop 100%x94% -resize 50% global_$out_file 
    convert global_$out_file zoom_$out_file -gravity NorthEast -geometry +10+10 -composite   $out_file
    rm  zoom_$out_file global_$out_file
}


export -f process_images
export -f convert_global_image
export -f convert_zoom_qp
export -f convert_with_inset


#main script starts

if [[ $# -lt 3 ]]
then
    echo "Usage: $0 folder_with_original_images  first_frame last_frame [movie_name]"
    echo "Modify script to change how the picture is produced from original one. Movie only made if its name is defined"
    exit
fi

folder=$1
first_frame=$2
last_frame=$3
movie_name=$4

#process_global_image out.png 400
rm frame*.png

label="With Hall Term"
geometry_options="-gravity South -crop 100%x94% -resize 50%"

parallel  process_images :::  $(ls ${folder}/*.png) ::: $first_frame ::: $last_frame ::: convert_with_inset

if [[ $# -eq 4 ]]
then
    avconv -threads auto -f image2 -framerate 10  -i frame_%d.png  -pix_fmt yuv420p  -b 4000k  $movie_name
    rm frame*.png
fi
