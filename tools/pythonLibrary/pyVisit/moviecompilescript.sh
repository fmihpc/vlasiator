#!/bin/bash

filePath=$1
fileName=$2
frameRate=$3
cd ${filePath}
avconv -r ${frameRate} -i ${fileName}%04d.png -r ${frameRate} -b 40M ${fileName}.mp4
cd -
