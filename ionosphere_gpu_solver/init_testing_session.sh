#!/bin/sh


/usr/bin/srun --interactive -n1 -c8 --mem=4G -t01:00:00 -Mukko --constraint=v100 -pgpu --pty bash


