import glob
import os

def get_sorted_file_names(name="*.vlsv"):
   fileNames=glob.glob(name)
   fileNames.sort()
   return fileNames
   
