import subprocess
import numpy
import time
import random
import multiprocessing
import os


frameStart = 1000
frameEnd = 3069
jobNumber = 2

frames=numpy.arange(frameStart,frameEnd+1)
frameLists=numpy.array_split(frames, jobNumber)

def parallel_worker(frameList):
   time.sleep((os.getpid()%jobNumber)*6)
   
   subprocess.call(["/home/kempf/visit/bin/visit", str(frameList[0]), str(frameList[-1]), "-lb-random", "-cli", "-nowin", "-debug", "1", "-s", "generate_frames.py"])
   
   return
   


if __name__ == '__main__':
   pool = multiprocessing.Pool(jobNumber)
   pool.map(parallel_worker, frameLists)
   
