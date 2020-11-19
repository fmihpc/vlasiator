import subprocess
import os.path
import sys

timeList = range(int(sys.argv[1]), int(sys.argv[2])+1)



def parallel_worker(timeValue):
   time = str(timeValue).rjust(7, '0')
   print time
   subprocess.call(["mkdir", time])
   subprocess.call(["ln", "-s", "/lustre/tmp/alfthan/2D/ACB_yk/bulk."+time+".vlsv", "."])
   
   subprocess.call(["/lustre/tmp/yann/visit/bin/visit", "-cli", "-nowin", "-s", "mosaic_worker.py", time])
   
   subprocess.call(["rm", "bulk."+time+".vlsv"])



from multiprocessing import Pool
if __name__ == '__main__':
   pool = Pool(6)
   pool.map(parallel_worker, timeList)






### serial:
#for cell in cellList:
   #parallel_worker(cell)

exit()

