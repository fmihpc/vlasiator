import re
import os
import argparse


parser = argparse.ArgumentParser(
    prog="Logparser",
    description="For pruning restarts from the logfile such that the latest successful timestep is left in the log",
)
parser.add_argument("logfile", type=str)
parser.add_argument("outputfile", type=str)
args = parser.parse_args()
input = args.logfile
output = args.outputfile

logdict = {}
with open(input, "r") as logfile:
    loglines = logfile.readlines()
    iterator = iter(loglines)
    n = 0
    line = iterator.__next__()
    logdict[-1] = []
    while n < len(loglines):
        try:
            regex = re.search(r"^-{10}\Wtstep\W=\W(\d+).+?-{10}", line)
            if regex:
                tstep = regex.group(1)

                logdict[tstep] = [line]
                line = iterator.__next__()
                n += 1
                while not re.search(r"^-{10}\Wtstep\W=\W(\d+).+?-{10}", line):
                    logdict[tstep].append(line)
                    # we dont want to do this at last step
                    try:
                        line = iterator.__next__()
                    except StopIteration:
                        break
                    n += 1
            else:
                logdict[-1].append(line)
                line = iterator.__next__()

        except StopIteration:
            break
    logfile.close()

with open(output, "w") as fileoutput:
    for lines in logdict.values():
        for line in lines:
            fileoutput.write(line)
