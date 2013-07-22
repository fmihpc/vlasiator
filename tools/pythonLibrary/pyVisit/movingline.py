#!/usr/bin/python -i

import visit as vis
import loadvisitsettings as visSettings
import subprocess
from useroptions import *

def make_moving_frame_of_reference_line_plot(  point1, point2, velocity, variable_name, input_directory, input_file_name, output_directory, output_file_name, start_frame=-1, end_frame=-1, frame_skip_dt=1.0 ):
   '''
   Function for making a line plot of some variable with a moving frame of reference
   :param point1              The starting point of the line (must be an array of size 3)
   :param point2              The ending point of the line (must be an array of size 3)
   :param velocity            The velocity vector of the frame of reference (must be an array of size 3)
   :param variable_name       Name of the variable (For ex \"rho\")
   :param input_directory     The path to the directory where the files are
   :param input_file_name     Name of the files (For ex \"bulk.*.silo\")
   :param output_directory    Directory where to output the movie
   :param output_file_name    Name of the outputted file (For ex \"RHOMOVIE\")
   :param start_frame         Starting frame for the movie (if -1, equals 0, -1 by default)
   :param end_frame           Ending frame for the movie (if -1, equals the last frame, -1 by default)
   :param frame_skip_dt       The number of seconds one skip in frame equals (1.0 by default) (Note: This may change depending on the run and should always be checked)
   '''
   if len(point1) != 3 or len(point2) != 3 or len(velocity) != 3:
      print "BAD INPUT IN make_moving_frame_of_reference_line_plot, POINT1, POINT2 AND VELOCITY MUST BE ARRAYS OF SIZE 3"
   
   # OPTIONS
   #################################################################
   # Input the boundary box for starting coordinates (Starting values)
   startX = point1[0] # The left x-boundary of the box
   endX = point2[0] # The right x-boundary of the box
   startY = point1[1] # The bottom y-boundary of the box
   endY = point2[1] # The upper y-boundary of the box
   startZ = poin1[2] # The left z-boundary of the box
   endZ = point2[2] # The right z-boundary of the box
   
   # Input frame properties
   startFrame = start_frame # Note: if startFrame is set to -1 the start frame gets set to 0
   endFrame = end_frame # Note: if endFrame is set to -1 the endFrame is automatically the number of frames in the database
   frameInSeconds = frame_skip_dt # Set how many seconds one frame skip is

   screenWidth = 3000
   screenHeight = 3000
   
   # Input speed in x and y direction
   speedX = velocity[0] # Meters per second
   speedY = velocity[1] # Meters per second
   speedZ = velocity[2] # Meters per second

   # Input variable name
   # Note: needs to have operators/Lineout/ for visit to recognize it as line plot. Additionally, visit does not accept any '/' in the variable name which is why they're removed. The curve definitions are in loadvisitsettings.py and in there the curve expressions are defined so that there's no '/' in the variable name
   variableName = "operators/Lineout/" + variable_name.replace("/","")
   
   # Input directory and file names
   outputDir = output_directory # Set the output directory (Where .png s are saved)
   outputFileName = output_file_name # The file names for the png files. These for ex. will be saved visit0000.png, visit0001.png, ..
   databaseName = "localhost:" + input_directory + input_file_name + " database" # For navigating to the silo files
   # visitBinDirectory = "/usr/lvariableNameocal/visit/bin" #Nevermind this
   # Note: a slice of the plot in z-axis is taken automatically
   #################################################################
   
   
   dx = speedX * frameInSeconds # Note: This is in meters per frame!
   dy = speedY * frameInSeconds # Note: This is in meters per frame!
   dz = speedZ * frameInSeconds # Note: This is in meters per frame!
   
   
   
   vis.OpenDatabase(databaseName, 0)
   #Load settings
   visSettings.load_visit_settings()

   vis.AddPlot("Curve", variableName, 1, 1)
   vis.LineoutAtts = vis.LineoutAttributes()
   vis.LineoutAtts.point1 = (startX, startY, 0)
   vis.LineoutAtts.point2 = (endX, endY, 0)
   vis.LineoutAtts.interactive = 0
   vis.LineoutAtts.ignoreGlobal = 0
   vis.LineoutAtts.samplingOn = 0
   vis.LineoutAtts.numberOfSamplePoints = 50
   vis.LineoutAtts.reflineLabels = 0
   vis.SetOperatorOptions(vis.LineoutAtts, 1)
   vis.CurveAtts = vis.CurveAttributes()
   vis.CurveAtts.showLines = 1
   vis.CurveAtts.lineStyle = vis.CurveAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.CurveAtts.lineWidth = 2
   vis.CurveAtts.showPoints = 1
   vis.CurveAtts.symbol = vis.CurveAtts.Point  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
   vis.CurveAtts.pointSize = 5
   vis.CurveAtts.pointFillMode = vis.CurveAtts.Static  # Static, Dynamic
   vis.CurveAtts.pointStride = 1
   vis.CurveAtts.symbolDensity = 50
   vis.CurveAtts.curveColorSource = vis.CurveAtts.Custom  # Cycle, Custom
   vis.CurveAtts.curveColor = (0, 0, 0, 255)
   vis.CurveAtts.showLegend = 1
   vis.CurveAtts.showLabels = 0
   vis.CurveAtts.designator = ""
   vis.CurveAtts.doBallTimeCue = 0
   vis.CurveAtts.ballTimeCueColor = (0, 0, 0, 255)
   vis.CurveAtts.timeCueBallSize = 0.01
   vis.CurveAtts.doLineTimeCue = 0
   vis.CurveAtts.lineTimeCueColor = (0, 0, 0, 255)
   vis.CurveAtts.lineTimeCueWidth = 0
   vis.CurveAtts.doCropTimeCue = 0
   vis.CurveAtts.timeForTimeCue = 0
   vis.SetPlotOptions(vis.CurveAtts)
   vis.DrawPlots()
   
   # Iterate through frames
   for i in xrange(startFrame, endFrame+1):
      vis.SetTimeSliderState(i)
      frame = i - startFrame
      vis.LineoutAtts = vis.LineoutAttributes()
      vis.LineoutAtts.point1 = (startX + frame*dx, startY + frame*dy, 0)
      vis.LineoutAtts.point2 = (endX + frame*dx, endY + frame*dy, 0)
      vis.LineoutAtts.interactive = 0
      vis.LineoutAtts.ignoreGlobal = 0
      vis.LineoutAtts.samplingOn = 0
      vis.LineoutAtts.numberOfSamplePoints = 50
      vis.LineoutAtts.reflineLabels = 0
      vis.SetOperatorOptions(vis.LineoutAtts, 1)
      vis.SaveWindowAtts = vis.SaveWindowAttributes()
      vis.SaveWindowAtts.outputToCurrentDirectory = 0
      vis.SaveWindowAtts.outputDirectory = outputDir
      vis.SaveWindowAtts.fileName = outputFileName
      vis.SaveWindowAtts.family = 1
      vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
      vis.SaveWindowAtts.width = screenWidth
      vis.SaveWindowAtts.height = screenHeight
      vis.SaveWindowAtts.screenCapture = 0
      vis.SaveWindowAtts.saveTiled = 0
      vis.SaveWindowAtts.quality = 80
      vis.SaveWindowAtts.progressive = 0
      vis.SaveWindowAtts.binary = 0
      vis.SaveWindowAtts.stereo = 0
      vis.SaveWindowAtts.compression = vis.SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
      vis.SaveWindowAtts.forceMerge = 0
      vis.SaveWindowAtts.resConstraint = vis.SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
      vis.SaveWindowAtts.advancedMultiWindowSave = 0
      vis.SetSaveWindowAttributes(vis.SaveWindowAtts)
      vis.SaveWindow()
   vis.DeleteActivePlots()
   vis.CloseDatabase(databaseName)
   # Make the movie:
   framerate = 5
   subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", outputDir, outputFileName, framerate])
