#!/usr/bin/python -i

import visit as vis
import loadvisitsettings as visSettings
import subprocess
from useroptions import *

def make_moving_frame_of_reference_movie( x_begin, x_end, y_begin, y_end, speed_x, speed_y, variable_name, minThreshold, maxThreshold, input_directory, input_file_name, output_directory, output_file_name, color_table="hot_desaturated", start_frame=-1, end_frame=-1, frame_skip_dt=1.0 ):
   '''
   Function for making a movie with a moving frame of reference.
   :param x_begin             The starting frame's beginning x-coordinate
   :param x_end               The starting frame's ending x-coordinate
   :param y_begin             The starting frame's beginning x-coordinate
   :param y_end               The starting frame's ending y-coordinate
   :param speed_x             The speed at which the frame moves in the x direction
   :param speed_y             The speed at which the frame moves in the y direction
   :param variable_name       Name of the variable (For ex \"rho\")
   :param minThreshold        Minimum threshold for the variable
   :param maxThreshold        Maximum threshold for the variable
   :param input_directory     The path to the directory where the files are
   :param input_file_name     Name of the files (For ex \"bulk.*.silo\")
   :param output_directory    Directory where to output the movie
   :param output_file_name    Name of the outputted file (For ex \"RHOMOVIE\")
   :param color_table         Name of the color table (\"hot_desaturated\" by default)
   :param start_frame         Starting frame for the movie (if -1, equals 0, -1 by default)
   :param end_frame           Ending frame for the movie (if -1, equals the last frame, -1 by default)
   :param frame_skip_dt       The number of seconds one skip in frame equals (1.0 by default) (Note: This may change depending on the run and should always be checked)
   '''
   # OPTIONS
   #################################################################
   # Input the boundary box for starting coordinates (Starting values)
   startX = x_begin # The left x-boundary of the box
   endX = x_end # The right x-boundary of the box
   startY = y_begin # The bottom y-boundary of the box
   endY = y_end # The upper y-boundary of the box
   
   # Input frame properties
   startFrame = start_frame # Note: if startFrame is set to -1 the start frame gets set to 0
   endFrame = end_frame # Note: if endFrame is set to -1 the endFrame is automatically the number of frames in the database
   frameInSeconds = frame_skip_dt # Set how many seconds one frame skip is
   
   # Input speed in x and y direction
   speedX = speed_x # Meters per second
   speedY = speed_y # Meters per second
   
   # Input variable
   variableName = variable_name
   minVariableValue = minThreshold
   maxVariableValue = maxThreshold
   colorTableName = color_table
   
   # Input directory and file names
   outputDir = output_directory # Set the output directory (Where .png s are saved)
   outputFileName = output_file_name # The file names for the png files. These for ex. will be saved visit0000.png, visit0001.png, .
   databaseName = "localhost:" + input_directory + input_file_name + " database" # For navigating to the silo files
   # visitBinDirectory = "/usr/local/visit/bin" #Nevermind this
   # Note: a slice of the plot in z-axis is taken automatically
   #################################################################
   
   
   # Launch visit
   # LaunchNowin(vdir=visitBinDirectory)
   dx = speedX * frameInSeconds # Note: This is in meters per frame!
   dy = speedY * frameInSeconds # Note: This is in meters per frame!
   #Set up window and annotations
   OpenDatabase(databaseName, 0)
   #Load settings
   visSettings.load_visit_settings()

   vis.AddPlot("Pseudocolor", variableName, 1, 1)
   vis.SetActivePlots(0)
   vis.PseudocolorAtts = vis.PseudocolorAttributes()
   vis.PseudocolorAtts.legendFlag = 1
   vis.PseudocolorAtts.lightingFlag = 1
   vis.PseudocolorAtts.minFlag = 1
   vis.PseudocolorAtts.maxFlag = 1
   vis.PseudocolorAtts.centering = vis.PseudocolorAtts.Natural  # Natural, Nodal, Zonal
   vis.PseudocolorAtts.scaling = vis.PseudocolorAtts.Linear  # Linear, Log, Skew
   vis.PseudocolorAtts.limitsMode = vis.PseudocolorAtts.CurrentPlot  # OriginalData, CurrentPlot
   vis.PseudocolorAtts.min = minVariableValue
   vis.PseudocolorAtts.max = maxVariableValue
   vis.PseudocolorAtts.pointSize = 0.05
   vis.PseudocolorAtts.pointType = vis.PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
   vis.PseudocolorAtts.skewFactor = 1
   vis.PseudocolorAtts.opacity = 1
   vis.PseudocolorAtts.colorTableName = vis.colorTableName
   vis.PseudocolorAtts.invertColorTable = 0
   vis.PseudocolorAtts.smoothingLevel = 0
   vis.PseudocolorAtts.pointSizeVarEnabled = 0
   vis.PseudocolorAtts.pointSizeVar = "default"
   vis.PseudocolorAtts.pointSizePixels = 2
   vis.PseudocolorAtts.lineStyle = vis.PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.PseudocolorAtts.lineWidth = 0
   vis.PseudocolorAtts.opacityType = vis.PseudocolorAtts.Explicit  # Explicit, ColorTable
   vis.SetPlotOptions(PseudocolorAtts)
   vis.SetActivePlots(0)
   vis.AddOperator("Slice", 1)
   vis.AddOperator("Threshold", 1)
   vis.ThresholdAtts = vis.ThresholdAttributes()
   vis.ThresholdAtts.outputMeshType = 0
   vis.ThresholdAtts.listedVarNames = ("Boundary_type")
   vis.ThresholdAtts.zonePortions = (1)
   vis.ThresholdAtts.lowerBounds = (1)
   vis.ThresholdAtts.upperBounds = (1)
   vis.ThresholdAtts.defaultVarName = variableName
   vis.ThresholdAtts.defaultVarIsScalar = 1
   vis.SetOperatorOptions(ThresholdAtts, 1)
   vis.ThresholdAtts = vis.ThresholdAttributes()
   vis.ThresholdAtts.outputMeshType = 0
   vis.ThresholdAtts.listedVarNames = ("Boundary_type")
   vis.ThresholdAtts.zonePortions = (1)
   vis.ThresholdAtts.lowerBounds = (1)
   vis.ThresholdAtts.upperBounds = (1)
   vis.ThresholdAtts.defaultVarName = variableName
   vis.ThresholdAtts.defaultVarIsScalar = 1
   vis.SetOperatorOptions(ThresholdAtts, 1)
   vis.SetActivePlots(0)
   vis.SliceAtts = vis.SliceAttributes()
   vis.SliceAtts.originType = vis.SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
   vis.SliceAtts.originPoint = (0, 0, 0)
   vis.SliceAtts.originIntercept = 0
   vis.SliceAtts.originPercent = 0
   vis.SliceAtts.originZone = 0
   vis.SliceAtts.originNode = 0
   vis.SliceAtts.normal = (0, 0, 1)
   vis.SliceAtts.axisType = vis.SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
   vis.SliceAtts.upAxis = (0, 1, 0)
   vis.SliceAtts.project2d = 1
   vis.SliceAtts.interactive = 1
   vis.SliceAtts.flip = 0
   vis.SliceAtts.originZoneDomain = 0
   vis.SliceAtts.originNodeDomain = 0
   vis.SliceAtts.meshName = "SpatialGrid"
   vis.SliceAtts.theta = 0
   vis.SliceAtts.phi = 90
   vis.SetOperatorOptions(SliceAtts, 1)
   vis.DrawPlots()
   
   if endFrame == -1:
      endFrame = vis.TimeSliderGetNStates() - 1
   
   if startFrame == -1:
      startFrame = 0
   
   # Iterate through frames
   for i in xrange(startFrame, endFrame+1):
      vis.SetTimeSliderState(i)
      frame = i - startFrame
      vis.View2DAtts = vis.View2DAttributes()
      vis.View2DAtts.windowCoords = (startX + frame*dx, endX + frame*dx, startY + frame*dy, endY + frame*dy)
      vis.View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
      vis.View2DAtts.fullFrameActivationMode = vis.View2DAtts.Auto  # On, Off, Auto
      vis.View2DAtts.fullFrameAutoThreshold = 100
      vis.View2DAtts.xScale = vis.View2DAtts.LINEAR  # LINEAR, LOG
      vis.View2DAtts.yScale = vis.View2DAtts.LINEAR  # LINEAR, LOG
      vis.View2DAtts.windowValid = 1
      vis.SetView2D(vis.View2DAtts)
      vis.SaveWindowAtts = vis.SaveWindowAttributes()
      vis.SaveWindowAtts.outputToCurrentDirectory = 0
      vis.SaveWindowAtts.outputDirectory = outputDir
      vis.SaveWindowAtts.fileName = outputFileName
      vis.SaveWindowAtts.family = 1
      vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
      vis.SaveWindowAtts.width = 1024
      vis.SaveWindowAtts.height = 1024
      vis.SaveWindowAtts.screenCapture = 0
      vis.SaveWindowAtts.saveTiled = 0
      vis.SaveWindowAtts.quality = 100
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
   framerate = 7
   subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", outputDir, outputFileName, framerate])
