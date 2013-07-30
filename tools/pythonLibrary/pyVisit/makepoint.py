import numpy as np
import os
import visit as vis
import loadvisitsettings as visSettings
import subprocess
from useroptions import *

def create_point_vtk( fileName, coordinates ):
   f = open(fileName,"w")
   f.truncate()
   f.write("# vtk DataFile Version 2.0\nWedges test\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 1 float\n")
   f.write(str(coordinates[0]) + " " + str(coordinates[1]) + " " + str(coordinates[2]) + "\n")
   f.write("CELLS 1 2\n1 0\nCELL_TYPES 1\n1\n")
   f.close()

def draw_point_picture( variableName, minValue, maxValue, inputDirectory, inputFileNames, coordinate, outputDirectory, outputFileName, colorTable="hot_desaturated"):
   '''
   Function for making a visit plot with a point
   
   Arguments:
   :param variableName                  Name of the variable
   :param minValue                      Minimum value of the variable
   :param maxValue                      Maximum value of the variable
   :param inputDirectory                Path to input vlsv/silo files
   :param inputFileNames                Name of the file, for example \"bulk.00000.silo\"
   :param coordinates                   Coordinates corresponding to the files so for example [ [[0,0,0], [0,1,0]], [[2,1,2], [2,1,4]] ]
   :param outputDirectory               Path to output directory
   :param outputFileName                Name of the output file
   :param colorTable="hot_desaturated"  Color table for the plots
   '''
   # OPTIONS
   #################################################################
   
   # Input variable
   _variableName = variableName
   minVariableValue = minValue
   maxVariableValue = maxValue
   colorTableName = colorTable
   
   # Input directory and file names
   _outputDir = outputDirectory
   _outputFileName = outputFileName # The file names for the png files.
   inputFileName = inputFileNames
   databaseName = "localhost:" + inputDirectory + inputFileName # For navigating to the silo files
   # Note: a slice of the plot in z-axis is taken automatically
   #################################################################

   inputFileName2 = "point.vtk"
   databaseName2 = "localhost:" + os.getcwd() + "/" + inputFileName2

   currentPlot = 0

   vis.OpenDatabase(databaseName, 0)
   #vis.ActiveDatabase("localhost:" + inputDirectory + inputFileName)
   #Load settings
   visSettings.load_visit_settings()
   vis.AddPlot("Pseudocolor", _variableName, 1, 1) #CONTINUE
   vis.SetActivePlots(currentPlot)
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
   vis.PseudocolorAtts.colorTableName = colorTableName
   vis.PseudocolorAtts.invertColorTable = 0
   vis.PseudocolorAtts.smoothingLevel = 0
   vis.PseudocolorAtts.pointSizeVarEnabled = 0
   vis.PseudocolorAtts.pointSizeVar = "default"
   vis.PseudocolorAtts.pointSizePixels = 2
   vis.PseudocolorAtts.lineStyle = vis.PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.PseudocolorAtts.lineWidth = 0
   vis.PseudocolorAtts.opacityType = vis.PseudocolorAtts.Explicit  # Explicit, ColorTable
   vis.SetPlotOptions(vis.PseudocolorAtts)

   
   vis.AddOperator("Slice", 1)
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
   vis.SetOperatorOptions(vis.SliceAtts, 1)
   vis.DrawPlots()

   create_point_vtk( fileName=inputFileName2, coordinates=coordinate )
   vis.OpenDatabase(databaseName2, 0)
   currentPlot = currentPlot + 1
   vis.SetActivePlots(currentPlot)
   vis.AddPlot("Mesh", "mesh", 1, 1)
   vis.MeshAtts = vis.MeshAttributes()
   vis.MeshAtts.legendFlag = 1
   vis.MeshAtts.lineStyle = vis.MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.MeshAtts.lineWidth = 0
   vis.MeshAtts.meshColor = (0, 0, 0, 255)
   vis.MeshAtts.outlineOnlyFlag = 0
   vis.MeshAtts.errorTolerance = 0.01
   vis.MeshAtts.meshColorSource = vis.MeshAtts.Foreground  # Foreground, MeshCustom
   vis.MeshAtts.opaqueColorSource = vis.MeshAtts.Background  # Background, OpaqueCustom
   vis.MeshAtts.opaqueMode = vis.MeshAtts.Auto  # Auto, On, Off
   vis.MeshAtts.pointSize = 0.05
   vis.MeshAtts.opaqueColor = (255, 255, 255, 255)
   vis.MeshAtts.smoothingLevel = vis.MeshAtts.None  # None, Fast, High
   vis.MeshAtts.pointSizeVarEnabled = 0
   vis.MeshAtts.pointSizeVar = "default"
   vis.MeshAtts.pointType = vis.MeshAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
   vis.MeshAtts.showInternal = 0
   vis.MeshAtts.pointSizePixels = 25
   vis.MeshAtts.opacity = 1
   vis.SetPlotOptions(vis.MeshAtts)
   vis.DrawPlots()
   vis.SaveWindowAtts = vis.SaveWindowAttributes()
   vis.SaveWindowAtts.outputToCurrentDirectory = 0
   vis.SaveWindowAtts.outputDirectory = _outputDir
   vis.SaveWindowAtts.fileName = _outputFileName
   vis.SaveWindowAtts.family = 1
   vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
   vis.SaveWindowAtts.width = 3000
   vis.SaveWindowAtts.height = 3000
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
   vis.SetActivePlots((0,1))
   vis.DeleteActivePlots()
   vis.CloseDatabase(databaseName2)
   vis.CloseDatabase(databaseName)

def create_visit_point_movie( variableName, minValue, maxValue, inputDirectory, inputFileNames, coordinates, outputDirectory, outputFileName, colorTable="hot_desaturated"):
   '''
   Function for making a movie
   
   Arguments:
   :param variableName                  Name of the variable
   :param minValue                      Minimum value of the variable
   :param maxValue                      Maximum value of the variable
   :param inputDirectory                Path to input vlsv/silo files
   :param inputFileNames                Name of the files for example [\"bulk.00000.silo\", \"bulk.00001.silo\"]
   :param coordinates                   Coordinates corresponding to the files so for example [ [[0,0,0], [0,1,0]], [[2,1,2], [2,1,4]] ]
   :param outputDirectory               Path to output directory
   :param outputFileName                Name of the output file
   :param colorTable="hot_desaturated"  Color table for the plots
   '''
   coordinates=[coordinates]

   for i in xrange(len(inputFileNames)):
      # OPTIONS
      #################################################################
      
      # Input variable
      _variableName = variableName
      minVariableValue = minValue
      maxVariableValue = maxValue
      colorTableName = colorTable
      
      # Input directory and file names
      #_outputDir = "/home/hannukse/MOVINGFRAME_MOVIES/AAJ_BZ_REMAKE/" # Set the output directory (Where .png s are saved)
      _outputDir = outputDirectory
      #_outputFileName = "BZ_FORESHOCK_2_" # The file names for the png files. These for ex. will be saved visit0000.png, visit0001.png, ..
      _outputFileName = outputFileName # The file names for the png files.
      #databaseName = "localhost:/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/bulk.*.silo database" # For navigating to the silo files
      inputFileName = inputFileNames[i]
      databaseName = "localhost:" + inputDirectory + inputFileName # For navigating to the silo files
      # Note: a slice of the plot in z-axis is taken automatically
      #################################################################
      # LaunchNowin(vdir=visitBinDirectory)
      #dx = speedX * frameInSeconds # Note: This is in meters per frame!
      #dy = speedY * frameInSeconds # Note: This is in meters per frame!
      #LaunchNowin(vdir="/usr/local/visit/bin")
      #Set up window and annotations
      #vis.LaunchNowin(vdir="/usr/local/visit/bin")

      inputFileName2 = "point.vtk"
      databaseName2 = "localhost:" + os.getcwd() + "/" + inputFileName2

      vis.OpenDatabase(databaseName, 0)
      #vis.ActiveDatabase("localhost:" + inputDirectory + inputFileName)
      #Load settings
      visSettings.load_visit_settings()
      vis.AddPlot("Pseudocolor", _variableName, 1, 1) #CONTINUE
      vis.SetActivePlots(1)
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
      vis.PseudocolorAtts.colorTableName = colorTableName
      vis.PseudocolorAtts.invertColorTable = 0
      vis.PseudocolorAtts.smoothingLevel = 0
      vis.PseudocolorAtts.pointSizeVarEnabled = 0
      vis.PseudocolorAtts.pointSizeVar = "default"
      vis.PseudocolorAtts.pointSizePixels = 2
      vis.PseudocolorAtts.lineStyle = vis.PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
      vis.PseudocolorAtts.lineWidth = 0
      vis.PseudocolorAtts.opacityType = vis.PseudocolorAtts.Explicit  # Explicit, ColorTable
      vis.SetPlotOptions(vis.PseudocolorAtts)
   
      
      vis.SetActivePlots(1)
      vis.AddOperator("Slice", 1)
      vis.SetActivePlots(1)
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
      vis.SetOperatorOptions(vis.SliceAtts, 1)
      vis.DrawPlots()
      vis.SetActivePlots(0)

      for coordinate in coordinates[i]:
         print str(coordinate)
         create_point_vtk( fileName=inputFileName2, coordinates=coordinate )
         vis.OpenDatabase(databaseName2, 0)
         vis.AddPlot("Mesh", "mesh", 1, 1)
         vis.SetActivePlots(vis.GetNumPlots())
         vis.MeshAtts = vis.MeshAttributes()
         vis.MeshAtts.legendFlag = 1
         vis.MeshAtts.lineStyle = vis.MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
         vis.MeshAtts.lineWidth = 0
         vis.MeshAtts.meshColor = (0, 0, 0, 255)
         vis.MeshAtts.outlineOnlyFlag = 0
         vis.MeshAtts.errorTolerance = 0.01
         vis.MeshAtts.meshColorSource = vis.MeshAtts.Foreground  # Foreground, MeshCustom
         vis.MeshAtts.opaqueColorSource = vis.MeshAtts.Background  # Background, OpaqueCustom
         vis.MeshAtts.opaqueMode = vis.MeshAtts.Auto  # Auto, On, Off
         vis.MeshAtts.pointSize = 0.05
         vis.MeshAtts.opaqueColor = (255, 255, 255, 255)
         vis.MeshAtts.smoothingLevel = vis.MeshAtts.None  # None, Fast, High
         vis.MeshAtts.pointSizeVarEnabled = 0
         vis.MeshAtts.pointSizeVar = "default"
         vis.MeshAtts.pointType = vis.MeshAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
         vis.MeshAtts.showInternal = 0
         vis.MeshAtts.pointSizePixels = 10
         vis.MeshAtts.opacity = 1
         vis.SetPlotOptions(vis.MeshAtts)
         vis.DrawPlots()
         # Iterate through frames
         vis.SaveWindowAtts = vis.SaveWindowAttributes()
         vis.SaveWindowAtts.outputToCurrentDirectory = 0
         vis.SaveWindowAtts.outputDirectory = _outputDir
         vis.SaveWindowAtts.fileName = _outputFileName
         vis.SaveWindowAtts.family = 1
         vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
         vis.SaveWindowAtts.width = 3000
         vis.SaveWindowAtts.height = 3000
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
      vis.DeleteActivePlots()
      vis.CloseDatabase(databaseName)
   vis.CloseDatabase(databaseName2)
   # Make the movie:
   #subprocess.call("./moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   pyVisitPath = "pyVisit/"
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh")
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   framerate = "10"
   subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", _outputDir, _outputFileName, framerate])
   # Delete the point vtk file:
   os.remove(os.getcwd() + "/" + inputFileName2)


