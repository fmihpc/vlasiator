#!/usr/bin/python -i

import visit as vis
import loadvisitsettings as visSettings
import subprocess
from useroptions import *

#cellids = [162870, 163500, 164115, 162885, 163515, 164130, 162900, 163530, 164145, 163440, 164055, 164685, 163455, 164070, 164700, 163470, 164085, 163485, 164100]

def make_distribution_movie(cellids, rotated, inputDirectory, outputDirectory, outputFileName, zoom=1.0, viewNormal=[0.488281, 0.382966, -0.784167], minThreshold=1e-18, maxThreshold=1e37):
   '''Makes a distribution movie of some given distribution data
      Example usage:
      make_distribution_movie(cellids=[18302, 19432, 19042], rotated=True, inputDirectory=\"/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/\", outputDirectory=\"/home/hannukse/MOVIES/\", outputFileName=\"testmovie\", zoom=0.8, viewNormal=[0.488281, 0.382966, -0.784167], minThreshold=1e-17, maxThreshold=1.2e37)
      Note: viewNormal determines the angle of view (straight from visit)
   '''
   if len(viewNormal) != 3:
      print "ERROR, INVALID VIEWNORMAL LENGTH, SHOULD BE 3"
      return
   for cell in sorted(cellids):
      # OPTIONS
      ###########################################################
      cellid = str(cell)
      #databaseName = "localhost:/home/hannukse/meteo/lustre/tmp/hannuksela/AAM/velgrid.rotated." + cellid + ".*.silo database"
      if rotated == True:
         rotateFix = "rotated."
      else:
         rotateFix = ""
      inputFileName = "velgrid." + rotateFix + cellid + ".*.silo"
      databaseName = "localhost:" + inputDirectory + inputFileName + " database"
      outputDir = outputDirectory
      fileName = outputFileName + "_" + cellid + "_"
      WIDTH = 3000
      HEIGHT = 3000
      # Threshold values:
      # TODO: USE VLSV READER TO AUTOMATE THIS
      minimumThreshold = minThreshold
      maximumThreshold = maxThreshold
      ###########################################################
      
      
      vis.OpenDatabase(databaseName, 0)
      #Load settings
      visSettings.load_visit_settings()
      #Make a plot
      vis.AddPlot("Pseudocolor", "avgs", 1, 1)
      vis.SetActivePlots(0)
      vis.AddOperator("Threshold", 1)
      vis.ThresholdAtts = vis.ThresholdAttributes()
      vis.ThresholdAtts.outputMeshType = 0
      vis.ThresholdAtts.listedVarNames = ("default")
      vis.ThresholdAtts.zonePortions = (1)
      vis.ThresholdAtts.lowerBounds = (minimumThreshold)
      vis.ThresholdAtts.upperBounds = (maximumThreshold)
      vis.ThresholdAtts.defaultVarName = "avgs"
      vis.ThresholdAtts.defaultVarIsScalar = 1
      vis.SetOperatorOptions(vis.ThresholdAtts, 1)
      vis.DrawPlots()
      # Begin spontaneous state
      vis.View3DAtts = vis.View3DAttributes()
      vis.View3DAtts.viewNormal = (viewNormal[0], viewNormal[1], viewNormal[2])
      vis.View3DAtts.focus = (-634.56, 91.3781, -13.7891)
      vis.View3DAtts.viewUp = (-0.102795, 0.917551, 0.3841)
      vis.View3DAtts.viewAngle = 30
      vis.View3DAtts.parallelScale = 1.45614e+06
      vis.View3DAtts.nearPlane = -2.91228e+06
      vis.View3DAtts.farPlane = 2.91228e+06
      vis.View3DAtts.imagePan = (0, 0)
      vis.View3DAtts.imageZoom = zoom
      vis.View3DAtts.perspective = 1
      vis.View3DAtts.eyeAngle = 2
      vis.View3DAtts.centerOfRotationSet = 0
      vis.View3DAtts.centerOfRotation = (-634.56, 91.3781, -13.7891)
      vis.View3DAtts.axis3DScaleFlag = 0
      vis.View3DAtts.axis3DScales = (1, 1, 1)
      vis.View3DAtts.shear = (0, 0, 1)
      vis.SetView3D(vis.View3DAtts)
      # End spontaneous state
      vis.ViewCurveAtts = vis.ViewCurveAttributes()
      vis.ViewCurveAtts.domainCoords = (0, 1)
      vis.ViewCurveAtts.rangeCoords = (0, 1)
      vis.ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
      vis.ViewCurveAtts.domainScale = vis.ViewCurveAtts.LINEAR  # LINEAR, LOG
      vis.ViewCurveAtts.rangeScale = vis.ViewCurveAtts.LINEAR  # LINEAR, LOG
      vis.SetViewCurve(vis.ViewCurveAtts)
      vis.View2DAtts = vis.View2DAttributes()
      vis.View2DAtts.windowCoords = (0, 1, 0, 1)
      vis.View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
      vis.View2DAtts.fullFrameActivationMode = vis.View2DAtts.Auto  # On, Off, Auto
      vis.View2DAtts.fullFrameAutoThreshold = 100
      vis.View2DAtts.xScale = vis.View2DAtts.LINEAR  # LINEAR, LOG
      vis.View2DAtts.yScale = vis.View2DAtts.LINEAR  # LINEAR, LOG
      vis.View2DAtts.windowValid = 0
      vis.SetView2D(vis.View2DAtts)
      vis.View3DAtts = vis.View3DAttributes()
      vis.View3DAtts.viewNormal = (viewNormal[0], viewNormal[1], viewNormal[2])
      vis.View3DAtts.focus = (-634.56, 91.3781, -13.7891)
      vis.View3DAtts.viewUp = (-0.102795, 0.917551, 0.3841)
      vis.View3DAtts.viewAngle = 30
      vis.View3DAtts.parallelScale = 1.45614e+06
      vis.View3DAtts.nearPlane = -2.91228e+06
      vis.View3DAtts.farPlane = 2.91228e+06
      vis.View3DAtts.imagePan = (0, 0)
      vis.View3DAtts.imageZoom = zoom
      vis.View3DAtts.perspective = 1
      vis.View3DAtts.eyeAngle = 2
      vis.View3DAtts.centerOfRotationSet = 0
      vis.View3DAtts.centerOfRotation = (-634.56, 91.3781, -13.7891)
      vis.View3DAtts.axis3DScaleFlag = 0
      vis.View3DAtts.axis3DScales = (1, 1, 1)
      vis.View3DAtts.shear = (0, 0, 1)
      vis.SetView3D(vis.View3DAtts)
      vis.ViewAxisArrayAtts = vis.ViewAxisArrayAttributes()
      vis.ViewAxisArrayAtts.domainCoords = (0, 1)
      vis.ViewAxisArrayAtts.rangeCoords = (0, 1)
      vis.ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
      vis.SetViewAxisArray(vis.ViewAxisArrayAtts)
      for i in range(0, vis.GetDatabaseNStates()):
         vis.SetTimeSliderState(i)
         vis.SaveWindowAtts = vis.SaveWindowAttributes()
         vis.SaveWindowAtts.outputToCurrentDirectory = 0
         vis.SaveWindowAtts.outputDirectory = outputDir
         vis.SaveWindowAtts.fileName = fileName
         vis.SaveWindowAtts.family = 1
         vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
         vis.SaveWindowAtts.width = WIDTH
         vis.SaveWindowAtts.height = HEIGHT
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
      framerate = 5
      subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", outputDir, fileName, framerate])
