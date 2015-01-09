import sys
import subprocess

time = sys.argv[1]

xCells = 666
yCells = 583

xStep = 50
yStep = 50

xList = range(xStep+1, xCells, xStep)
yList = range(yStep+1, yCells, yStep)

cellList = []

for x in xList:
   for y in yList:
      cellList.append((y-1)*xCells + x)

print len(cellList)

for n in cellList:
   number = str(n)
   siloFile = "velgrid."+number+"."+time+".silo"
   fileToExtract = "bulk."+time+".vlsv"
   subprocess.call(["aprun", "/home/users/alfthan/opt/bin/vlsvextract_DP", fileToExtract, "--cellid", number])
   
   if os.path.isfile(siloFile) == False:
      subprocess.call(["ln", "-s", "velgrid.1.0000000.silo", siloFile])
   
   subprocess.call(["chmod", "a+r", siloFile])
   
   fileLocation = "/stornext/field/vlasiator/visualizations/ACB/mosaics_tmp3/velgrid."+number+"."+time+".silo"
   outputLocation = "/stornext/field/vlasiator/visualizations/ACB/mosaics_tmp3/"
   fileName = "ACB_3dv_"+time+"_"+number.rjust(7, '0')

   OpenDatabase(fileLocation, 0)

   AddPlot("Contour", "avgs", 1, 1)
   ContourAtts = ContourAttributes()
   ContourAtts.defaultPalette.GetControlPoints(0).colors = (255, 0, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(0).position = 0
   ContourAtts.defaultPalette.GetControlPoints(1).colors = (0, 255, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(1).position = 0.034
   ContourAtts.defaultPalette.GetControlPoints(2).colors = (0, 0, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(2).position = 0.069
   ContourAtts.defaultPalette.GetControlPoints(3).colors = (0, 255, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(3).position = 0.103
   ContourAtts.defaultPalette.GetControlPoints(4).colors = (255, 0, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(4).position = 0.138
   ContourAtts.defaultPalette.GetControlPoints(5).colors = (255, 255, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(5).position = 0.172
   ContourAtts.defaultPalette.GetControlPoints(6).colors = (255, 135, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(6).position = 0.207
   ContourAtts.defaultPalette.GetControlPoints(7).colors = (255, 0, 135, 255)
   ContourAtts.defaultPalette.GetControlPoints(7).position = 0.241
   ContourAtts.defaultPalette.GetControlPoints(8).colors = (168, 168, 168, 255)
   ContourAtts.defaultPalette.GetControlPoints(8).position = 0.276
   ContourAtts.defaultPalette.GetControlPoints(9).colors = (255, 68, 68, 255)
   ContourAtts.defaultPalette.GetControlPoints(9).position = 0.31
   ContourAtts.defaultPalette.GetControlPoints(10).colors = (99, 255, 99, 255)
   ContourAtts.defaultPalette.GetControlPoints(10).position = 0.345
   ContourAtts.defaultPalette.GetControlPoints(11).colors = (99, 99, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(11).position = 0.379
   ContourAtts.defaultPalette.GetControlPoints(12).colors = (40, 165, 165, 255)
   ContourAtts.defaultPalette.GetControlPoints(12).position = 0.414
   ContourAtts.defaultPalette.GetControlPoints(13).colors = (255, 99, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(13).position = 0.448
   ContourAtts.defaultPalette.GetControlPoints(14).colors = (255, 255, 99, 255)
   ContourAtts.defaultPalette.GetControlPoints(14).position = 0.483
   ContourAtts.defaultPalette.GetControlPoints(15).colors = (255, 170, 99, 255)
   ContourAtts.defaultPalette.GetControlPoints(15).position = 0.517
   ContourAtts.defaultPalette.GetControlPoints(16).colors = (170, 79, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(16).position = 0.552
   ContourAtts.defaultPalette.GetControlPoints(17).colors = (150, 0, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(17).position = 0.586
   ContourAtts.defaultPalette.GetControlPoints(18).colors = (0, 150, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(18).position = 0.621
   ContourAtts.defaultPalette.GetControlPoints(19).colors = (0, 0, 150, 255)
   ContourAtts.defaultPalette.GetControlPoints(19).position = 0.655
   ContourAtts.defaultPalette.GetControlPoints(20).colors = (0, 109, 109, 255)
   ContourAtts.defaultPalette.GetControlPoints(20).position = 0.69
   ContourAtts.defaultPalette.GetControlPoints(21).colors = (150, 0, 150, 255)
   ContourAtts.defaultPalette.GetControlPoints(21).position = 0.724
   ContourAtts.defaultPalette.GetControlPoints(22).colors = (150, 150, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(22).position = 0.759
   ContourAtts.defaultPalette.GetControlPoints(23).colors = (150, 84, 0, 255)
   ContourAtts.defaultPalette.GetControlPoints(23).position = 0.793
   ContourAtts.defaultPalette.GetControlPoints(24).colors = (160, 0, 79, 255)
   ContourAtts.defaultPalette.GetControlPoints(24).position = 0.828
   ContourAtts.defaultPalette.GetControlPoints(25).colors = (255, 104, 28, 255)
   ContourAtts.defaultPalette.GetControlPoints(25).position = 0.862
   ContourAtts.defaultPalette.GetControlPoints(26).colors = (0, 170, 81, 255)
   ContourAtts.defaultPalette.GetControlPoints(26).position = 0.897
   ContourAtts.defaultPalette.GetControlPoints(27).colors = (68, 255, 124, 255)
   ContourAtts.defaultPalette.GetControlPoints(27).position = 0.931
   ContourAtts.defaultPalette.GetControlPoints(28).colors = (0, 130, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(28).position = 0.966
   ContourAtts.defaultPalette.GetControlPoints(29).colors = (130, 0, 255, 255)
   ContourAtts.defaultPalette.GetControlPoints(29).position = 1
   ContourAtts.defaultPalette.smoothing = ContourAtts.defaultPalette.None  # None, Linear, CubicSpline
   ContourAtts.defaultPalette.equalSpacingFlag = 1
   ContourAtts.defaultPalette.discreteFlag = 1
   #ContourAtts.defaultPalette.externalFlag = 0
   ContourAtts.changedColors = (0)
   ContourAtts.colorType = ContourAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
   ContourAtts.colorTableName = "Default"
   ContourAtts.invertColorTable = 0
   ContourAtts.legendFlag = 1
   ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   ContourAtts.lineWidth = 0
   ContourAtts.singleColor = (255, 0, 0, 255)
   ContourAtts.SetMultiColor(0, (255, 0, 0, 130))
   ContourAtts.SetMultiColor(1, (0, 255, 0, 255))
   ContourAtts.SetMultiColor(2, (0, 0, 255, 255))
   ContourAtts.SetMultiColor(3, (0, 255, 255, 255))
   ContourAtts.SetMultiColor(4, (255, 0, 255, 255))
   ContourAtts.SetMultiColor(5, (255, 255, 0, 255))
   ContourAtts.SetMultiColor(6, (255, 135, 0, 255))
   ContourAtts.SetMultiColor(7, (255, 0, 135, 255))
   ContourAtts.SetMultiColor(8, (168, 168, 168, 255))
   ContourAtts.SetMultiColor(9, (255, 68, 68, 255))
   ContourAtts.contourNLevels = 10
   ContourAtts.contourValue = (1e-15)
   ContourAtts.contourPercent = ()
   ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
   ContourAtts.minFlag = 0
   ContourAtts.maxFlag = 0
   ContourAtts.min = 0
   ContourAtts.max = 1
   ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
   ContourAtts.wireframe = 0
   SetPlotOptions(ContourAtts)
   DrawPlots()


   View3DAtts = View3DAttributes()
   View3DAtts.viewNormal = (0, -1, 1)
   View3DAtts.focus = (0, 0, 0)
   View3DAtts.viewUp = (0, 0, 1)
   View3DAtts.viewAngle = 30
   View3DAtts.parallelScale = 1.63228e+06
   View3DAtts.nearPlane = -3.54897e+06
   View3DAtts.farPlane = 3.54897e+06
   View3DAtts.imagePan = (0, 0)
   View3DAtts.imageZoom = 1
   View3DAtts.perspective = 1
   View3DAtts.eyeAngle = 2
   View3DAtts.centerOfRotationSet = 0
   View3DAtts.centerOfRotation = (0, 0, 0)
   View3DAtts.axis3DScaleFlag = 0
   View3DAtts.axis3DScales = (1, 1, 1)
   View3DAtts.shear = (0, 0, 1)
   SetView3D(View3DAtts)

   ViewAxisArrayAtts = ViewAxisArrayAttributes()
   ViewAxisArrayAtts.domainCoords = (0, 1)
   ViewAxisArrayAtts.rangeCoords = (0, 1)
   ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
   SetViewAxisArray(ViewAxisArrayAtts)

   AnnotationAtts = AnnotationAttributes()
   AnnotationAtts.axes3D.visible = 1
   #AnnotationAtts.axes3D.autoSetTicks = 1
   #AnnotationAtts.axes3D.autoSetScaling = 1
   AnnotationAtts.axes3D.lineWidth = 5
   AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside    # Inside, Outside, Both
   AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.OutsideEdges  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
   AnnotationAtts.axes3D.triadFlag = 0
   AnnotationAtts.axes3D.bboxFlag = 1
   AnnotationAtts.axes3D.xAxis.title.visible = 0
   #AnnotationAtts.axes3D.xAxis.title.font.font = 
   #AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.xAxis.title.font.scale = 2
   #AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.xAxis.title.font.bold = 0
   #AnnotationAtts.axes3D.xAxis.title.font.italic = 0
   #AnnotationAtts.axes3D.xAxis.title.userTitle = 0
   #AnnotationAtts.axes3D.xAxis.title.userUnits = 0
   #AnnotationAtts.axes3D.xAxis.title.title = " "
   #AnnotationAtts.axes3D.xAxis.title.units = ""
   AnnotationAtts.axes3D.xAxis.label.visible = 0
   #AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.xAxis.label.font.scale = 3
   #AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.xAxis.label.font.bold = 0
   #AnnotationAtts.axes3D.xAxis.label.font.italic = 0
   #AnnotationAtts.axes3D.xAxis.label.scaling = 0
   AnnotationAtts.axes3D.xAxis.tickMarks.visible = 0
   #AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
   #AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
   #AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
   #AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
   #AnnotationAtts.axes3D.xAxis.grid = 0
   AnnotationAtts.axes3D.yAxis.title.visible = 0
   #AnnotationAtts.axes3D.yAxis.title.font.font = 
   #AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.yAxis.title.font.scale = 1
   #AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.yAxis.title.font.bold = 0
   #AnnotationAtts.axes3D.yAxis.title.font.italic = 0
   #AnnotationAtts.axes3D.yAxis.title.userTitle = 0
   #AnnotationAtts.axes3D.yAxis.title.userUnits = 0
   #AnnotationAtts.axes3D.yAxis.title.title = " "
   #AnnotationAtts.axes3D.yAxis.title.units = ""
   AnnotationAtts.axes3D.yAxis.label.visible = 0
   #AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.yAxis.label.font.scale = 3
   #AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.yAxis.label.font.bold = 0
   #AnnotationAtts.axes3D.yAxis.label.font.italic = 0
   #AnnotationAtts.axes3D.yAxis.label.scaling = 0
   AnnotationAtts.axes3D.yAxis.tickMarks.visible = 0
   #AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
   #AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
   #AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
   #AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
   #AnnotationAtts.axes3D.yAxis.grid = 0
   AnnotationAtts.axes3D.zAxis.title.visible = 0
   #AnnotationAtts.axes3D.zAxis.title.font.font = 
   #AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.zAxis.title.font.scale = 1
   #AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.zAxis.title.font.bold = 0
   #AnnotationAtts.axes3D.zAxis.title.font.italic = 0
   #AnnotationAtts.axes3D.zAxis.title.userTitle = 0
   #AnnotationAtts.axes3D.zAxis.title.userUnits = 0
   #AnnotationAtts.axes3D.zAxis.title.title = " "
   #AnnotationAtts.axes3D.zAxis.title.units = ""
   AnnotationAtts.axes3D.zAxis.label.visible = 0
   #AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
   #AnnotationAtts.axes3D.zAxis.label.font.scale = 3
   #AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
   #AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
   #AnnotationAtts.axes3D.zAxis.label.font.bold = 0
   #AnnotationAtts.axes3D.zAxis.label.font.italic = 0
   #AnnotationAtts.axes3D.zAxis.label.scaling = 0
   AnnotationAtts.axes3D.zAxis.tickMarks.visible = 0
   #AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
   #AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
   #AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
   #AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
   #AnnotationAtts.axes3D.zAxis.grid = 0
   AnnotationAtts.axes3D.setBBoxLocation = 1
   AnnotationAtts.axes3D.bboxLocation = (-1.0e6, 1.0e6, -1.0e6, 1.0e6, -1.0e6, 1.0e6)

   AnnotationAtts.userInfoFlag = 0
   AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.userInfoFont.scale = 1
   AnnotationAtts.userInfoFont.useForegroundColor = 1
   AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
   AnnotationAtts.userInfoFont.bold = 0
   AnnotationAtts.userInfoFont.italic = 0
   AnnotationAtts.databaseInfoFlag = 0
   AnnotationAtts.timeInfoFlag = 1
   AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.databaseInfoFont.scale = 1
   AnnotationAtts.databaseInfoFont.useForegroundColor = 1
   AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
   AnnotationAtts.databaseInfoFont.bold = 0
   AnnotationAtts.databaseInfoFont.italic = 0
   AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
   AnnotationAtts.databaseInfoTimeScale = 1
   AnnotationAtts.databaseInfoTimeOffset = 0
   AnnotationAtts.legendInfoFlag = 0
   SetAnnotationAttributes(AnnotationAtts)



   # MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
   #SaveSession("/zhome/academic/HLRS/pri/iprykemp/.visit/crash_recovery.session")
   # MAINTENANCE ISSUE: SetSuppressMessagesRPC is not handled in Logging.C. Please contact a VisIt developer.
   SaveWindowAtts = SaveWindowAttributes()
   SaveWindowAtts.outputToCurrentDirectory = 0
   SaveWindowAtts.outputDirectory = outputLocation
   SaveWindowAtts.fileName = fileName
   SaveWindowAtts.family = 0
   SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
   SaveWindowAtts.width = 1024
   SaveWindowAtts.height = 1024
   SaveWindowAtts.screenCapture = 0
   SaveWindowAtts.saveTiled = 0
   SaveWindowAtts.quality = 80
   SaveWindowAtts.progressive = 0
   SaveWindowAtts.binary = 0
   SaveWindowAtts.stereo = 0
   SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
   SaveWindowAtts.forceMerge = 0
   SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
   SaveWindowAtts.advancedMultiWindowSave = 0
   SetSaveWindowAttributes(SaveWindowAtts)
   SaveWindow()
   
   DeleteAllPlots()
   CloseDatabase(fileLocation)
   
   subprocess.call(["mv", fileName+".png", str(time).rjust(7, '0')+"/"])
   subprocess.call(["rm", siloFile])
   



exit()
