# In order to use this script, you need to at least edit the lines with the following comment:
# EDIT

# For more information on the usage of this script:
# https://github.com/fmihpc/vlasiator/wiki/Basic-movie-making-script-documentation,-using-VisIt-on-voima

import sys

fileName = "BCB_rho_" # EDIT
outputLocation = "/lustre/tmp/alfthan/2D/BCB/visualizations/movies/rho/" # EDIT


frameList = range(int(sys.argv[1]),int(sys.argv[2])+1)

for entry in frameList:
   number=str(entry).rjust(7, '0')
   fileLocation="/lustre/tmp/alfthan/2D/BCB/bulk."+number+".vlsv" # EDIT
   OpenDatabase(fileLocation)
   
   AddPlot("Pseudocolor", "rho", 0, 1) # EDIT
   AddOperator("Slice")
   AddOperator("Threshold")
   PseudocolorAtts = PseudocolorAttributes()
   PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
   PseudocolorAtts.skewFactor = 1
   PseudocolorAtts.limitsMode = PseudocolorAtts.CurrentPlot  # OriginalData, CurrentPlot
   PseudocolorAtts.minFlag = 1
   PseudocolorAtts.min = 100000 # EDIT
   PseudocolorAtts.maxFlag = 1
   PseudocolorAtts.max = 1e+07 # EDIT
   PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
   PseudocolorAtts.colorTableName = "hot_desaturated"
   PseudocolorAtts.invertColorTable = 0
   PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
   PseudocolorAtts.opacityVariable = ""
   PseudocolorAtts.opacity = 1
   PseudocolorAtts.opacityVarMin = 0
   PseudocolorAtts.opacityVarMax = 1
   PseudocolorAtts.opacityVarMinFlag = 0
   PseudocolorAtts.opacityVarMaxFlag = 0
   PseudocolorAtts.pointSize = 0.05
   PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
   PseudocolorAtts.pointSizeVarEnabled = 0
   PseudocolorAtts.pointSizeVar = "default"
   PseudocolorAtts.pointSizePixels = 2
   PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
   PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   PseudocolorAtts.lineWidth = 0
   PseudocolorAtts.tubeDisplayDensity = 10
   PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
   PseudocolorAtts.tubeRadiusAbsolute = 0.125
   PseudocolorAtts.tubeRadiusBBox = 0.005
   PseudocolorAtts.varyTubeRadius = 0
   PseudocolorAtts.varyTubeRadiusVariable = ""
   PseudocolorAtts.varyTubeRadiusFactor = 10
   PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
   PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
   PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
   PseudocolorAtts.endPointRadiusAbsolute = 1
   PseudocolorAtts.endPointRadiusBBox = 0.005
   PseudocolorAtts.endPointRatio = 2
   PseudocolorAtts.renderSurfaces = 1
   PseudocolorAtts.renderWireframe = 0
   PseudocolorAtts.renderPoints = 0
   PseudocolorAtts.smoothingLevel = 0
   PseudocolorAtts.legendFlag = 0
   PseudocolorAtts.lightingFlag = 1
   SetPlotOptions(PseudocolorAtts)

   SliceAtts = SliceAttributes()
   SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
   SliceAtts.originPoint = (0, 0, 0)
   SliceAtts.originIntercept = 0
   SliceAtts.originPercent = 0
   SliceAtts.originZone = 0
   SliceAtts.originNode = 0
   SliceAtts.normal = (0, 1, 0)
   SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
   SliceAtts.upAxis = (0, 0, 1)
   SliceAtts.project2d = 1
   SliceAtts.interactive = 1
   SliceAtts.flip = 0
   SliceAtts.originZoneDomain = 0
   SliceAtts.originNodeDomain = 0
   SliceAtts.meshName = "SpatialGrid"
   SliceAtts.theta = 0
   SliceAtts.phi = 90
   SetOperatorOptions(SliceAtts)

   ThresholdAtts = ThresholdAttributes()
   ThresholdAtts.outputMeshType = 0
   ThresholdAtts.listedVarNames = ("Boundary_type")
   ThresholdAtts.zonePortions = (1)
   ThresholdAtts.lowerBounds = (1)
   ThresholdAtts.upperBounds = (1)
   ThresholdAtts.defaultVarName = "rho"
   ThresholdAtts.defaultVarIsScalar = 1
   SetOperatorOptions(ThresholdAtts)

   # Logging for SetAnnotationObjectOptions is not implemented yet.
   AnnotationAtts = AnnotationAttributes()
   AnnotationAtts.axes2D.visible = 1
   AnnotationAtts.axes2D.autoSetTicks = 1
   AnnotationAtts.axes2D.autoSetScaling = 1
   AnnotationAtts.axes2D.lineWidth = 0
   AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
   AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
   AnnotationAtts.axes2D.xAxis.title.visible = 1
   AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
   AnnotationAtts.axes2D.xAxis.title.font.scale = 1
   AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
   AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes2D.xAxis.title.font.bold = 1
   AnnotationAtts.axes2D.xAxis.title.font.italic = 1
   AnnotationAtts.axes2D.xAxis.title.userTitle = 0
   AnnotationAtts.axes2D.xAxis.title.userUnits = 0
   AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
   AnnotationAtts.axes2D.xAxis.title.units = ""
   AnnotationAtts.axes2D.xAxis.label.visible = 1
   AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
   AnnotationAtts.axes2D.xAxis.label.font.scale = 1
   AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
   AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes2D.xAxis.label.font.bold = 1
   AnnotationAtts.axes2D.xAxis.label.font.italic = 1
   AnnotationAtts.axes2D.xAxis.label.scaling = 0
   AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
   AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
   AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
   AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
   AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
   AnnotationAtts.axes2D.xAxis.grid = 0
   AnnotationAtts.axes2D.yAxis.title.visible = 1
   AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
   AnnotationAtts.axes2D.yAxis.title.font.scale = 1
   AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
   AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes2D.yAxis.title.font.bold = 1
   AnnotationAtts.axes2D.yAxis.title.font.italic = 1
   AnnotationAtts.axes2D.yAxis.title.userTitle = 0
   AnnotationAtts.axes2D.yAxis.title.userUnits = 0
   AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
   AnnotationAtts.axes2D.yAxis.title.units = ""
   AnnotationAtts.axes2D.yAxis.label.visible = 1
   AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
   AnnotationAtts.axes2D.yAxis.label.font.scale = 1
   AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
   AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axes2D.yAxis.label.font.bold = 1
   AnnotationAtts.axes2D.yAxis.label.font.italic = 1
   AnnotationAtts.axes2D.yAxis.label.scaling = 0
   AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
   AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
   AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
   AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
   AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
   AnnotationAtts.axes2D.yAxis.grid = 0

   AnnotationAtts.userInfoFlag = 0
   AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.userInfoFont.scale = 1
   AnnotationAtts.userInfoFont.useForegroundColor = 1
   AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
   AnnotationAtts.userInfoFont.bold = 0
   AnnotationAtts.userInfoFont.italic = 0
   AnnotationAtts.databaseInfoFlag = 1
   AnnotationAtts.timeInfoFlag = 0
   AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
   AnnotationAtts.databaseInfoFont.scale = 1
   AnnotationAtts.databaseInfoFont.useForegroundColor = 0
   AnnotationAtts.databaseInfoFont.color = (255, 255, 255, 255)
   AnnotationAtts.databaseInfoFont.bold = 0
   AnnotationAtts.databaseInfoFont.italic = 0
   AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
   AnnotationAtts.databaseInfoTimeScale = 1
   AnnotationAtts.databaseInfoTimeOffset = 0
   AnnotationAtts.legendInfoFlag = 0
   AnnotationAtts.backgroundColor = (255, 255, 255, 255)
   AnnotationAtts.foregroundColor = (0, 0, 0, 255)
   AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
   AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
   AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
   AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
   AnnotationAtts.backgroundImage = ""
   AnnotationAtts.imageRepeatX = 1
   AnnotationAtts.imageRepeatY = 1

   AnnotationAtts.axesArray.visible = 1
   AnnotationAtts.axesArray.ticksVisible = 1
   AnnotationAtts.axesArray.autoSetTicks = 1
   AnnotationAtts.axesArray.autoSetScaling = 1
   AnnotationAtts.axesArray.lineWidth = 0
   AnnotationAtts.axesArray.axes.title.visible = 1
   AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
   AnnotationAtts.axesArray.axes.title.font.scale = 1
   AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
   AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
   AnnotationAtts.axesArray.axes.title.font.bold = 0
   AnnotationAtts.axesArray.axes.title.font.italic = 0
   AnnotationAtts.axesArray.axes.title.userTitle = 0
   AnnotationAtts.axesArray.axes.title.userUnits = 0
   AnnotationAtts.axesArray.axes.title.title = ""
   AnnotationAtts.axesArray.axes.title.units = ""
   AnnotationAtts.axesArray.axes.label.visible = 1
   AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
   AnnotationAtts.axesArray.axes.label.font.scale = 1
   AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
   AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
   AnnotationAtts.axesArray.axes.label.font.bold = 0
   AnnotationAtts.axesArray.axes.label.font.italic = 0
   AnnotationAtts.axesArray.axes.label.scaling = 0
   AnnotationAtts.axesArray.axes.tickMarks.visible = 1
   AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
   AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
   AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
   AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
   AnnotationAtts.axesArray.axes.grid = 0
   SetAnnotationAttributes(AnnotationAtts)

   View2DAtts = View2DAttributes()
   View2DAtts.windowCoords = (-3e+08, 3e+08, -3e+08, 3e+08)
   View2DAtts.viewportCoords = (0, 1, 0, 1)
   View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
   View2DAtts.fullFrameAutoThreshold = 100
   View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
   View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
   View2DAtts.windowValid = 1
   SetView2D(View2DAtts)
      
   DrawPlots()

   SaveWindowAtts = SaveWindowAttributes()
   SaveWindowAtts.outputToCurrentDirectory = 0
   SaveWindowAtts.outputDirectory = outputLocation
   SaveWindowAtts.fileName = fileName+number
   SaveWindowAtts.family = 0
   SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
   SaveWindowAtts.width = 2000
   SaveWindowAtts.height = 2000
   SaveWindowAtts.screenCapture = 0
   SaveWindowAtts.saveTiled = 0
   SaveWindowAtts.quality = 80
   SaveWindowAtts.progressive = 0
   SaveWindowAtts.binary = 0
   SaveWindowAtts.stereo = 0
   SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
   SaveWindowAtts.forceMerge = 0
   SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
   SaveWindowAtts.advancedMultiWindowSave = 0
   SetSaveWindowAttributes(SaveWindowAtts)
   SaveWindow()
   
   DeleteAllPlots()
   #ResetView()
   CloseDatabase(fileLocation)

exit()
