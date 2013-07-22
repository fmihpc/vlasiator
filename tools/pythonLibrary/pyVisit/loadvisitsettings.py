#!/usr/bin/python -i

import visit as vis

def load_visit_settings():
   # Logging for SetAnnotationObjectOptions is not implemented yet.
   vis.AnnotationAtts = vis.AnnotationAttributes()
   vis.AnnotationAtts.axes2D.visible = 1
   vis.AnnotationAtts.axes2D.autoSetTicks = 1
   vis.AnnotationAtts.axes2D.autoSetScaling = 1
   vis.AnnotationAtts.axes2D.lineWidth = 0
   vis.AnnotationAtts.axes2D.tickLocation = vis.AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
   vis.AnnotationAtts.axes2D.tickAxes = vis.AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
   vis.AnnotationAtts.axes2D.xAxis.title.visible = 1
   vis.AnnotationAtts.axes2D.xAxis.title.font.font = vis.AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
   vis.AnnotationAtts.axes2D.xAxis.title.font.scale = 1
   vis.AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes2D.xAxis.title.font.bold = 1
   vis.AnnotationAtts.axes2D.xAxis.title.font.italic = 1
   vis.AnnotationAtts.axes2D.xAxis.title.userTitle = 0
   vis.AnnotationAtts.axes2D.xAxis.title.userUnits = 0
   vis.AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
   vis.AnnotationAtts.axes2D.xAxis.title.units = ""
   vis.AnnotationAtts.axes2D.xAxis.label.visible = 1
   vis.AnnotationAtts.axes2D.xAxis.label.font.font = vis.AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
   vis.AnnotationAtts.axes2D.xAxis.label.font.scale = 1
   vis.AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes2D.xAxis.label.font.bold = 1
   vis.AnnotationAtts.axes2D.xAxis.label.font.italic = 1
   vis.AnnotationAtts.axes2D.xAxis.label.scaling = 0
   vis.AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
   vis.AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axes2D.xAxis.grid = 0
   vis.AnnotationAtts.axes2D.yAxis.title.visible = 1
   vis.AnnotationAtts.axes2D.yAxis.title.font.font = vis.AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
   vis.AnnotationAtts.axes2D.yAxis.title.font.scale = 1
   vis.AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes2D.yAxis.title.font.bold = 1
   vis.AnnotationAtts.axes2D.yAxis.title.font.italic = 1
   vis.AnnotationAtts.axes2D.yAxis.title.userTitle = 0
   vis.AnnotationAtts.axes2D.yAxis.title.userUnits = 0
   vis.AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
   vis.AnnotationAtts.axes2D.yAxis.title.units = ""
   vis.AnnotationAtts.axes2D.yAxis.label.visible = 1
   vis.AnnotationAtts.axes2D.yAxis.label.font.font = vis.AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
   vis.AnnotationAtts.axes2D.yAxis.label.font.scale = 1
   vis.AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes2D.yAxis.label.font.bold = 1
   vis.AnnotationAtts.axes2D.yAxis.label.font.italic = 1
   vis.AnnotationAtts.axes2D.yAxis.label.scaling = 0
   vis.AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
   vis.AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axes2D.yAxis.grid = 0
   vis.AnnotationAtts.axes3D.visible = 1
   vis.AnnotationAtts.axes3D.autoSetTicks = 1
   vis.AnnotationAtts.axes3D.autoSetScaling = 1
   vis.AnnotationAtts.axes3D.lineWidth = 0
   vis.AnnotationAtts.axes3D.tickLocation = vis.AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
   vis.AnnotationAtts.axes3D.axesType = vis.AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
   vis.AnnotationAtts.axes3D.triadFlag = 1
   vis.AnnotationAtts.axes3D.bboxFlag = 1
   vis.AnnotationAtts.axes3D.xAxis.title.visible = 1
   vis.AnnotationAtts.axes3D.xAxis.title.font.font = vis.AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.xAxis.title.font.scale = 1
   vis.AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.xAxis.title.font.bold = 0
   vis.AnnotationAtts.axes3D.xAxis.title.font.italic = 0
   vis.AnnotationAtts.axes3D.xAxis.title.userTitle = 0
   vis.AnnotationAtts.axes3D.xAxis.title.userUnits = 0
   vis.AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
   vis.AnnotationAtts.axes3D.xAxis.title.units = ""
   vis.AnnotationAtts.axes3D.xAxis.label.visible = 1
   vis.AnnotationAtts.axes3D.xAxis.label.font.font = vis.AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.xAxis.label.font.scale = 1
   vis.AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.xAxis.label.font.bold = 0
   vis.AnnotationAtts.axes3D.xAxis.label.font.italic = 0
   vis.AnnotationAtts.axes3D.xAxis.label.scaling = 0
   vis.AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
   vis.AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axes3D.xAxis.grid = 0
   vis.AnnotationAtts.axes3D.yAxis.title.visible = 1
   vis.AnnotationAtts.axes3D.yAxis.title.font.font = vis.AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.yAxis.title.font.scale = 1
   vis.AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.yAxis.title.font.bold = 0
   vis.AnnotationAtts.axes3D.yAxis.title.font.italic = 0
   vis.AnnotationAtts.axes3D.yAxis.title.userTitle = 0
   vis.AnnotationAtts.axes3D.yAxis.title.userUnits = 0
   vis.AnnotationAtts.axes3D.yAxis.title.title = "Y-Axis"
   vis.AnnotationAtts.axes3D.yAxis.title.units = ""
   vis.AnnotationAtts.axes3D.yAxis.label.visible = 1
   vis.AnnotationAtts.axes3D.yAxis.label.font.font = vis.AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.yAxis.label.font.scale = 1
   vis.AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.yAxis.label.font.bold = 0
   vis.AnnotationAtts.axes3D.yAxis.label.font.italic = 0
   vis.AnnotationAtts.axes3D.yAxis.label.scaling = 0
   vis.AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
   vis.AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axes3D.yAxis.grid = 0
   vis.AnnotationAtts.axes3D.zAxis.title.visible = 1
   vis.AnnotationAtts.axes3D.zAxis.title.font.font = vis.AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.zAxis.title.font.scale = 1
   vis.AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.zAxis.title.font.bold = 0
   vis.AnnotationAtts.axes3D.zAxis.title.font.italic = 0
   vis.AnnotationAtts.axes3D.zAxis.title.userTitle = 0
   vis.AnnotationAtts.axes3D.zAxis.title.userUnits = 0
   vis.AnnotationAtts.axes3D.zAxis.title.title = "Z-Axis"
   vis.AnnotationAtts.axes3D.zAxis.title.units = ""
   vis.AnnotationAtts.axes3D.zAxis.label.visible = 1
   vis.AnnotationAtts.axes3D.zAxis.label.font.font = vis.AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axes3D.zAxis.label.font.scale = 1
   vis.AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axes3D.zAxis.label.font.bold = 0
   vis.AnnotationAtts.axes3D.zAxis.label.font.italic = 0
   vis.AnnotationAtts.axes3D.zAxis.label.scaling = 0
   vis.AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
   vis.AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axes3D.zAxis.grid = 0
   vis.AnnotationAtts.axes3D.setBBoxLocation = 0
   vis.AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
   vis.AnnotationAtts.userInfoFlag = 0
   vis.AnnotationAtts.userInfoFont.font = vis.AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.userInfoFont.scale = 1
   vis.AnnotationAtts.userInfoFont.useForegroundColor = 1
   vis.AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
   vis.AnnotationAtts.userInfoFont.bold = 0
   vis.AnnotationAtts.userInfoFont.italic = 0
   vis.AnnotationAtts.databaseInfoFlag = 1
   vis.AnnotationAtts.timeInfoFlag = 0
   vis.AnnotationAtts.databaseInfoFont.font = vis.AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.databaseInfoFont.scale = 0.4
   vis.AnnotationAtts.databaseInfoFont.useForegroundColor = 1
   vis.AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
   vis.AnnotationAtts.databaseInfoFont.bold = 0
   vis.AnnotationAtts.databaseInfoFont.italic = 0
   vis.AnnotationAtts.databaseInfoExpansionMode = vis.AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
   vis.AnnotationAtts.databaseInfoTimeScale = 1
   vis.AnnotationAtts.databaseInfoTimeOffset = 0
   vis.AnnotationAtts.legendInfoFlag = 1
   vis.AnnotationAtts.backgroundColor = (255, 255, 255, 255)
   vis.AnnotationAtts.foregroundColor = (0, 0, 0, 255)
   vis.AnnotationAtts.gradientBackgroundStyle = vis.AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
   vis.AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
   vis.AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
   vis.AnnotationAtts.backgroundMode = vis.AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
   vis.AnnotationAtts.backgroundImage = ""
   vis.AnnotationAtts.imageRepeatX = 1
   vis.AnnotationAtts.imageRepeatY = 1
   vis.AnnotationAtts.axesArray.visible = 1
   vis.AnnotationAtts.axesArray.ticksVisible = 1
   vis.AnnotationAtts.axesArray.autoSetTicks = 1
   vis.AnnotationAtts.axesArray.autoSetScaling = 1
   vis.AnnotationAtts.axesArray.lineWidth = 0
   vis.AnnotationAtts.axesArray.axes.title.visible = 1
   vis.AnnotationAtts.axesArray.axes.title.font.font = vis.AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axesArray.axes.title.font.scale = 1
   vis.AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
   vis.AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axesArray.axes.title.font.bold = 0
   vis.AnnotationAtts.axesArray.axes.title.font.italic = 0
   vis.AnnotationAtts.axesArray.axes.title.userTitle = 0
   vis.AnnotationAtts.axesArray.axes.title.userUnits = 0
   vis.AnnotationAtts.axesArray.axes.title.title = ""
   vis.AnnotationAtts.axesArray.axes.title.units = ""
   vis.AnnotationAtts.axesArray.axes.label.visible = 1
   vis.AnnotationAtts.axesArray.axes.label.font.font = vis.AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
   vis.AnnotationAtts.axesArray.axes.label.font.scale = 1
   vis.AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
   vis.AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
   vis.AnnotationAtts.axesArray.axes.label.font.bold = 0
   vis.AnnotationAtts.axesArray.axes.label.font.italic = 0
   vis.AnnotationAtts.axesArray.axes.label.scaling = 0
   vis.AnnotationAtts.axesArray.axes.tickMarks.visible = 1
   vis.AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
   vis.AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
   vis.AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
   vis.AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
   vis.AnnotationAtts.axesArray.axes.grid = 0
   vis.SetAnnotationAttributes(vis.AnnotationAtts)

   # Load expressions:
   vis.DefineScalarExpression("volume/Vx", "rho_v[0] / rho")
   vis.DefineScalarExpression("volume/Vy", "rho_v[1]/rho")
   vis.DefineScalarExpression("volume/Vz", "rho_v[2]/rho")
   vis.DefineScalarExpression("edge/Ex", "E[0]")
   vis.DefineScalarExpression("edge/Ey", "E[1]")
   vis.DefineScalarExpression("edge/Ez", "E[2]")
   vis.DefineScalarExpression("face/Bx", "B[0]")
   vis.DefineScalarExpression("face/By", "B[1]")
   vis.DefineScalarExpression("face/Bz", "B[2]")
   vis.DefineScalarExpression("volume/Ex", "E_vol[0]")
   vis.DefineScalarExpression("volume/Ey", "E_vol[1]")
   vis.DefineScalarExpression("volume/Ez", "E_vol[2]")
   vis.DefineScalarExpression("volume/Bx", "B_vol[0]")
   vis.DefineScalarExpression("volume/By", "B_vol[1]")
   vis.DefineScalarExpression("volume/Bz", "B_vol[2]")
   vis.DefineScalarExpression("Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineScalarExpression("Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineScalarExpression("Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineScalarExpression("Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineScalarExpression("Thermal energy", "1.5 * Pressure / rho")
   vis.DefineScalarExpression("divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineScalarExpression("divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineTensorExpression("PTensor", "{{PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]},{PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]},{PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]}}")
   vis.DefineVectorExpression("V", "{<volume/Vx>, <volume/Vy>, <volume/Vz>}")
   vis.DefineScalarExpression("face/perBx", "perturbed_B[0]")
   vis.DefineScalarExpression("face/perBy", "perturbed_B[1]")
   vis.DefineScalarExpression("face/perBz", "perturbed_B[2]")
   vis.DefineScalarExpression("V_magnitude", "magnitude(V)")
   vis.DefineScalarExpression("face/bgBx", "background_B[0]")
   vis.DefineScalarExpression("face/bgBy", "background_B[1]")
   vis.DefineScalarExpression("face/bgBz", "background_B[2]")
   vis.DefineVectorExpression("J", "(curl({<face/perBx>, <face/perBy>, <face/perBz>})) / 1.25663706144e-6")
   vis.DefineScalarExpression("J_magnitude", "magnitude(J)")
   vis.DefineVectorExpression("Hall", "cross(J,B)/(1.602177e-19 * (rho+0.1) )")
   vis.DefineScalarExpression("Hall_magnitude", "magnitude(Hall)")
   vis.DefineVectorExpression("V_parallel", "dot( V,B ) * B / (B_magnitude * B_magnitude )")
   vis.DefineVectorExpression("V_perp", "V-V_parallel")
   vis.DefineScalarExpression("V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineScalarExpression("V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineScalarExpression("alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineScalarExpression("Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineScalarExpression("gyrotropicity/theta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineVectorExpression("gyrotropicity/u", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {0.0, 0.0, 0.0},   {<volume/By>, -<volume/Bx>, 0.0} / (B_vol_magnitude * sin(<gyrotropicity/theta>)))")
   vis.DefineScalarExpression("gyrotropicity/ux", "<gyrotropicity/u>[0]")
   vis.DefineScalarExpression("gyrotropicity/uy", "<gyrotropicity/u>[1]")
   vis.DefineScalarExpression("gyrotropicity/uz", "<gyrotropicity/u>[2]")
   vis.DefineTensorExpression("gyrotropicity/R", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {{1.0, 0.0, 0.0},    {0.0, 1.0, 0.0},    {0.0, 0.0, 1.0}},   {      {cos(<gyrotropicity/theta>) + <gyrotropicity/ux>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/ux> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      <gyrotropicity/ux> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uy> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uy> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uy>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/uy> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/ux> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uz> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uy> * sin(<gyrotropicity/theta>),      <gyrotropicity/uz> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/ux> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uz>^2 * (1.0 - cos(<gyrotropicity/theta>))}   })")
   vis.DefineTensorExpression("PTensor_rotated", "<gyrotropicity/R>*PTensor*transpose(<gyrotropicity/R>)")
   vis.DefineVectorExpression("gyrotropicity/V_to_B", "<gyrotropicity/R>*B_vol")
   vis.DefineScalarExpression("gyrotropicity/minDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineTensorExpression("gyrotropicity/invR", "inverse(<gyrotropicity/R>)")
   vis.DefineScalarExpression("gyrotropicity/maxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineScalarExpression("Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineTensorExpression("TemperatureTensor", "{{PTensor[0][0]/rho, PTensor[0][1]/rho, PTensor[0][2]/rho},{PTensor[1][0]/rho, PTensor[1][1]/rho, PTensor[1][2]/rho},{PTensor[2][0]/rho, PTensor[2][1]/rho, PTensor[2][2]/rho}}")
   vis.DefineScalarExpression("TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineScalarExpression("volume/Vx", "rho_v[0] / rho")
   vis.DefineScalarExpression("volume/Vy", "rho_v[1]/rho")
   vis.DefineScalarExpression("volume/Vz", "rho_v[2]/rho")
   vis.DefineScalarExpression("edge/Ex", "E[0]")
   vis.DefineScalarExpression("edge/Ey", "E[1]")
   vis.DefineScalarExpression("edge/Ez", "E[2]")
   vis.DefineScalarExpression("face/Bx", "B[0]")
   vis.DefineScalarExpression("face/By", "B[1]")
   vis.DefineScalarExpression("face/Bz", "B[2]")
   vis.DefineScalarExpression("volume/Ex", "E_vol[0]")
   vis.DefineScalarExpression("volume/Ey", "E_vol[1]")
   vis.DefineScalarExpression("volume/Ez", "E_vol[2]")
   vis.DefineScalarExpression("volume/Bx", "B_vol[0]")
   vis.DefineScalarExpression("volume/By", "B_vol[1]")
   vis.DefineScalarExpression("volume/Bz", "B_vol[2]")
   vis.DefineScalarExpression("Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineScalarExpression("Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineScalarExpression("Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineScalarExpression("Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineScalarExpression("Thermal energy", "1.5 * Pressure / rho")
   vis.DefineScalarExpression("divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineScalarExpression("divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineTensorExpression("PTensor", "{{PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]},{PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]},{PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]}}")
   vis.DefineVectorExpression("V", "{<volume/Vx>, <volume/Vy>, <volume/Vz>}")
   vis.DefineScalarExpression("face/perBx", "perturbed_B[0]")
   vis.DefineScalarExpression("face/perBy", "perturbed_B[1]")
   vis.DefineScalarExpression("face/perBz", "perturbed_B[2]")
   vis.DefineScalarExpression("V_magnitude", "magnitude(V)")
   vis.DefineScalarExpression("face/bgBx", "background_B[0]")
   vis.DefineScalarExpression("face/bgBy", "background_B[1]")
   vis.DefineScalarExpression("face/bgBz", "background_B[2]")
   vis.DefineVectorExpression("J", "(curl({<face/perBx>, <face/perBy>, <face/perBz>})) / 1.25663706144e-6")
   vis.DefineScalarExpression("J_magnitude", "magnitude(J)")
   vis.DefineVectorExpression("Hall", "cross(J,B)/(1.602177e-19 * (rho+0.1) )")
   vis.DefineScalarExpression("Hall_magnitude", "magnitude(Hall)")
   vis.DefineVectorExpression("V_parallel", "dot( V,B ) * B / (B_magnitude * B_magnitude )")
   vis.DefineVectorExpression("V_perp", "V-V_parallel")
   vis.DefineScalarExpression("V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineScalarExpression("V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineScalarExpression("alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineScalarExpression("Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineScalarExpression("gyrotropicity/theta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineVectorExpression("gyrotropicity/u", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {0.0, 0.0, 0.0},   {<volume/By>, -<volume/Bx>, 0.0} / (B_vol_magnitude * sin(<gyrotropicity/theta>)))")
   vis.DefineScalarExpression("gyrotropicity/ux", "<gyrotropicity/u>[0]")
   vis.DefineScalarExpression("gyrotropicity/uy", "<gyrotropicity/u>[1]")
   vis.DefineScalarExpression("gyrotropicity/uz", "<gyrotropicity/u>[2]")
   vis.DefineTensorExpression("gyrotropicity/R", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {{1.0, 0.0, 0.0},    {0.0, 1.0, 0.0},    {0.0, 0.0, 1.0}},   {      {cos(<gyrotropicity/theta>) + <gyrotropicity/ux>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/ux> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      <gyrotropicity/ux> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uy> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uy> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uy>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/uy> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/ux> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uz> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uy> * sin(<gyrotropicity/theta>),      <gyrotropicity/uz> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/ux> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uz>^2 * (1.0 - cos(<gyrotropicity/theta>))}   })")
   vis.DefineTensorExpression("PTensor_rotated", "<gyrotropicity/R>*PTensor*transpose(<gyrotropicity/R>)")
   vis.DefineVectorExpression("gyrotropicity/V_to_B", "<gyrotropicity/R>*B_vol")
   vis.DefineScalarExpression("gyrotropicity/minDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineTensorExpression("gyrotropicity/invR", "inverse(<gyrotropicity/R>)")
   vis.DefineScalarExpression("gyrotropicity/maxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineScalarExpression("Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineTensorExpression("TemperatureTensor", "{{PTensor[0][0]/rho, PTensor[0][1]/rho, PTensor[0][2]/rho},{PTensor[1][0]/rho, PTensor[1][1]/rho, PTensor[1][2]/rho},{PTensor[2][0]/rho, PTensor[2][1]/rho, PTensor[2][2]/rho}}")
   vis.DefineScalarExpression("TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineScalarExpression("volume/Vx", "rho_v[0] / rho")
   vis.DefineScalarExpression("volume/Vy", "rho_v[1]/rho")
   vis.DefineScalarExpression("volume/Vz", "rho_v[2]/rho")
   vis.DefineScalarExpression("edge/Ex", "E[0]")
   vis.DefineScalarExpression("edge/Ey", "E[1]")
   vis.DefineScalarExpression("edge/Ez", "E[2]")
   vis.DefineScalarExpression("face/Bx", "B[0]")
   vis.DefineScalarExpression("face/By", "B[1]")
   vis.DefineScalarExpression("face/Bz", "B[2]")
   vis.DefineScalarExpression("volume/Ex", "E_vol[0]")
   vis.DefineScalarExpression("volume/Ey", "E_vol[1]")
   vis.DefineScalarExpression("volume/Ez", "E_vol[2]")
   vis.DefineScalarExpression("volume/Bx", "B_vol[0]")
   vis.DefineScalarExpression("volume/By", "B_vol[1]")
   vis.DefineScalarExpression("volume/Bz", "B_vol[2]")
   vis.DefineScalarExpression("Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineScalarExpression("Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineScalarExpression("Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineScalarExpression("Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineScalarExpression("Thermal energy", "1.5 * Pressure / rho")
   vis.DefineScalarExpression("divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineScalarExpression("divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineTensorExpression("PTensor", "{{PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]},{PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]},{PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]}}")
   vis.DefineVectorExpression("V", "{<volume/Vx>, <volume/Vy>, <volume/Vz>}")
   vis.DefineScalarExpression("face/perBx", "perturbed_B[0]")
   vis.DefineScalarExpression("face/perBy", "perturbed_B[1]")
   vis.DefineScalarExpression("face/perBz", "perturbed_B[2]")
   vis.DefineScalarExpression("V_magnitude", "magnitude(V)")
   vis.DefineScalarExpression("face/bgBx", "background_B[0]")
   vis.DefineScalarExpression("face/bgBy", "background_B[1]")
   vis.DefineScalarExpression("face/bgBz", "background_B[2]")
   vis.DefineVectorExpression("J", "(curl({<face/perBx>, <face/perBy>, <face/perBz>})) / 1.25663706144e-6")
   vis.DefineScalarExpression("J_magnitude", "magnitude(J)")
   vis.DefineVectorExpression("Hall", "cross(J,B)/(1.602177e-19 * (rho+0.1) )")
   vis.DefineScalarExpression("Hall_magnitude", "magnitude(Hall)")
   vis.DefineVectorExpression("V_parallel", "dot( V,B ) * B / (B_magnitude * B_magnitude )")
   vis.DefineVectorExpression("V_perp", "V-V_parallel")
   vis.DefineScalarExpression("V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineScalarExpression("V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineScalarExpression("alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineScalarExpression("Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineScalarExpression("gyrotropicity/theta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineVectorExpression("gyrotropicity/u", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {0.0, 0.0, 0.0},   {<volume/By>, -<volume/Bx>, 0.0} / (B_vol_magnitude * sin(<gyrotropicity/theta>)))")
   vis.DefineScalarExpression("gyrotropicity/ux", "<gyrotropicity/u>[0]")
   vis.DefineScalarExpression("gyrotropicity/uy", "<gyrotropicity/u>[1]")
   vis.DefineScalarExpression("gyrotropicity/uz", "<gyrotropicity/u>[2]")
   vis.DefineTensorExpression("gyrotropicity/R", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {{1.0, 0.0, 0.0},    {0.0, 1.0, 0.0},    {0.0, 0.0, 1.0}},   {      {cos(<gyrotropicity/theta>) + <gyrotropicity/ux>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/ux> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      <gyrotropicity/ux> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uy> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uy> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uy>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/uy> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/ux> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uz> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uy> * sin(<gyrotropicity/theta>),      <gyrotropicity/uz> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/ux> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uz>^2 * (1.0 - cos(<gyrotropicity/theta>))}   })")
   vis.DefineTensorExpression("PTensor_rotated", "<gyrotropicity/R>*PTensor*transpose(<gyrotropicity/R>)")
   vis.DefineVectorExpression("gyrotropicity/V_to_B", "<gyrotropicity/R>*B_vol")
   vis.DefineScalarExpression("gyrotropicity/minDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineTensorExpression("gyrotropicity/invR", "inverse(<gyrotropicity/R>)")
   vis.DefineScalarExpression("gyrotropicity/maxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineScalarExpression("Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineTensorExpression("TemperatureTensor", "{{PTensor[0][0]/rho, PTensor[0][1]/rho, PTensor[0][2]/rho},{PTensor[1][0]/rho, PTensor[1][1]/rho, PTensor[1][2]/rho},{PTensor[2][0]/rho, PTensor[2][1]/rho, PTensor[2][2]/rho}}")
   vis.DefineScalarExpression("TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineScalarExpression("volume/Vx", "rho_v[0] / rho")
   vis.DefineScalarExpression("volume/Vy", "rho_v[1]/rho")
   vis.DefineScalarExpression("volume/Vz", "rho_v[2]/rho")
   vis.DefineScalarExpression("edge/Ex", "E[0]")
   vis.DefineScalarExpression("edge/Ey", "E[1]")
   vis.DefineScalarExpression("edge/Ez", "E[2]")
   vis.DefineScalarExpression("face/Bx", "B[0]")
   vis.DefineScalarExpression("face/By", "B[1]")
   vis.DefineScalarExpression("face/Bz", "B[2]")
   vis.DefineScalarExpression("volume/Ex", "E_vol[0]")
   vis.DefineScalarExpression("volume/Ey", "E_vol[1]")
   vis.DefineScalarExpression("volume/Ez", "E_vol[2]")
   vis.DefineScalarExpression("volume/Bx", "B_vol[0]")
   vis.DefineScalarExpression("volume/By", "B_vol[1]")
   vis.DefineScalarExpression("volume/Bz", "B_vol[2]")
   vis.DefineScalarExpression("Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineScalarExpression("Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineScalarExpression("Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineScalarExpression("Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineScalarExpression("Thermal energy", "1.5 * Pressure / rho")
   vis.DefineScalarExpression("divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineScalarExpression("divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineTensorExpression("PTensor", "{{PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]},{PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]},{PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]}}")
   vis.DefineVectorExpression("V", "{<volume/Vx>, <volume/Vy>, <volume/Vz>}")
   vis.DefineScalarExpression("face/perBx", "perturbed_B[0]")
   vis.DefineScalarExpression("face/perBy", "perturbed_B[1]")
   vis.DefineScalarExpression("face/perBz", "perturbed_B[2]")
   vis.DefineScalarExpression("V_magnitude", "magnitude(V)")
   vis.DefineScalarExpression("face/bgBx", "background_B[0]")
   vis.DefineScalarExpression("face/bgBy", "background_B[1]")
   vis.DefineScalarExpression("face/bgBz", "background_B[2]")
   vis.DefineVectorExpression("J", "(curl({<face/perBx>, <face/perBy>, <face/perBz>})) / 1.25663706144e-6")
   vis.DefineScalarExpression("J_magnitude", "magnitude(J)")
   vis.DefineVectorExpression("Hall", "cross(J,B)/(1.602177e-19 * (rho+0.1) )")
   vis.DefineScalarExpression("Hall_magnitude", "magnitude(Hall)")
   vis.DefineVectorExpression("V_parallel", "dot( V,B ) * B / (B_magnitude * B_magnitude )")
   vis.DefineVectorExpression("V_perp", "V-V_parallel")
   vis.DefineScalarExpression("V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineScalarExpression("V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineScalarExpression("alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineScalarExpression("Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineScalarExpression("gyrotropicity/theta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineVectorExpression("gyrotropicity/u", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {0.0, 0.0, 0.0},   {<volume/By>, -<volume/Bx>, 0.0} / (B_vol_magnitude * sin(<gyrotropicity/theta>)))")
   vis.DefineScalarExpression("gyrotropicity/ux", "<gyrotropicity/u>[0]")
   vis.DefineScalarExpression("gyrotropicity/uy", "<gyrotropicity/u>[1]")
   vis.DefineScalarExpression("gyrotropicity/uz", "<gyrotropicity/u>[2]")
   vis.DefineTensorExpression("gyrotropicity/R", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   {{1.0, 0.0, 0.0},    {0.0, 1.0, 0.0},    {0.0, 0.0, 1.0}},   {      {cos(<gyrotropicity/theta>) + <gyrotropicity/ux>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/ux> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      <gyrotropicity/ux> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uy> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uy> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/uz> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uy>^2 * (1.0 - cos(<gyrotropicity/theta>)),      <gyrotropicity/uy> * <gyrotropicity/uz> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/ux> * sin(<gyrotropicity/theta>)},      {<gyrotropicity/uz> * <gyrotropicity/ux> * (1.0 - cos(<gyrotropicity/theta>)) - <gyrotropicity/uy> * sin(<gyrotropicity/theta>),      <gyrotropicity/uz> * <gyrotropicity/uy> * (1.0 - cos(<gyrotropicity/theta>)) + <gyrotropicity/ux> * sin(<gyrotropicity/theta>),      cos(<gyrotropicity/theta>) + <gyrotropicity/uz>^2 * (1.0 - cos(<gyrotropicity/theta>))}   })")
   vis.DefineTensorExpression("PTensor_rotated", "<gyrotropicity/R>*PTensor*transpose(<gyrotropicity/R>)")
   vis.DefineVectorExpression("gyrotropicity/V_to_B", "<gyrotropicity/R>*B_vol")
   vis.DefineScalarExpression("gyrotropicity/minDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineTensorExpression("gyrotropicity/invR", "inverse(<gyrotropicity/R>)")
   vis.DefineScalarExpression("gyrotropicity/maxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineScalarExpression("Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineTensorExpressiourveExpression("operators/Lineout/volumeVx", "rho_v[0] / rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVy", "rho_v[1]/rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVz", "rho_v[2]/rho")
   vis.DefineCurveExpression("operators/Lineout/edgeEx", "E[0]")
   vis.DefineCurveExpression("operators/Lineout/edgeEy", "E[1]")
   vis.DefineCurveExpression("operators/Lineout/edgeEz", "E[2]")
   vis.DefineCurveExpression("operators/Lineout/faceBx", "B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceBy", "B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceBz", "B[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeEx", "E_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeEy", "E_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeEz", "E_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeBx", "B_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeBy", "B_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeBz", "B_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineCurveExpression("operators/Lineout/Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineCurveExpression("operators/Lineout/Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineCurveExpression("operators/Lineout/Thermal energy", "1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineCurveExpression("operators/Lineout/divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/faceperBx", "perturbed_B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceperBy", "perturbed_B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceperBz", "perturbed_B[2]")
   vis.DefineCurveExpression("operators/Lineout/V_magnitude", "magnitude(V)")
   vis.DefineCurveExpression("operators/Lineout/facebgBx", "background_B[0]")
   vis.DefineCurveExpression("operators/Lineout/facebgBy", "background_B[1]")
   vis.DefineCurveExpression("operators/Lineout/facebgBz", "background_B[2]")
   vis.DefineCurveExpression("operators/Lineout/J_magnitude", "magnitude(J)")
   vis.DefineCurveExpression("operators/Lineout/Hall_magnitude", "magnitude(Hall)")
   vis.DefineCurveExpression("operators/Lineout/V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineCurveExpression("operators/Lineout/V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineCurveExpression("operators/Lineout/alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitytheta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityux", "<gyrotropicity/u>[0]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuy", "<gyrotropicity/u>[1]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuz", "<gyrotropicity/u>[2]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityminDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitymaxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineCurveExpression("operators/Lineout/TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineCurveExpression("operators/Lineout/volumeVx", "rho_v[0] / rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVy", "rho_v[1]/rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVz", "rho_v[2]/rho")
   vis.DefineCurveExpression("operators/Lineout/edgeEx", "E[0]")
   vis.DefineCurveExpression("operators/Lineout/edgeEy", "E[1]")
   vis.DefineCurveExpression("operators/Lineout/edgeEz", "E[2]")
   vis.DefineCurveExpression("operators/Lineout/faceBx", "B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceBy", "B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceBz", "B[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeEx", "E_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeEy", "E_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeEz", "E_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeBx", "B_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeBy", "B_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeBz", "B_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineCurveExpression("operators/Lineout/Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineCurveExpression("operators/Lineout/Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineCurveExpression("operators/Lineout/Thermal energy", "1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineCurveExpression("operators/Lineout/divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/faceperBx", "perturbed_B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceperBy", "perturbed_B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceperBz", "perturbed_B[2]")
   vis.DefineCurveExpression("operators/Lineout/V_magnitude", "magnitude(V)")
   vis.DefineCurveExpression("operators/Lineout/facebgBx", "background_B[0]")
   vis.DefineCurveExpression("operators/Lineout/facebgBy", "background_B[1]")
   vis.DefineCurveExpression("operators/Lineout/facebgBz", "background_B[2]")
   vis.DefineCurveExpression("operators/Lineout/J_magnitude", "magnitude(J)")
   vis.DefineCurveExpression("operators/Lineout/Hall_magnitude", "magnitude(Hall)")
   vis.DefineCurveExpression("operators/Lineout/V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineCurveExpression("operators/Lineout/V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineCurveExpression("operators/Lineout/alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitytheta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityux", "<gyrotropicity/u>[0]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuy", "<gyrotropicity/u>[1]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuz", "<gyrotropicity/u>[2]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityminDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitymaxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineCurveExpression("operators/Lineout/TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineCurveExpression("operators/Lineout/volumeVx", "rho_v[0] / rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVy", "rho_v[1]/rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVz", "rho_v[2]/rho")
   vis.DefineCurveExpression("operators/Lineout/edgeEx", "E[0]")
   vis.DefineCurveExpression("operators/Lineout/edgeEy", "E[1]")
   vis.DefineCurveExpression("operators/Lineout/edgeEz", "E[2]")
   vis.DefineCurveExpression("operators/Lineout/faceBx", "B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceBy", "B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceBz", "B[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeEx", "E_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeEy", "E_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeEz", "E_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeBx", "B_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeBy", "B_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeBz", "B_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineCurveExpression("operators/Lineout/Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineCurveExpression("operators/Lineout/Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineCurveExpression("operators/Lineout/Thermal energy", "1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineCurveExpression("operators/Lineout/divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/faceperBx", "perturbed_B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceperBy", "perturbed_B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceperBz", "perturbed_B[2]")
   vis.DefineCurveExpression("operators/Lineout/V_magnitude", "magnitude(V)")
   vis.DefineCurveExpression("operators/Lineout/facebgBx", "background_B[0]")
   vis.DefineCurveExpression("operators/Lineout/facebgBy", "background_B[1]")
   vis.DefineCurveExpression("operators/Lineout/facebgBz", "background_B[2]")
   vis.DefineCurveExpression("operators/Lineout/J_magnitude", "magnitude(J)")
   vis.DefineCurveExpression("operators/Lineout/Hall_magnitude", "magnitude(Hall)")
   vis.DefineCurveExpression("operators/Lineout/V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineCurveExpression("operators/Lineout/V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineCurveExpression("operators/Lineout/alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitytheta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityux", "<gyrotropicity/u>[0]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuy", "<gyrotropicity/u>[1]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuz", "<gyrotropicity/u>[2]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityminDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitymaxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineCurveExpression("operators/Lineout/TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")
   vis.DefineCurveExpression("operators/Lineout/volumeVx", "rho_v[0] / rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVy", "rho_v[1]/rho")
   vis.DefineCurveExpression("operators/Lineout/volumeVz", "rho_v[2]/rho")
   vis.DefineCurveExpression("operators/Lineout/edgeEx", "E[0]")
   vis.DefineCurveExpression("operators/Lineout/edgeEy", "E[1]")
   vis.DefineCurveExpression("operators/Lineout/edgeEz", "E[2]")
   vis.DefineCurveExpression("operators/Lineout/faceBx", "B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceBy", "B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceBz", "B[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeEx", "E_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeEy", "E_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeEz", "E_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/volumeBx", "B_vol[0]")
   vis.DefineCurveExpression("operators/Lineout/volumeBy", "B_vol[1]")
   vis.DefineCurveExpression("operators/Lineout/volumeBz", "B_vol[2]")
   vis.DefineCurveExpression("operators/Lineout/Total energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>) + (<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6) + 0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>) + 1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/Electric energy", "0.5 * 8.85418781762e-12 * (<volume/Ex>*<volume/Ex> + <volume/Ey>*<volume/Ey> + <volume/Ez>*<volume/Ez>)")
   vis.DefineCurveExpression("operators/Lineout/Magnetic energy", "(<volume/Bx>*<volume/Bx> + <volume/By>*<volume/By> + <volume/Bz>*<volume/Bz>) / (2 * 1.25663706144e-6)")
   vis.DefineCurveExpression("operators/Lineout/Kinetic energy", "0.5 * 1.6726217e-27 * rho * (<volume/Vx>*<volume/Vx> + <volume/Vy>*<volume/Vy> + <volume/Vz>*<volume/Vz>)")
   vis.DefineCurveExpression("operators/Lineout/Thermal energy", "1.5 * Pressure / rho")
   vis.DefineCurveExpression("operators/Lineout/divB", "divergence({<face/Bx>, <face/By>, <face/Bz>})")
   vis.DefineCurveExpression("operators/Lineout/divB_scaled", "relative_size(SpatialGrid) * divergence({<face/Bx>, <face/By>, <face/Bz>}) / B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/faceperBx", "perturbed_B[0]")
   vis.DefineCurveExpression("operators/Lineout/faceperBy", "perturbed_B[1]")
   vis.DefineCurveExpression("operators/Lineout/faceperBz", "perturbed_B[2]")
   vis.DefineCurveExpression("operators/Lineout/V_magnitude", "magnitude(V)")
   vis.DefineCurveExpression("operators/Lineout/facebgBx", "background_B[0]")
   vis.DefineCurveExpression("operators/Lineout/facebgBy", "background_B[1]")
   vis.DefineCurveExpression("operators/Lineout/facebgBz", "background_B[2]")
   vis.DefineCurveExpression("operators/Lineout/J_magnitude", "magnitude(J)")
   vis.DefineCurveExpression("operators/Lineout/Hall_magnitude", "magnitude(Hall)")
   vis.DefineCurveExpression("operators/Lineout/V_parallel_magnitude", "magnitude(V_parallel)")
   vis.DefineCurveExpression("operators/Lineout/V_perp_magnitude", "magnitude(V_perp)")
   vis.DefineCurveExpression("operators/Lineout/alfvenicMachNumber", "V_magnitude*sqrt(1.256637e-6*1.6726e-27*rho)/B_magnitude")
   vis.DefineCurveExpression("operators/Lineout/Heli's criterion", "rho*<volume/Vx>*<volume/Vx>/ (1.0e6*25.0e10)")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitytheta", "if(   and(      eq(<volume/By>, 0.0),      eq(<volume/Bx>, 0.0)   ),   0.0,   acos(<volume/Bz> / B_vol_magnitude))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityux", "<gyrotropicity/u>[0]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuy", "<gyrotropicity/u>[1]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityuz", "<gyrotropicity/u>[2]")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicityminDIagonalPressureDifference", "min(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/gyrotropicitymaxDIagonalPressureDifference", "max(   2.0 * abs(PTensor_rotated[0][0] - PTensor_rotated[1][1]) / (PTensor_rotated[0][0] + PTensor_rotated[1][1]),   2.0 * abs(PTensor_rotated[1][1] - PTensor_rotated[2][2]) / (PTensor_rotated[1][1] + PTensor_rotated[2][2]))")
   vis.DefineCurveExpression("operators/Lineout/Temperature", "Pressure / ((rho + 1) * 1.38065e-23)")
   vis.DefineCurveExpression("operators/Lineout/TemperatureTensor_rotated", "<gyrotropicity/R>*TemperatureTensor*transpose(<gyrotropicity/R>)")

