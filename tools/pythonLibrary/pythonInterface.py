# Import files
visitImported = False
import sys
import subprocess
from numpy import *
from termcolor import colored
if visitImported == False:
   from visit import *
   visitImported = True
from useroptions import *
sys.path.insert(0, pythonLibDirectoryPath + pyVlsvPath)
sys.path.insert(0, pythonLibDirectoryPath + pyVisitPath)
sys.path.insert(0, pythonLibDirectoryPath + pyMiscPath)
from vlsvreader import *
from makeamovie import *
from makeamovieauto import *
from spacecraft import *
from distributionplot import *
from movingframeofreference import *
from movingline import *
from vlsvplots import *
from plotfunctions import *
from filenames import *
from getvariables import *
from fourier2d import *
from plots import *
from pitchangle import *
from rotation import *

visitLaunched = False

def launch_visit(noWindow=True):
   '''
   Launches visit
   :param noWindow=True   Determines whether window is shown or not
   '''
   global visitLaunched
   if visitLaunched == True:
      print "Visit already launched"
      return
   if noWindow == True:
      LaunchNowin(vdir=pathToVisit)
   else:
      Launch(vdir=pathToVisit)
   visitLaunched = True
   print "Visit launched"
   return

def list_functions():
   print "Class: " + colored("VlsvFile", "red")
   print "   Functions:"
   print "      vlsvReader = VlsvFile(\"distribution.vlsv\") -- opens a vlsv file for reading"
   print "      vlsvReader.list() -- Gives a list of attributes in the vlsv file"
   print "      vlsvReader.read(name=\"rho\", tag=\"VARIABLE\", mesh=\"SpatialGrid\", read_single_cellid=-1) -- Reads array with the name 'rho', tag 'VARIABLE', mesh 'SpatialGrid'. if read_single_cellid is specified then the reader reads only the given cell id, if specified as -1 it reads the whole array."
   print "      vlsvReader.read_variables(\"rho\") -- Reads values of rho in the form of array"
   print "      vlsvReader.read_variable(\"rho\", 16) -- Reads the 16th cell id's value of rho and returns it"
   print "      vlsvReader.read_blocks(16) -- Reads the raw blocks of cell 16\n"
   print "      vlsvReader.get_cell_coordinates(cellid=21421) -- Returns a given cell's coordinates"
   print "      vlsvReader.get_cellid(coordinates=[2e6, 120e6, 0]) -- Returns the cell at given coordinates"
   print "Function: " + colored("launch_visit", "red")
   print "Function: " + colored("make_movie", "red")
   print "Function: " + colored("make_movie_auto", "red")
   print "Function: " + colored("make_moving_frame_of_reference_movie", "red")
   print "Function: " + colored("make_moving_frame_of_reference_line_plot", "red")
   print "Function: " + colored("make_distribution_movie", "red")
   print "Function: " + colored("make_moving_spacecraft_plot", "red")
   print "Function: " + colored("draw_plots_by_cellid", "red")
   print "Function: " + colored("draw_point_picture", "red")
   print "Function: " + colored("draw_plots_by_cellid", "red")
   print "Function: " + colored("draw_by_cellid_array", "red")
   print "Function: " + colored("plot_time_evolution", "red")
   print "Function: " + colored("take_cut_through", "red")
   print "Function: " + colored("take_cut_through_array", "red")
   print "Function: " + colored("time_evolution_array", "red")
   print "Function: " + colored("list_fit_functions", "red")
   print "Function: " + colored("get_cell_coordinates_from_file", "red")
   print "Function: " + colored("get_cell_id_from_file", "red")



def old_list_functions():
   print "Class: " + colored("VlsvFile", "red")
   print "   Example usage:"
   print "      vlsvReader = VlsvFile(\"distribution.vlsv\") -- opens a vlsv file for reading"
   print "      vlsvReader.list() -- Gives a list of attributes in the vlsv file"
   print "      vlsvReader.read(name=\"rho\", tag=\"VARIABLE\", mesh=\"SpatialGrid\", read_single_cellid=-1) -- Reads array with the name 'rho', tag 'VARIABLE', mesh 'SpatialGrid'. if read_single_cellid is specified then the reader reads only the given cell id, if specified as -1 it reads the whole array."
   print "      vlsvReader.read_variables(\"rho\") -- Reads values of rho in the form of array"
   print "      vlsvReader.read_variable(\"rho\", 16) -- Reads the 16th cell id's value of rho and returns it"
   print "      vlsvReader.read_blocks(16) -- Reads the raw blocks of cell 16\n"
   print "Function: " + colored("launch_visit(noWindow=True)", "red")
   print "   Example usage:"
   print "      launch_visit(noWindow=False) NOTE: this must be launched before doing anything regarding visit\n"
   print "Function: " + colored("make_movie( variableName, minValue, maxValue, inputDirectory, inputFileName, outputDirectory, outputFileName, colorTable=\"hot_desaturated\", startFrame=-1, endFrame=-1 )", "red")
   print "   Example usage:"
   print "      make_movie(variableName=\"rho\", minValue=1.0e6, maxValue=5.0e6, inputDirectory=\"/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/\", inputFileName=\"bulk.*.silo\", outputDirectory=\"/home/hannukse/MOVINGFRAME_MOVIES/AAJ_BZ_REMAKE/\", outputFileName=\"RHO_FORESHOCK_\", colorTable=\"hot_desaturated\", startFrame=30, endFrame=120)\n"
   print "Function: " + colored("make_movie_auto(BBOX)", "red")
   print "   Example usage:"
   print "      make_movie_auto( variableName=\"rho\", boundaryBox=array([1.14322e8, 2.51245e8, -3.38029e7, 3.89375e7, 0, 0]), vlsvFileName=\"bulk.0000970.vlsv\", inputDirectory=\"/home/hannukse/meteo/lustre/tmp/hannuksela/pythonTests/\", inputFileName=\"bulk.*.silo\", outputDirectory=\"/home/hannukse/testFile/\", outputFileName=\"RHO_TEST\", colorTableName=\"hot_desaturated\", startFrame=-1, endFrame=-1 )\n"
   print "Function: " + colored("make_moving_frame_of_reference_movie( x_begin, x_end, y_begin, y_end, speed_x, speed_y, variable_name, minThreshold, maxThreshold, input_directory, input_file_name, output_directory, output_file_name, color_table=\"hot_desaturated\", start_frame=-1, end_frame=-1, frame_skip_dt=1.0 )", "red")
   print "   Example usage:"
   print "      make_moving_frame_of_reference_movie( x_begin=50e6, x_end=150e6, y_begin=-150e6, y_end=-50e6, speed_x=-500000, speed_y=1000, variable_name=\"rho\", minThreshold=1.0e5, maxThreshold=1.0e6, input_directory=\"/stornext/field/vlasiator/AAJ/distributions/\", input_file_name=\"bulk.*.silo\", output_directory=\"/stornext/field/vlasiator/visualizations/\", output_file_name=\"RHO_MOVIE\", color_table=\"hot_desaturated\", start_frame=0, end_frame=85, frame_skip_dt=2.0 )\n"
   print "Function: " + colored("make_moving_frame_of_reference_line_plot(  point1, point2, velocity, variable_name, input_directory, input_file_name, output_directory, output_file_name, start_frame=-1, end_frame=-1, frame_skip_dt=1.0 )", "red")
   print "   Example usage:"
   print "      make_moving_frame_of_reference_line_plot( point1=[40e6, 140e6, 0], point2=[12e6, 523e6, 0], velocity=[-500000, 0, 0], variable_name=\"rho\", input_directory=\"/stornext/field/vlasiator/2D/AAJ/silo_files/\", input_file_name=\"bulk.*.silo\", output_directory=\"/stornext/field/vlasiator/visualizations/movies/\", output_file_name=\"RHO_MOVIE_\", start_frame=10, end_frame=50, frame_skip_dt=1.5 )\n"
   print "Function: " + colored("make_distribution_movie(cellids, rotated, inputDirectory, outputDirectory, outputFileName, zoom=1.0, viewNormal=[0.488281, 0.382966, -0.784167], minThreshold=1e-18, maxThreshold=1e37)", "red")
   print "   Example usage:"
   print "      make_distribution_movie(cellids=[18302, 19432, 19042], rotated=True, inputDirectory=\"/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/\", outputDirectory=\"/home/hannukse/MOVIES/\", outputFileName=\"testmovie\", zoom=0.8, viewNormal=[0.488281, 0.382966, -0.784167], minThreshold=1e-17, maxThreshold=1.2e37)"
   print "   Note: viewNormal determines the angle of view (straight from visit)\n"
   print "Function: " + colored("make_moving_spacecraft_plot()", "red")
   print "   Example usage:"
   print "      make_moving_spacecraft_plot()\n"
   print "Function: " + colored("draw_plots_by_cellid(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[])", "red")
   print "Function: " + colored("list_fit_functions()", "red")
   print "   Example usage:"
   print "      list_fit_functions()"
   print "Function: " + colored("take_cut_through( fileName, variables1, variables2, point1, point2 )", "red")
   print "   Example usage:"
   print "      Note: this would create plots of (x, rho), (x, B[1]), (x, B[2]), (y, rho), (y, B[1]), (y, B[2]) and try to fit a*x^2 + b*x + c curve into the data"
   print "      variables2 = {}"
   print "      variables2[\"rho\"] = 0"
   print "      variables2[\"B\"] = [1, 2]"
   print "      variables1 = {}"
   print "      variables1[\"coordinates\"] = [0, 1]"
   print "      take_cut_through( fileName=\"bulk.00001.vlsv\", variables1=variables1, variables2=variables2, point1=[0,0,0], point2=[150e6, 0, 0], fitFunction=\"polynomialsecond\" )\n"
   print "Function: " + colored("take_cut_through_array(fileName, variables1, variables2, point1, point2 )", "red")
   print "   Example usage:"
   print "      Note: Same as take_cut_through but instead of plotting, returns array with variables"
   print "      variables2 = {}"
   print "      variables2[\"rho\"] = 0"
   print "      variables2[\"B\"] = [1, 2]"
   print "      variables1 = {}"
   print "      variables1[\"coordinates\"] = [0, 1]"
   print "      take_cut_through( fileName=\"bulk.00001.vlsv\", variables1=variables1, variables2=variables2, point1=[0,0,0], point2=[150e6, 0, 0] )\n"
   print "Function: " + colored("draw_plots_by_cellid(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[], fitFunction=nullfit)", "red")
   print "   Example usage:"
   print "      Note: Draws a plot of given variables for the given cell ids and fits a function if specified"
   print "      variables2 = {}"
   print "      variables2[\"rho\"] = 0"
   print "      variables2[\"B\"] = [1, 2]"
   print "      variables1 = {}"
   print "      variables1[\"coordinates\"] = [0, 1]"
   print "      draw_plots_by_cellid(vlsvReader=VlsvFile(\"bulk.00001.vlsv\"), variables1=variables1, variables2=variables2, cellids=[15600, 1200, 2000], fitFunction=\"polynomialsecond\")\n"
   print "Function: " + colored("draw_by_cellid_array(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[], fitFunction=nullfit)", "red")
   print "   Example usage:"
   print "      Note: Returns variables for the given cell ids as well as fitted function if specified Format: np.array([variables1list, variables2list, variables1listfitted, variables2listfitted])"
   print "      variables2 = {}"
   print "      variables2[\"rho\"] = 0"
   print "      variables2[\"B\"] = [1, 2]"
   print "      variables1 = {}"
   print "      variables1[\"coordinates\"] = [0, 1]"
   print "      draw_by_cellid_array(vlsvReader=VlsvFile(\"bulk.00001.vlsv\"), variables1=variables1, variables2=variables2, cellids=[15600, 1200, 2000], fitFunction=\"polynomialsecond\")\n"
   print "Function: " + colored("plot_time_evolution( fileNames, dt, cellid, variables, forceConstantAmplitude=False, fitFunction=nullfit )", "red")
   print "   Example usage:"
   print "      Note: Plots the time evolution of some cell and fits a fourier series in the plot"
   print "      variables = {}"
   print "      variables[\"rho\"] = 0"
   print "         plot_time_evolution( fileNames=[\"bulk.00000.vlsv\", \"bulk.00001.vlsv\", \"bulk.00002.vlsv\"], dt=1, cellid=15000, variables=variables, forceConstantAmplitude=True, fitFunction=\"polynomialsecond\" )\n"
   print "Function: " + colored("time_evolution_array( fileNames, dt, cellid, variables, forceConstantAmplitude=False, fitFunction=nullfit )", "red")
   print "   Example usage:"
   print "      Note: Same as plot_time_evolution but instead of plotting returns an array"
   print "      variables = {}"
   print "      variables[\"rho\"] = 0"
   print "      result = time_evolution_array( fileNames=[\"bulk.00000.vlsv\", \"bulk.00001.vlsv\", \"bulk.00002.vlsv\"], dt=1, cellid=15000, variables=variables, forceConstantAmplitude=True, fitFunction=\"polynomialsecond\" )"
   print "      fourier_variables = result[0]"
   print "      y = result[1]"
   print "      t = np.arange(3) * dt\n"
   print "Function: " + colored("draw_point_picture( variableName, minValue, maxValue, inputDirectory, inputFileNames, coordinates, outputDirectory, outputFileName, colorTable=\"hot_desaturated\")", "red")
   print "   Example usage:"
   print "      draw_point_picture( variableName=\"rho\", minValue=9.2e5, maxValue=1.2e6, inputDirectory=\"/stornext/field/vlasiator/2D/AAJ/silo_files/\", inputFileNames=\"bulk.*.silo\", coordinates=[[200e6, 0, 0], [150e6, -50e6, 0], [120e6, -90e6, 0]], outputDirectory, outputFileName, colorTable=\"hot_desaturated\")\n"





