# This function includes helper functions for fitting fourier series and plotting frequencies/amplitudes

import numpy as np
import pylab as pl


def saveplotdata( xcoordinates, ycoordinates, fileName ):
   xcoordinates = np.atleast_1d(xcoordinates)
   ycoordinates = np.atleast_1d(ycoordinates)
   if len(xcoordinates) != len(ycoordinates):
      print "BAD LENGTHS"
      return False
   if len(np.atleast_1d(xcoordinates[0])) != 1:
      data = np.array([[xcoordinates[i], ycoordinates[i]] for i in xrange(len(xcoordinates))])
   else:
      data = np.array([[xcoordinates, ycoordinates]])
   np.save(fileName, data)
   return True

def loadplotdata( fileName, showplot=True ):
   data = np.load(fileName)
   xcoordinates = np.array([data[i][0] for i in xrange(len(data))])
   ycoordinates = np.array([data[i][1] for i in xrange(len(data))])
   maxPlots = len(data)
   subplotnums = [maxPlots, 1, 1]
   pl.figure()
   for i in xrange(maxPlots):
      pl.subplot(subplotnums[0], subplotnums[1], subplotnums[2])
      pl.plot(xcoordinates[i], ycoordinates[i])
      subplotnums[2] = subplotnums[2] + 1
   if showplot == True:
      pl.ion()
      pl.show()


def fourier_array(t, y, kaiserwindowparameter=0):
   ''' Function for returning fourier series of some given arrays t and y
       :param t         time variable (along the x-axis)
       :param y         variable along the y-axis
       :returns         numpy array of the fourier series [t2, y2] and frequency [freq]: np.array([t2,y2], freq)
   '''
   # First check the t array whether it has a constant dt
   dt = t[1] - t[0]
   for i in xrange(len(t)-1):
      if dt != t[i+1] - t[i]:
         print "Gave bad timestep to plot_fourier, the time step in array t must be constant (for now)"
   # Use kaiser window on y
   y = y*np.kaiser(len(y), kaiserwindowparameter)
   # Do FFT on the data
   fourier=np.fft.fft(y)*(1/(float)(len(t)))
   # Get frequencies of the fourier
   freq=np.fft.fftfreq(len(fourier), d=dt)
   # Declare t2 (Note: This is the same as t but we want the steps to be thicker so the data looks smoother
   dt2=dt*0.01
   t2=np.arange(len(t)*100)*dt2
   # Declare y2
   y2=np.array([np.sum(fourier*np.exp(complex(0,1)*2*np.pi*freq*T)) for T in t2])
   return np.array([[t2,y2], freq])

def plot_fourier(t, y, subplotnums=[[2,1,1],[2,1,2]], savedata="none", kaiserwindowparameter=0):
   ''' Function for plotting fourier series of some given arrays t and y
       :param t         time variable (plotted along the x-axis)
       :param y         variable to be plotted along the y-axis
       :subplotnums     subplotnum for plotting multiple plots in one window
   '''
   if len(subplotnums) - 1 != len(np.atleast_1d(kaiserwindowparameter)):
      print "Bad subplotnum vs kaiserwindowparameter"
      return
   # First check the t array whether it has a constant dt
   dt = t[1] - t[0]
   for i in xrange(len(t)-1):
      if dt != t[i+1] - t[i]:
         print "Gave bad timestep to plot_fourier, the time step in array t must be constant (for now)"
   # Use kaiser window on y
   _y = y*np.kaiser(len(y), np.atleast_1d(kaiserwindowparameter)[0])
   # Plot the raw data as well as the fitted data (fourier series)
   pl.subplot(subplotnums[0][0],subplotnums[0][1],subplotnums[0][2])
   pl.plot(t,_y,'.',color='r')
   # Do FFT on the data
   fourier=np.fft.fft(_y)*(1/(float)(len(t)))
   # Get frequencies of the fourier
   freq=np.fft.fftfreq(len(fourier), d=dt)
   # Declare t2 (Note: This is the same as t but we want the steps to be thicker so the data looks smoother
   dt2=dt*0.01
   t2=np.arange(len(t)*100)*dt2
   # Declare y2
   y2=np.array([np.sum(fourier*np.exp(complex(0,1)*2*np.pi*freq*T)) for T in t2])
   pl.plot(t2,y2,'-',color='b')
   # Plot the frequency spectrums
   j = 0
   for i in subplotnums[1:]:
      pl.subplot(i[0],i[1],i[2])
      # Use window function on the data
      _y = y*np.kaiser(len(y), np.atleast_1d(kaiserwindowparameter)[j])
      # Do FFT on the data
      fourier=np.fft.fft(_y)*(1/(float)(len(t)))
      # Get frequencies of the fourier
      freq=np.fft.fftfreq(len(fourier), d=dt)
      # Get the indexes:
      toIndex = (int)((len(freq)/2)/2.0 + 1)
      #maxValue = np.max(np.abs(fourier[1:len(fourier)/2]))
      #while True:
      #   if toIndex <= 2 or np.abs(fourier[toIndex-1]) > 0.02*maxValue:
      #      break
      #   toIndex = toIndex - 1
      #pl.plot(freq[1:len(freq)/2],2*np.abs(fourier[1:len(fourier)/2]), marker='.', linestyle='-', linewidth=0.5)
      if savedata != "none":
         saveplotdata( freq[1:toIndex], np.log(2*np.abs(fourier[1:toIndex])), savedata + "_" + str(j) + ".npy" )
      pl.plot(freq[1:toIndex],2*np.abs(fourier[1:toIndex]), marker='.', linestyle='-', linewidth=0.5)
      pl.ylim([0,1.05*max(2*np.abs(fourier[1:len(fourier)/2]))])
      #pl.xlim([0, 0.13])
      #xTicks = np.arange(14)/100.0
      #pl.xticks(xTicks)
      # Put interactive mode on and show the plot:
      #pl.ion()
      #pl.show()
      j = j + 1
   #pl.show()
