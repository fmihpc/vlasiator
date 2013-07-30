#Doing fourier series here
import numpy as np
import pylab as pl

dt=1/10.0
t=np.arange(100)*dt
y=np.cos(2*np.pi*t) + np.sin(2*np.pi*3*t) + np.cos(2*np.pi*5*t)
fourier=np.fft.fft(y)*(1/(float)(len(t)))
freq=np.fft.fftfreq(len(fourier), d=dt)
dt2=1/1000.0
t2=np.arange(10000)*dt2
y2=np.array([np.sum(fourier*np.exp(complex(0,1)*2*np.pi*freq*T)) for T in t2])
pl.subplot(2,1,1)
pl.plot(t,y,'.',color='r')
pl.plot(t2,y2,'-',color='b')
#pl.xlim([-0.05,1.05])
#pl.ylim([min(y),max(y)])
pl.subplot(2,1,2)
pl.plot(freq[1:len(freq)/2],2*np.abs(fourier[1:len(fourier)/2]), '.')
pl.ylim([0,1.05*max(2*np.abs(fourier[1:len(fourier)/2]))])
pl.show()
