import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like, size 

x = numpy.loadtxt('../output/x.dat');
 

u_init = numpy.loadtxt('../output/u_initial.dat');
u_comp = numpy.loadtxt('../output/u_final.dat');

pyplot.figure();

pyplot.plot(x,u_init,'--k',label='Initial sol');
pyplot.plot(x,u_comp,'-or',label='Numerical sol');

pyplot.show();
