import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like 

x = numpy.loadtxt('x.dat');  

u_init = numpy.loadtxt('u_initial.dat');
u_comp = numpy.loadtxt('u_final.dat');


pyplot.figure();

pyplot.plot(x,u_init,'--b',label='Initial sol');
pyplot.plot(x,u_comp,'-or',label='Numerical sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('x');
pyplot.ylabel('u');



pyplot.show()





