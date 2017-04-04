import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
<<<<<<< HEAD
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
#pyplot.title('comparison of all solutions');


pyplot.show()






=======
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like, size 

x = numpy.loadtxt('../output/x.dat');
 

u_init = numpy.loadtxt('../output/u_initial.dat');
u_comp = numpy.loadtxt('../output/u_final.dat');

pyplot.figure();

pyplot.plot(x,u_init,'--k',label='Initial sol');
pyplot.plot(x,u_comp,'-or',label='Numerical sol');

pyplot.show();
>>>>>>> master
