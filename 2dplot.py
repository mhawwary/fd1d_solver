import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like 

data_comp = numpy.loadtxt('./Results/FD2nd_RK3_test/nodal/u_nodal_N20_CFL1.73_150000.0T.dat');
x_comp = data_comp[:,0];
u_comp = data_comp[:,1];
u_exact = data_exact = numpy.loadtxt('./Results/FD2nd_RK3_test/nodal/u_nodal_exact_N20_CFL1.73_150000.0T.dat');
x_exact = data_exact[:,0];
u_exact = data_exact[:,1];

pyplot.figure();

pyplot.plot(x_exact,u_exact,'--k',label='Exact sol');
pyplot.plot(x_comp,u_comp,'-or',label='Numerical sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('x');
pyplot.ylabel('u');

pyplot.show()






