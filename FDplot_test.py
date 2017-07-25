from matplotlib import pyplot ,ticker   #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size, loadtxt, arange , log
import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_FD_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='mode':
            mode=str(row[1]);
        elif row[0]=='CFL':    
            CFL=Decimal(row[1]);
        elif row[0]=='FDOA':
            FD=str(row[1]);
        elif row[0] == 'RK':    
            RK=str(row[1]);
        elif row[0] == 'Nelem':    
            Nelem=str(int(row[1]));
        elif row[0] == 'Nexact':    
            Nexact=str(int(row[1]));
        elif row[0]=='dt':
            dt_=str(row[1]);
        elif row[0] == 'T':     
            T=Decimal(row[1]);
        elif row[0] =='dir':    
            dir1 = str(row[1]);
        elif row[0] == 'exact': 
            sol_exact = str(row[1]);
        elif row[0] == 'numerical': 
            sol_comp = str(row[1]);

CFL=Decimal(CFL.quantize(Decimal('.001')));
T=Decimal(T.quantize(Decimal('.001')));
dt = float(dt_);

fname_exact = dir1+sol_exact+"_"+str(T)+str("T.dat");
data_exact= loadtxt(fname_exact);

#if mode=='CFL_const':
#	fname_comp = dir1+sol_comp+str("_N")+Nelem\
#	+str("_CFL")+str(CFL)+str("_")+str(T)+str("T.dat");
#	data_num= loadtxt(fname_comp);

if mode=='dt_const':
	fname_comp = dir1+sol_comp+str("_N")+Nelem\
	+str("_dt")+dt_+str("_")+str(T)+str("T.dat");
	data_num= loadtxt(fname_comp);
else:
	fname_comp = dir1+sol_comp+str("_N")+Nelem\
	+str("_CFL")+str(CFL)+str("_")+str(T)+str("T.dat");
	data_num= loadtxt(fname_comp);

x_exact = data_exact[:,0]; 
u_exact = data_exact[:,1]; 

x_num = data_num[:,0];
u_num = data_num[:,1];

pyplot.rc('legend',**{'loc':'upper left'});
pyplot.rcParams[u'legend.fontsize'] = 15
pyplot.rcParams[u'legend.edgecolor']='white'
pyplot.rcParams[u'font.weight']='normal'
#pyplot.rcParams['font.serif']='false'
pyplot.rcParams[u'xtick.labelsize']=14
pyplot.rcParams[u'ytick.labelsize']=14
pyplot.rcParams[u'axes.titlesize']=18
pyplot.rcParams[u'axes.labelsize']=15
pyplot.rcParams[u'axes.spines.right']='false';
pyplot.rcParams[u'axes.spines.top']='false';
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;

fig = pyplot.figure();

ylim_0 = list();
ylim_1 = list();
ylim_0.append(min(u_exact));
ylim_1.append(max(u_exact));
ylim_0.append(min(u_num));
ylim_1.append(max(u_num));

pyplot.plot(x_exact,u_exact,'--k',label='Exact sol'); 
pyplot.plot(x_num,u_num,'-r',label='Numerical sol'); 

pyplot.legend();

if int(FD)==1:
    FDOA = str("1st order upwind");
elif int(FD)==2:
    FDOA = str("2nd order central");
elif int(FD)==3:
    FDOA = str("3rd order biased upwind");
elif int(FD)==4:
    FDOA = str("4th order central");
elif int(FD)==6:
    FDOA = str("6th order central"); 

if mode=='dt_const':
	title_a = str("FD")+ FDOA +  " RK"+ RK\
	+" with dt="+ '{:1.2e}'.format(dt)+ ' at t/T='+str(Decimal(T.quantize(Decimal('.01'))));
else:
	title_a = str("FD ")+ FDOA + ", RK"+ RK \
	+ " with CFL="+ str(CFL)+ ' and at t/T='+str(Decimal(T.quantize(Decimal('.01'))));

pyplot.title(title_a);
pyplot.xlabel('X',labelpad=10);
pyplot.ylabel('u(x)',labelpad=10);

pyplot.xlim(min(x_exact),max(x_exact));
pyplot.ylim(min(ylim_0)*1.05,max(ylim_1)*1.05);

xlabels=linspace(min(x_exact),max(x_exact),5);
xlocs=xlabels;
pyplot.xticks(xlocs, xlabels);

ylabels=arange(0,1.2,0.2);
ylocs=arange(0,1.2,0.2);
pyplot.yticks(ylocs,ylabels);

fig.tight_layout()
pyplot.show()










