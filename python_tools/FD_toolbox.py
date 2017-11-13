from matplotlib import pyplot as plt    #and the useful plotting library
from matplotlib import ticker

import numpy as np
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp,sqrt\
    , shape, empty_like , size, loadtxt, savetxt, arange , log , hstack, vstack, trapz, append, savez
    
import argparse
from decimal import Decimal
import csv

from subprocess import call, run, PIPE, Popen
from sys_cmd_toolbox import system_process    # locally defined file

plt.rc('legend',**{'loc':'upper right'});
plt.rcParams[u'legend.fontsize'] = 15
plt.rcParams[u'legend.edgecolor']='white'
plt.rcParams[u'font.weight']='normal'
#plt.rcParams['font.serif']='false'
plt.rcParams[u'xtick.labelsize']=14
plt.rcParams[u'ytick.labelsize']=14
plt.rcParams[u'axes.titlesize']=18
plt.rcParams[u'axes.labelsize']=15
plt.rcParams[u'axes.spines.right']='false';
plt.rcParams[u'axes.spines.top']='false';
plt.rcParams[u'lines.linewidth'] = 1.5;
plt.rcParams[u'lines.markersize'] = 8;



def FD_input_reader(input_file_):
    with open(input_file_) as file: 
        reader=csv.reader(file, delimiter=':');
        
        filter_flag = 0   # by default no filter is used
        for row in reader:
            if row[0] =='dir':    
                dir_result = str(row[1])
            elif row[0] =='errors':   
                errors_dir = str(row[1])
            elif row[0] == 'exact':
                exact_sol_dir = str(row[1])
            elif row[0] == 'numerical':
                num_sol_dir = str(row[1])
            elif row[0] =='unsteady_num':   
                sol_dir = str(row[1])
            elif row[0]=='Eqn_set':
                eqn_set=str(row[1]);
            elif row[0]=='Eqn_type':
                eqn_type=str(row[1]);
            elif row[0]=='wave_form':
                wave_form=str(row[1]);
            elif row[0]=='mode':
                mode=str(row[1]);
            elif row[0]=='scheme_type':
                scheme_type=str(row[1]);
            elif row[0]=='FDOA':
                FD=str(row[1])
                if int(FD)==2:
                    FDOA = FD + str('nd')
                elif int(FD)==3:
                    FDOA = FD + str('rd')
                else:
                    FDOA = FD + str('th')
            elif row[0]=='filter_flag':
                filter_flag=int(row[1]);
            elif row[0]=='filter_type':
                filter_type=str(row[1]);
            elif row[0]=='filterOA':
                filOA=str(row[1]);
            elif row[0]=='filter_alpha':
                alpha_filter=str(row[1]);
            elif row[0] == 'RK':    
                RK=str(row[1]);
            elif row[0] == 'Nelem':    
                Nelem=str(int(row[1]));
            elif row[0] == 'Nexact':    
                Nexact=str(int(row[1]));
            elif row[0]=='CFL':
                CFL=Decimal(row[1])
            elif row[0] == 'dt':
                dt_=str(row[1]);
            elif row[0] == 't_end':
                t_end=Decimal(row[1]);
            elif row[0] == 'T':     
                T=Decimal(row[1]);
          
    if(scheme_type=='implicit'):
        if(filter_flag==1):
            scheme_name = 'C' + FD + 'F'+ filOA + r'$^{'+str(alpha_filter)+'}$' + ' and RK'+ RK
        elif(filter_flag==0):
            scheme_name = 'CD' + FD + '-RK'+ RK
    elif(scheme_type=='explicit'):
        if(filter_flag==1):
            scheme_name = 'F' + FD + 'F'+ filOA + r'$^{'+str(alpha_filter)+'}$' + ' and RK'+ RK
        elif(filter_flag==0):
            scheme_name = 'FD'+FD + '-RK'+ RK

    CFL = Decimal(CFL.quantize(Decimal('0.0001')))
    t_end = Decimal(t_end.quantize(Decimal('0.001')))
    
    cmd=['mkdir',dir_result+'tempfig/']
    cmd_out,cmd_err=system_process(cmd,1000)
                
    return dir_result, mode, eqn_type, scheme_name, FD, RK, Nelem, CFL, dt_, sol_dir, exact_sol_dir
  
def load_sol(fname):
    data= loadtxt(fname);
    x_ = data[:,0];
    u_ = data[:,1]; 
    
    return x_,u_   

def decay_burg_turb_temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, tt_, sol_dir):

    if mode=='dt_const':
        sim_name = str('_dt') + dt_
        sim_name1 = str('dt=') + dt_
    else:
        sim_name = str('_CFL') + str(CFL)
        sim_name1 = str('CFL=') + str(CFL)
	
    fname = dir_res+sol_dir+str('_N')+Nelem+sim_name+str('_')+str(tt_)+str('t.dat')
    x_num, u_num = load_sol(fname);

    fig = plt.figure();

    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_num));
    ylim_1.append(max(u_num));

    if mode=='dt_const':
	    label_ = scheme_name\
	    +" with dt="+ '{:1.2e}'.format(dt)+ ', t='+str(tt_);
    else:
	    label_ = scheme_name \
	    + " with CFL="+ str(CFL)+ ', t ='+str(tt_);
	    
    plt.plot(x_num,u_num,'-r',label=label_); 

    plt.legend();

    plt.xlabel('x',labelpad=10);
    plt.ylabel('u(x)',labelpad=10);

    plt.grid()
    plt.legend()
    fig.tight_layout()
    figname = dir_res + 'tempfig/' + 'contsol_N'+Nelem+sim_name+str('_')+str(tt_)+str('t.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname)
    plt.show()
    
    return
    
def temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, tt_, sol_dir, exact_sol_dir):

    if mode=='dt_const':
        sim_name = str('_dt') + dt_
        sim_name1 = str('dt=') + dt_
        label_ = scheme_name+" with dt="+ '{:1.2e}'.format(dt)+ ', t='+str(tt_);
    else:
        sim_name = str('_CFL') + str(CFL)
        sim_name1 = str('CFL=') + str(CFL)
        label_ = scheme_name+ " with CFL="+ str(CFL)+ ', t ='+str(tt_);

    fname = dir_res+sol_dir+str('_N')+Nelem+sim_name+str('_')+str(tt_)+str('t.dat')
    x_num, u_num = load_sol(fname);
    
    fname = dir_res+'time_data/u_exact_'+str(tt_)+str('t.dat')
    x_exact, u_exact = load_sol(fname);

    fig = plt.figure();
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_num));
    ylim_1.append(max(u_num));
    print('u_max: ',max(u_num));
    plt.plot(x_num,u_num,'-or',label=label_); 
    plt.legend();
    plt.xlabel('x',labelpad=10);
    plt.ylabel('u(x)',labelpad=10);
    plt.grid()
    plt.legend()
    fig.tight_layout()
    figname = dir_res + 'tempfig/' + 'contsol_N'\
    +Nelem+sim_name+str('_')+str(tt_)+str('t.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname)
    plt.show()
    
    fig = plt.figure();
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_exact));
    ylim_1.append(max(u_exact));
    plt.plot(x_num,u_num,'-or',label=label_); 
    plt.plot(x_exact,u_exact,'--k',label='exact solution'); 
    plt.legend();
    plt.xlabel('x',labelpad=10);
    plt.ylabel('u(x)',labelpad=10);
    plt.xlim(min(x_num), max(x_num))
    plt.ylim(min(ylim_0), max(ylim_1)*1.3)
    plt.grid()
    plt.legend()
    fig.tight_layout()
    figname = dir_res + 'tempfig/' + 'contsol_N'\
    +Nelem+sim_name+str('_')+str(tt_)+str('t_comp.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname)
    plt.show()
    
    return

