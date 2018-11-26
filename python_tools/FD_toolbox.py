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
from fft_toolbox_python import load_data, compute_fft, compute_Etotal, compute_exact_init_Energy

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
    elif(scheme_type=='DRP4s7'):
        if(filter_flag==1):
            scheme_name = 'DRP4s7'+ 'F'+ filOA + r'$^{'+str(alpha_filter)+'}$' + ' and RK'+ RK
        elif(filter_flag==0):
            scheme_name = 'DRP4s7'+ '-RK'+ RK
    elif(scheme_type=='Rem2s7'):
        if(filter_flag==1):
            scheme_name = 'Rem2s7'+ 'F'+ filOA + r'$^{'+str(alpha_filter)+'}$' + ' and RK'+ RK
        elif(filter_flag==0):
            scheme_name = 'Rem2s7'+ '-RK'+ RK

    CFL = Decimal(CFL.quantize(Decimal('0.0001')))
    t_end = Decimal(t_end.quantize(Decimal('0.001')))
    T = Decimal(T.quantize(Decimal('0.001')))
    
    cmd=['mkdir',dir_result+'tempfig/']
    cmd_out,cmd_err=system_process(cmd,1000)
                
    return dir_result, mode, eqn_type, scheme_name, FD, RK, Nelem, CFL, dt_, sol_dir, exact_sol_dir,T
  
def load_sol(fname):
    data= loadtxt(fname);
    x_ = data[:,0];
    u_ = data[:,1]; 
    
    return x_,u_   

def decay_burg_turb_temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, tt_, sol_dir):

    #from fft_toolbox_python import load_data, compute_fft

    dt =float(dt_);
    
    if mode=='dt_const':
        sim_name = str('_dt') + dt_
        sim_name1 = str('dt=') + dt_
    else:
        sim_name = str('_CFL') + str(CFL)
        sim_name1 = str('CFL=') + str(CFL)
	
    fname = dir_res+sol_dir+str('_N')+Nelem+sim_name+str('_')+str(tt_)+str('t.dat')
    x_num, u_num = load_sol(fname);
    
    # compute fft:
    k_freq, u_amp, KEnerg = compute_fft(u_num);
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()

    del fname

    ###################### Plot the solution ##########################
    fig, ax = plt.subplots(frameon='True')

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
    #fig.tight_layout()
    
    fig.set_size_inches(13.0, 9.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0,rect=(0.0,0.0,1.0,0.985))
    
    temp_name = 'sol_vs_x_p' + FD + 'RK' + RK +'_Nn' + str(Nelem) + '_'+ sim_name \
              + '_t'+ str(tt_);
           
    figname = dir_res + str('tempfig/') + temp_name+'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_res + str('tempfig/') + temp_name+'.png'
    plt.savefig(figname,format='png',bbox='tight')
    
    #figname = dir_res + 'tempfig/' + 'contsol_N'+Nelem+sim_name+str('_')+str(tt_)+str('t.png')
    #fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname)
    
    ###################### Plot the fft computed Energy spectrum ##########################
    fig, ax = plt.subplots(frameon='True')
    ax.plot(k_freq, KEnerg,'--b', label=r'E$_{num}$',lw=0.7)
    ax.plot(k_init, E_init, ':k', label=r'E$_{init}$',lw=0.7)
    ax.plot(k_ex, E_ex*0.08, '-k', label=r'E $\propto$ k$^{-2}$',lw=0.7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=19,edgecolor='black')
    ax.set_facecolor('white')
    ax.set_xlabel(r'$k$', labelpad=2,fontsize=24);
    ax.set_ylabel(r'$E$', labelpad=2, fontsize=24);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-10, 10 ** -1)
    ax.tick_params(axis='both', which='both', labelsize=24)
    
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    
    fig.set_size_inches(13.0, 9.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0,rect=(0.0,0.0,1.0,0.985))
    
    temp_name = 'KE_FD' + FD + 'RK' + RK +'_Nn' + str(Nelem) +'_'+ sim_name \
              + '_t'+ str(tt_);
    figname = dir_res + str('tempfig/') + temp_name+'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_res + str('tempfig/') + temp_name+'.png'
    plt.savefig(figname,format='png',bbox='tight')
    
    plt.show()
    
    return
    
def temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, tt_, sol_dir, exact_sol_dir,T_period):

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
    k_freq, u_amp, KE = compute_fft(u_num)
    E_tot = compute_Etotal(k_freq,KE)

    print('u_max: ',max(u_num));
    print('error: ',abs(1-max(u_num)));
    
    fname = dir_res+'time_data/u_exact_'+str(tt_)+str('t.dat')
    x_exact, u_exact = load_sol(fname);
    k_freq_exact, u_amp_exact, KE_exact = compute_fft(u_exact)
    E_ex = compute_Etotal(k_freq_exact,KE_exact)

    fig = plt.figure();
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(1.02*min(u_num));
    ylim_1.append(1.02*max(u_num));
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
    #plt.show()
    
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
    # for gaussian
    #if min(ylim_0)<=1e-5:
     # ylim_0.append(-0.052632);
    #plt.ylim(0.95*min(ylim_0), max(ylim_1)*1.2)
    #for sine wave
    plt.ylim(1.05*min(ylim_0), max(ylim_1)*1.05)
    plt.grid()
    plt.legend()
    fig.tight_layout()
    figname = dir_res + 'tempfig/' + 'contsol_N'\
    +Nelem+sim_name+str('_')+str(tt_)+str('t_comp.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname)
    
    # Read DG data:
    res_dir = './Results/'
    fname = dir_res + 'errors/errors_N'+str(Nelem)+'_CFL'+str(CFL)+'_'+str(T_period)+'T.dat'
    print('fname, ',fname)
    data = loadtxt(fname);  # continuous exact nodal solution
    time  = data[:, 0];
    L1err = data[:, 1];
    L2err = data[:, 2];

    fig = plt.figure();
    
    plt.plot(time,L1err,'-.b',label=r'L$_{1}$')
    plt.plot(time, L2err,'-ok',label=r'L$_{2}$')
    plt.xlabel('time')
    plt.ylabel('Errors')
    plt.legend(loc='upper left')
    
    figname = dir_res + 'tempfig/'+ 'errors_N'+Nelem\
    +sim_name+str('_')+str(T_period)+str('T.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir_res + 'tempfig/'+ 'errors_N'+Nelem\
    +sim_name+str('_')+str(T_period)+str('T.eps')
    plt.savefig(figname,format='eps')
    
    #================== Plot FFT ==========================#
    k_max_ = int(k_freq[-1]);
    print('k_max: ',k_max_, '  k_max_ex: ',int(k_freq_exact[-1]))
    print('exact Nelem:',size(x_exact))
    fig, ax = plt.subplots()
    
    plt.plot(2*pi*k_freq_exact/int(Nelem), u_amp_exact, '-.k',markevery=1\
        , label=r'exact, E_${tot}$= '+str(np.round(E_ex,4)))
    plt.plot(2*pi*k_freq/int(Nelem), u_amp, '-or',markevery=1, \
        label=label_+r', E$_{tot}$='+str(np.round(E_tot,4)))
    
    xlabels = ['0',r'$\pi$/8',r'$\pi$/4',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/8,pi/4,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    plt.legend();
    plt.xlabel('K/(P+1)', labelpad=2);
    plt.ylabel('|u|', labelpad=2);
    plt.grid()

    fig.set_size_inches(13.0, 10.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir_res + 'tempfig/'+ 'fft_N'+Nelem\
    +sim_name+str('_')+str(tt_)+str('t_comp.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir_res + 'tempfig/'+ 'fft_N'+Nelem\
    +sim_name+str('_')+str(tt_)+str('t_comp.eps')
    plt.savefig(figname,format='eps')
    
    plt.show()
    
    return

