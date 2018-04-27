
import argparse
from decimal import Decimal
import csv
from time import sleep
from subprocess import call, run, PIPE, Popen
import fileinput
from numpy import size

from sys_cmd_toolbox import system_process

parser = argparse.ArgumentParser(description='FD_argparsing');
parser.add_argument('-f', type=str, dest='inputfname')
parser.add_argument('-t', type=str, dest='time_for_plot');
args = parser.parse_args();

elem_num = [255]  

fname = args.inputfname

# Computing and Running the first no. of elements :
#--------------------------------------------------

for lines in fileinput.input(fname, inplace=True):
    a = lines.split('=')
    b = a[0].strip()
    if b == 'spectrum_restart_flag':
        print('   {} = {}'.format(a[0].strip(),0))
    elif b == 'num_elements':
        print(' {} = {}'.format(a[0].strip(),elem_num[0]))
    elif b == 'title':
        sim_title = a[1].strip()
        print(' {} = {}'.format(b,a[1].strip()))
    elif b == 'order':
        OA = a[1].strip()
        print('  {} = {}'.format(b,a[1].strip()))
    elif b == 'RK_order':
        Prk = a[1].strip()
        print('   {} = {}'.format(b,a[1].strip()))
    else:
       print(lines.rstrip())

log_name_dir = './Results_AdvecDiffus/'+ str(sim_title)+'/FD' + str(OA) + 'RK' + str(Prk) +'/case'
print('input fname:  ',fname,' .............')
print('title: ',sim_title,'    .............')
print('log_name_dir: ',log_name_dir,' .............')

print('\n==============================================\nnum of elements: ',elem_num[0],'\n==============================================')  
        
for i in range(1,65):
    case_index = i   
    if case_index < 10:
        case_index_s = str('0')+str(i)
    else:
        case_index_s = str(i)
     
    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'case_no':
            print(' {} = {}'.format(a[0].strip(),case_index_s))
        else:
            print(lines.rstrip())

    cmd = ['./bin/FD1DFlow.exe', fname]
    process = Popen(cmd, stdin=PIPE, stdout=PIPE,start_new_session=True)

    try:
        outs, errs = process.communicate(timeout=5000)
        if not(outs is None):
            output = outs.decode("utf-8")
            #print(output)

    except TimeoutExpired:
        process.kill()
        outs, errs = process.communicate()
        if not(outs is None):
            output = outs.decode("utf-8")
          #  print(output)
           
    process.terminate()

    log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
    if not(outs is None):
        log = open(log_name,'a+')
        log.writelines(output)
    
    print('case_no: ',case_index_s)
