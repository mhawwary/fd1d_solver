import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_FD_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input')
parser.add_argument('-t', type=str, dest='time_for_plot');

args = parser.parse_args();


from FD_toolbox import FD_input_reader, decay_burg_turb_temp_sol_plot, temp_sol_plot

dir_res, mode, eqn_type, scheme_name, FD, RK, Nelem, CFL, dt_, sol_dir, exact_sol_dir = FD_input_reader(args.python_input)

t_plot = Decimal(args.time_for_plot)
t_plot =  Decimal(t_plot.quantize(Decimal('.001')))

if eqn_type=='linear_advec':
    temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, t_plot, sol_dir, exact_sol_dir)
else:
    decay_burg_turb_temp_sol_plot(dir_res, mode, scheme_name, FD, RK, Nelem, CFL, dt_, t_plot, sol_dir)
