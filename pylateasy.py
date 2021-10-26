from tqdm.auto import tqdm
from time import sleep
import os
import subprocess

LE_DIR = os.path.join(os.getcwd(), 'LATTICEEASY') # latticeeasy directory; by default LATTICEEASY in the same path
PLOTTER = "plotly"  # which plotting engine to use

class plotting:
    """ Change plotting parameters """

    @staticmethod
    def use(plotter):

        global PLOTTER
        PLOTTER = plotter

def write_param_file(params, DIR, NDIMS):
    """ updates the parameter values in parameters.h with values in params (given as a dictionary) """

    param_path = os.path.join(DIR, "parameters.h")

    # read the parameters.h file 
    with open(param_path, "r") as f:
        lines = f.readlines()

    # update parameters    
    for i,line in enumerate(lines):
        if "const int" in line:
            param_name = line[10:].split("=")[0]
            comment_line = line.split(";")[1]
            if param_name in params.keys():
                lines[i] = "const int " + param_name + "=" + str(params[param_name]) + ";" + comment_line
        elif "const float" in line:
            param_name = line[12:].split("=")[0]
            comment_line = line.split(";")[1]
            if param_name in params.keys():
                lines[i] = "const float " + param_name + "=" + str(params[param_name]) + ";" + comment_line
        elif "#define NDIMS" in line:
            lines[i] = "#define NDIMS " + str(NDIMS) + '\n'
    
    # write the new parameters.h file
    f = open(param_path, "w")
    for line in lines:
        f.write(line)
    f.close()


def run(params, DIR=LE_DIR, NDIMS=3, DATADIR='', save_data=True, output_to_txt=False):
    """ Runs latticeeasy; params= dictionary of parameter values"""

    pwd = os.getcwd()
        
    # directory where to put data files; by assumption ./LE_DATA in the same directory
    if DATADIR == '':
        DATADIR = os.path.join(pwd, 'LE_DATA')
        
    write_param_file(params, DIR, NDIMS)
    
    # Go to the LATTICEEASY directory, compile and run LATTICEEASY
    output_to_txt = " > out.txt" if output_to_txt else ""
    os.chdir(DIR)
    subprocess.run("make cleaner" + output_to_txt, shell=True)
    subprocess.run("make all" + output_to_txt, shell=True)
    run_statement = f"Running LATTICEEASY in {NDIMS}D with "
    for key in sorted(params.keys()):
        run_statement += key + "=" + str(params[key]) + ", "
    print(run_statement + "...")
    subprocess.run(os.path.join(os.path.curdir, "latticeeasy" + output_to_txt), shell=True)
    
    # move data into a new directory
    if save_data:
        print("Moving data...")
        param_string = param_string= ''.join(["_%s=%s" % (key, params[key]) for key in params.keys()])
        data_path = os.path.join(DATADIR, f'LE{NDIMS}D' + param_string)
        if not os.path.isdir(DATADIR):
            os.mkdir(DATADIR)
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        os.system('cp ./*.dat ' + data_path)
    
    # return to the original directory
    os.chdir(pwd)
    print('Done!')

def progress(tf, update_interval=5):
    """ Get progress by periodically reading the output file """
    
    with tqdm(total=100, desc="Progress") as pbar:
        
        file = os.path.join(LE_DIR, "output.txt")
        while True:
            with open(file, "r") as f:
                last = f.readlines()[-1]
                if "LATTICEEASY" in last:
                    pbar.n = 100
                    pbar.refresh()
                    break
                pbar.n = int(100*float(last)/tf)
                pbar.refresh()
                sleep(update_interval)

class RunTemplate:
    """ Generic blueprint for a run  object """
    pass

class Run(RunTemplate):

    """ Run object """

    def __init__(self, parameters):

        pass