from tqdm.auto import tqdm
from time import sleep
from threading import Thread
import pandas as pd
import numpy as np
import subprocess
import os

# plotting with plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import sample_colorscale

LE_DIR = os.path.join(os.getcwd(), 'LATTICEEASY') # latticeeasy directory; by default LATTICEEASY in the same path
PLOTTER = "plotly"  # which plotting engine to use
MASS_UNIT = "Mpl"  # whether to use reduced/non-reduced Planck mass (Mpl/mpl)
mpl = np.sqrt(8.*np.pi)


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
            param_name = line[10:].split("=")[0].strip()
            comment_line = line.split(";")[1]
            if param_name in params.keys():
                lines[i] = "const int " + param_name + "=" + str(params[param_name]) + ";" + comment_line
        elif "const float" in line:
            param_name = line[12:].split("=")[0].strip()
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
    subprocess.run("rm *.txt" + output_to_txt, shell=True)
    subprocess.run("make all" + output_to_txt, shell=True)
    run_statement = f"Running LATTICEEASY in {NDIMS}D with "
    for key in sorted(params.keys()):
        run_statement += key + "=" + str(params[key]) + ", "
    print(run_statement + "...")
    subprocess.run(os.path.join(os.path.curdir, "latticeeasy" + output_to_txt), shell=True)
    
    # copy data into a new directory
    if save_data:
        param_string = param_string= ''.join(["_%s=%s" % (key, params[key]) for key in params.keys()])
        data_path = os.path.join(DATADIR, f'LE{NDIMS}D' + param_string)
        if not os.path.isdir(DATADIR):
            os.mkdir(DATADIR)
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        os.system('cp ./*.dat ' + data_path)
    
    # return to the original directory
    os.chdir(pwd)

def progress(tf, update_interval=2):
    """ Get progress by periodically reading the output file """
    
    with tqdm(total=100, desc="Progress") as pbar:
        
        file = os.path.join(LE_DIR, "output.txt")
        while True:

            # check that output file exists
            if not os.path.isfile(file):
                print("No file yet. Waiting...")
                sleep(update_interval)
                continue

            with open(file, "r") as f:
                lines = f.readlines()
                if lines:
                    last = lines[-1]
                    if "LATTICEEASY" in last:
                        pbar.n = 100
                        pbar.refresh()
                        break
                    pbar.n = int(100*float(last)/tf)
                    pbar.refresh()
                sleep(update_interval)

def run_with_progress(params, NDIMS=3, save_data=True):
    """ Run latticeeasy in a separate thread while monitoring progress in the main """

    tf = params["tf"]
    le_thread = Thread(target=run,args=(params,), kwargs={"output_to_txt": True, 'NDIMS': NDIMS, "save_data": save_data})
    le_thread.start()
    tf = params["tf"]
    sleep(2)
    progress(tf)
    #le_thread.join()

def get_rescalings(DIR=LE_DIR):
    """ Read rescaling values from info.dat """

    file_path = os.path.join(DIR, "info_0.dat")
    with open(file_path) as f:
        lines = f.readlines()
    
    rescalings = {}
    for line in lines:
        if "rescale_" in line:
            line = line.split('=')
            rescalings[line[0].strip()] = float(line[1])

    return rescalings

def get_fields(DIR=LE_DIR, rescale=1):
    """ Import fields and variances from corresponding .dat file. Provide rescaling for physical units; defaults to program units. """

    means = pd.read_csv(os.path.join(DIR,[name for name in os.listdir(DIR) if 'means' in name][0]), names=["t_pr","phi","chi"], sep=' ')
    means["phi"] = means["phi"]/rescale
    means["chi"] = means["chi"]*rescale

    variances = pd.read_csv(os.path.join(DIR,[name for name in os.listdir(DIR) if 'variance' in name][0]), names=["t_pr","phi2","chi2"], sep=' ')
    variances["phi2"] = variances["phi2"]/rescale**2
    variances["chi2"] = variances["chi2"]/rescale**2

    return (means, variances)

def plot_energies(energies, yscale='log', show=True, format=None):
    """ Plot energy densities """

    fig = go.Figure()
    time = energies["t_pr"]
    
    for component in energies.columns[1:]:
        fig.add_trace(go.Scatter(x=time, y=energies[component], name=component))
    
    fig.update_yaxes(exponentformat="power", type=yscale, title="energy density")
    fig.update_xaxes(title="t_pr")
    
    if show:
        fig.show(format)

    return fig

class Energies(pd.DataFrame):
    """ Object containing energy densities. Extends pandas dataframe with specific plotting. """

    def plot(self, format=None):
        plot_energies(self, format=format)

def get_energies(DIR=LE_DIR, names=''):
    """ Import energies """

    if not names:
        names = ["t_pr", "phi_kinetic", "chi_kinetic", "phi_gradient", "chi_gradient", "potential_1", "potential_2"]

    # import energy data if exists
    energy_files = [name for name in os.listdir(DIR) if 'energy' in name]
    if energy_files:
        energy  = pd.read_csv(os.path.join(DIR, energy_files[0]), sep=" ", names=names)
    else:
        energy = pd.DataFrame()

    return Energies(energy)

def plot_spectra(nk_data, t_min='all', t_max='all', time_skip=1, colorscale='balance', xscale='log', yscale='auto', title='', format=None):
    """ Plot spectra given a povot table of data"""
    
    # choose the initial and final times for the spectra
    t_min = nk_data.columns[0] if t_min=='all' else t_min
    t_max = nk_data.columns[-1] if t_max=='all' else t_max

    # create a curresponding subset of the spectra
    nk_slice = nk_data[nk_data.columns[t_min<=nk_data.columns]]
    nk_slice = nk_slice[nk_slice.columns[t_max>=nk_slice.columns]]
    if time_skip>1:
        nk_slice = nk_slice[nk_slice.columns[::time_skip]]
        
    # determine range and type of plot
    vmax = 1.5*nk_slice.max().max()
    vmin = 0
    if vmax > 100:
        yscale = 'log' if yscale == 'auto' else yscale
    else:
        yscale = 'linear' if yscale == 'auto' else yscale
    if yscale == 'log':
        vmax = np.log10(vmax)
        vmin = -.5


    if type(colorscale)==int:
            colorscale=["balance", "icefire", "jet"][colorscale]
            
    colors = sample_colorscale(colorscale,nk_slice.shape[1])
    fig = go.Figure()

    for i, tau in enumerate(nk_slice.columns):
        spectrum = nk_slice[tau]
        fig.add_trace(go.Scatter(x=spectrum.index, y=spectrum, line=dict(color=colors[i])))

    fig.update_yaxes(exponentformat="power", title="n_k", type=yscale, range=(vmin,vmax))
    fig.update_xaxes(type=xscale, title="k_pr")
    fig.update_layout(showlegend=False, title=title)
    fig.show(format)


class PowerSpectrum(pd.DataFrame):
    
    """ Spectrum object. Extends DataFrame by making it callable. Normally returns a pivot table of nk(tau, k).
        If called with no arguments returns the final spectrum. If called with a tau value, returns the spectrum n_k(k) for the nearest tau. If called with k value,
        returns the evolution of the occupation number n_k(tau) for the nearest k. If called with both,
        returns the closest value in both dimensions."""
    
    def __call__(self, t="end", k="full"):
        
        # by default return the final spectrum
        if t=="end" and k=="full":
            return self.iloc[:,-1]
        
        # if called with a tau value return nearest
        elif t != "end" and k=="full":
            
            for tvar in self.columns:
                
                if t <= tvar:
                    return self[tvar]
        
        # if called with k value 
        elif t=="end" and k != "full":
            
            for kvar in self.index:
                
                if kvar >= k:
                    return self.loc[kvar]

class Spectra:
    """ Object containing the spectra for both fields for plotting"""

    def __init__(self, spectra):

        self.nk_phi, self.nk_chi = spectra
    
    def plot(self, format=None, t_min='all', t_max='all', time_skip=1, colorscale='balance', xscale='log', yscale='auto'):
        
        plot_spectra(self.nk_phi, format=format, t_min='all', t_max='all', time_skip=1, colorscale='balance', xscale='log', yscale='auto', title="phi spectrum")
        plot_spectra(self.nk_chi, format=format, t_min='all', t_max='all', time_skip=1, colorscale='balance', xscale='log', yscale='auto', title="chi spectrum")



def get_spectra(DIR=LE_DIR):

    # import spectra data if exists
    spectra_files = sorted([file for file in os.listdir(DIR) if "spectra" in file])
    if spectra_files:
        phi_spectra = pd.read_csv(os.path.join(DIR,spectra_files[0]), sep=" ", names = ["k", "kpoints", "omega_k", "|f_k|^2", "|df_k|^2", "n_k", "rho_k"])
        chi_spectra = pd.read_csv(os.path.join(DIR,spectra_files[1]), sep=" ", names = ["k", "kpoints", "omega_k", "|f_k|^2", "|df_k|^2", "n_k", "rho_k"])
        spectra_times = np.loadtxt(os.path.join(DIR,spectra_files[2]))
        spectr_num = spectra_times.shape[0]
        k_num = chi_spectra.shape[0]//spectr_num
        time_vals = np.concatenate([[tpr]*k_num for tpr in spectra_times])
        chi_spectra.insert(1, "t_pr", time_vals)
        phi_spectra.insert(1, "t_pr", time_vals)
        nk_phi = phi_spectra.pivot_table(index="k", columns="t_pr", values="n_k")
        nk_chi = chi_spectra.pivot_table(index="k", columns="t_pr", values="n_k")
        nk_phi = PowerSpectrum(nk_phi)
        nk_chi = PowerSpectrum(nk_chi)

    else:
        nk_phi = pd.DataFrame([])
        nk_chi = pd.DataFrame([])        

    return (nk_phi, nk_chi)

class Slices:

    def __init__(self, DIR=LE_DIR):
        
        self.phi_slices, self.chi_slices, self.slicetimes = get_slices(DIR=DIR)

    def plot(self, t=-1, format=None):

        fig = make_subplots(1, 2, horizontal_spacing=0.2)
        fig.add_trace(go.Heatmap(z=self.phi_slices[t], colorbar_x=0.45, colorscale='viridis'), 1, 1)
        fig.add_trace(go.Heatmap(z=self.chi_slices[t], colorscale= 'RdBu'), 1, 2)
        if format:
            fig.update_layout(width=1000, height=500)
        fig.show(format)

def get_slices(DIR=LE_DIR):

    # import slice data for the fields
    phi_slices = np.loadtxt(os.path.join(DIR, "slices0_0.dat"))
    chi_slices = np.loadtxt(os.path.join(DIR, "slices1_0.dat"))

    # import slice time data
    slicetimes = np.loadtxt(os.path.join(DIR, "slicetimes_0.dat"))
    side = int(np.sqrt(len(phi_slices)/len(slicetimes)))

    # reshape
    phi_slices = phi_slices.reshape(len(slicetimes), side, side)
    chi_slices = chi_slices.reshape(len(slicetimes), side, side)
    
    return phi_slices, chi_slices, slicetimes

def compare_models(models, legend='', legend_title='', title='', returns=False, format=None):
    
    if legend == '':
        legend = [f"model {i}" for i in range(len(models))]
    
    fig1 = go.Figure()
    fig2 = go.Figure()
    
    # plot mean phi and <chi^2> for all models
    for i, model in enumerate(models):
    
        fig1.add_trace(go.Scatter(x=model.means['t_pr'], y=model.means['phi'], name=legend[i]))
        fig2.add_trace(go.Scatter(x=model.variances['t_pr'], y=model.variances['chi2'], name=legend[i]))

    fig1.update_yaxes(exponentformat="power")
    fig1.update_layout(legend_title_text=legend_title, title=title, xaxis_title="t_pr", yaxis_title="phi")
    fig1.show(format)
    
    fig2.update_layout(legend_title_text=legend_title, title=title, xaxis_title="t_pr", yaxis_title="<chi^2>")
    fig2.update_yaxes(exponentformat="power", type="log")
    fig2.show(format) 

    if returns:
        return [fig1, fig2]


class Run:

    """ Run object """

    # get the lattice dimension
    def __init__(self, parameters, DIR=LE_DIR, save_LE_data=False, spectra=True, slices=False):

        # if lattice dimension not given default to two
        if 'NDIMS' not in parameters:
            parameters['NDIMS'] =2

        parameters["sspectra"] = int(spectra)
        parameters["sslices"] = int(slices)

        # save parameters
        self.parameters = parameters
        for param, value in parameters.items():
            self.__dict__[param] = value


        # run latticeeasy
        run_with_progress(parameters, NDIMS=self.NDIMS, save_data=save_LE_data)

        # import rescalings
        rescalings = get_rescalings(DIR=DIR)
        self.rescale_A = rescalings["rescale_A"]
        self.rescale_B = rescalings["rescale_B"]
        self.rescale_r = rescalings["rescale_r"]
        self.rescale_s = rescalings["rescale_s"]

        # import field data
        self.means, self.variances = get_fields(DIR=DIR, rescale=self.rescale_A/mpl)

        # import energies
        self.energies = get_energies(DIR=DIR)

        # import spectra
        if spectra:
            self.spectra = Spectra(get_spectra(DIR=DIR))

        # import slices
        if slices:
            self.slices = Slices(DIR=DIR)

    def plot(self,format=None,**kwargs):

        compare_models([self], format=format, **kwargs)

    
