import numpy as np
import pandas as pd
import sqlite3 as sql
import sys
import matplotlib

#matplotlib.use('pdf')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import warnings

from astropy import units as u
from astropy.modeling import models

import pickle as pk
from pathlib import Path
import argparse

from atm import modifyErrors
from atm import multiFit
from atm.models import NEATM
from atm.obs import WISE, WISE_ACT
from atm.analysis import calcMagChi2
from atm.analysis import calcMagReducedChi2
from atm.functions import calcFluxLambdaSED, calcFluxLambdaAtObsWithSunlight
from atm.plotting import plotObservations
from atm.plotting import plotSED
from atm.helpers import __handleParameters as handleParameters
from atm.functions.hg import calcQ

#Surpresses a warming having to do with assignments around line 185
pd.options.mode.chained_assignment = None







def get_desig(id_num):
    home = str(Path.home())
    with open(home+'/dev/minorplanets/asteroids.pk', 'rb') as f:
        df = pk.load(f)
        name = df['name'][id_num]
        desig = df['designation'][id_num]
        semimajor = df['semimajor'][id_num]
    return desig, name, semimajor

class asteroid():
    def __init__(self, asteroid_name, desig, act_flux_dict, semimajor = False, show_plot = False, obs = WISE(), atm_dir = '/home/r/rbond/jorlo/dev/atm/', 
            save_path = '/scratch/r/rbond/jorlo/actxminorplanets/sigurd/', plot_path = '/scratch/r/rbond/jorlo/actxminorplanets/sigurd/plots/'):
        self.name = asteroid_name
        self.desig = desig
        self.act_flux_dict = act_flux_dict
        self.path = atm_dir
        self.save_path = save_path
        self.plot_path = plot_path
        self.show_plot = show_plot
        self.semimajor = semimajor

        #Put your atm_data dir here
        con = sql.connect(atm_dir +"atm_data/paper1/sample.db")
        self.observations = pd.read_sql("""SELECT * FROM observations""", con)
        self.additional = pd.read_sql("""SELECT * FROM additional""", con)

        # Only keep clipped observations
        self.observations = self.observations[self.observations["keep"] == 1]
        self.additional = self.additional[self.additional["obs_id"].isin(self.observations["obs_id"].values)]

        # Remove missing H value, G value objects... 
        self.observations = self.observations[~self.observations["designation"].isin(['2010 AJ104', '2010 BM69', '2010 DZ64', '2010 EL27', '2010 EW144',
           '2010 FE82', '2010 FJ48', '2010 HK10', '2010 LE80'])]

        # Convert phase angle to radians
        self.observations["alpha_rad"] = np.radians(self.observations["alpha_deg"])
        self.ran_override = False

        # Initialize observatory 
        self.obs = obs

        #Initialize ACT obs with NaNs so that we can do modifyErrors        
        if self.obs.acronym == 'WISE_ACT':
            self.observations['mag_090PA5_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['mag_090PA6_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['mag_150PA4_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['mag_150PA5_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['mag_150PA6_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['mag_220PA4_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            
            self.observations['magErr_090PA5_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['magErr_090PA6_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['magErr_150PA4_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['magErr_150PA5_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['magErr_150PA6_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            self.observations['magErr_220PA4_si'] = np.ones(len(self.observations['obs_id']))*np.nan
            
            self.columnMapping = {
                "designation" : "designation",
                "obs_id": "obs_id",
                "exp_mjd": "mjd",
                "r_au": "r_au",
                "delta_au": "delta_au",
                "alpha_rad": "alpha_rad",
                "G": "G",
                "logD": "logD",
                "logT1" : "logT1",
                "eta": "eta",
                "eps": "eps",
                "flux_si": ["flux_W1_si", "flux_W2_si", "flux_W3_si", "flux_W4_si", "flux_220PA4_si", "flux_150PA4_si", "flux_150PA5_si", "flux_150PA6_si", "flux_090PA5_si", "flux_090PA6_si"],
                "fluxErr_si": ["fluxErr_W1_si", "fluxErr_W2_si", "fluxErr_W3_si", "fluxErr_W4_si", "fluxErr_220PA4_si", "fluxErr_150PA4_si", 
                               "fluxErr_150PA5_si", "fluxErr_150PA6_si", "fluxErr_090PA5_si", "fluxErr_090PA6_si"], 
                "mag" : ["mag_W1", "mag_W2", "mag_W3", "mag_W4", "mag_220PA4_si", "mag_150PA4_si", "mag_150PA5_si", "mag_150PA6_si", "mag_090PA5_si", "mag_090PA6_si"],
                "magErr" : ["magErr_W1", "magErr_W2", "magErr_W3", "magErr_W4", "magErr_220PA4_si", "magErr_150PA4_si", "magErr_150PA5_si", "magErr_150PA6_si", "magErr_090PA5_si", "magErr_090PA6_si"]
                }
        else:
            self.columnMapping = {
                        "obs_id" : "obs_id",
                        "designation" : "designation",
                        "exp_mjd" : "mjd",
                        "r_au" : "r_au",
                        "delta_au" : "delta_au",
                        "alpha_rad" : "alpha_rad",
                        "eps" : None,
                        "p" : None,
                        "G" : "G",
                        "logT1" : None,
                        "logD" : None,
                        "flux_si" : ["flux_W1_si", "flux_W2_si", "flux_W3_si", "flux_W4_si"],
                        "fluxErr_si" : ["fluxErr_W1_si", "fluxErr_W2_si", "fluxErr_W3_si", "fluxErr_W4_si"],
                        "mag" : ["mag_W1", "mag_W2", "mag_W3", "mag_W4"],
                        "magErr" : ["magErr_W1", "magErr_W2", "magErr_W3", "magErr_W4"]
                }

        if self.ran_override == False:
            self.observations = modifyErrors(self.observations, self.obs, sigma=0.15, columnMapping = self.columnMapping)
            self.ran_override = True
        else:
            print("No need to run this again!")

        self.DPI = 600
        self.SAVE_DIR = "../plots/"
        self.FIG_FORMAT = "png"

        self.SAVE_FIGS = False

        #Get WISE data
        self.data = self.observations[self.observations["designation"].isin([self.desig])]
        self.data_additional = self.additional[self.additional["obs_id"].isin(self.data["obs_id"].values)] # Contains published magnitudes

        if self.data.empty:
            sys.exit('No Wise data for {}'.format(self.name))

        #Add ACT data
        if True: #self.obs.acronym == 'WISE_ACT':
            flux_dict = self.act_flux_dict
            flux_pa4_150, err_pa4_150 = flux_dict['night']['pa4']['150']['flux'], flux_dict['night']['pa4']['150']['var']
            flux_pa4_220, err_pa4_220 = flux_dict['night']['pa4']['220']['flux'], flux_dict['night']['pa4']['220']['var']
            
            flux_pa5_090, err_pa5_090 = flux_dict['night']['pa5']['090']['flux'], flux_dict['night']['pa5']['090']['var']
            flux_pa5_150, err_pa5_150 = flux_dict['night']['pa5']['150']['flux'], flux_dict['night']['pa5']['150']['var']
            
            flux_pa6_090, err_pa6_090 = flux_dict['night']['pa6']['090']['flux'], flux_dict['night']['pa6']['090']['var']
            flux_pa6_150, err_pa6_150 = flux_dict['night']['pa6']['150']['flux'], flux_dict['night']['pa6']['150']['var']
            
            act_fluxes = [flux_pa5_090, flux_pa6_090, flux_pa4_150, flux_pa5_150, flux_pa6_150, flux_pa4_220]*u.mJy
            act_errs = [err_pa5_090, err_pa6_090, err_pa4_150, err_pa5_150, err_pa6_150, err_pa4_220]*u.mJy
            act_freqs = np.array([0.00333103, 0.00333103, 0.00199862, 0.00199862, 0.00199862, 0.00136269])*u.m
            
            act_flux_units = np.zeros(len(act_fluxes))
            act_errs_units = np.zeros(len(act_fluxes))
            for i in range(len(act_fluxes)):
                #Convert 
                #Note factor of 1e6 converts from micron to m
                act_flux_units[i] = ((2.99792458e+14 * act_fluxes[i].to(u.W*u.m**-2*u.Hz**-1)).value / (act_freqs[i].to(u.um)**2).value)*1e6
                act_errs_units[i] = ((2.99792458e+14 * act_errs[i].to(u.W*u.m**-2*u.Hz**-1)).value / (act_freqs[i].to(u.um)**2).value)*1e6
            flux_090PA5_si = np.ones(len(self.data['obs_id']))*act_flux_units[0]
            flux_090PA6_si = np.ones(len(self.data['obs_id']))*act_flux_units[1]
            flux_150PA4_si = np.ones(len(self.data['obs_id']))*act_flux_units[2]
            flux_150PA5_si = np.ones(len(self.data['obs_id']))*act_flux_units[3]
            flux_150PA6_si = np.ones(len(self.data['obs_id']))*act_flux_units[4]
            flux_220PA4_si = np.ones(len(self.data['obs_id']))*act_flux_units[5]
            
            fluxErr_090PA5_si = np.ones(len(self.data['obs_id']))*act_errs_units[0]
            fluxErr_090PA6_si = np.ones(len(self.data['obs_id']))*act_errs_units[1]
            fluxErr_150PA4_si = np.ones(len(self.data['obs_id']))*act_errs_units[2]
            fluxErr_150PA5_si = np.ones(len(self.data['obs_id']))*act_errs_units[3]
            fluxErr_150PA6_si = np.ones(len(self.data['obs_id']))*act_errs_units[4]
            fluxErr_220PA4_si = np.ones(len(self.data['obs_id']))*act_errs_units[5]
            
            self.data['flux_090PA5_si'] = flux_090PA5_si
            self.data['flux_090PA6_si'] = flux_090PA6_si
            self.data['flux_150PA4_si'] = flux_150PA4_si
            self.data['flux_150PA5_si'] = flux_150PA5_si
            self.data['flux_150PA6_si'] = flux_150PA6_si
            self.data['flux_220PA4_si'] = flux_220PA4_si
            
            self.data['fluxErr_090PA5_si'] = fluxErr_090PA5_si
            self.data['fluxErr_090PA6_si'] = fluxErr_090PA6_si
            self.data['fluxErr_150PA4_si'] = fluxErr_150PA4_si
            self.data['fluxErr_150PA5_si'] = fluxErr_150PA5_si
            self.data['fluxErr_150PA6_si'] = fluxErr_150PA6_si
            self.data['fluxErr_220PA4_si'] = fluxErr_220PA4_si
        #if semimajor:
        #    factor = self.semimajor**2 * self.semimajor**2 /(self.data['delta_au']**2 + self.data['r_au']**2)

        #    for key in self.data.keys():
        #        if 'flux' in key and 'PA' not in key:
        #            self.data[key] = self.data[key] * factor


        self.defaultRunConfig = {
            "fitParameters" : ["logT1", "logD", "eps"],
            "emissivitySpecification" : None,
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : self.columnMapping.copy() 
            }

        self.runDict = {}
        self.dataDict = {}

 
    def inv_var(self, flux_dict, freq):

        
        num = 0
        denom = 0
        
        for i in range(len(flux_dict[freq]['flux'])):
            num += flux_dict[freq]['flux'][i]/flux_dict[freq]['var'][i]**2
            denom += 1/flux_dict[freq]['var'][i]**2
            
        return num/denom, denom

    def normalize_act_fluxes(self, verbose = False, write = True):
        #Adjusts ACT fluxes from 1AU, 0 alpha to the r_sun, r_earth, and alpha of the asteroid
        r_sun = np.mean(self.data['r_au'].to_numpy())
        r_earth = np.mean(self.data['delta_au'].to_numpy())
        alpha = np.mean(self.data['alpha_deg'].to_numpy())
        factor = np.sqrt(r_sun)*r_earth**2*10**(0.004*alpha)
        if verbose: print('normalize act: ', factor)
        for pa_key in self.act_flux_dict['night'].keys():
            for freq_key in self.act_flux_dict['night'][pa_key].keys():
                if verbose:
                    print('Before normalization, PA {} Freq {} is {}'.format(pa_key, freq_key,self.act_flux_dict['night'][pa_key][freq_key]['flux']))
                self.act_flux_dict['night'][pa_key][freq_key]['flux'] = self.act_flux_dict['night'][pa_key][freq_key]['flux'] / factor

                self.act_flux_dict['night'][pa_key][freq_key]['var'] = self.act_flux_dict['night'][pa_key][freq_key]['var'] / factor

                if verbose:
                    print('After normalization, PA {} Freq {} is {}'.format(pa_key, freq_key,self.act_flux_dict['night'][pa_key][freq_key]['flux']))
        if write:  
            with open('/scratch/r/rbond/jorlo/actxminorplanets/sigurd/fluxes/{}_normalized_flux_dict.pk'.format(self.name), 'wb') as f:
                pk.dump(self.act_flux_dict, f)


    def add_run(self, key, eps = None, p = None, runData = None, runConfig = None):
    
        self.dataDict[key] = runData

        #It doesn't like making a none comparison to a dataframe, so we keep it if it's not none and set otherwise
        if runData is None:  
            runData = self.data.copy()
    
        if runConfig == None:
            runConfig = self.defaultRunConfig
 
        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = eps
        curColumnMapping['p'] = p

        runConfig['columnMapping'] = curColumnMapping

        self.runDict[str(key)] = runConfig

        self.dataDict[str(key)] = runData

    def atm_fit(self, fitConfig = None): 

        if fitConfig == None: 
            fitConfig = {
                    "chains" : 20,
                    "samples" : 3000,
                    "burnInSamples": 500,
                    "threads": 1,
                    "scaling": 0.01,
                    "plotTrace" : True,
                    "plotCorner" : True,
                    "progressBar" : True,
                    "figKwargs" : {"dpi" : self.DPI}
                } 

        obs = self.obs
        self.model = NEATM(verbose=False, tableDir=self.path+'atm/models/tables/neatm/')

        self.summary, self.model_observations = multiFit(self.model, self.obs, self.dataDict, self.runDict, fitConfig, 
                                               saveDir=self.save_path + 'atm_fits/' + self.name + '_' + self.obs.acronym)
        

    def post_process(self):
        observations_pp, model_observations_pp, observed_stats, model_stats = postProcess(obs, 
                                                                                      data,
                                                                                      model_observations, 
                                                                                      summary)
    def make_SED(self, lambdaNum=250, write = False):
        lambdaRange=[1.0e-6, 30e-5]
        lambdaNum=lambdaNum
        lambdaEdges=[3.9e-6, 6.5e-6, 18.5e-6]

        SEDs = {}
        for code in self.runDict.keys():                
                
 
            SEDs[code] = calcFluxLambdaSED(self.model, self.obs, self.dataDict[code], 
                                           summary=self.summary[self.summary["code"] == code], 
                                           fitParameters=self.runDict[code]["fitParameters"],
                                           emissivitySpecification=self.runDict[code]["emissivitySpecification"],
                                           albedoSpecification=self.runDict[code]["albedoSpecification"],
                                           columnMapping=self.runDict[code]["columnMapping"],
                                           lambdaRange=lambdaRange,
                                           lambdaNum=lambdaNum,
                                           lambdaEdges=lambdaEdges,
                                           linearInterpolation=False)
        
        lambdaRange=[30e-5, 30e-2]
        lambdaNum=lambdaNum
        if self.obs.acronym == 'WISE_ACT': 
            lambdaEdges=[3.9e-6, 6.5e-6, 18.5e-6, 1e-4, 5e-4, 2e-3, 2.1e-3,1.3e-2]
        else:
            lambdaEdges=[3.9e-6, 6.5e-6, 18.5e-6]
        SEDs_ACT = {}
        for code in self.runDict.keys():
            #I have forgotten what this does but it's only active if you've done joint ACT Temperature fitting
            if 'logT1ACT' in self.runDict[code]['fitParameters']:
                temp_summary = self.summary[self.summary["code"] == code]
                array = temp_summary['median'].to_numpy()
                array[0] = array[1]
                
                new_df = pd.DataFrame({'median': array})
                temp_summary.update(new_df)
            else:
                temp_summary = self.summary[self.summary["code"] == code]
            
            temp_data = self.dataDict[code].copy()
            #if self.semimajor:
                #If semimajor axis is specified, then we want to compute the ACT flux at that distance
            #    temp_data['r_au'] = self.semimajor
            #    temp_data['delta_au'] = self.semimajor

            
            SEDs_ACT[code] = calcFluxLambdaSED(self.model, self.obs, temp_data,#self.dataDict[code], 
                                           summary=temp_summary,
                                           fitParameters=self.runDict[code]["fitParameters"],
                                           emissivitySpecification=self.runDict[code]["emissivitySpecification"],
                                           albedoSpecification=self.runDict[code]["albedoSpecification"],
                                           columnMapping=self.runDict[code]["columnMapping"],
                                           lambdaRange=lambdaRange,
                                           lambdaNum=lambdaNum,
                                           lambdaEdges=lambdaEdges,
                                           linearInterpolation=False)
        self.WISE_SEDs = SEDs
        self.ACT_SEDs = SEDs_ACT
        
        for code in SEDs.keys():
            SEDs[code] = SEDs[code].append(SEDs_ACT[code])
            
        self.FULL_SEDs = SEDs
        if write:
            with open(self.save_path+'SEDs/{}_sed.pk'.format(self.name), 'wb') as f:
                pk.dump(SEDs, f)

    def plt_SEDs(self, spectral = None, verbose = False, subplot = False, for_pub = False):
        fig, ax = plt.subplots(1, 1, dpi=self.DPI)

        if spectral != 'bbody':
            plotObservations(self.obs, self.data, spectral = spectral, columnMapping = self.columnMapping, ax=ax, plotMedian=True)
            for code in self.dataDict.keys():
                freqs_for_plot, SED_for_plot = plotSED(self.FULL_SEDs[code], ax=ax, spectral = spectral, plotKwargs={"label": code})

        if self.obs.acronym != 'WISE_ACT':
            flux_pa4_150, err_pa4_150 = self.act_flux_dict['night']['pa4']['150']['flux'], self.act_flux_dict['night']['pa4']['150']['var']
            flux_pa4_220, err_pa4_220 = self.act_flux_dict['night']['pa4']['220']['flux'], self.act_flux_dict['night']['pa4']['220']['var']
            
            flux_pa5_090, err_pa5_090 = self.act_flux_dict['night']['pa5']['090']['flux'], self.act_flux_dict['night']['pa5']['090']['var']
            flux_pa5_150, err_pa5_150 = self.act_flux_dict['night']['pa5']['150']['flux'], self.act_flux_dict['night']['pa5']['150']['var']
            
            flux_pa6_090, err_pa6_090 = self.act_flux_dict['night']['pa6']['090']['flux'], self.act_flux_dict['night']['pa6']['090']['var']
            flux_pa6_150, err_pa6_150 = self.act_flux_dict['night']['pa6']['150']['flux'], self.act_flux_dict['night']['pa6']['150']['var']
            
            act_fluxes = [flux_pa5_090, flux_pa6_090, flux_pa4_150, flux_pa5_150, flux_pa6_150, flux_pa4_220]*u.mJy
            act_errs = [err_pa5_090, err_pa6_090, err_pa4_150, err_pa5_150, err_pa6_150, err_pa4_220]*u.mJy
            act_freqs = np.array([0.00333103, 0.00333103, 0.00199862, 0.00199862, 0.00199862, 0.00136269])*u.m

            act_flux_units = np.zeros(len(act_fluxes))
            act_errs_units = np.zeros(len(act_fluxes))
            for i in range(len(act_fluxes)):
                act_flux_units[i] = ((2.99792458e+14 * act_fluxes[i].to(u.W*u.m**-2*u.Hz**-1)).value / (act_freqs[i].to(u.um)**2).value)
                act_errs_units[i] = ((2.99792458e+14 * act_errs[i].to(u.W*u.m**-2*u.Hz**-1)).value / (act_freqs[i].to(u.um)**2).value)
                if verbose: 
                    print('Converted s/n: ', act_flux_units[i]/act_errs_units[i]) 
                    print('Original s/n: ', act_fluxes[i]/act_errs[i])
            if spectral == 'lambda4':
                act_flux_units *= (act_freqs.to(u.um).value)**4
                act_errs_units *= (act_freqs.to(u.um).value)**4 
            
            if spectral == 'bbody':
                #hard coded run5a
                logT1 = self.summary[self.summary['parameter'] == 'logT1']['median'].values[-2]
                r = np.median(self.data[self.columnMapping["r_au"]].values)
                T_ss = 10**logT1 / np.sqrt(r)
                bbody_model = models.BlackBody(temperature = T_ss*u.K, scale = 1*u.mJy/u.sr)(act_freqs.to(u.GHz, equivalencies=u.spectral())).value

                for i in range(len(act_flux_units)):
                    act_flux_units[i] /= bbody_model[i]
                    act_errs_units[i] /= bbody_model[i]

                plotObservations(self.obs, self.data, spectral = spectral, T_ss = T_ss, columnMapping = self.columnMapping, ax=ax, plotMedian=True)
                for code in self.dataDict.keys():
                    ignore_fig, ignore_ax, freqs_for_plot, SED_for_plot = plotSED(self.FULL_SEDs[code], ax=ax, spectral = spectral, T_ss = T_ss, plotKwargs={"label": code})


            #print(freqs, SED_for_plot)
            ax.errorbar(act_freqs.to(u.um).value, 
                                    act_flux_units, 
                                    act_errs_units, 
                                    fmt='o',
                                    c="k",
                                    ms=2,
                                    capsize=2,
                                    elinewidth=2)
        
        if not for_pub: ax.set_title(str(self.name)) 
        if for_pub:
            a = plt.axes([0.45, 0.35, .4, .3])
            a.errorbar(act_freqs.to(u.um).value,
                                    act_flux_units,
                                    0.03*act_flux_units,#act_errs_units,
                                    fmt='none',
                                    c="k",
                                    ms=2,
                                    capsize=2,
                                    elinewidth=2)
            a.plot(freqs_for_plot, SED_for_plot)
 
            a.set_xlim(1e3,5e3)
            a.set_ylim(1e-6, 5e-6)
            a.set_xscale('log')
            a.set_yscale('log')
            a.set_xticks([1e3, 2e3, 4e3], [1e3, 2e3, 4e3])
            locs_x = np.array([act_freqs.to(u.um).value[0], act_freqs.to(u.um).value[1], act_freqs.to(u.um).value[-1]])
            locs_y = np.array([act_flux_units[0], act_flux_units[1], act_flux_units[2]])
            #plt.annotate(['f090', 'f150', 'f220'], (locs_x, locs_y),(locs_x, locs_y), 
            #         arrowprops=dict(facecolor='black', shrink=0.05))
            #plt.yticks([])

            mark_inset(ax, a, loc1=2, loc2=1, fc="none", ec="0.5")
            #a.annotate(['f090'], ([act_freqs.to(u.um).value[0]], [act_flux_units[0]]))
        #plt.legend() 
        if not subplot:
            if for_pub : plot_name = 'pub_{}_atm_{}_plot'.format(self.name, self.obs.acronym)
            else: plot_name = '{}_atm_{}_plot'.format(self.name, self.obs.acronym)
            ax.set_xlim(1, 1e5)
            plt.savefig(self.plot_path + plot_name +'.png')
            plt.savefig(self.plot_path + plot_name +'.pdf') 
        else: 
            #reset scaling
            ax.relim()
            ax.autoscale_view()
            if subplot == 'ACT': 
                ax.set_xlim(1e3, 1e4)
                flags = np.where((1e3 <= (self.FULL_SEDs['run5a']["lambda"]*1e6)) & ((self.FULL_SEDs['run5a']["lambda"]*1e6)<= 1e4))[0]
                
                adj_fluxes = (self.FULL_SEDs['run5a']["lambda"] * 1e6)**4 * 1/1e6*self.FULL_SEDs['run5a']["flux"]
                adj_fluxes = adj_fluxes.to_numpy()
                
                y_lim_lower, y_lim_upper = min(adj_fluxes[flags]), max(adj_fluxes[flags])
                ax.set_ylim(y_lim_lower*0.05, y_lim_upper*10.)
                ax.legend(loc='lower center', ncol=5, prop={'size': 8})  
            elif subplot == 'WISE':  
                ax.set_xlim(1, 1e2)

                flags = np.where((1e0 <= (self.FULL_SEDs['run5a']["lambda"]*1e6)) & ((self.FULL_SEDs['run5a']["lambda"]*1e6)<= 1e2))[0]

                adj_fluxes =  1/1e6*self.FULL_SEDs['run5a']["flux"]
                adj_fluxes = adj_fluxes.to_numpy()

                y_lim_lower, y_lim_upper = min(adj_fluxes[flags]), max(adj_fluxes[flags])

                ax.set_ylim(y_lim_lower*0.05, y_lim_upper*10.)
                ax.legend(loc='lower center', ncol=5, prop={'size': 8})

            else: 
                print('Error, invalid subplot type')
                return 1        

            plt.savefig(self.plot_path + 'atmplots/{}_atm_{}_{}_plot.pdf'.format(self.name, self.obs.acronym, subplot))
            plt.savefig(self.plot_path + 'atmplots/{}_atm_{}_{}_plot.png'.format(self.name, self.obs.acronym, subplot))

        if self.show_plot:
            plt.show()
        
        #self.act_fluxes = act_fluxes
        #self.act_flux_units = act_flux_units

    def eval_SED(self, r_au, r_delt, lambd, code, linearInterpolation=True, threads = 1, normalize = False, verbose = False):
        
        if type(lambd) != list: lambd = [lambd]
        if type(lambd) != np.ndarray: lambd = np.array(lambd)

        summary = self.summary[self.summary["code"] == code]

        fitParameters=self.runDict[code]["fitParameters"]
        emissivitySpecification=self.runDict[code]["emissivitySpecification"]
        albedoSpecification=self.runDict[code]["albedoSpecification"]
        columnMapping=self.runDict[code]["columnMapping"]
        data = self.dataDict[code] 

        obs = self.obs
 
        fitParametersSet, parametersSet, emissivityParameters, albedoParameters, dataParametersToIgnoreSet = handleParameters(
            self.obs,
            fitParameters,
            data.columns.tolist(),
            emissivitySpecification=emissivitySpecification,
            albedoSpecification=albedoSpecification,
            columnMapping=columnMapping)

     
        if "logT1" in parametersSet and "logT1" not in fitParametersSet:
            logT1 = data[columnMapping["logT1"]].values[0] * np.ones(len(lambd))
        else:
            logT1 = summary[summary["parameter"] == "logT1"]["median"].values[0] * np.ones(len(lambd))
    
        T_ss = 10**logT1 / np.sqrt(r_au)

        if "logD" in parametersSet and "logD" not in fitParametersSet:
            logD = data[columnMapping["logD"]].values[0] * np.ones(len(lambd))
        else:
            logD = summary[summary["parameter"] == "logD"]["median"].values[0] * np.ones(len(lambd))

        if "alpha_rad" in parametersSet and "alpha_rad" not in fitParametersSet:
            alpha = np.median(data[columnMapping["alpha_rad"]].values) * np.ones(len(lambd))
        else:
            alpha = summary[summary["parameter"].isin(["alpha_rad__{}".format(i) for i in range(0, len(data))])]["median"].values * np.ones(len(lambd))

        if "G" in parametersSet and "G" not in fitParametersSet:
            G = data[columnMapping["G"]].values[0] * np.ones(len(lambd))
        else:
            G = summary[summary["parameter"] == "G"]["median"].values[0] * np.ones(len(lambd))

        if emissivityParameters == "eps" and emissivitySpecification != "auto":
            if "eps" in parametersSet and "eps" not in fitParametersSet:
                eps = data[columnMapping["eps"]].values[0] * np.ones(len(lambd))
            else:
                eps = summary[summary["parameter"] == "eps"]["median"].values[0] * np.ones(len(lambd))
    
            if albedoSpecification == "auto":
                p = (1 - eps) / calcQ(G)
    
        if type(emissivityParameters) is list and emissivitySpecification != "auto":
            eps_values = np.zeros_like(emissivityParameters, dtype=float)
            for i, parameter in enumerate(emissivityParameters):
                if parameter in parametersSet and parameter not in fitParametersSet:
                    eps_values[i] = data[parameter].values[0]
                else:
                    eps_values[i] = summary[summary["parameter"] == parameter]["median"].values[0]
    
            eps = np.zeros_like(lambd)
            if linearInterpolation is True:
                eps = np.interp(lambd, obs.filterEffectiveLambdas, eps_values)
            else:
                for i, (edge_start, edge_end) in enumerate(zip([lambdaRange[0]] + lambdaEdges,
                                                       lambdaEdges + [lambdaRange[-1]])):
                    eps = np.where((lambd >= edge_start) & (lambd <= edge_end), eps_values[i], eps)
    
            if albedoSpecification == "auto":
                p = (1 - eps) / calcQ(G)
    
        if albedoParameters == "p" and albedoSpecification != "auto":
            if "p" in parametersSet and "p" not in fitParametersSet:
                p = data[columnMapping["p"]].values[0] * np.ones(len(lambd))
            else:
                p = summary[summary["parameter"] == "p"]["median"].values[0] * np.ones(len(lambd))
    
            if emissivitySpecification == "auto":
                eps = 1 - p * calcQ(G)
    
        if type(albedoParameters) is list and albedoSpecification != "auto":
            p_values = np.zeros_like(albedoParameters, dtype=float)
            for i, parameter in enumerate(albedoParameters):
                if parameter in parametersSet and parameter not in fitParametersSet:
                    p_values[i] = data[parameter].values[0]
                else:
                    p_values[i] = summary[summary["parameter"] == parameter]["median"].values[0]
    
            p = np.zeros_like(lambd)
            if linearInterpolation is True:
                p = np.interp(lambd, obs.filterEffectiveLambdas, p_values)
            else:
                for i, (edge_start, edge_end) in enumerate(zip([lambdaRange[0]] + lambdaEdges,
                                                       lambdaEdges + [lambdaRange[-1]])):
                    p = np.where((lambd >= edge_start) & (lambd <= edge_end), p_values[i], p)
    
            if emissivitySpecification == "auto":
                eps = 1 - p * calcQ(G)

        if "G" in parametersSet and "G" not in fitParametersSet:
    	        G = data[columnMapping["G"]].values[0] * np.ones(len(lambd))
        else:
            G = summary[summary["parameter"] == "G"]["median"].values[0] * np.ones(len(lambd))
        
        
       
        fluxes_W = calcFluxLambdaAtObsWithSunlight(self.model, r_au, r_delt, lambd, T_ss, 10**logD, alpha, eps, p, G, threads=threads) * (u.W*u.m**-2*u.Hz**-1)
        #print(fluxes_W*(lambd*1e6)**4/1e6)
        if normalize:
            r_sun = np.mean(self.data['r_au'].to_numpy())
            r_earth = np.mean(self.data['delta_au'].to_numpy())
            alpha_my = np.mean(self.data['alpha_deg'].to_numpy())
            
            weight = np.sqrt(r_sun)*r_earth**2*10**(0.004*alpha_my) 
            if verbose: print("Make SED weight: ", weight)
            fluxes_W *= weight 
        lambd = lambd * u.m
        fluxes_mJy = np.zeros(len(lambd))

        for i in range(len(fluxes_mJy)):
            #Inverse conversion as above. Note the factor of 1e-6 which is cause calcFluxLambda gives in units of W /m^2 /Hz /um
            fluxes_mJy[i] = fluxes_W[i].to(u.mJy).value * 1e-6 / 2.99792458e+14 * (lambd[i].to(u.um)**2).value
        
        return fluxes_mJy

    def load_SEDs(self):
        with open(self.save_path+'SEDs/{}_sed.pk'.format(self.name), 'rb') as f:
            self.FULL_SEDs = pk.load(f)

if __name__ == '__main__':
        dataDict["run0"] = self.data.copy()

        #Need to modify base column mapping to account for what combo of params are free
        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = None
        curColumnMapping['p'] = None

        runDict["run0"] = {
            "fitParameters" : ["logT1", "logD", "eps"],
            "emissivitySpecification" : None,
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }
        dataDict["run1"] = self.data.copy()
        dataDict["run1"]["eps_W3W4"] = np.ones(len(self.data)) * 0.9

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ["eps_W3W4"]
        curColumnMapping['p'] = None

        runDict["run1"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1W2"],
            "emissivitySpecification" : {
                        "eps_W1W2" : ["W1","W2"],
                        "eps_W3W4" : ["W3","W4"]},
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" :self.columnMapping
        }
    
        dataDict["run2a"] = self.data.copy()
        dataDict["run2a"]["eps_W3"] = np.ones(len(self.data)) * 0.70
        dataDict["run2a"]["eps_W4"] = np.ones(len(self.data)) * 0.86

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ["eps_W3", "eps_W4"]
        curColumnMapping['p'] = None

        runDict["run2a"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1W2"],
            "emissivitySpecification" : {
                        "eps_W1W2" : ["W1","W2"],
                        "eps_W3" : ["W3"],
                        "eps_W4" : ["W4"]},
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }
        dataDict["run2b"] = self.data.copy()
        dataDict["run2b"]["eps_W3"] = np.ones(len(self.data)) * 0.70
        dataDict["run2b"]["eps_W4"] = np.ones(len(self.data)) * 0.86

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ["eps_W3", "eps_W4"]
        curColumnMapping['p'] = None

        runDict["run2b"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1", "eps_W2"],
            "emissivitySpecification" : "perBand",
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run3a"] = self.data.copy()
        dataDict["run3a"]["eps_W3"] = np.ones(len(self.data)) * 0.76
        dataDict["run3a"]["eps_W4"] = np.ones(len(self.data)) * 0.93

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ['eps_W3', 'eps_W4']
        curColumnMapping['p'] = None

        runDict["run3a"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1W2"],
            "emissivitySpecification" : {
                        "eps_W1W2" : ["W1","W2"],
                        "eps_W3" : ["W3"],
                        "eps_W4" : ["W4"]},
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run3b"] = self.data.copy()
        dataDict["run3b"]["eps_W3"] = np.ones(len(self.data)) * 0.76
        dataDict["run3b"]["eps_W4"] = np.ones(len(self.data)) * 0.93

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ['eps_W3', 'eps_W4']
        curColumnMapping['p'] = None

        runDict["run3b"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1", "eps_W2"],
            "emissivitySpecification" : "perBand",
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run4a"] = self.data.copy()
        dataDict["run4a"]["eps_W3"] = np.ones(len(self.data)) * 0.80
        dataDict["run4a"]["eps_W4"] = np.ones(len(self.data)) * 0.98

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ['eps_W3', 'eps_W4']
        curColumnMapping['p'] = None

        runDict["run4a"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1W2"],
            "emissivitySpecification" : {
                        "eps_W1W2" : ["W1","W2"],
                        "eps_W3" : ["W3"],
                        "eps_W4" : ["W4"]},
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run4b"] = self.data.copy()
        dataDict["run4b"]["eps_W3"] = np.ones(len(self.data)) * 0.80
        dataDict["run4b"]["eps_W4"] = np.ones(len(self.data)) * 0.98

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = ['eps_W3', 'eps_W4']
        curColumnMapping['p'] = None

        runDict["run4b"] = {
            "fitParameters" : ["logT1", "logD", "eps_W1", "eps_W2"],
            "emissivitySpecification" : "perBand",
            "albedoSpecification": "auto",
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run5a"] = self.data.copy()
        dataDict["run5a"]["eps"] = np.ones(len(self.data)) * 0.9
        dataDict["run5a"]["p_W3"] = np.ones(len(self.data)) * 0.0
        dataDict["run5a"]["p_W4"] = np.ones(len(self.data)) * 0.0

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = 'eps'
        curColumnMapping['p'] = ['p_W3', 'p_W4']

        runDict["run5a"] = {
            "fitParameters" : ["logT1", "logD", "p_W1W2"],
            "emissivitySpecification" : None,
            "albedoSpecification": {
                        "p_W1W2" : ["W1", "W2"],
                        "p_W3" : ["W3"],
                        "p_W4" : ["W4"]},
            "fitFilters" : "all",
            "columnMapping" : curColumnMapping
        }

        dataDict["run5b"] = self.data.copy()
        dataDict["run5b"]["eps"] = np.ones(len(self.data)) * 0.9
        dataDict["run5b"]["p_W3"] = np.ones(len(self.data)) * 0.0
        dataDict["run5b"]["p_W4"] = np.ones(len(self.data)) * 0.0

        curColumnMapping = self.columnMapping.copy()
        curColumnMapping['eps'] = 'eps'
        curColumnMapping['p'] = ['p_W3', 'p_W4']

        runDict["run5b"] = {
            "fitParameters" : ["logT1", "logD", "p_W1W2"],
            "emissivitySpecification" : None,
            "albedoSpecification": {
                        "p_W1W2" : ["W1", "W2"],
                        "p_W3" : ["W3"],
                        "p_W4" : ["W4"]},
            "fitFilters" : ["W1", "W3", "W4"],
            "columnMapping" : curColumnMapping
        }


