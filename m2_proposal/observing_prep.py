#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 18:11:19 2020

@author: Emily
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd

from astropy.coordinates import AltAz, EarthLocation, get_sun, SkyCoord
from astropy.table import Table,vstack
from astropy.time import Time
import astropy.units as u

from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()
# https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
# -------------------------------------------------------------------------------------------- #
# locations of observatories
alma = EarthLocation.of_site('alma')
gbt = EarthLocation.of_site('GBT')
geminiN = EarthLocation.of_site('Gemini North')
geminiS = EarthLocation.of_site('Gemini South')
lofar = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=4080*u.m)
palomar = EarthLocation.of_site('Palomar')
sma = EarthLocation(lat=19.8243*u.deg, lon=-155.478*u.deg, height=4080*u.m)
vla = EarthLocation.of_site('vla')
telescope_locations = {'ALMA':alma, 'GBT':gbt, 'GeminiN':geminiN, 'GeminiS':geminiS, 'lofar':lofar, 'Palomar':palomar, 'SMA':sma, 'VLA':vla}
# -------------------------------------------------------------------------------------------- #
# elevation limits
alma_elev_lim = 10*u.deg
gbt_elev_lim = 25.*u.deg
gemini_elev_lim = 20.*u.deg
lofar_elev_lim = 15.*u.deg
palomar_elev_lim = 25*u.deg
sma_elev_lim = 15*u.deg
vla_elev_lim = 15*u.deg
elev_lims = {'ALMA':alma_elev_lim, 'GBT':gbt_elev_lim, 'GeminiN':gemini_elev_lim, 'GeminiS':gemini_elev_lim, 'lofar':lofar_elev_lim, 'Palomar':palomar_elev_lim, 'SMA':sma_elev_lim, 'VLA':vla_elev_lim}
# -------------------------------------------------------------------------------------------- #
# other
colors = ['k','r','orange','lightcoral','gold','limegreen','g','cornflowerblue','deepskyblue','b','mediumslateblue','magenta','mediumorchid','plum','hotpink','peru']
day_colors = ['blue','royalblue','cornflowerblue','dodgerblue','deepskyblue','darkturquoise','aqua','aquamarine','springgreen']
location = '/Users/Emily/Documents/Research/Observing/'
# -------------------------------------------------------------------------------------------- #
def read_catalog(filepath):
    '''Read and return the catalog.'''
    catalog = Table.read(filepath, format="fits")    
    return catalog
    
def visibilities_loop_days_and_objects(table_path,
                                        table,
                                        telescope,
                                        time_zone,
                                        dates,
                                        semester,
                                        file_save_location,
                                        plot_suffix):
    '''Loops through a list of days and objects and for each day, plots the visibilities for a list of objects.'''
    # read in table
    if table_path != 'N/A':
        t = pd.read_csv(table_path)
    else:
        t = table
    # create coordinate arrays
    clu_coord = SkyCoord(ra=t['RA'].values*u.deg,dec=t['Dec'].values*u.deg)
    # replace _ in MOO names with spaces
    t['Name'] = t['Name'].str.replace('_', ' J')

    for dt in dates:
        print(dt)

        midnight = Time(dt)
        delta_midnight = np.linspace(0, 24, 1000)*u.hour
        times_day = midnight + delta_midnight
        frame_day = AltAz(obstime=times_day, location=telescope_locations[telescope])
        # sun (necessary for the definition of night)
        sunaltazs = get_sun(times_day).transform_to(frame_day)

        ### PLOT ###
        fig = plt.figure()
        ax = plt.subplot(111)

        # data
        for i,c in enumerate(clu_coord):
            clualtazs = c.transform_to(frame_day)
            ax.plot(delta_midnight, clualtazs.alt,
                    label=t['Name'][i],
                    color=colors[i])

        # limits
        ax.set_xlim(0*u.hour, 24*u.hour)
        ax.set_xticks(np.arange(25))
        ax.set_ylim(0*u.deg, 90*u.deg)

        # define night
        ax.fill_between(delta_midnight, 
                        0*u.deg, 
                        90*u.deg,
                        sunaltazs.alt < -0*u.deg, 
                        color='0.5', 
                        zorder=0)
                        
        if telescope == 'SMA':
            # observing period of SMA (20:00 - 8:00)
            ax.fill_between(x=delta_midnight, 
                            y1=0*u.deg, 
                            y2=90*u.deg,
                            where= (delta_midnight<8*u.hour) | (delta_midnight > 20*u.hour), 
                            color='b', 
                            zorder=0,
                            alpha=0.5)

        # show elevation limit
        ax.axhspan(0*u.deg,elev_lims[telescope], alpha=0.4, color='k',zorder=20)
        
        # plot labels
        day = dt.strftime("%Y-%m-%d")
        ax.set_title(semester+' '+ day)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel(time_zone)
        ax.set_ylabel('Altitude [deg]')
        
        plt.tight_layout()
        plt.savefig(file_save_location+telescope+'_'+semester+plot_suffix+day+'.pdf',
                    bbox_inches='tight',
                    dpi=300,
                    format='pdf')
        #file_save_location+semester+'/visibilities/'+semester+plot_suffix+day+'.png',
        #plt.show()

def visibilities_loop_days_for_one_object(table_path,
                                        table,
                                        telescope,
                                        time_zone,
                                        dates,
                                        semester,
                                        file_save_location,
                                        plot_suffix):
    '''For each object, loops through a list of days.'''
    # read in table
    if table_path != 'N/A':
        t = pd.read_csv(table_path)
    else:
        t = table
    clu_coord = SkyCoord(ra=t['RA'].values*u.deg,dec=t['Dec'].values*u.deg)
    # replace _ in MOO names with spaces
    t['Name'] = t['Name'].str.replace('_', ' J')

    rise_times = []
    set_times = []

    for i,c in enumerate(t['Name']):
        print(c)
        ### PLOT ###
        fig = plt.figure()
        ax = plt.subplot(111)
        # limits
        ax.set_xlim(0*u.hour, 24*u.hour)
        ax.set_xticks(np.arange(25))
        ax.set_ylim(0*u.deg, 90*u.deg)
        # show elevation limit for M2
        ax.axhspan(0*u.deg,elev_lims[telescope], alpha=0.4, color='k',zorder=20)

        # plot labels
        ax.set_title(c+ ' ' +semester)
        ax.set_xlabel(time_zone)
        ax.set_ylabel('Altitude [deg]')

        # plot elevation curve for each day
        for d,dt in enumerate(dates):
            print(dt)

            midnight = Time(dt)
            delta_midnight = np.linspace(0, 24, 1000)*u.hour
            times_day = midnight + delta_midnight
            frame_day = AltAz(obstime=times_day, location=telescope_locations[telescope])

            # sun (necessary for the definition of night)
            sunaltazs = get_sun(times_day).transform_to(frame_day)
            sunup = sunaltazs.alt > -0*u.deg
            rise_times.append(delta_midnight[sunup][0].value)
            set_times.append(delta_midnight[sunup][-1].value)

            # data
            clualtaz = clu_coord[i].transform_to(frame_day)
            day = dt.strftime("%Y-%m-%d")
            ax.plot(delta_midnight, clualtaz.alt,
                    label=day,
                    color=day_colors[d])

        # plot average night-time
        ax.axvspan(0, np.average(rise_times), color='0.5',zorder=0)
        ax.axvspan(np.average(set_times), 24, color='0.5',zorder=0)
        if telescope == 'SMA':
            # observing period of SMA (20:00 - 8:00)
            ax.fill_between(x=delta_midnight, 
                            y1=0*u.deg, 
                            y2=90*u.deg,
                            where= (delta_midnight<8*u.hour) | (delta_midnight > 20*u.hour), 
                            color='b', 
                            zorder=0,
                            alpha=0.5)
                            
        # date based legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # show and save plot
        plt.tight_layout()
        plt.savefig(file_save_location+semester+'/'+c+'_'+semester+plot_suffix+'.png',
                    bbox_inches='tight',
                    dpi=300,
                    format='png')

        plt.show()
