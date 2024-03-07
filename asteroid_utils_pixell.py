
from pixell import enmap, utils, reproject, bunch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import os,sys,argparse 
from scipy.interpolate import interp1d
import math
import pandas as pd
import pickle as pk
import h5py
import time, datetime
from pathlib import Path

from astropy import wcs
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
import astropy.table as tb 
from astropy.time import Time 
#from astroquery.jplhorizons import Horizons

import ephem

import re
from numba import jit

def compute_alpha(ra_sun, dec_sun, d_earth_sun, ra_ast, dec_ast, d_earth_ast):
    sun_earth_vec = utils.ang2rect([ra_sun, dec_sun])*d_earth_sun
    earth_ast_vec = utils.ang2rect([ra_ast, dec_ast])*d_earth_ast
    sun_ast_vec = earth_ast_vec - sun_earth_vec
    angle = utils.vec_angdist(-earth_ast_vec, -earth_ast_vec+sun_earth_vec)

    return angle

def get_desig(id_num):
    home = str(Path.home())
    with open(home+'/dev/minorplanets/asteroids.pk', 'rb') as f:
        df = pk.load(f)
        name = df['name'][id_num]
        desig = df['designation'][id_num]
        semimajor = df['semimajor'][id_num]
    return desig, name, semimajor


  
def get_index(name):
    '''
    Inputs:
      name, type: string, name of object

    Output:
      desig, type: integer, index of object in asteroids.pk file
    '''
  
    with open('/home/r/rbond/ricco/minorplanets/asteroids.pk', 'rb') as f:
        df = pk.load(f)    
        idx = np.where((df['name'] == name))[0]    
        desig = df['designation'][idx]

    string = desig.to_string()

    num_string = ''
  
    for s in string:
        if s == ' ':
            break
        else:
            num_string += s

    try:
        indx = int(num_string)
        return indx
    except ValueError:
        print('Object not in current data set')

def inv_var(data, variances):
    '''
    Inputs
      data, type: array of ints/floats, data to be weighted
      variances, type: array of ints/floats, variances of data
    Output
      Inverse-variance weighted average
    '''
    ave = 0
    var = 0
    for i in range(len(data)):
        ave += data[i]/variances[i]
        var += 1/variances[i]
    return ave/var, 1/var    

def tnoStamp(ra, dec, imap, width = 0.5):
    #Takes an individual stamp of a map at the requested ra/dec. We form the ratio later  
    #frhs is a matched filter map and kmap is the inverse variance per pixel
    #both maps are at a ~3 day cadence

    #Inputs: ra, dec in degrees, j2000
    #kmap or frhs, described above. They must have the same wcs
    #but if you're using this code and sigurd's maps they will
    #width, the desired width of the stamp in degrees

    #Output: a stamp, centered on ra/dec, with width width, where
    #each pixel in the map records the S/N. This stamp is for one
    #object, one 3 day map. These must then be stacked for each object
    #and then the objects can be stacked together

    #Find the pixel 
    coords = np.deg2rad(np.array((dec,ra)))
    ypix,xpix = enmap.sky2pix(imap.shape,imap.wcs,coords)
        
    #nans are formed when try to form S/n for pixels with no hits
    #I just set the S/N to 0 which I think is safe enough
    imap[~np.isfinite(imap)] = 0

    #Reproject will attempt to take a stamp at the desired location: if it can't
    #for whatever reason just return None. We don't want it throwing errors
    #while on the wall, however the skips should probably be checked after
    #try:
    stamp = reproject.thumbnails(imap, [coords])
    #except:
    #    return None
    
    return stamp

class OrbitInterpolator:
    '''
    Constructs a class that can predict, using an interpolation scheme, the location of an object given an identifier and time
    '''
    def __init__(self, table):
        '''
        Requires a table generated from astroquery, querying JPL HORIZONS. The 
        interpolation is done automatically, but, of course, only works
        if the table is sampled densely enough and in the time range for which
        the positions were queried
        '''
        self.table = table
        self.targets = np.unique(table['targetname'])

        self._construct_dictionary()

    def _interpolate_radec(self, target):
        '''
        "Hidden" function that constructs the interpolations for each target
        '''        
        table = self.table[self.table['targetname'] == target]
        zero = np.min(table['datetime_jd'])
        ra_interp = interp1d(table['datetime_jd'] - zero, table['RA'])
        dec_interp = interp1d(table['datetime_jd'] - zero, table['DEC'])
        delta_interp = interp1d(table['datetime_jd'] - zero, table['delta'])
        r_interp = interp1d(table['datetime_jd'] - zero, table['r'])

        return zero, ra_interp, dec_interp, delta_interp, r_interp

    def _construct_dictionary(self):
        '''
        "Hidden" function that creates the look-up dictionary of targets for simplicity of usage
        '''
        self.obj_dic = {}
            
        for j,i in enumerate(self.targets):
            z, ra, dec, delta, r = self._interpolate_radec(i)
            self.obj_dic[i] = {}
            self.obj_dic[i]['zero'] = z
            self.obj_dic[i]['RA'] = ra
            self.obj_dic[i]['DEC'] = dec
            self.obj_dic[i]['delta'] = delta
            self.obj_dic[i]['r'] = r

    def get_radec_dist(self, target, time):
        '''
        Specifying a target name (see self.obj_dic.keys() for a list of targets) and a time (in JD), finds
        the interpolated RA and Dec for the objects
        '''
        time = time + 2400000.5

        t_intep = time - self.obj_dic[target]['zero']
        
        ra = self.obj_dic[target]['RA'](t_intep)
        dec = self.obj_dic[target]['DEC'](t_intep)
        dist = self.obj_dic[target]['delta'](t_intep)
        r = self.obj_dic[target]['r'](t_intep)

        return ra, dec, dist, r

class QueryHorizons:
    '''
    Constructs a class that can query a bunch of positions for objects from JPL-Horizons given a set of times and an observatory location
    '''

    def __init__(self, time_start, time_end, observer_location, step = '1d'):
        '''
        Initialization function

        Arguments:
        - time_start: start time for the query, should be in MJD
        - time_end: end time for the query, should be in MJD
        - observer_location: location of the observer
        - step: time step for the ephemeris
        The simples way to get the observer location variable is via the
        list of IAU observatory codes: https://en.wikipedia.org/wiki/List_of_observatory_codes
        
        Custom locations are also accepted by JPL
        '''
        if type(time_start) == str: 
            t_st = Time(time_start, format='isot', scale='utc')
            
        else: 
            t_st = Time(time_start, format='mjd')
        self.time_start = t_st.utc.iso  
        if type(time_end) == str: 
            t_en = Time(time_end, format='isot', scale='utc')
            
        else: 
            t_en = Time(time_end, format='mjd')
            
        self.time_end = t_en.utc.iso 
        
        self.observer = observer_location

        self.step = step 


    def queryObjects(self, objects):
        '''
        Returns a table (and saves it on the object as well) for the provided list of objects
        '''
        self.table = [] 

        for i in objects:
            query = Horizons(id = i, location = self.observer, 
                epochs = {'start' : self.time_start, 'stop' : self.time_end, 'step' : self.step})

            eph = query.ephemerides()
            
            self.table.append(eph['RA', 'DEC', 'datetime_jd', 'targetname', 'delta', 'r', 'alpha'])
        
        self.table = tb.vstack(self.table)

        return self.table

class minorplanet():
    '''
    Constructs a class that includes everything we want to do with a minor planet
    '''
    def __init__(self, name, semimajor, path = '/scratch/r/rbond/jorlo/actxminorplanets/sigurd/asteroids/', obs = 'W99', 
                 t_start ='2010-01-01T00:00:00', t_end='2020-01-01T00:00:00'):
        '''
        Requires a name for the object, as well as an observatory code which can be found at
        https://en.wikipedia.org/wiki/List_of_observatory_codes
        See the QueryHorizons class for full details of observatory code options
        '''
        if name[0] not in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            name.capitalize()
        self.name = name
        self.path = path
        self.obs = obs
        self.semimajor = semimajor
        
        self.sun = ephem.Sun() 

        #Initialize orbit and flux stuff
        self.t_start = t_start
        self.t_end = t_end
        
        self.map_dict = {'day':{'pa4':{'150':{'flux':0, 'snr':0}, '220':{'flux':0, 'snr':0}},
                                'pa5':{'090':{'flux':0, 'snr':0}, '150':{'flux':0, 'snr':0}},
                                'pa6':{'090':{'flux':0, 'snr':0}, '150':{'flux':0, 'snr':0}}
                               },
                         'night':{'pa4':{'150':{'flux':0, 'snr':0}, '220':{'flux':0, 'snr':0}},
                                'pa5':{'090':{'flux':0, 'snr':0}, '150':{'flux':0, 'snr':0}},
                                'pa6':{'090':{'flux':0, 'snr':0}, '150':{'flux':0, 'snr':0}}
                                 }
                         }                        
 
        self.flux_dict = {'day':{'pa4':{'150':{'flux':0, 'var':0}, '220':{'flux':0, 'var':0}},
                                 'pa5':{'090':{'flux':0, 'var':0}, '150':{'flux':0, 'var':0}},
                                 'pa6':{'090':{'flux':0, 'var':0}, '150':{'flux':0, 'var':0}}
                                },
                          'night':{'pa4':{'150':{'flux':0, 'var':0}, '220':{'flux':0, 'var':0}},
                                   'pa5':{'090':{'flux':0, 'var':0}, '150':{'flux':0, 'var':0}},
                                   'pa6':{'090':{'flux':0, 'var':0}, '150':{'flux':0, 'var':0}}
                                  }
                          }                              
        #We interpolate the orbit for the minor planet at initialization as it's fairly fast
        #Stacking is a method as it is much slower so we want to call it only when we actually care
        #self.eph_table = QueryHorizons(time_start=self.t_start, time_end=self.t_end, 
        #                               observer_location=obs, step = '1d').queryObjects([str(self.name)])
        
        #self.interp_orbit = OrbitInterpolator(self.eph_table)
        if self.name in ['Jupiter', 'Mars', 'Mercury', 'Moon', 'Neptune', 'Saturn', 'Sun', 'Uranus', 'Venus']:
            #We use planets as calibration checks, this just loads them from the right place if they're a planet
            info = np.load('/home/r/rbond/sigurdkn/project/actpol/ephemerides/objects_tmp/{}.npy'.format(self.name)).view(np.recarray)
        else: 
            info = np.load('/gpfs/fs0/project/r/rbond/sigurdkn/actpol/ephemerides/objects/{}.npy'.format(self.name)).view(np.recarray)
        self.orbit   = interp1d(info.ctime, [utils.unwind(info.ra*utils.degree), 
                                                    info.dec*utils.degree, info.r, info.rsun, info.ang*utils.arcsec], kind=3)
        

    def show_orbit(self, t_orb_start = '2013-01-01T00:00:00', t_orb_end = '2020-01-01T00:00:00', 
                   directory = None):
        '''
        Function that plots the orbit of the minor planet from t_orb_start to t_orb_end
        '''
        
        t_start = Time(t_orb_start, format='isot', scale='utc')
        t_start = t_start.mjd

        t_en = Time(t_orb_end, format='isot', scale='utc')
        t_en = t_en.mjd
        mjds = np.linspace(t_start, t_en, 1000)

        ras, decs, delts, rs = self.interp_orbit.get_radec_dist(self.eph_table['targetname'][0], mjds)

        plt.plot(ras, decs)
        plt.xlabel('ra')
        plt.ylabel('dec')
        plt.title('Orbital Plot of {}'.format(self.eph_table['targetname'][0]))
        if directory is not None:
            plt.savefig(directory+'{}_orbit.pdf'.format(str(self.eph_table['targetname'][0]).replace(' ', '_').replace('(','').replace(')','')))
        #plt.show()
        plt.close()
        
    def make_stack(self, pa, freq = 'f150', pol='T', restrict_time = False, weight_type=None, plot = False, verbose = False, weight_debug = False, time = 'night', lightcurve=False, freq_adjust = None, extern_calib = None, time_cut = None):
        pol_dict = {'T':0, 'Q':1, 'U':2}
        pol = pol_dict[pol]

        if not os.path.isdir(self.path+'/'+self.name+'/'):
            print('Depth 1 stamps not made, please make them')
            return
        rho_stack = 0
        kappa_stack = 0
        
        rho_un = 0 
        kappa_un = 0
        
        rho_weight = 0
        kappa_weight = 0 
 
        if verbose or weight_debug: 
            hours = [] 
            weights = []
            alphas = []
            fluxes = []
            r_suns = []
            r_earths = []
            isos = []
            sns = []
        if lightcurve:
            fluxes_lc = []
            errs_lc = []
            times_lc = []
            Fs_lc = []           
            rho_lc = []
            kappa_lc = [] 
        for i, dirname in enumerate(os.listdir(path=self.path+'/'+self.name+'/')):
            
            if dirname.find('kappa.fits') == -1: continue #Just look at the kappa maps
            if dirname.find(pa) == -1: continue
            if dirname.find('f'+str(freq)) == -1: continue #Check each file to see if it's the specified freq
            if verbose: print('In map: ', dirname)
         
            rhofile    = utils.replace(dirname, "kappa.fits", "rho.fits")
            infofile = utils.replace(dirname, "kappa.fits", "info.hdf")
            kappa = enmap.read_map(self.path+'/'+self.name+'/' + dirname)
            rho = enmap.read_map(self.path+'/'+self.name+'/' + rhofile)
            info = bunch.read(self.path+self.name+'/' +infofile)
     
            kappa = np.maximum(kappa, np.max(kappa)*1e-2) 
           
            ctime0 = info.ctime_ast
            
            tol = 1e-2
            r = 5 
            mask = kappa > np.max(kappa)*tol
            mask = mask.distance_transform(rmax=r) >= r
            rho   *= mask
            kappa *= mask
            
            if kappa[pol,:,:].at([0,0]) <= 1e-9:  
                continue #checks if in masked area

            if time_cut is not None:
                flag = False
                for j in range(len(time_cut)):
                    if (time_cut[j][0] < ctime0) and (ctime0 < time_cut[j][1]): 
                        flag = True
                        break
                if flag:  
                    continue
            

            if restrict_time and (ctime0 < restrict_time[0] or restrict_time[1] < ctime0): continue #restrict to only within a certain time range
            hour = (ctime0/3600)%24

            if time == 'night' and (11<hour<23): continue
            elif time == 'day' and (hour<11 or hour >23): continue
            #else: print('Something wrong with your time: {} and split {}'.format(hour, time))

            if verbose:
                ts = int(ctime0)  
                print(ctime0) 
                print(datetime.datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))

            adata = self.orbit(ctime0)
            ra_ast, dec_ast, delta_earth, delta_sun, ignore_ang = adata
            
            if weight_type == 'flux':
                cur_time = utils.ctime2djd(ctime0) 

                self.sun.compute(cur_time)

                alpha = compute_alpha(self.sun.ra, self.sun.dec, self.sun.earth_distance, ra_ast*np.pi/180, dec_ast*np.pi/180, delta_earth)
                alpha *= (180/np.pi)

                weight = 1/(delta_sun**(-2) * delta_earth**(-2) *10**(-0.004*alpha))
                if verbose: print('Weight = ', weight)
                
                
                
            elif weight_type == 'invflux':
                d_sun_0, d_earth_0 = self.semimajor,self.semimajor
                weight = 1 / (d_sun_0**2 * d_earth_0**2/(delta_earth**2*delta_sun**2))
                if verbose: print('Weight = ', weight)
            elif weight_type == 'paper_weight': 
           
                cur_time = utils.ctime2djd(ctime0) 
                
                self.sun.compute(cur_time)
                
                alpha = compute_alpha(self.sun.ra, self.sun.dec, self.sun.earth_distance, ra_ast, dec_ast, delta_earth)
                alpha *= (180/np.pi)
                
               
                weight = 1/(delta_sun**(-1/2) * delta_earth**(-2) *10**(-0.004*alpha)) 
                if verbose: 
                    print('Weight = ', weight) 
            elif weight_type == 'earth_only':
                weight = 1/(delta_earth**(-2))
                if verbose: print('Weight = ', weight)
            
            else:
                weight = 1
            adjust = 1 
            if freq_adjust:
                if int(freq) == 90: old_freq = 98e9
                elif int(freq) == 150: old_freq = 150e9
                elif int(freq) == 220: old_freq = 220e9
                adjust *= utils.dplanck(freq_adjust)/utils.dplanck(old_freq)
         
            if extern_calib:
                adjust *= extern_calib
            
            sn_map = (rho[pol,:,:]/weight)/np.sqrt( kappa[pol,:,:]/weight**2)
            #flux_map = (rho[pol,:,:]/weight) /  (kappa[pol,:,:]/weight**2)
            #if np.amax(np.abs(flux_map)) >10000: continue
            if np.median(np.abs(sn_map))>1.4 and weight_type != 'earth_only': 
                print('Bad Tile: S/n')
                continue            
           
            rho_stack += rho[pol,:,:]/weight
            kappa_stack += kappa[pol,:,:]/weight**2
            if verbose or weight_debug:
                flux = rho[pol,:,:].at([0,0]) /  kappa[pol,:,:].at([0,0])
                sn = rho[pol,:,:].at([0,0]) / np.sqrt(kappa[pol,:,:].at([0,0]))
                if verbose:
                    print('Rho \t Kappa \t Flux \t SN\t Scaled Flux\n')
                    print(rho[pol,:,:].at([0,0]), '\t', kappa[pol,:,:].at([0,0]), '\t', flux*adjust, '\t', sn, '\t', weight*flux)
                    print("LC comp: ", rho[pol,:,:].at([0,0]) /  kappa[pol,:,:].at([0,0]) * adjust)
                    if pol >0 and sn >3:
                        print("High S/n")
            rho_un += rho[pol,:,:]
            kappa_un+= kappa[pol,:,:]
            
            rho_weight += weight
            kappa_weight += weight**2  

            if plot and verbose:
                flux_map = (rho[pol,:,:]/weight) /  (kappa[pol,:,:]/weight**2)
             
                plt.scatter(40,40, marker = '+', color = 'r')
                print('Stack weight: ', weight)
                print('D_earth: ', delta_earth)
                print('D_sun: ', delta_sun)
                plt.title(datetime.datetime.utcfromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
                plt.imshow(flux_map)
                plt.colorbar()
                path = '/scratch/r/rbond/jorlo/actxminorplanets/sigurd/plots/stamps/{}/'.format(self.name)
                if not os.path.exists(path):
                    os.makedirs(path)
                plt.savefig(path + '{}_{}_{}_{}.pdf'.format(int((ctime0/3600)%24), pa, freq, int(ctime0)))
                plt.savefig(path + '{}_{}_{}_{}.png'.format(int((ctime0/3600)%24), pa, freq, int(ctime0))) 
                
                print("Kappa: ", np.sqrt( kappa[pol,:,:]/weight**2).at([0,0]))
                print('S/n: ', sn_map.at([0,0]))        
                print("Median S/n: ", np.median(np.abs(sn_map)))
                print("Mean Flux ", np.mean(flux_map))
                print("Mean abs Flux ", np.mean(np.abs(flux_map)))
                print("Max flux: ", np.amax(np.abs(flux_map)))
                plt.show()
                plt.close()

                plt.scatter(40,40, marker = '+', color = 'r')
                plt.title('Kappa')
                plt.imshow(kappa[pol,:,:]/weight**2)
                plt.colorbar()
                plt.show()
                plt.close()
                
                plt.scatter(40,40, marker = '+', color = 'r')
                plt.title('SNR')
                plt.imshow(sn_map)
                plt.colorbar()
                plt.show()
                plt.close()


                hours.append(int((ctime0/3600)%24)) 
                
            if lightcurve:
                fluxes_lc.append(rho[pol,:,:].at([0,0]) /  kappa[pol,:,:].at([0,0]) * adjust)
                errs_lc.append(1/kappa[pol,:,:].at([0,0])**(1/2) * adjust) 
                times_lc.append(ctime0) 
                Fs_lc.append(weight)
                rho_lc.append(rho[pol, :, :].at([0,0]))
                kappa_lc.append(kappa[pol, :, :].at([0,0]))

        try: 
            stack = (rho_stack/kappa_stack) #/ (kappa_stack/kappa_weight)
        except ZeroDivisionError:
            print("Zero Division Error")  
            stack = np.zeros((81,81))
            self.map_dict[time][pa][freq]['flux'] = stack
            self.map_dict[time][pa][freq]['snr'] = stack

            self.flux_dict[time][pa][freq]['flux'] = 0
            self.flux_dict[time][pa][freq]['var'] = 0
            return 

        
        plt.scatter(40,40, marker = '+', color = 'r') 
        plt.imshow(stack)
        plt.colorbar()
        plt.title('{} Flux of {} from PA {} at Freq {}, Flux {}'.format(time, self.name, pa, freq,stack.at([0,0])*adjust)) 
        plt.savefig('/scratch/r/rbond/jorlo/actxminorplanets/sigurd/plots/stacks/{}_{}_{}_{}.pdf'.format(time, self.name, pa, freq)) 
        if plot: 
            print('{} Flux of {} from PA {} at Freq {}, Flux {}'.format(time, self.name, pa, freq, stack.at([0,0])*adjust))
            plt.show()
        plt.close()
        
        if plot and verbose:
            plt.scatter(40,40, marker = '+', color = 'r')
            plt.imshow(kappa_stack)
            plt.colorbar()
            plt.title('Var of {} from PA {} at Freq {}'.format(self.name, pa, freq))
            plt.show()
            plt.close()


            fig, ax = plt.subplots(1,1)
            plt.hist(hours)
            plt.xlabel('Hour of Obs')
            plt.ylabel('# Obs')
            ax.axvline(11, color='black')
            ax.axvline(23, color='black')
            plt.title('Hours of Observation for {}'.format(pa))
            plt.savefig('./plots/hour_hist_{}.pdf'.format(pa))
            plt.show()
         
        self.map_dict[time][pa][freq]['flux'] = stack * adjust
        self.map_dict[time][pa][freq]['snr'] = rho_stack/np.sqrt(kappa_stack)         
        
        self.flux_dict[time][pa][freq]['flux'] = stack.at([0,0]) * adjust
        self.flux_dict[time][pa][freq]['var'] = 1/np.sqrt(np.mean(kappa_stack)) * adjust
        
        plt.imshow(kappa_stack)
        plt.colorbar()
        #plt.show()
        plt.close()

        if weight_debug:
            debug_dict = {'sn':sns,'iso':isos, 'hour':hours, 'weights':weights, 'alphas':alphas, 'r_sun':r_suns, 'r_earth':r_earths, 'fluxes': fluxes}
            with open('/home/r/rbond/jorlo/dev/minorplanets/pks/{}_debug_{}_{}_{}.pk'.format(self.name, time, pa, freq), 'wb') as f: 
                pk.dump(debug_dict, f)           
        if lightcurve: 
            lc_dict = {'flux': fluxes_lc, 'err':errs_lc, 'time':times_lc, 'F':Fs_lc, "rho":rho_lc, "kappa":kappa_lc}
            if weight_type == 1:
                name = '{}_unweighted_lc_{}_{}_{}'.format(self.name, time, pa, freq)
            else:
                cur_pol = [i for i in pol_dict if pol_dict[i]==pol] #Gets key from pol index
                name = '{}_{}_{}_lc_{}_{}_{}'.format(pol, self.name, weight_type, time, pa, freq)
            with open('/scratch/r/rbond/jorlo/actxminorplanets/sigurd/lightcurves/'+name+".pk", 'wb') as f:
                pk.dump(lc_dict, f)

            c1 = fits.Column(name = "Time", array = np.array(times_lc), format = "D", unit = "Unix Timestamp")
            c2 = fits.Column(name = "Flux", array = np.array(fluxes_lc), format = "D", unit = "mJy")
            c3 = fits.Column(name = "FluxUncertainty", array = np.array(errs_lc), format = "D", unit = "mJy")
            c4 = fits.Column(name = "Weight", array = np.array(Fs_lc), format = "D", unit = "unitless")
            #cols = fits.ColDefs([c1, c2, c3])
            hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
            
            with open('/scratch/r/rbond/jorlo/actxminorplanets/sigurd/lightcurves/'+name+".fits", 'wb') as f: 
                
                hdu.writeto(f, overwrite = True)


    def make_all_stacks(self, weight_type = 1, pol='T', directory = '/scratch/r/rbond/jorlo/actxminorplanets/sigurd/', plot = False, verbose = False, weight_debug = False, lightcurve = False, freq_adjust = None, extern_calib = None, time_cut = None):
        for pa_key in self.map_dict['night'].keys(): 
            for freq_key in self.map_dict['night'][pa_key].keys():
                
                if type(freq_adjust) == dict: #If freq_adjust not specified, set it to None
                    cur_freq_adjust = freq_adjust[pa_key][freq_key]
                else:
                    cur_freq_adjust = None

                if type(extern_calib) == dict:
                    cur_extern_calib = extern_calib[pa_key][freq_key][0]
                elif type(extern_calib) == float or type(extern_calib) == int:
                    cur_extern_calib = extern_calib
                else:
                    cur_extern_calib = None

                self.make_stack(pa = pa_key, freq = freq_key, pol=pol, plot = plot, verbose = verbose, weight_type = weight_type, time = 'night', weight_debug = weight_debug, 
                                lightcurve = lightcurve, freq_adjust = cur_freq_adjust, extern_calib = cur_extern_calib, time_cut = time_cut)
                self.make_stack(pa = pa_key, freq = freq_key, pol=pol, plot = plot, verbose = verbose, weight_type = weight_type, time = 'day', weight_debug = weight_debug,
                                lightcurve = lightcurve, freq_adjust = cur_freq_adjust, extern_calib = cur_extern_calib, time_cut = time_cut)
        if pol == 'T':
            if weight_type == 1:
                fname = "fluxes/{}_unweight_flux_dict.pk".format(self.name)
            else:
                fname = "fluxes/{}_{}_weight_flux_dict.pk".format(self.name, weight_type)
            with open(directory+fname, 'wb') as f:
                print('Flux dict written for ', self.name)
                pk.dump(self.flux_dict, f)
        else:
            with open(directory+'fluxes/{}_pol_{}_dict.pk'.format(self.name, pol), 'wb') as f:
                #print(f)
                pk.dump(self.flux_dict, f)

    def save_stack(self, directory):
        aster_dict = {'flux_stack':self.flux_stack, 'fstack':self.fstack, 'kstack':self.kstack, 'flux_scale':self.flux_scale}
        with open(directory+'{}_stamp.pk'.format(str(self.eph_table['targetname'][0]).replace(' ', '_').replace('(','').replace(')','')), 'wb') as f:
            pk.dump(aster_dict, f)
            
    
            
    def plot_stack(self, directory = None, scale = None):
        plt.imshow(self.flux_stack)
        plt.xlabel('ra')
        plt.ylabel('dec')
        
        plt.title('Plot of {} Stack'.format(self.eph_table['targetname'][0]))
        if directory is not None:
            plt.savefig(directory + '{}_stamp.pdf'.format(str(self.eph_table['targetname'][0]).replace(' ', '_').replace('(','').replace(')','')))
        #plt.show()
        
    def dif_pa4_pa5(self, weight_type = None, time = 'night', pol = 'T', verbose = False):
        pol_dict = {'T':0, 'Q':1, 'U':2}
        pol = pol_dict[pol]


        pas = ['pa4', 'pa5']

        pa_dict = {'pa4':{'flux':[], 'time':[], 'kappa':[]},
                   'pa5':{'flux':[], 'time':[], 'kappa':[]}}
        if weight_type:
            pa_dict['pa4']['weight'] = []
            pa_dict['pa5']['weight'] = []

            pa_dict['pa4']['weighted_flux'] = []
            pa_dict['pa5']['weighted_flux'] = []

 
        for i, dirname in enumerate(os.listdir(path=self.path+'/'+self.name+'/')):
            
            for pa in pas:
                if dirname.find('kappa.fits') == -1: continue #Just look at the kappa maps
                if dirname.find(pa) == -1: continue
                if dirname.find('f150') == -1: continue #Check each file to see if it's the specified freq
                if verbose: print('In map: ', dirname)
                #hdu = fits.open(path+'/'+name+'/' + dirname)
                #data = hdu[0].data[pol]
                rhofile    = utils.replace(dirname, "kappa.fits", "rho.fits")
                infofile = utils.replace(dirname, "kappa.fits", "info.hdf")

                kappa = enmap.read_map(self.path+'/'+self.name+'/' + dirname)
                rho = enmap.read_map(self.path+'/'+self.name+'/' + rhofile)
                info = bunch.read(self.path+'/'+self.name+'/' +infofile)
                ctime0   = np.mean(info.period)

                hour = (ctime0/3600)%24

                if time == 'night' and (11<hour<23): continue
                elif time == 'day' and (hour<11 or hour >23): continue

                adata = self.orbit(ctime0)
                ra_ast, dec_ast, delta_earth, delta_sun, ignore_ang = adata
                flux = rho[pol,:,:].at([0,0]) /  kappa[pol,:,:].at([0,0])
                pa_dict[pa]['flux'].append(flux)
                pa_dict[pa]['time'].append(ctime0)
                pa_dict[pa]['kappa'].append(kappa[pol,:,:].at([0,0]))
                if weight_type == 'spt':
                    adata = self.orbit(ctime0)
                    ra_ast, dec_ast, delta_earth, delta_sun, ignore_ang = adata

                    cur_time = utils.ctime2djd(ctime0)

                    self.sun.compute(cur_time)

                    alpha = compute_alpha(self.sun.ra, self.sun.dec, self.sun.earth_distance, ra_ast, dec_ast, delta_earth)
                    alpha *= (180/np.pi)

          
                    weight = 1/(delta_sun**(-1/2) * delta_earth**(-2) *10**(-0.004*alpha))
           
                    pa_dict[pa]['weight'].append(weight)
             
                    pa_dict[pa]['weighted_flux'].append((rho[pol,:,:]/weight).at([0,0])/(kappa[pol,:,:]/weight**2).at([0,0]))

        return pa_dict

#if args.index == -1:
#        name     = args.name 

#else:
#	desig, name, semimajor = get_desig(args.index)


#ast = minorplanet(name)
#ast.make_all_stacks(weight_type='flux')





