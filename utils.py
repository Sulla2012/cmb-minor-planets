import numpy as np
import ephem

from pixell import utils, enmap, bunch, reproject, mpi
from scipy import interpolate, optimize
import glob

from pathlib import Path

import os

import pickle as pk

def in_box(box, point): #checks if points are inside box or not
        box   = np.asarray(box)
        point = np.asarray(point)
        point[1] = utils.rewind(point[1], box[0,1])
        # This assumes reverse RA axis in box
        return point[0] > box[0,0] and point[0] < box[1,0] and point[1] > box[1,1] and point[1] < box[0,1]

def make_box(point, rad): #making box
        box     = np.array([point - rad, point + rad])
        box[:,1]= box[::-1,1] # reverse ra
        return box

def filter_map(map, lknee=3000, alpha=-3, beam=0): #filtering map somehow (FFT)
        fmap  = enmap.fft(map)
        l     = np.maximum(0.5, map.modlmap())
        filter= (1+(l/lknee)**alpha)**-1
        if beam:
                filter *= np.exp(-0.5*l**2*beam**2)
        fmap *= filter
        omap  = enmap.ifft(fmap).real
        return omap

def calc_obs_ctime(orbit, tmap, ctime0):
        def calc_chisq(x):
                ctime = ctime0+x
                try:
                        adata = orbit(ctime)
                        mtime = tmap.at(adata[1::-1], order=1)
                except ValueError:
                        mtime = 0
                return (mtime-ctime)**2
        ctime = optimize.fmin_powell(calc_chisq, 0, disp=False)+ctime0
        err   = calc_chisq(ctime-ctime0)**0.5
        return ctime, err

def calc_sidereal_time(lon, ctime):
        obs      = ephem.Observer()
        obs.lon  = lon
        obs.date = utils.ctime2djd(ctime)
        return obs.sidereal_time()

def geocentric_to_site(pos, dist, site_pos, site_alt, ctime):
        """Given a geocentric position pos[{ra,dec},...] and distance dist [...]
        in m, transform it to coordinates relative to the given site, with
        position pos[{lon,lat}] and altitude alt_site in m, returns the
        position observed from the site, as well as the distance from the site in m"""
        # This function isn't properly debugged. I should check the sign of
        # the sidereal time shift. But anyway, this function is a 0.2 arcmin
        # effect in the asteroid belt.
        sidtime    = calc_sidereal_time(site_pos[0], ctime)
        site_radec = np.array([site_pos[0]+sidtime*15*utils.degree,site_pos[1]])
        vec_site   = utils.ang2rect(site_radec)*(utils.R_earth + site_alt)
        vec_obj    = utils.ang2rect(pos)*dist
        vec_rel  = vec_obj-vec_site
        dist_rel = np.sum(vec_rel**2,0)**0.5
        pos_rel  = utils.rect2ang(vec_rel)
        return pos_rel, dist_rel

def get_index(name, verbose = False):
    '''
 
    Given an asteroid name, returns the ACT internal asteroid index, which is a rough proxy for 200GHz flux
    Parameters
    ---------- 
    
    name: str
          name of object

    Returns
    -------
    
    desig: int 
      index of object in asteroids.pk file
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
        if verbose: print('Object not in current data set')
        return 999999

def get_orbits(ast_dir, max_idx = 500):
    """

    This is not the ideal way to do this. I'd prefer to make an x-array of times and a y-array of [ra, dec, r, ang] and call interp on that, the issue is that the asteroids are not all sampled with the same ammount of points, so that the x and y arrays end up ragged, and numpy doesn't like interpolating that. There may be a work around. For now just returning a list of interp obejcts

    """
    orbits = []
    
    for asteroid in os.listdir(ast_dir):
        if get_index(asteroid.strip('.npy')) > max_idx: continue #If the index of the asteroid is greater than the specified max index, continue since the asteroid is too faint to matter

        info = np.load(ast_dir + asteroid).view(np.recarray)
        orbit = interpolate.interp1d(info.ctime, 
            [utils.unwind(info.ra*utils.degree), info.dec*utils.degree, info.r, info.ang*utils.arcsec],
            kind=3)
        orbits.append(orbit)

    return orbits

def check_map_ast(infofile, ast_dir):
    """
    Checks a depth 1 map to see if any asteroids are present in it

    Parameters
    ----------
    infofile: pixel.bunch.Bunch
        pixel bunch file containing at minimum the 'box' attribute defining the depth-1 map bounding box
    ast_dir: str
        string of directory containing asteroid ephemerides to check #TODO probably a better way to specify which asteroids to check, we only need to look at bright ones

    Returns
    -------
    
    """

    for asteroid in os.listdir(ast_dir):
        info = np.load(ast_dir + asteroid)
        orbit = interpolate.interp1d(info.ctime, [utils.unwind(info.ra*utils.degree), info.dec*utils.degree, info.r, info.ang*utils.arcsec], kind=3) 
            
def check_loc_ast(ra, dec, time, ast_dir, tol = 2*utils.arcmin, max_idx = 500):
    """
    Checks a specific location to see if there is an asteroid within tol at the specified time

    Parameters
    ----------
    ra: float
        ra of interest in radians
    dec: float
        dec of interest in radians
    time: float
        time of interest in unix time
    ast_dir: str
        string of directory containing asteroid ephemerides to check #TODO probably a better way to specify which asteroids to check, we only need to look at bright ones
    tol: float
        minimum allowed distance between asteroid and point of interest
    max_idx: int
        maximum ACT interal asteroid index to consider. Roughly corresponds to 150GHz flux

    Returns
    -------
    near_ast: bool
        True if there is an asteroid within tol at time, else false #TODO may want to include more detailed info here, i.e. ast by ast breakdown
    """

    ast_orbits = get_orbits(ast_dir, max_idx = max_idx)
    ast_locs = np.zeros((len(ast_orbits), 4)) #an individual interp returns [ra, dec, r, ang] so second dim is 4 
    for i, orbit in enumerate(ast_orbits):
        ast_locs[i] = orbit(time)
    near_ast = np.any(np.array([(np.sqrt((ra-ast_locs[...,0])**2 + (dec-ast_locs[...,1]**2))<tol)]))

    return near_ast 

    
