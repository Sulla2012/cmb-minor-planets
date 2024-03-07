import argparse, os
parser = argparse.ArgumentParser()

parser.add_argument("ifiles", nargs="+") 
parser.add_argument("odir")
parser.add_argument("-A", "--astinfo",       type=str,   default=None)
parser.add_argument("-N", "--name",  type=str,   default=None)
parser.add_argument("-r", "--rad",   type=float, default=10.0)
parser.add_argument("-p", "--pad",   type=float, default=30.0)
parser.add_argument("-t", "--tol",   type=float, default=0.1)
parser.add_argument("-l", "--lknee", type=float, default=1500)
parser.add_argument("-a", "--alpha", type=float, default=3.5)
parser.add_argument("-b", "--beam",  type=float, default=0)
parser.add_argument("-v", "--verbose", default=1, action="count")
parser.add_argument("-e", "--quiet",   default=0, action="count")
parser.add_argument("-i", "--index", type=int, default=-1)
parser.add_argument("-f", "--freq", type=str, default=None)
args = parser.parse_args()
import numpy as np, ephem
from numpy.lib import recfunctions
from pixell import utils, enmap, bunch, reproject, colors, coordinates#, mpi
from scipy import interpolate, optimize
import glob
import matplotlib.pyplot as plt
import matplotlib
import pickle as pk

matplotlib.use('pdf')

from pathlib import Path
from asteroid_utils_pixell import minorplanet, get_desig

try:
    import mpi4py.rc

    mpi4py.rc.threads = False
    from mpi4py import MPI

    print("mpi4py imported")
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nproc = comm.Get_size()
    print("nproc:, ", nproc)
    if nproc > 1:
        have_mpi = True
    else:
        have_mpi = False
except:
    MPI = None
    comm = None
    have_mpi = False
    myrank = 0
    nproc = 1

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

comm    = MPI.COMM_WORLD
verbose = args.verbose - args.quiet

# Radius of output thumbnails
r_thumb = args.rad*utils.arcmin
# ... and temporary radius to use before we know the true position
r_full  = r_thumb + args.pad*utils.arcmin

if args.index == -1:
	name     = args.name or utils.replace(os.path.basename(args.astinfo), ".npy", "").lower()

else:
	desig, name, semimajor = get_desig(args.index)

lknee    = args.lknee
alpha    = -args.alpha
beam     = args.beam*utils.fwhm*utils.arcmin
site     = coordinates.default_site
site_pos = np.array([site.lon,site.lat])*utils.degree
time_tol  = 600
time_sane = 3600*24

# Expand any globs in the input file names

ifiles  = sum([sorted(utils.glob(ifile)) for ifile in args.ifiles],[]) 
# Read in and spline the orbit
#print(ifiles)
if args.index == -1:
	info    = np.load(args.astinfo).view(np.recarray)
else:
	info    = np.load('/gpfs/fs0/project/r/rbond/sigurdkn/actpol/ephemerides/objects/{}.npy'.format(name.capitalize())).view(np.recarray)

orbit   = interpolate.interp1d(info.ctime, [utils.unwind(info.ra*utils.degree), info.dec*utils.degree, info.r, info.ang*utils.arcsec], kind=3)
utils.mkdir(args.odir)

print('Starting {}'.format(name))
print('Index {}'.format(args.index))

for fi in range(myrank, len(ifiles), nproc):
    ifile    = ifiles[fi]
    infofile = utils.replace(ifile, "map.fits", "info.hdf")
    tfile    = utils.replace(ifile, "map.fits", "time.fits")
    varfile    = utils.replace(ifile, "map.fits", "ivar.fits")
    rhofile    = utils.replace(ifile, "map.fits", "rho.fits")
    kfile    = utils.replace(ifile, "map.fits", "kappa.fits")
    
    if args.freq:
        if args.freq not in kfile: continue
    
    ofname   = "%s/%s/%s_%s" % (args.odir, name.capitalize(), name, os.path.basename(ifile))
    varname   = "%s/%s/%s_%s" % (args.odir, name.capitalize(), name, os.path.basename(varfile))
    rhoname   = "%s/%s/%s_%s" % (args.odir, name.capitalize(), name, os.path.basename(rhofile))
    kname   = "%s/%s/%s_%s" % (args.odir, name.capitalize(), name, os.path.basename(kfile))
    infoname = "%s/%s/%s_%s" % (args.odir, name.capitalize(), name, os.path.basename(infofile))
    
    info     = bunch.read(infofile)
    # Get the asteroid coordinates
    ctime0   = np.mean(info.period)
    adata0   = orbit(ctime0)
    ast_pos0 = utils.rewind(adata0[1::-1])
    message  = "%.0f  %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f" % (info.t, ast_pos0[1]/utils.degree, ast_pos0[0]/utils.degree, info.box[0,1]/utils.degree, info.box[1,1]/utils.degree, info.box[0,0]/utils.degree, info.box[1,0]/utils.degree)
    # Check if we're in bounds
    if not in_box(info.box, ast_pos0):
        if verbose >= 3:
    	    print(colors.lgray + message + " outside" + colors.reset)
    	    continue
    #Ok, should try to read in this map. Decide on
    # bounding box to read in
    
    
    full_box  = make_box(ast_pos0, r_full)
    # Read it and check if we have enough hits in the area we want to use
    try:
        tmap = enmap.read_map(tfile, box=full_box)
        tmap[tmap!=0] += info.t
    except (TypeError, FileNotFoundError):
        print("Error reading %s. Skipping" % ifile)
        continue
    # Break out early if nothing is hit
    if np.all(tmap == 0): continue
    	#if verbose >= 2:
    	#	print(colors.white + message + " unhit" + colors.reset)
    	#continue
    # Figure out what time the asteroid was actually observed
    ctime, err = calc_obs_ctime(orbit, tmap, ctime0)
    if err > time_tol or abs(ctime-ctime0) > time_sane:
        if verbose >= 2:
    	    print(colors.white + message + " time" + colors.reset)
    		#print(err)
    		#print(ifile)
        continue
    # Now that we have the proper time, get the asteroids actual position
    adata    = orbit(ctime)
    ast_pos  = utils.rewind(adata[1::-1]).reshape((2,))
    # optionally transform to topocentric here. ~0.1 arcmin effect
    thumb_box = make_box(ast_pos, r_thumb)
    
    
    # Read the actual data
    try:
        imap = enmap.read_map(ifile, box=full_box)
        var = enmap.read_map(varfile, box=full_box)
        rho = enmap.read_map(rhofile, box=full_box)
        kap = enmap.read_map(kfile, box=full_box)
    except Exception as e:
        print(colors.red + str(e) + colors.reset)
        continue
    
    if np.mean(imap.submap(thumb_box) == 0) > args.tol:
        if verbose >= 2:
            print(colors.white + message + " unhit" + colors.reset)
            continue
    # Filter the map
    wmap     = filter_map(imap, lknee=lknee, alpha=alpha, beam=beam)
    # And reproject it
    omap = reproject.thumbnails(wmap, ast_pos, r=r_thumb)
    enmap.write_map(ofname, omap)
    
    vmap = reproject.thumbnails(var, ast_pos, r=r_thumb)
    enmap.write_map(varname, vmap)
    
    rmap = reproject.thumbnails(rho, ast_pos, r=r_thumb)
    enmap.write_map(rhoname, rmap)
    
    kmap = reproject.thumbnails(kap, ast_pos, r=r_thumb)
    enmap.write_map(kname, kmap)
    
    info.ctime_ast = ctime[0]
    info.ctime_err = err 
    bunch.write(infoname, info)
    print('info: ', infoname)
    if verbose >= 1:
        print(colors.lgreen + message + " ok" + colors.reset)
    	#print("time: ", ctime[0])
    
		
print('Finished')

if myrank == 0:
    with open('/home/r/rbond/jorlo/dev/minorplanets/calibs_dict_ACTxACT-dr4.pkl', 'rb') as f:
        dr6_calib = pk.load(f)
    calib_dict = {'pa4':{'150':dr6_calib['dr6_pa4_f150']['calibs'], '220':dr6_calib['dr6_pa4_f220']['calibs']},
                  'pa5':{'090':dr6_calib['dr6_pa5_f090']['calibs'], '150':dr6_calib['dr6_pa5_f150']['calibs']},
                  'pa6':{'090':dr6_calib['dr6_pa6_f090']['calibs'], '150':dr6_calib['dr6_pa6_f150']['calibs']}
             }
    freq_dict = {'pa4':{'150':148.7e9, '220':227.2e9},
                 'pa5':{'090':96.7e9,  '150':149.5e9},
                 'pa6':{'090':95.5e9,  '150':148.1e9}}#Effective centers from Hasslefield paper 2022 in prep
    ast = minorplanet(name, semimajor)
    ast.make_all_stacks(weight_type='paper_weight', lightcurve = True, freq_adjust = freq_dict, extern_calib = calib_dict)
    print("Made T for ", name)
    ast.make_all_stacks(weight_type='paper_weight', pol="U", lightcurve = True, freq_adjust = freq_dict, extern_calib = calib_dict)
    print("Made U for", name)
    ast.make_all_stacks(weight_type='paper_weight', pol="Q", lightcurve = True, freq_adjust = freq_dict, extern_calib = calib_dict)
    print("All done for", name)

