#!/usr/bin/env python


"""
Script to calculate the array factor and determine the tied-array beam of the MWA.

Author: Bradley Meyers
Date: 2016-12-8
"""



# numerical and maths modules
import numpy as np
from astropy.coordinates import SkyCoord,EarthLocation,AltAz
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from astropy.constants import c,k_B

#utility and processing modules
import os,sys
from mpi4py import MPI
import argparse
from mwapy import ephem_utils,metadata
from mwapy.pb import primary_beam as pb
import urllib
import urllib2
import json



def getmeta(service='obs', params=None):
    """
    getmeta(service='obs', params=None)
    
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.
    """
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''
    #Validate the service name
    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return
    #Get the data
    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return
    #Return the result dictionary
    return result


def get_delay_steps(obs):
    #for line in txtfile:
    beam_meta_data = getmeta(service='obs', params={'obs_id':obs})
    obn = beam_meta_data[u'obsname']
    ra = beam_meta_data[u'metadata'][u'ra_pointing']
    dec = beam_meta_data[u'metadata'][u'dec_pointing']
    dura = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
    mode = beam_meta_data[u'mode']
    Tsky = beam_meta_data[u'metadata'][u'sky_temp']
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']

    minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)
    #print "{0}  {1:<15}  {2:3.2f}   {3:4.2f}    {4:>4}      {5}   {6}".format(ob,obn,ra,dec,dura,centrefreq,Tsky)
    #print xdelays

    return [obs,ra,dec,dura,xdelays,centrefreq,channels]


def getTileLocations(obsid,flags=[],fdir="."):
    """
    Function grab the MWA tile locations for any given observation ID. Downloads the relevant metafits file from the database, saves it as <obdID>_metafits_ppds.fits.
    
    Input:
      obsid - the GPS observation ID
      flags - RTS tile flags (i.e. the first entry in the metafits correspond to "tile 0", irrespective of what the antenna name is) 

    Return:
      a list of lists containing the following:
        list[0] = a list of tile positions East of the array centre
        list[1] = a list of tile positions North of the array centre
        list[2] = a list of tile heights about sea-level 
    """

    f = fits.open('{0}/{1}_metafits_ppds.fits'.format(fdir,obsid))        
    #names = f[1].data['Antenna'][::2]
    east = f[1].data['East'][::2]
    north = f[1].data['North'][::2]
    height = f[1].data['Height'][::2] # height above sea-level
    # MWA array centre height above sea-level
    mwacentre_h = 377.827
    height = height - mwacentre_h

    #flag the tiles from the x,y,z positions
    east = np.delete(east,flags)
    north = np.delete(north,flags)
    height = np.delete(height,flags)
            
    return east,north,height


def get_obstime_duration(obsid,fdir="."):
    """
    Funciton to grab the recorded start-time and duration of the observation
    
    Input:
      obsid - the GPS observation ID

    Return:
      a list containing the folloing two items (in order):
        list[0] = observation starting time in UTC
        list[1] = observation duration in seconds
    """
    # metafits file will already have been downloaded
    f = fits.open('{0}/{1}_metafits_ppds.fits'.format(obsid,fdir))
    
    return [f[0].header['DATE-OBS'],f[0].header['EXPOSURE']]


def getTargetAZZA(ra,dec,time,lat=-26.7033,lon=116.671,height=377.827):
    """
    Function to get the target position in alt/az at a given EarthLocation and Time.
    
    Default lat,lon,height is the centre of  MWA.

    Input:
      ra - target right ascension in astropy-readable format
      dec - target declination in astropy-readable format
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target azimuth in radians
        list[1] = target zenith angle in radians
        list[2] = target azimuth in degrees
        list[3] = target zenith angle in degrees
    """
    #print "Creating EarthLocation for: lat = {0} deg, lon = {1} deg".format(lat,lon)
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    
    #print "Creating SkyCoord for target at (RA,DEC) = ({0},{1})".format(ra,dec)
    coord = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    #print "Target at: ({0},{1}) deg".format(coord.ra.deg,coord.dec.deg)
    
    obstime = Time(time)
    #print "Observation time: {0}".format(obstime.iso)
    
    #print "Converting to Alt,Az for given time and observatory location..."
    altaz = coord.transform_to(AltAz(obstime=obstime,location=location))
    #print "Target (Alt,Az) = ({0},{1}) deg".format(altaz.alt.deg,altaz.az.deg)
    
    #print "Converting to (Az,ZA)"
    az = altaz.az.rad 
    azdeg = altaz.az.deg
     
    za = np.pi/2 - altaz.alt.rad
    zadeg = 90 - altaz.alt.deg
    
    #print "Target (Az,ZA) = ({0},{1}) deg".format(azdeg,zadeg)

    return [az,za,azdeg,zadeg]


def calcWaveNumbers(wl,p,t):
    """
    Function to calculate the 3D wavenumbers for a given wavelength and az/za grid.

    Input:
      wl - central wavelength for the observation
      p - azimuth/phi (either a scalar or an array)
      t - zenith angle/theta (either a scalar or an array)
    
    Return:
      [kx,ky,kz] - the 3D wavenumbers 
    """
    # the standard equations are:
    #    a = 2 * pi / lambda
    #    kx = a * sin(theta) * cos(phi)
    #    ky = a * sin(theta) * sin(phi)
    #    kz = a * cos(theta)
    # this is assuming that theta,phi are in the convention from Sutinjo et al. 2015
    #    i.e. phi = pi/2 - az
    kx = (2*np.pi/wl) * np.sin(t) * np.cos(p)
    ky = (2*np.pi/wl) * np.sin(t) * np.sin(p)
    kz = (2*np.pi/wl) * np.cos(t)

    return [kx,ky,kz]


# Make generator functions for Azimuth and Zenith so we don't have to create the whole
# meshgrid at once and then store it.
# This is far more memory efficient that storing the entire AZ and ZA planes in memory, 
# especially with the required resolution.
def genAZZA(start,stop,step,end=False):
    """
    Generator function to use for iterating over angles (both ZA and Azimuth).

    Input:
      start - angle to start iteration from 
      stop - angle to finish iteration before
      step - step size between iterations
      end - return the "stop" parameter to make range [start,stop] rather than [start,stop)

    Return:
      None - is a generator
    """
    i = 0
    num = int((stop-start)/step)
    while i<num:
        yield start+i*step
        i += 1
    if end:
        yield stop
    return
    

def createArrayFactor(targetAZ,targetZA,obsid,delays,time,obsfreq,eff,flagged_tiles,theta_res,phi_res,
                        coplanar,start,end,azstart,azstop,write,radius,strip):
    """
    Primary function to calcaulte the array factor with the given information.

    Input:
      targetRA,targetDEC - the desired target Right Ascension and Declination
      obsid - the observation ID for which to create the phased array beam
      delays - the beamformer delays required for point tile beam (should be a set of 16 numbers)
      obstime -  the time at which to evaluate the target position (as we need Azimuth and Zenith angle, which are time dependent)
      obsfreq - the centre observing frequency for the observation
      eff - the array efficiency (frequency and pointing dependent, currently require engineers to calculate for us...)
      flagged_tiles - the flagged tiles from the calibration solution (in the RTS format)
      theta_res - the zenith angle resolution in degrees
      phi_res - the azimuth resolution in degrees
      coplanar - whether or not to assume that the z-components of antenna locations is 0m
      start - the starting angle for the ZA band to calcualted across
      end - the end angle for the ZA band to be calculated across
	  radius - radius of square to be calculated, If None will be all sky

    Return:
      results -  a list of lists cotaining [ZA, Az, beam power], for all ZA and Az in the given band
    """    
    # convert frequency to wavelength
    obswl = obsfreq/c.value
    
    # get the target azimuth and zenith anglpython /group/mwaops/nswainston/pabeam/grid.py --plot -o 1163853320 -p "00:24:30.00 -72:04:30.0"e in radians and degrees
    # these are define in the normal sense: za = 90 - elevation, az = angle east of North (i.e. E=90)
    #targetAZ,targetZA,targetAZdeg,targetZAdeg = getTargetAZZA(targetRA,targetDEC,time)

    # calculate the target (kx,ky,kz)
    #target_kx,target_ky,target_kz = calcWaveNumbers(obswl,np.pi/2-targetAZ,targetZA)
    targetAZ = np.pi/2-targetAZ
    target_kx,target_ky,target_kz = calcWaveNumbers(obswl,targetAZ,targetZA)


    results = []

    # is this the last process?
    lastrank = rank == size-1

    # we will also calculate the beam area contribution of this part of the sky
    omega_A = 0
    array_factor_max = -1

    # for each ZA "pixel", 90deg inclusive
    for za in genAZZA(start,end,np.radians(theta_res),end=lastrank):
        # for each Az "pixel", 360deg not included
        for az in genAZZA(azstart,azstop,np.radians(phi_res)):
            
            # calculate the relevent wavenumber for (theta,phi)
            #kx,ky,kz = calcWaveNumbers(obswl,np.pi/2-az,za)
            kx,ky,kz = calcWaveNumbers(obswl,np.pi/2-az,za)
            array_factor = 0+0.j
            # determine the interference pattern seen for each tile
            for x,y,z in zip(xpos,ypos,zpos):
                ph = kx*x+ky*y+kz*z
                ph_target = target_kx*x+target_ky*y+target_kz*z
                array_factor += np.cos(ph-ph_target) + 1.j*np.sin(ph-ph_target)
            # normalise to unity at pointing position
            array_factor /= len(xpos)
            
            # keep track of maximum value calculated
            if array_factor > array_factor_max: array_factor_max = array_factor
    
            # calculate the tile beam at the given Az,ZA pixel
            
            print "za ",za 
            print "az ",az
            print "freq ",obsfreq
            print "delays ",delays
            
            # keep track of maximum value calculated
            if np.abs(array_factor)**2 > array_factor_max: array_factor_max = np.abs(array_factor)**2

            # calculate the tile beam at the given Az,ZA pixel
            #tile_xpol,tile_ypol = pb.MWA_Tile_full_EE([[za]],[[az]],freq=obsfreq,delays=[delays,delays],power=True,zenithnorm=True,interp=False)
            tile_xpol,tile_ypol = pb.MWA_Tile_analytic(za,az,freq=obsfreq,delays=delays,
                                                        power=True,zenithnorm=True)
            tile_pattern = (tile_xpol+tile_ypol)/2.0

            # calculate the phased array power pattern 
            phased_array_pattern = tile_pattern * np.abs(array_factor)**2          
            #phased_array_pattern = tile_pattern[0][0] * np.abs(array_factor)**2 # indexing due to tile_pattern now being a 2-D array           
        
            # append results to a reference list for later    
            results.append([np.degrees(za),np.degrees(az),phased_array_pattern])
        
            # add this contribution to the beam solid angle
            omega_A += np.sin(za) * np.abs(array_factor)**2 * np.radians(theta_res) * np.radians(phi_res)

    print "worker {0}, array factor maximum = {1}".format(rank,array_factor_max)

    # write a file based on rank of the process being used
    if write:
        print "writing file from worker {0}".format(rank)
        with open(oname.replace(".dat",".{0}.dat".format(rank)),'w') as f:
            if rank == 0:
                # if master process, write the header information first and then the data
                f.write("##File Type: Far field\n##File Format: 3\n##Source: mwa_tiedarray\n##Date: {0}\n".format(time.iso))
                f.write("** File exported by FEKO kernel version 7.0.1-482\n\n")
                f.write("#Request Name: FarField\n#Frequency: {0}\n".format(obsfreq))
                f.write("#Coordinate System: Spherical\n#No. of Theta Samples: {0}\n#No. of Phi Samples: {1}\n".format(ntheta,nphi))
                f.write("#Result Type: Gain\n#No. of Header Lines: 1\n")
                f.write('#\t"Theta"\t"Phi"\t"Re(Etheta)"\t"Im(Etheta)"\t"Re(Ephi)"\t"Im(Ephi)"\t"Gain(Theta)"\t"Gain(Phi)"\t"Gain(Total)"\n')

            for res in results:
                # write each line of the data
                # we actually need to rotate out phi values by: phi = pi/2 - az because that's what FEKO expects.
                    # values are calculated using that convetion, so we need to represent that here
                #f.write("{0}\t{1}\t0\t0\t0\t0\t0\t0\t{2}\n".format(res[0],res[1],res[2][0][0]))
                f.write("{0}\t{1}\t0\t0\t0\t0\t0\t0\t{2}\n".format(res[0],res[1],res[2]))
        print "worker {0} pattern maximum ".format(rank)    
    else:
        print "worker {0} not writing".format(rank)
        
    return omega_A
    
def two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    return values 
    

#####################
##  OPTION PARSING ##
#####################
parser = argparse.ArgumentParser(description="""Script to calculate the array factor required to model the tied-array beam for the MWA. 
                        This is an MPI-based simulation code and will use all available processes when run 
                        (i.e. there is no user choice in how many to use)""",\
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-o","--obsid",type=str,action='store',metavar="obsID",\
        help="""Observation ID (used to figure out spatial configuration of array). 
            Also used to retrieve observation start-time and duration.""",default=None)

parser.add_argument("--flagged_tiles",type=str,nargs='+',action='store',metavar="tile",\
        help="The tiles flagged as in when running the RTS. Must be a list of space separated tile numbers, e.g. 0 1 2 5",default=None)

parser.add_argument("-p","--target",type=two_floats,action='store',\
        help="The RA and DEC of the target pointing (i.e the desired phase-centre). Should be formtted like: hh:mm:ss.sss dd\"mm\'ss.ssss",)
        #default=("00:00:00.0000","00:00:00.0000")
        #metavar=("ra","dec")

parser.add_argument("-t","--time",type=float,action='store',metavar="time",\
        help="""The GPS time desired for beam evaluation. This will override whatever is read from the metafits.""",default=None)

parser.add_argument("-f","--freq",type=float,action='store',metavar="frequency",help="The centre observing frequency for the observation (in Hz!)",default=184.96e6)

parser.add_argument("-e","--efficiency",type=float,action='store',metavar="eta",help="Frequency and pointing dependent array efficiency",default=1)

parser.add_argument("--grid_res",type=float,action='store',nargs=2,metavar=("theta_res","phi_res"),
        help="""Resolution of the Azimuth (Az) and Zenith Angle (ZA) grid to be created in degrees. 
            Be warned: setting these too small will result in a MemoryError and the job will die.""",default=(0.1,0.1))

parser.add_argument("--coplanar",action='store_true',help="Assume the array is co-planar (i.e. height above array centre is 0 for all tiles)")

parser.add_argument("--zenith",action='store_true',help="Assume zenith pointing (i.e  delays are 0), ZA = 0 and AZ = 0")

parser.add_argument("--out_dir",type=str,action='store',help="Location (full path) to write the output data files",default=".")

parser.add_argument("--write",action='store_true',
        help="""Write the beam pattern to disk when done calculating. 
            If this option is not passed, you will just get a '.stats' files containing basic information about simulation parameters and the calculated gain.""")
    
#Nicks extra options        
parser.add_argument("--plot_centre",type=two_floats,action='store',help="Plots around this point (used with -r)")

parser.add_argument("--out_name",type=str,action='store',help="Name of the output file.")
        
parser.add_argument("-r","--radius",type=float,action='store',metavar="radius",\
        help="""The radius of the square to calculate in degrees. Will calculate only around the pointing. Default is none which will search the whole sky.""",default=None)
        
        
parser.add_argument("-s","--strip",action='store_true',\
        help="""Will calculate a for one value of zenith. This is used to get an approxiamte FBHM.""")

# parse the arguments
args = parser.parse_args()

###############
## Setup MPI ##
###############
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# small calculations and data gethering from arguments is fine and won't run into trouble by having multiple processes do it simultaneously
ra,dec = args.target
if dec.startswith("\""):
    dec = dec[1:]
if dec.endswith("\""):
    dec = dec[:-1]
    
tres,pres = args.grid_res
ntheta,nphi = 90/tres,360/pres

if args.flagged_tiles == None:
    flags = []
else:
    flags = args.flagged_tiles

# For jobs involving downloads, reading files and accessing databases, ONLY let the master node gather and then broadcast

# if the observation ID metafits file doesn't exist, get the master node to download it
if rank == 0:
    print "will use {0} processes".format(size)
    print "gathering required data"
    os.system('wget -O {0}/{1}_metafits_ppds.fits mwa-metadata01.pawsey.org.au/metadata/fits?obs_id={1}'.format(args.out_dir,args.obsid))

    # for delays, which requires reading the metafits file, only let master node do it and then broadcast to workers
    if args.zenith:
        delays = [0]*16
    else:
        delays = get_delay_steps(args.obsid)[4]

    # same for obs time, master reads and then distributes
    if args.time is None:
        time = Time(get_obstime_duration(args.obsid)[0],format='gps',fdir=args.out_dir)
    else:
        time = Time(args.time,format='gps')

    # get the tile locations from the metafits
    xpos,ypos,zpos = getTileLocations(args.obsid,flags,fdir=args.out_dir)    
    if args.coplanar:
        zpos = np.zeros(len(xpos))
    
# wait for the master node to gather the data
comm.barrier()

if rank == 0:
    # create the object that contains all the data from master
    data = (delays,time,xpos,ypos,zpos)
elif rank != 0:
    # create the data object, but leave it empty (will be filled by bcast call)
    data = None

# now broadcast the data from the master to the slaves
data = comm.bcast(data,root=0)
if data:
    print "broadcast received by worker {0} successfully".format(rank)
    delays,time,xpos,ypos,zpos = data
else:
    print "broadcast failed to worker {0}".format(rank)
    print "!!! ABORTING !!!"
    comm.Abort(errorcode=1)

# wait for all processes to have recieved the data
comm.barrier()    

# set the base output file name (will be editted based on the worker rank)
# TODO: make this more generic - maybe as an option for what directory?
if args.out_name:
    oname="{0}/{1}.dat".format(args.out_dir, args.out_name)
else:
    oname = "{0}/{1}_{2}_{3}MHz_{4}_{5}.dat".format(args.out_dir,args.obsid,time.gps,args.freq/1e6,ra,dec)

#if no deg they're in radians
if args.zenith:
    pointingAZ,pointingZA,pointingAZdeg,pointingZAdeg = 0,0,0,0
else:
    pointingAZ,pointingZA,pointingAZdeg,pointingZAdeg = getTargetAZZA(ra,dec,time) 

if args.plot_centre:
    plot_ra,plot_dec = args.plot_centre
    plot_AZ,plot_ZA,plot_AZdeg,plot_ZAdeg = getTargetAZZA(plot_ra,plot_dec,time)
    
#Calc za start and end (limits of where to plot)
if args.radius:
    #Find the target azumith and zenith then make a box to calculate inside of it
    
    if args.strip:
        totalZAangle = np.radians(tres)
        ZAstart = pointingZA
        ZAend = pointingZA
    elif args.plot_centre:
        totalZAangle = 2.0*np.radians(args.radius)
        ZAstart = plot_ZA - np.radians(args.radius)
        ZAend = plot_ZA + np.radians(args.radius)
    else:
        totalZAangle = 2.0*np.radians(args.radius)
        ZAstart = pointingZA - np.radians(args.radius)
        ZAend = pointingZA + np.radians(args.radius)
    
    totalAZangle = 2.0*np.radians(args.radius)
    
else:#No radius so calculate for all sky
    totalZAangle = np.pi/2
    ZAstart = 0.0
    ZAend = np.pi/2
    
    totalAZangle = 2*np.pi

totalZAevals = totalZAangle/np.radians(tres) #total number of ZA evaluations required
assert totalZAevals >= size, "Total number of ZA evalutions must be >= the number of processors available"

# figure out how many chunks to split up ZA into
if args.radius:
    if args.plot_centre:
        AZstart = plot_AZ - np.radians(args.radius)
        AZstop = plot_AZ + np.radians(args.radius)
    else:
        AZstart = pointingAZ - np.radians(args.radius)
        AZstop = pointingAZ + np.radians(args.radius)
    
else:#No radius so calculate for all sky
    AZstart = 0.0
    AZstop = 2*np.pi

# iterate through the ZA range given the process rank
ZA_rank_start = ZAstart + rank * np.radians(tres) * (totalZAevals//size)
ZA_rank_end = ZAstart + (rank+1) * np.radians(tres) * (totalZAevals//size)
if rank == (size-1):
    ZA_rank_end = ZAend # give the last process anything that's left

# calculate how many calculations this worker will do
numZA_calcs = (ZA_rank_end-ZA_rank_start)/np.radians(tres)
numAZ_calcs = totalAZangle/np.radians(pres)
numWorker_calcs = numZA_calcs * numAZ_calcs
print "worker:",rank,"total calcs:",numWorker_calcs,"start ZA:",np.degrees(ZA_rank_start),"end ZA",np.degrees(ZA_rank_end)

# create array factor for given ZA band and write to file
beam_area = createArrayFactor(pointingAZ,pointingZA,args.obsid,delays,time,args.freq,args.efficiency,flags,
        tres,pres,args.coplanar,ZA_rank_start,ZA_rank_end,AZstart,AZstop,
        args.write,args.radius,args.strip)
        
#beam_area = createArrayFactor(ra,dec,args.obsid,delays,time,args.freq,args.efficiency,flags,tres,pres,args.coplanar,args.zenith,start,end,args.write)

# collect results for the beam area calculation
if rank != 0:
    comm.send(beam_area,dest=0)
elif rank == 0:
    result = beam_area
    for i in range(1,size):
        result += comm.recv(source=i)

# wait for everything to be collected (not really sure if this is necessary...)
comm.barrier()

# calculate the gain for that pointing and frequency and write a "stats" file
if rank == 0:
    eff_area = args.efficiency*((c.value/args.freq)**2/result)
    gain = (1e-26)*eff_area/(2*k_B.value)
    
    if args.write:
        nfiles = size
    else:
        nfiles = 0

    # write the stats file
    with open(oname.replace(".dat",".stats"),"w") as f:
        f.write("#observation ID\n{0}\n".format(args.obsid))
        f.write("#pointing (ra,dec)\n{0} {1}\n".format(ra,dec))
        f.write("#time (gps seconds)\n{0}\n".format(time.gps))
        f.write("#frequency (MHz)\n{0}\n".format(args.freq/1e6))
        f.write("#telescope efficiency\n{0}\n".format(args.efficiency))
        f.write("#flagged tiles\n{0}\n".format(flags))
        f.write("#beam solid angle (sr)\n{0}\n".format(result))
        f.write("#effective area (m^2)\n{0}\n".format(eff_area))
        f.write("#gain (K/Jy)\n{0}\n".format(gain))
        f.write("#ZA and Az resolution (degrees per pixel)\n({0} , {1})\n".format(tres,pres))
        f.write("#number of data files written\n{0}".format(nfiles))
    
    print "Gain for {0} at {1} MHz pointed at {2} {3} is: {4:.6f} K/Jy".format(args.obsid,args.freq,ra,dec,gain)
