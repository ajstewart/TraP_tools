#!/usr/bin/env python

# GetThumbnails.py
#
# A code to genreate thumbnails of all extracted sources for a given dataset.
#
# Import all the dependencies and generic setup
import scipy as sp
import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship
import tkp.db
import logging
logging.basicConfig(level=logging.INFO)
query_loglevel = logging.WARNING  # Set to INFO to see queries, otherwise WARNING
import sys
sys.path.append('../')
from dblogin import * # This file contains all the variables required to connect to the database
from databaseTools import dbtools
from tools import tools
# from plotting import plot_varib_params as pltvp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import multiprocessing
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
from astropy.wcs import WCS
from functools import partial
import os
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from astropy.visualization import ZScaleInterval,ImageNormalize

import warnings

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)

os.nice(5)

def load_image(imgpath):
    with fits.open(imgpath) as thefits:
        hdu = thefits[0]
        wcs = WCS(hdu.header, naxis=2)
        full_header = hdu.header
        try:
            img_data = hdu.data[0,0,:,:]
        except:
            img_data = hdu.data[:,:]
    
    return hdu, wcs, full_header, img_data

def show_scale(ax, scale_length = 1.*u.arcmin, scale_text = None, scale_loc = 'lower right'):
   ''' Adds a scalebar to a given 2D plot. 
   
   :param ax: The axes to which to add the sale.
   :param scale_length: The scale length in astropy.units
   :param scale_text: the string to print as the scalebar, or None
   :param scale_loc: string [default: 'top left']
                 The location of the scale bar. e.g 'top right', 'bottom left', etc ...
                 
   '''
   
   # raise Exception('Ouch! The scalebar is not working yet ...!')
   
   # Make a default scale text in case none is specified
   if scale_text is None:
      
      the_unit = scale_length.unit.to_string()
      
      if the_unit == 'arcsec':
         the_unit = r'$^{\prime\prime}$'
      elif the_unit == 'arcmin':
         the_unit = r'$^{\prime}$'
      else:
         the_unit = ' ' + the_unit
                
      scale_text = '%.1f%s' % (scale_length.value, the_unit)
   
   # TODO:
   # Looks like I still need to divide the length by the "pixel scale" ? 
   # Unless I need to divide by the Dec ?
   # But how do I extract this from the ax ???
      
   # Create the scalebar
   bar = AnchoredSizeBar(ax.get_transform('world'), scale_length.to(u.degree).value, 
                         scale_text, scale_loc, 
                         sep = 5, pad = 0.5, borderpad = 0.3)
    
   # Adjust the bar thickness
   bar.size_bar.get_children()[0].set_linewidth(2)
   # Add it to the plot
   return bar

def create_thumbnails(source_ids, size, full_wcs, img_data, rms, coords, ellipses, outdir):
    # src_coord = SkyCoord(ra, decl, unit=(u.deg, u.deg))
    # inverval = ZScaleInterval()
    # vmin, vmax = inverval.get_limits(img_data)
    norm = ImageNormalize(img_data, interval=ZScaleInterval(contrast=0.2))
    # vmin = -2.5*rms
    # vmax = 4*rms
    fig = plt.figure(figsize=(3.0,3.0))
  
    for source_id in source_ids:
        # print "Creating cutout for ID:{}".format(source_id)
        cutout = Cutout2D(img_data, position=coords[source_id], size=size, wcs=full_wcs)
        ax = fig.add_subplot(111,projection=cutout.wcs)
        im_data = cutout.data
        # std_dev = np.std(im_data)
        # sigma_clipped = im_data[im_data < 3*std_dev]
  
        # rms = np.std(sigma_clipped)
        # rms =
  

        
        im = ax.imshow(im_data,norm=norm,cmap='gray_r')
        # ax.add_collection(ellipses, transform=ax.transAxes)
    # ax.set_xlabel('Right Ascenscion')
    # ax.set_ylabel('Declination')
        
        ellipse = Ellipse((ellipses[source_id][0],ellipses[source_id][1]),width=ellipses[source_id][3], height=ellipses[source_id][2], angle=ellipses[source_id][4],
             transform=ax.get_transform('world'), fill=False,color='C{}'.format(ellipses[source_id][5]),linewidth=0.5, 
                linestyle="--", alpha=0.8, label="Extract Type: {}".format(ellipses[source_id][5]))
  
        ax.add_patch(ellipse)
        ax.legend(prop={'size': 6})
  
        ax.invert_yaxis()
        
        lon = ax.coords[0]
        lat = ax.coords[1]
        
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)
        
        bar = show_scale(ax)
        ax.add_artist(bar)
        
        # ax.colorbar(im)
        # fig.colorbar(im)
        
        fig.savefig(os.path.join(outdir,"{}.png".format(source_id)), dpi=300, bbox_inches="tight")
        fig.clf()
    #plt.show()
    plt.close(fig)
    return



# The input database, dataset and thresholds
dataset_id = 4
database = 's190814bv'
# sigma1 = 2 # Threshold on the reduced weighted chi^2
# sigma2 = 2 # Threshold on the variability parameter
# websiteURL = 'http://banana.transientskp.org/r4/vlo_'+database+'/runningcatalog/'

outdir = os.path.join("thumbnails",engine+"_"+database,str(dataset_id))

if os.path.isdir(outdir):
    print "Directory already exists!"
    sys.exit()
else:
    os.makedirs(outdir)

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
# VarParams = dbtools.GetVarParams(session,dataset_id)
images = dbtools.GetImages(session,dataset_id)
extracted_sources = dbtools.GetExtractedSources(session,dataset_id)
session.close()
n_workers = 18
workers = multiprocessing.Pool(processes=n_workers)
# now create thumbnails of every single source for each image
for i,row in images.iterrows():
    this_url = row["url"]
    this_rms = row["rms_qc"]
    print "Image: {}".format(this_url.split("/")[-1])
    full_hdu, full_wcs, full_header, img_data = load_image(this_url)
  
    sources = extracted_sources[extracted_sources["image"]==row["id"]].reset_index(drop=True)
    sources_dict = {row["id"]:SkyCoord(row["ra"], row["decl"], unit=(u.deg, u.deg)) for i, row in sources.iterrows()}
    ellipse_dict = {row["id"]:[row["ra"],row["decl"],2.*row["semimajor"]/3600.,2.*row["semiminor"]/3600., 
        row["pa"], row["extract_type"]] for i, row in sources.iterrows()}
    # ellipses=[Ellipse((row["ra"],row["decl"]),width=2.*row["semimajor"]/3600., height=2.*row["semiminor"]/3600.,
        # angle=row["pa"], fill=False,color='C0',linewidth=1.5) for i, row in sources.iloc[0:100].iterrows()]
    # e = PatchCollection(ellipses)
    num_sources = len(sources_dict)
    master_loops = {i:[] for i in range(n_workers)}
    i=0
    for j,val in enumerate(sources_dict):
        master_loops[i].append(val)
        if i==max(master_loops.keys()):
            i=0
        else:
            i+=1
    
    to_loop = [master_loops[i] for i in master_loops]
        
    create_thumbnails_multi = partial(create_thumbnails, size=Angle(5, unit=u.arcmin), full_wcs = full_wcs, img_data=img_data, rms=this_rms, 
        coords=sources_dict, ellipses=ellipse_dict, outdir=outdir)
    print "Launching thumbnail production with {} processes. {} to do, be patient!".format(n_workers, len(to_loop))
    workers.map(create_thumbnails_multi, to_loop, n_workers)
    
workers.close()
    