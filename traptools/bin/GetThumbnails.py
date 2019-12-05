#!/usr/bin/env python

# GetThumbnails.py
#
# A code to genreate thumbnails of all extracted sources for a given dataset.
#
# Import all the dependencies and generic setup
import numpy as np
import logging
import sys
# from plotting import plot_varib_params as pltvp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
import pandas as pd
import os
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from astropy.visualization import ZScaleInterval,ImageNormalize

from traptools.dbtools import access
from traptools.dbtools import GetExtractedSources
from traptools.dbtools import GetRunningCatalogs
from traptools.dbtools import GetImages
from traptools.dbtools import GetNewsources

import argparse
import getpass

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

def filter_selavy_components(catalog, catalog_coords, src_coord, imsize):
    #Filter out selavy components outside field of image
    seps = src_coord.separation(catalog_coords)
    mask = seps <= imsize/1.4 #I think cutout2d angle means the width of the image, not a radius hence /2
    #drop the ones we don't need
    return catalog[mask].reset_index(drop=True)

def create_thumbnails(source_ids, size, full_wcs, img_data, norm, rms, coords, ellipses, outdir, legend_elements):
    # src_coord = SkyCoord(ra, decl, unit=(u.deg, u.deg))
    # inverval = ZScaleInterval()
    # vmin, vmax = inverval.get_limits(img_data)
    # vmin = -2.5*rms
    # vmax = 4*rms
    fig = plt.figure(figsize=(3.0,3.0))

    for source_id in source_ids:
        # print "Creating cutout for ID:{}".format(source_id)
        cutout = Cutout2D(img_data, position=coords[source_id], size=size, wcs=full_wcs)
        ax = fig.add_subplot(111,projection=cutout.wcs)
        # ax.set_autoscale_on(False)
        im_data = cutout.data
        # std_dev = np.std(im_data)
        # sigma_clipped = im_data[im_data < 3*std_dev]
  
        # rms = np.std(sigma_clipped)
        # rms =
  

        
        im = ax.imshow(im_data,norm=norm,cmap='gray_r')
        # ax.add_collection(ellipses, transform=ax.transAxes)
    # ax.set_xlabel('Right Ascenscion')
    # ax.set_ylabel('Declination')
        collection = PatchCollection(ellipses[source_id][0], facecolors="None", edgecolors=ellipses[source_id][1], linestyle="-", lw=1, transform=ax.get_transform('world'))
        ax.add_collection(collection, autolim=False)
        # ellipse = Ellipse((,ellipses[source_id][1]),width=ellipses[source_id][3], height=ellipses[source_id][2], angle=ellipses[source_id][4],
        #      transform=ax.get_transform('world'), fill=False,color='C{}'.format(ellipses[source_id][5]),linewidth=0.5,
        #         linestyle="--", alpha=0.8, label="Extract Type: {}".format(ellipses[source_id][5]))
        #
        # ax.add_patch(ellipse)
        ax.legend(handles=legend_elements, prop={'size': 6})
  
        # ax.invert_yaxis()
        
        lon = ax.coords[0]
        lat = ax.coords[1]
        
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)
        
        # bar = show_scale(ax)
        # ax.add_artist(bar)
        
        # ax.colorbar(im)
        # fig.colorbar(im)
        
        fig.savefig(os.path.join(outdir,"{}.png".format(source_id)), dpi=300, bbox_inches="tight")
        fig.clf()
    #plt.show()
    plt.close(fig)
    return

parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('database', type=str, help='Database name.')
parser.add_argument('dataset_id', type=int, help='Dataset id.')
parser.add_argument('--angle-size', type=float, help='Angle size of thumbnail (in arcmin).', default=3.0)
parser.add_argument('--monitor-only', action="store_true", help='Only generate stamps for monitor sources.')
parser.add_argument('--new-source-only', action="store_true", help='Only generate stamps for new sources.')
parser.add_argument('--user-sources', action="store_true", help='Only generate stamps for new sources.')
parser.add_argument('--user-sources-file', type=str, help='File containing runcat ids to query.', default="sources.csv")
parser.add_argument('--db-engine', type=str, help='Database engine', default="postgresql")
parser.add_argument('--db-user', type=str, help='Database user', default=getpass.getuser())
parser.add_argument('--db-password', type=str, help='Database password', default=getpass.getuser())
parser.add_argument('--db-host', type=str, help='Database host', default="localhost")
parser.add_argument('--db-port', type=int, help='Database port', default=5432)
parser.add_argument('--nworkers', type=int, help='Number of workers', default=1)
parser.add_argument('--selavy-angle-corr', action="store_true", help='Apply the correction to the selavy position angle.')

args=parser.parse_args()

# The input database, dataset and thresholds
dataset_id = args.dataset_id
# database = 'dev_testing'
database = args.database
engine = args.db_engine
user = args.db_user
password = args.db_password
host = args.db_host
port = args.db_port
monitor_only = args.monitor_only
newsource_only = args.new_source_only
#run cat ids
user_choice = args.user_sources
user_choice_file = args.user_sources_file

n_workers = args.nworkers
angle_size = args.angle_size

theangle = Angle(args.angle_size, unit=u.arcmin)

selavy_corr = args.selavy_angle_corr

if monitor_only and newsource_only:
    sys.exit()
# sigma1 = 2 # Threshold on the reduced weighted chi^2
# sigma2 = 2 # Threshold on the variability parameter
# websiteURL = 'http://banana.transientskp.org/r4/vlo_'+database+'/runningcatalog/'

outdir = os.path.join("thumbnails",engine.replace("postgresql", "postgres")+"_"+database,str(dataset_id))

if os.path.isdir(outdir):
    print "Directory already exists!"
    sys.exit()
else:
    os.makedirs(outdir)

# Connect to the database and run the queries
session = access(engine,host,port,user,password,database)

images = GetImages(session,dataset_id)
extracted_sources = GetExtractedSources(session,dataset_id)
session.close()
if newsource_only:
    newsources = GetNewsources(session, dataset_id)
    newsource_xtrsrc_id = newsources["trigger_xtrsrc"].tolist()
# now create thumbnails of every single source for each image
for i,row in images.iterrows():
    this_url = row["url"]
    this_rms = row["rms_qc"]
    print "Image: {}".format(this_url.split("/")[-1])
    full_hdu, full_wcs, full_header, img_data = load_image(this_url)
    norm = ImageNormalize(img_data, interval=ZScaleInterval(contrast=0.15))
    all_sources = extracted_sources[extracted_sources["image"]==row["id"]].reset_index(drop=True)
    all_sources_coords = SkyCoord(all_sources["ra"]*u.deg, all_sources["decl"]*u.deg)
    if monitor_only:
        sources = all_sources[all_sources["extract_type"]==2].reset_index(drop=True)
    elif user_choice:
        user_file = pd.read_csv(user_choice_file)
        runcat_ids = user_file["runcat"].tolist()
        sources = all_sources[all_sources["runcat"].isin(runcat_ids)].reset_index(drop=True)
    elif newsource_only:
        sources = all_sources[all_sources["runcat"].isin(newsources.runcat)].reset_index(drop=True)
        if len(sources.index)==0:
            print "Skipping {} no new sources.".format(this_url.split("/")[-1])
            continue
    else:
        sources = all_sources
    # print "Building ellipse collection..."
    # ww = 2.*all_sources["semimajor"].astype(float)/3600.
    # hh = 2.*all_sources["semiminor"].astype(float)/3600.
    # aa = all_sources["pa"].astype(float)
    # x = all_sources["ra"].astype(float)
    # y = all_sources["decl"].astype(float)
    # patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, 90.+(180.-aa[i])) for i in range(len(x))]
    # ellipses = [patches, ["C{}".format(i) for i in sources.extract_type.tolist()]]
    # print "Done!"

    # sources = sources.iloc[0:50]
    sources_dict = {row["id"]:SkyCoord(row["ra"], row["decl"], unit=(u.deg, u.deg)) for i, row in sources.iterrows()}
    ellipse_dict = {}
    num_sources = len(sources_dict)
    
    if num_sources > 0:
        if num_sources >= n_workers:
            actual_workers = n_workers
            workers = multiprocessing.Pool(processes=n_workers) 
        else:
            actual_workers = num_sources
            workers = multiprocessing.Pool(processes=num_sources) 
        print "Building ellipse collection..."
        for i, row in sources.iterrows():
            target = SkyCoord(row.ra*u.deg, row.decl*u.deg)
            sources_trimmed = filter_selavy_components(all_sources, all_sources_coords, target, theangle)
            ww = 2.*sources_trimmed["semimajor"].astype(float)/3600.
            hh = 2.*sources_trimmed["semiminor"].astype(float)/3600.
            aa = sources_trimmed["pa"].astype(float)
            x = sources_trimmed["ra"].astype(float)
            y = sources_trimmed["decl"].astype(float)
            if selavy_corr:
                patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, 90.+(180.-aa[i])) for i in range(len(x))]
            else:
                patches = [Ellipse((x[i], y[i]), ww[i]*1.1, hh[i]*1.1, aa[i]) for i in range(len(x))]
            ellipse_dict[row["id"]]=[patches, ["C{}".format(i) for i in sources_trimmed.extract_type.tolist()]]
        print "Done!"
        # ellipses=[Ellipse((row["ra"],row["decl"]),width=2.*row["semimajor"]/3600., height=2.*row["semiminor"]/3600.,
            # angle=row["pa"], fill=False,color='C0',linewidth=1.5) for i, row in sources.iloc[0:100].iterrows()]
        # e = PatchCollection(ellipses)
        master_loops = {i:[] for i in range(n_workers)}
        i=0
        for j,val in enumerate(sources_dict):
            master_loops[i].append(val)
            if i==max(master_loops.keys()):
                i=0
            else:
                i+=1
        legend_elements = [Line2D([0], [0], markeredgecolor='C0', markerfacecolor="none", marker="o", label='Blind', markersize=5, linestyle="none"),
                           Line2D([0], [0], markeredgecolor='C1', markerfacecolor="none", marker="o", label='Forced', markersize=5, linestyle="none"),
                           Line2D([0], [0], markeredgecolor='C2', markerfacecolor="none", marker="o", label='Monitor', markersize=5, linestyle="none")
                              ]
        to_loop = [master_loops[i] for i in master_loops]
        create_thumbnails_multi = partial(create_thumbnails, size=theangle, full_wcs = full_wcs, img_data=img_data, norm=norm, rms=this_rms, 
            coords=sources_dict, ellipses=ellipse_dict, outdir=outdir, legend_elements=legend_elements)

        print "Launching thumbnail production with {} processes. {} to do, be patient!".format(actual_workers, num_sources)
        workers.map(create_thumbnails_multi, to_loop)
        workers.close()
    else:
        print "No sources to do, going to next image."
        continue
    