import tkp.db
from tkp.db.model import Varmetric
from tkp.db.model import Runningcatalog
from tkp.db.model import RunningcatalogFlux
from tkp.db.model import Newsource
from tkp.db.model import Extractedsource
from tkp.db.model import Image
from tkp.db.model import Assocxtrsource
from tkp.db.model import Skyregion
from tkp.db.model import Frequencyband
from sqlalchemy import *
from sqlalchemy.orm import relationship
import pandas as pd


def access(engine,host,port,user,password,database):
    """ Access the database using sqlalchemy"""
    # make db global in order to be used in GetPandaExtracted
    global db
    db = tkp.db.Database(engine=engine, host=host, port=port,
                     user=user, password=password, database=database)
    db.connect()
    session = db.Session()
    print 'connected!'
    return session

def GetExtractedSources(session, dataset_id):
    extracted_sources_query = session.query(Extractedsource,
                                            Image.taustart_ts,
                                            Assocxtrsource.runcat_id, 
                                            Assocxtrsource.distance_arcsec, 
                                            Assocxtrsource.v_int, 
                                            Assocxtrsource.eta_int).select_from(join(Extractedsource,Image).join(Assocxtrsource)).filter(Image.dataset_id == dataset_id)
    extracted_sources = pd.read_sql(extracted_sources_query.statement, session.bind).sort_values(by="id")
    return extracted_sources

def GetRunningCatalogs(session, dataset_id):
    running_catalogs_query = session.query(Runningcatalog.id,
                                           Runningcatalog.dataset_id,
                                           Runningcatalog.wm_ra,
                                           Runningcatalog.wm_decl,
                                           Runningcatalog.wm_uncertainty_ew,
                                           Runningcatalog.wm_uncertainty_ns,
                                           Runningcatalog.avg_ra_err,
                                           Runningcatalog.avg_decl_err,
                                           Runningcatalog.datapoints,
                                           Runningcatalog.forcedfits_count,
                                           RunningcatalogFlux.band_id,
                                           RunningcatalogFlux.avg_f_peak,
                                           RunningcatalogFlux.avg_f_int,
                                           Runningcatalog.mon_src,
                                           Varmetric.sigma_rms_max,
                                           Varmetric.sigma_rms_min,
                                           Varmetric.lightcurve_max,
                                           Varmetric.lightcurve_avg,
                                           Varmetric.lightcurve_median,
                                           Varmetric.v_int,
                                           Varmetric.eta_int,
                                           Frequencyband.freq_central
                                          ).select_from(join(Runningcatalog,RunningcatalogFlux).join(Varmetric).join(Frequencyband)).filter(Runningcatalog.dataset_id == dataset_id)
    running_catalogs = pd.read_sql(running_catalogs_query.statement, session.bind).sort_values(by="id")
    return running_catalogs

def GetImages(session, dataset_id):
    images_query = session.query(Image, Skyregion.centre_ra, Skyregion.centre_decl).select_from(join(Image,Skyregion)).filter(Image.dataset_id == dataset_id)
    images = pd.read_sql(images_query.statement, session.bind).sort_values(by="id")
    return images

def GetNewsources(session, dataset_id):
    newsources_query = session.query(Newsource).select_from(join(Newsource,Runningcatalog)).filter(Runningcatalog.dataset_id == dataset_id)
    newsources = pd.read_sql(newsources_query.statement, session.bind).sort_values(by="id")
    return newsources