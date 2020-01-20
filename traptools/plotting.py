import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from scipy.stats import norm
from astroML import density_estimation
import numpy as np

from astropy.wcs import WCS
from astropy.io import fits
from astropy.visualization import ZScaleInterval,ImageNormalize,PercentileInterval
from astropy.visualization import LinearStretch
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

def make_cmap(frequencies):
    cm = matplotlib.cm.get_cmap('hsv')
    col = [cm(1.*i/len(frequencies)) for i in range(len(frequencies))]
    return col

def gaussian_fit(data,param):
    range_data=np.linspace(min(data),max(data),1000)
    fit=norm.pdf(range_data,loc=param[0],scale=param[1])
    return range_data,fit

def make_bins(x):
    new_bins = density_estimation.bayesian_blocks(x)
    binsx = [new_bins[a] for a in range(len(new_bins)-1) if abs((new_bins[a+1]-new_bins[a])/new_bins[a])>0.05]
    binsx = binsx + [new_bins[-1]]
    return binsx

def diagnostic_plot(df, sigcutx, sigcuty, peak=False):
    
    if peak:
        v_int = "v_int_peak"
        eta_int = "eta_int_peak"
    else:
        v_int = "v_int"
        eta_int = "eta_int"
    
    frequencies = df.freq_central.unique()
    freq_labels=["{:.0f} MHz".format(f/1.e6) for f in frequencies]
    nullfmt   = NullFormatter()         # no labels
    
    fig = plt.figure(1,figsize=(12,12))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    fontP = FontProperties()
    fontP.set_size('large')
    fig.subplots_adjust(hspace = .001, wspace = 0.001)
    ax1.set_ylabel(r'$\eta_\nu$', fontsize=28)
    ax3.set_ylabel(r'$V_\nu$', fontsize=28)
    ax3.set_xlabel('Max Flux (Jy)', fontsize=24)
    ax4.set_xlabel('Max Flux / Median Flux', fontsize=24)

    for i in range(len(frequencies)):
        dfTMP=df.loc[(df['freq_central']==frequencies[i])]
        xdata_ax3=dfTMP['lightcurve_max']
        xdata_ax4=dfTMP['lightcurve_max']/dfTMP['lightcurve_avg']
        ydata_ax1=dfTMP[eta_int]
        ydata_ax3=dfTMP[v_int]
        ax1.scatter(xdata_ax3, ydata_ax1, s=10., zorder=5)
        ax2.scatter(xdata_ax4, ydata_ax1, s=10., zorder=6)
        ax3.scatter(xdata_ax3, ydata_ax3, s=10., zorder=7)
        ax4.scatter(xdata_ax4, ydata_ax3, s=10., zorder=8)
        ax4.legend(freq_labels, loc=4, prop=fontP)

    Xax3=df['lightcurve_max']
    Xax4=df['lightcurve_max']/dfTMP['lightcurve_avg']
    Yax1=df[eta_int]
    Yax3=df[v_int]
    
    if sigcutx != 0 or sigcuty != 0:
        ax1.axhline(y=10.**sigcutx, linewidth=2, color='k', linestyle='--')
        ax2.axhline(y=10.**sigcutx, linewidth=2, color='k', linestyle='--')
        ax3.axhline(y=10.**sigcuty, linewidth=2, color='k', linestyle='--')
        ax4.axhline(y=10.**sigcuty, linewidth=2, color='k', linestyle='--')

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax4.set_yscale('log')
    xmin_ax3=10.**(int(np.log10(min(Xax3))-1.1))
    xmax_ax3=10.**(int(np.log10(max(Xax3))+1.2))
    xmin_ax4=0.8
    xmax_ax4=int(max(xdata_ax4)+0.5)
    ymin_ax1=10.**(int(np.log10(min(Yax1))-1.1))
    ymax_ax1=10.**(int(np.log10(max(Yax1))+1.2))
    ymin_ax3=10.**(int(np.log10(min(Yax3))-1.1))
    ymax_ax3=10.**(int(np.log10(max(Yax3))+1.2))
    ax1.set_ylim(ymin_ax1,ymax_ax1)
    ax3.set_ylim(ymin_ax3,ymax_ax3)
    ax3.set_xlim(xmin_ax3,xmax_ax3)
    ax4.set_xlim(xmin_ax4,xmax_ax4)
    ax1.set_xlim( ax3.get_xlim() )
    ax4.set_ylim( ax3.get_ylim() )
    ax2.set_xlim( ax4.get_xlim() )
    ax2.set_ylim( ax1.get_ylim() )
    ax1.xaxis.set_major_formatter(nullfmt)
    ax4.yaxis.set_major_formatter(nullfmt)
    ax2.xaxis.set_major_formatter(nullfmt)
    ax2.yaxis.set_major_formatter(nullfmt)
    
    plt.show()
    plt.close()
    
def eta_vs_v_plot(df, sigcutx, sigcuty, paramx, paramy, peak=False):

    if peak:
        v_int = "v_int_peak"
        eta_int = "eta_int_peak"
    else:
        v_int = "v_int"
        eta_int = "eta_int"
        
    frequencies = df.freq_central.unique()

    nullfmt   = NullFormatter()         # no labels
    fontP = FontProperties()
    fontP.set_size('large')
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    fig = plt.figure(figsize=(12,12))
    axScatter = fig.add_subplot(223, position=rect_scatter)
    plt.xlabel(r'$\eta_{\nu}$', fontsize=28)
    plt.ylabel(r'$V_{\nu}$', fontsize=28)
    axHistx=fig.add_subplot(221, position=rect_histx)
    axHisty=fig.add_subplot(224, position=rect_histy)
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    axHistx.axes.yaxis.set_ticklabels([])
    axHisty.axes.xaxis.set_ticklabels([])

    for i in range(len(frequencies)):
        dfTMP=df.loc[(df['freq_central']==frequencies[i])]
        xdata_var=np.log10(dfTMP[eta_int])
        ydata_var=np.log10(dfTMP[v_int])
        axScatter.scatter(xdata_var, ydata_var, s=10., zorder=5, color='C1')

    x = np.log10(df[eta_int])
    y = np.log10(df[v_int])

    axHistx.hist(x, bins=make_bins(x), normed=1, histtype='stepfilled', color='C0')
    axHisty.hist(y, bins=make_bins(y), normed=1, histtype='stepfilled', orientation='horizontal', color='C0')

    freq_labels=["{:.0f} MHz".format(f/1.e6) for f in frequencies]
    axScatter.legend(freq_labels,loc=4, prop=fontP)
    xmin=int(min(x)-1.1)
    xmax=int(max(x)+1.1)
    ymin=int(min(y)-1.1)
    ymax=int(max(y)+1.1)
    xvals=range(xmin,xmax)
    xtxts=[r'$10^{'+str(a)+'}$' for a in xvals]
    yvals=range(ymin,ymax)
    ytxts=[r'$10^{'+str(a)+'}$' for a in yvals]
    axScatter.set_xlim([xmin,xmax])
    axScatter.set_ylim([ymin,ymax])
    axScatter.set_xticks(xvals)
    axScatter.set_xticklabels(xtxts, fontsize=20)
    axScatter.set_yticks(yvals)
    axScatter.set_yticklabels(ytxts, fontsize=20)
    axHistx.set_xlim( axScatter.get_xlim())
    axHisty.set_ylim( axScatter.get_ylim())

    if sigcutx != 0 or sigcuty != 0:
        axHistx.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')
        axHisty.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')
        axScatter.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')
        axScatter.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')

    range_x,fitx = gaussian_fit(x,paramx)
    axHistx.plot(range_x,fitx, 'k:', linewidth=2)
    range_y,fity = gaussian_fit(y,paramy)
    axHisty.plot(fity,range_y, 'k:', linewidth=2)
    plt.show()
    plt.close()
    
def plot_lightcurve(xtr_sources, title="Lightcurve", peak_flux=False, save=False, plot_limits=False,
                    plot_limits_dates=[], plot_limits_values=[], plot_limits_sigma=5.):
    fig = plt.figure(figsize=(20,5))
    ax = fig.add_subplot(111)
    if peak_flux:
        flux_to_use = "f_peak"
        flux_err_to_use = "f_peak_err"
        flux_label = "Peak Flux (mJy/beam)"
    else:
        flux_to_use = "f_int"
        flux_err_to_use = "f_int_err"
        flux_label = "Int. Flux (mJy)"
    ax.plot(xtr_sources["taustart_ts"], xtr_sources[flux_to_use]*1.e3, marker="None")
    blind = xtr_sources[xtr_sources["extract_type"]==0]
    ax.errorbar(blind["taustart_ts"], blind[flux_to_use]*1.e3, yerr=blind[flux_err_to_use]*1.e3, marker="s", color="C0", linestyle="None")
    forced = xtr_sources[xtr_sources["extract_type"]==1]
    ax.errorbar(forced["taustart_ts"], forced[flux_to_use]*1.e3, yerr=forced[flux_err_to_use]*1.e3, marker="v", color="C0", linestyle="None")
    mon = xtr_sources[xtr_sources["extract_type"]==2]
    ax.errorbar(mon["taustart_ts"], mon[flux_to_use]*1.e3, yerr=mon[flux_err_to_use]*1.e3, marker="d", color="C0", linestyle="None")
    if plot_limits:
        ax.plot(plot_limits_dates, plot_limits_values*1.e3*plot_limits_sigma, marker="v", color="C1")
    ax.set_title(title)
    ax.set_xlabel("Date")
    ax.set_ylabel(flux_label)
    ax.grid(True)
    plt.show()
    if save:
        fig.savefig("{}_lc.png".format(title), bbox_inches="tight")
    
def create_img(target, extracted_sources, img_data, img_wcs, images, skyrgn, img_skyrgns, imgsize = 2., title="ASKAP", skip_first=False, max_cols=4, 
    all_images=False, percentile=99.9, zscale=False, zscale_contrast=0.2, save=False, scatter_x=[], scatter_y=[], flat_norm=True):
    num_subplots = len(extracted_sources.index)
    imgsize = imgsize * u.arcmin
    extracted_sources = extracted_sources.sort_values(by="taustart_ts")
    if skip_first and len(extracted_sources.index) > 1:
        num_subplots -= 1
    rows = num_subplots/max_cols
    if num_subplots%max_cols != 0:
        rows+=1
    key=1
    norm_done = False
    if rows/max_cols >= 2:
        fig = plt.figure(figsize=(15,15*int(rows/max_cols)))
    else:
        fig = plt.figure(figsize=(20,10))
    all_images_ids = []
    for i in sorted(img_data.keys()):
        if img_skyrgns[i] not in skyrgn:
            continue
        else:
            all_images_ids.append(i)
    if all_images:
        for i,val in enumerate(all_images_ids):
            
            if i==0 and skip_first:
                
                continue
            cutout = Cutout2D(img_data[val], target, imgsize, wcs=img_wcs[val])
            if not norm_done:
                if zscale:
                    norm = ImageNormalize(cutout.data, interval=ZScaleInterval(contrast=zscale_contrast))
                #Use the same scaling throughout (hope that the first image is good!)
                else:
                    norm = ImageNormalize(cutout.data, interval=PercentileInterval(percentile), stretch=LinearStretch())
                if flat_norm == True:
                    norm_done = True
            ax = fig.add_subplot(rows,max_cols,key, projection=cutout.wcs)
            ax.imshow(cutout.data, norm=norm, cmap="gray_r")
            ax.scatter([target.ra.deg], [target.dec.deg], transform=ax.get_transform('world'), marker="o", facecolors="none", edgecolors="r", zorder=10, s=300)
            if len(scatter_x) > 0:
                plt.autoscale(False)
                ax.scatter(scatter_x, scatter_y, transform=ax.get_transform('world'), marker="x", color="C0", zorder=10, s=300)
            thistitle = "{} {}".format(title, images[images["id"]==val].iloc[0]["taustart_ts"].strftime("%Y-%m-%d %H:%M:%S"))
            ax.set_title(thistitle)
            lon = ax.coords[0]
            lat = ax.coords[1]
            lon.set_ticks_visible(False)
            lon.set_ticklabel_visible(False)
            lat.set_ticks_visible(False)
            lat.set_ticklabel_visible(False)
            # lon.set_major_formatter('hh:mm:ss.s')
            key+=1
    else:
        for i,row in extracted_sources.iterrows():
            if i==0 and skip_first:
                continue
            cutout = Cutout2D(img_data[row.image], target, imgsize, wcs=img_wcs[row.image])
            # norm = ImageNormalize(cutout.data, interval=ZScaleInterval(contrast=0.2))
            if not norm_done:
                #Use the same scaling throughout (hope that the first image is good!)
                if zscale:
                    norm = ImageNormalize(cutout.data, interval=ZScaleInterval(contrast=zscale_contrast))
                else:
                    norm = ImageNormalize(cutout.data, interval=PercentileInterval(percentile), stretch=LinearStretch())
                if flat_norm == True:
                    norm_done = True
            ax = fig.add_subplot(rows,max_cols,key, projection=cutout.wcs)
            ax.imshow(cutout.data, norm=norm, cmap="gray_r")
            ax.scatter([target.ra.deg], [target.dec.deg], transform=ax.get_transform('world'), marker="o", facecolors="none", edgecolors="r", zorder=10, s=300)
            if len(scatter_x) > 0:
                plt.autoscale(False)
                ax.scatter(scatter_x, scatter_y, transform=ax.get_transform('world'), marker="x", color="C0", zorder=10, s=300)
            thistitle = "{} {}".format(title, row.taustart_ts.strftime("%Y-%m-%d %H:%M:%S"))
            ax.set_title(thistitle)
            lon = ax.coords[0]
            lat = ax.coords[1]
            lon.set_ticks_visible(False)
            lon.set_ticklabel_visible(False)
            lat.set_ticks_visible(False)
            lat.set_ticklabel_visible(False)
            # lon.set_major_formatter('hh:mm:ss.s')
            key+=1
    plt.show()
    if save:
        fig.savefig("{}_stamps.png".format(title), bbox_inches="tight")