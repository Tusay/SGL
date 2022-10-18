import sys
import turbo_seti
import glob
import os
import matplotlib
import matplotlib.pyplot as plt
from turbo_seti.find_event.plot_dat import make_plot
from blimpy import Waterfall
import numpy as np
import pandas as pd
import blimpy as bl
from blimpy.utils import rebin
from turbo_seti import find_event
from turbo_seti.find_event import plot_event
import gc
from astropy.time import Time
from turbo_seti.find_event.plot_event import plot_waterfall
from turbo_seti.find_event.plot_event import overlay_drift


fontsize=16
font = {'family' : 'DejaVu Sans','size' : fontsize}
MAX_IMSHOW_POINTS = (4096, 1268)


def get_source_name(fil):
    return fil.split("guppi_")[1][12:].split("_")[0]


def one_wf_plot(folder,fil_file_list,f_start=0,f_stop=0,drift_rate=0,source_name_list=[]):

    # prepare for plotting
    matplotlib.rc('font', **font)

    # set up the sub-plots
    n_plots = len(fil_file_list)
    fig = plt.subplots(n_plots, sharex=True, sharey=True,figsize=(10, 2*n_plots))

    # read in data for the first panel
    max_load = bl.calcload.calc_max_load(fil_file_list[0])
    if not f_start:
        wf=Waterfall(fil_file_list[0],load_data=False)
        f_start=wf.container.f_start
        del wf
        gc.collect()
    if not f_stop:
        wf=Waterfall(fil_file_list[0],load_data=False)
        f_stop=wf.container.f_stop
        del wf
        gc.collect()
    wf1 = bl.Waterfall(fil_file_list[0], f_start=f_start, f_stop=f_stop, max_load=max_load)
    t0 = wf1.header['tstart']
    plot_f1, plot_data1 = wf1.grab_data()

    # rebin data to plot correctly with fewer points
    dec_fac_x, dec_fac_y = 1, 1
    if plot_data1.shape[0] > MAX_IMSHOW_POINTS[0]:
        dec_fac_x = plot_data1.shape[0] / MAX_IMSHOW_POINTS[0]
    if plot_data1.shape[1] > MAX_IMSHOW_POINTS[1]:
        dec_fac_y =  int(np.ceil(plot_data1.shape[1] /  MAX_IMSHOW_POINTS[1]))
    plot_data1 = rebin(plot_data1, dec_fac_x, dec_fac_y)

    # define more plot parameters
    mid_f = np.abs(f_start+f_stop)/2.

    subplots = []
    del wf1, plot_f1, plot_data1
    gc.collect()

    if not source_name_list:
        source_name_list=[get_source_name(f) for f in fil_file_list]

    on_source_name=get_source_name(fil_file_list[0])

    # Fill in each subplot for the full plot
    for ii, filename in enumerate(fil_file_list):
        # identify panel
        subplot = plt.subplot(n_plots, 1, ii + 1)
        subplots.append(subplot)

        # read in data
        max_load = bl.calcload.calc_max_load(filename)
        wf = bl.Waterfall(filename, f_start=f_start, f_stop=f_stop, max_load=max_load)

        # make plot with plot_waterfall
        source_name = source_name_list[ii]
        this_plot = plot_waterfall(wf,
                                    source_name,
                                    f_start=f_start,
                                    f_stop=f_stop,
                                    )

        if drift_rate:
            # calculate parameters for estimated drift line
            t_elapsed = Time(wf.header['tstart'], format='mjd').unix - Time(t0, format='mjd').unix
            t_duration = (wf.n_ints_in_file - 1) * wf.header['tsamp']
            f_event = mid_f + drift_rate / 1e6 * t_elapsed

            # plot estimated drift line
            overlay_drift(f_event, f_start, f_stop, drift_rate, t_duration, offset=0)

        # Title the full plot
        if ii == 0:
            if drift_rate:
                plot_title = "%s \n $\\dot{\\nu}$ = %2.3f Hz/s, MJD:%5.5f" % (on_source_name, drift_rate, t0)
                plt.title(plot_title)
            else:
                plot_title = "%s \n MJD:%5.5f" % (on_source_name, t0)
                plt.title(plot_title)
        # Format full plot
        if ii < len(fil_file_list)-1:
            plt.xticks(np.linspace(f_start, f_stop, num=4), ['','','',''])

        del wf
        gc.collect()

    # More overall plot formatting, axis labelling
    factor = 1e6
    units = 'Hz'

    xloc = np.linspace(f_start, f_stop, 5)
    xticks = [round(loc_freq) for loc_freq in (xloc - mid_f)*factor]
    if np.max(xticks) > 1000:
        xticks = [xt/1000 for xt in xticks]
        units = 'kHz'
    plt.xticks(xloc, xticks)
    plt.xlabel("Relative Frequency [%s] from %f MHz"%(units,mid_f),fontdict=font)

    # Add colorbar
    cax = fig[0].add_axes([0.94, 0.11, 0.03, 0.77])
    fig[0].colorbar(this_plot,cax=cax,label='Normalized Power (Arbitrary Units)')

    # Adjust plots
    plt.subplots_adjust(hspace=0,wspace=0)

    # save the figures
    dirpath=folder+'special_plots/'
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)
    path_png = dirpath + 'all_hits_' + on_source_name + '_dr_' + "{:0.2f}".format(drift_rate) + '_freq_' "{:0.6f}".format(mid_f) + ".png"
    plt.savefig(path_png, bbox_inches='tight')

    # plt.show()
    return None

def get_arrays(indir):
    data = pd.read_csv(indir+'concat_dats.csv')
    searchfile = open(indir+"info/drift_rates.txt","r").readlines()
    n_cadence=len(glob.glob(indir+'*0.h5'))
    BaryDriftRates=[]
    for line in searchfile:
        BaryDriftRates.append(float(line.split('Rate: ')[1].split(' Hz/s')[0]))
    # diffs=[(data["freq_start "][x]-data["freq_end "][x])/2*10**6 for x in range(len(data))]
    BaryDrifts = [-d for d in BaryDriftRates]
    max_drift = min(BaryDrifts)
    min_drift = max(BaryDrifts)
    # avg=sum(diffs)/len(data)
    data=data[(data["Drift_Rate "] < min_drift) & (data["Drift_Rate "]> max_drift)].reset_index(drop=True)
    Z=[(data["freq_start "][x]+data["freq_end "][x])/2 for x in range(len(data))]
    data["Freq"]=Z
    data=data.sort_values(by="Freq")
    f_starts=[]
    f_stops=[]
    DRs=[]
    f_mids=[]
    for x in range(len(data)):
        DR=data["Drift_Rate "].iloc[x]
        max_diff=abs(DR*300*n_cadence)/10**6
        f_mid=data["Freq"].iloc[x]
        f_start=round(f_mid-max_diff,6)
        f_stop=round(f_mid+max_diff,6)
        f_mid=np.abs(f_start+f_stop)/2.
        f_mids.append("{:0.6f}".format(f_mid))
        f_starts.append(f_start)
        f_stops.append(f_stop)
        DRs.append(DR)
    return f_starts, f_stops, DRs, f_mids


def main():
    bands=['C_band','L_band','S_band']
    for band in bands:
        folder=f'/gpfs/group/jtw13/default/gbt_2020/2021/{band}/splice_no_overlap_ordered/'
        f_starts,f_stops,DRs,f_mids=get_arrays(folder)
        data={'f_starts':f_starts,'f_stops':f_stops,'DRs':DRs,'f_mids':f_mids}
        df=pd.DataFrame(data=data)
        print(f'{len(df)} hits found above the SNR, within the DR range in any of the ON source pointings.')
        plots=sorted(glob.glob(folder+'special_plots/*.png'))
        if plots:
            for plot in plots:
                df=df[df["f_mids"]!=plot.split("freq_")[1].split(".png")[0]]
            print(f'{len(plots)} plots already made. Making {len(df)} more plots.')
        if len(df)==0:
            continue
        fils=sorted(glob.glob(folder+'*0.h5'))
        for x in range(len(df)):
            print(f'Plotting frequency {f_mids[x]} MHz over range: {f_starts[x]} to {f_stops[x]}')
            f_start=df.f_starts.iloc[x]
            f_stop=df.f_stops.iloc[x]
            drift_rate=df.DRs.iloc[x]
            one_wf_plot(folder,fils,f_start=f_start,f_stop=f_stop,drift_rate=drift_rate)
    return None
# run it!
if __name__ == "__main__":
    main()