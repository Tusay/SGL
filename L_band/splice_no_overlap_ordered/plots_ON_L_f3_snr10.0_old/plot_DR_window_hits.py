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
            f_mid = (f_start+f_stop)/2
            f_event = f_mid + drift_rate / 1e6 * t_elapsed

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
    path_png = dirpath + 'all_hits_' + on_source_name + '_dr_' + "{:0.2f}".format(drift_rate) + '_freq_' "{:0.6f}".format(f_start) + ".png"
    plt.savefig(path_png, bbox_inches='tight')

    # plt.show()
    return None

def main():
    f_starts=[8437.9743432356,3870.116560212,3840.8196390538,3820.3117723938,3811.522695175,3799.803902015,3750.0514187119998,
    8437.974300735601,8312.4645893938,3749.9372403938,3687.4507163938,6943.358741272599,5821.2882522356,4040.0383651750003,
    4028.3195862355997,3966.7961943938,3937.5036329932,8334.960211053802,3925.7805577726,3905.272667175,3902.3429913938,3896.4835500149998,
    3873.046217212,3864.2571197725997,3858.397776712,3840.8195593938,3764.6477335538,3750.0188113332,8437.9872489932,8437.6872077726,
    8437.3871999932,4702.147781993201,4693.3587189932005,3925.7805161750002,3914.0617808938,3823.2415165538,3791.0149108938,
    3750.009829015,3749.9367929932,3625.0058058938,8393.5539777356,8367.186770515002,8358.3978404932,8343.749408493199,8329.1008527356,
    8317.3821870538,3916.9914467356,3893.5540390538004,3878.9056217726,3846.678984175,3840.8196984932,3793.9445635149996,3761.7180898332003]
    f_stops=[8437.9760307644,3870.117937788,3840.8211199462,3820.3133566062,3811.524313825,3799.8056239850002,3750.0527962879996,
    8437.9759882644,8312.4661736062,3749.9388246062003,3687.4523006062004,6943.3601877273995,5821.2899397644005,4040.039983825,
    4028.3212737644,3966.7977786062,3937.5050450068,8334.9616919462,3925.7820042274,3905.274285825,3902.3445756062,3896.485271985,
    3873.047594788,3864.2585662274,3858.399154288,3840.8211436062,3764.6492144462,3750.0203266668,8437.9886610068,8437.6886542274,
    8437.3886120068,4702.149194006801,4693.360131006801,3925.782134825,3914.0633651062003,3823.2429974461998,3791.0164951062,
    3750.011550985,3749.9382050068,3625.0073901062,8393.555665264399,8367.188492485,8358.3992525068,8343.750820506799,8329.102540264399,
    8317.383667946198,3916.9931342644004,3893.5555199462,3878.9070682274005,3846.680602825,3840.8211105068003,3793.946285485,3761.7196051668]
    DRs=[-0.468758,-0.38266,-0.411359,-0.440059,-0.449625,-0.478325,-0.38266,-0.468758,-0.440059,-0.440059,-0.440059,-0.401793,-0.468758,
    -0.449625,-0.468758,-0.440059,-0.392226,-0.411359,-0.401793,-0.449625,-0.440059,-0.478325,-0.38266,-0.401793,-0.38266,-0.440059,-0.411359,
    -0.420926,-0.392226,-0.401793,-0.392226,-0.392226,-0.392226,-0.449625,-0.440059,-0.411359,-0.440059,-0.478325,-0.392226,-0.440059,
    -0.468758,-0.478325,-0.392226,-0.392226,-0.468758,-0.411359,-0.468758,-0.411359,-0.401793,-0.449625,-0.392226,-0.478325,-0.420926]
    folder=f'/gpfs/group/jtw13/default/gbt_2020/2021/C_band/splice_no_overlap_ordered/'
    os.chdir(folder)
    fils=sorted(glob.glob('*0.h5'))
    for x in range(len(DRs)):
        print(f'Frequency Range: {f_starts[x]} to {f_stops[x]}')
        f_start=f_starts[x]
        f_stop=f_stops[x]
        drift_rate=DRs[x]
        one_wf_plot(folder,fils,f_start=f_start,f_stop=f_stop,drift_rate=drift_rate)
        return None
# run it!
if __name__ == "__main__":
    main()