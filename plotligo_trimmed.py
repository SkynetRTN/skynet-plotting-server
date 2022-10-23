# Standard python numerical analysis imports:
import numpy as np
import base64
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt
from gwpy.timeseries import TimeSeries
import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')

# LIGO-specific readligo.py 
import readligo as rl

# done importing now

# function to whiten data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    norm = 1. / np.sqrt(1. / (dt * 2))
    white_hf = hf / np.sqrt(interp_psd(freqs)) * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

#does all the work
def get_data_from_file(file_name, whiten_data=0, plot_spectrogram=0):

    # frequency band for bandpassing signal
    fband = [43.0, 400.0]
    # fband = [35.0, 350.0]
    fs = 16384
 #   fs = 4096
    
    
    #----------------------------------------------------------------
    # Load LIGO data from a single file.
    # FIRST, define the filenames fn_H1 and fn_L1, above.
    #----------------------------------------------------------------
    try:
        # read in data from H1 and L1, if available:
        strain, time, chan_dict, duration, file_fs, gpsStart, gpsEnd = rl.loaddata(file_name)
    except Exception as e:
        print(e)
        print("Cannot find data files! - - " + file_name)
        print("Quitting.")
        quit()

    if (duration == 4096):
        raise Exception("Please upload a 32s file")
    if (file_fs != 16384):
        raise Exception("Please upload a 16Khz file")
    #----------------------------------------------------------------
    # Grab time info and print some stuff
    #----------------------------------------------------------------
    
    # the time sample interval (uniformly sampled!)
    dt = time[1] - time[0]
    


    if whiten_data:
        # number of sample for the fast fourier transform:
        NFFT = 4*fs
        Pxx, freqs = mlab.psd(strain, Fs = fs, NFFT = NFFT)
    
        # We will use interpolations of the ASDs computed above for whitening:
        psd = interp1d(freqs, Pxx)

        # now whiten the data from H1 and L1, and the template (use H1 PSD):
        strain_whiten = whiten(strain,psd,dt)
        
        # We need to suppress the high frequency noise (no signal!) with some bandpassing:
        bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
        normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
        
        # final product here
        strain_whitenbp = filtfilt(bb, ab, strain_whiten) / normalization
        return np.concatenate((time.reshape(-1,1), strain_whitenbp.reshape(-1,1)), axis=1)

    # Finds the spectrogram and plots it. Returns the figure and the spectrogram object
    if plot_spectrogram:
        timewindow = 0.05  # hardcoded to plot 5 seconds centered on merger - can change

        ## Using GWPY to make graph
        strain_timeseries = TimeSeries(strain, times=time)

        midpoint = (gpsStart + gpsEnd)/2
        window = (gpsEnd-gpsStart)*timewindow

        hq = strain_timeseries.q_transform(outseg=(midpoint-window, midpoint+window), norm=False)
        fig = hq.plot()
        # vmin=0, vmax=25
        ax = fig.gca()
        fig.colorbar(label="Energy")
        ax.set(xlim=(midpoint-window, midpoint+window))
        ax.grid(False)
        ax.set_yscale('log')


        return fig, hq

def find_merger(spec_array, THRESHOLD=5 * 10**-14):
    row_event_count = []
    col_event_count = []
    # print(spec_array[0].shape)
    col_maxval_count = np.zeros(spec_array.shape[1])
    row_maxval_count = np.zeros(spec_array.shape[0])
    for row in spec_array:
        boolrow = row >= THRESHOLD
        row_event_count.append(np.count_nonzero(boolrow))

        argmax_row = np.argmax(row)
        col_maxval_count[argmax_row] += 1
    for col in np.transpose(spec_array):
        boolcol = col >= THRESHOLD
        col_event_count.append(np.count_nonzero(boolcol))

        argmax_col = np.argmax(col)
        row_maxval_count[argmax_col] += 1

    return col_maxval_count, row_maxval_count, col_event_count, row_event_count


#150914
file_name = "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"
#190814
# file_name = "L-L1_GWOSC_16KHZ_R1-1249852241-32.hdf5"

figor, hq = get_data_from_file(file_name, plot_spectrogram=1)
print(hq)

col_maxval_count, row_maxval_count, col_event_count, row_event_count = find_merger(hq)

figor.savefig("specplot_test.png")
fig = matplotlib.pyplot.figure(5)
ax = fig.add_subplot()
ax.plot(hq.yindex, col_maxval_count)
ax.set(xscale="log")
fig.savefig("ydim_maxvalcount.png")

fig = matplotlib.pyplot.figure(6)
ax = fig.add_subplot()
ax.plot(hq.xindex, row_maxval_count)
fig.savefig("xdim_maxvalcount.png")

fig = matplotlib.pyplot.figure(7)
ax = fig.add_subplot()
ax.plot(hq.xindex, row_event_count)
fig.savefig("xdim_eventcount.png")

fig = matplotlib.pyplot.figure(8)
ax = fig.add_subplot()
ax.plot(hq.yindex, col_event_count)
ax.set(xscale="log")
fig.savefig("ydim_eventcount.png")



