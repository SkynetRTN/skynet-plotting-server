# Standard python numerical analysis imports:
import numpy as np
import base64
import os
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt
from scipy.signal.windows import tukey
from gwpy.timeseries import TimeSeries
import matplotlib.mlab as mlab
import matplotlib
from scipy.special import erfinv

matplotlib.use('Agg')

# LIGO-specific readligo.py 
import readligo as rl

# done importing now

# Separate bandpassing function for use on whitened data

def bandpassData(freqHigh, freqLow, strain_whiten, time):
    # We need to suppress the high frequency noise (no signal!) with some bandpassing:
    fband = [freqLow, freqHigh]
    fs = 16384
    bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
        
    # final product here
    strain_whitenbp = filtfilt(bb, ab, strain_whiten) / normalization
    
    return np.concatenate((time.reshape(-1,1), strain_whitenbp.reshape(-1,1)), axis=1)


# function to whiten data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    # pronorm = np.sqrt(interp_psd(freqs)).max() ----- > trying to do normalizations the same as models
    norm = 1. / np.sqrt(1. / (dt * 2))
    white_hf = hf / np.sqrt(interp_psd(freqs)) * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

# Finds a resonable window for the event given a spectrogram
# spectro - gwpy spectrogram object
# mad_tresh - percentile to threshold the normalized medians of valid points
# buffer sclaers - add percentage of total spec size to each axis
def chauvenet_windowing(spectro, mad_thresh=68.0, buffer_scalers=[0.05, 0.05]):
    # Step 1: Set mu = 0

    spectro_np = np.asarray(spectro)
    
    # Step 2: Calculate sigma = stdev = rms
    rms = np.sqrt(np.mean(np.square(spectro_np)))
    sigma = rms
    
    # Step 3: Apply Chauvenet rejection
    N = spectro_np.shape[0] * spectro_np.shape[1]
    z = np.abs(spectro_np - 0) / sigma

    threshold_chauvenet = np.sqrt(2) * sigma * erfinv(1 - 0.5 / N)
    is_outlier = z > threshold_chauvenet
    
    # Step 4: Apply 68% median filtering
    
    valid_region = np.where(is_outlier, spectro_np, 0)
    med_spec = np.median(valid_region)
    abs_deviation = np.abs(valid_region - med_spec)
    mad_spec = np.median(abs_deviation)
    threshold = np.percentile(mad_spec, mad_thresh)
    is_significant = valid_region > threshold
    # Step 5: Define the window
    time_range = np.where(np.any(is_significant, axis=1))[0]
    if len(time_range) == 0:
        # No significant time range found
        print("No significant time range found")
        return (0,0)
    else:
        time_start_indx, time_end_indx = time_range[0], time_range[-1] + 1

    freq_range = np.where(np.any(is_significant, axis=0))[0]
    if len(freq_range) == 0:
        # No significant time range found
        print("No significant frequency range found")
        return (0,0)
    else:
        freq_start_indx, freq_end_indx = freq_range[0], freq_range[-1] + 1
        
        
    dx = spectro.dx
    x0 = spectro.x0
    dy = 0.5 #always 0.5 but for some reason doesn't come through metadata sometimes
    y0 = spectro.y0
    
    x_buffer = np.floor(spectro_np.shape[0] * buffer_scalers[0]) if spectro_np.shape[0] > len(time_range) else 0
    y_buffer = np.floor(spectro_np.shape[1] * buffer_scalers[1]) if spectro_np.shape[1] > len(freq_range) else 0
    time_lb = (time_start_indx - x_buffer) * dx + x0
    time_ub = (time_end_indx + x_buffer) * dx + x0
    freq_lb = (freq_start_indx - y_buffer) * dy + y0.value if freq_start_indx > y_buffer else y0.value
    freq_ub = (freq_end_indx + y_buffer) * dy + y0.value if freq_end_indx < spectro_np.shape[1] - y_buffer else spectro_np.shape[1] * dy + y0.value
    
    return (time_lb.value, time_ub.value), (freq_lb, freq_ub)

#does all the work
## Fundamental change - - only outputs whitened unbandpassed data when uploading files
def get_data_from_file(file_name, whiten_data=0, plot_spectrogram=0):
    # frequency band for bandpassing signal
    # fband = [freqLow, freqHigh]
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
        
        import h5py

        hdf5_file = h5py.File(file_name, 'r')
        timeOfRecord = hdf5_file['meta']['UTCstart'][()]
        timeOfRecord = timeOfRecord.decode('utf-8')
        timeZero =  hdf5_file['meta']['GPSstart'][()]
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
    strain_timeseries = TimeSeries(strain, times=time)


    if whiten_data:
        # number of sample for the fast fourier transform:
        NFFT = 4*fs
        psd_window = tukey(NFFT, alpha=1./4)
        Pxx, freqs = mlab.psd(strain, Fs = fs, NFFT = NFFT, window=psd_window)
    
        # We will use interpolations of the ASDs computed above for whitening:
        psd = interp1d(freqs, Pxx)

        # now whiten the data from H1 and L1, and the template (use H1 PSD):

        ## HUGE REMINDER i have no clue if you have to scale before or after whitening, here we will find out!!
        # for i in range(len(strain)):
        #     strain[i] = strain[i] * 26
        strain_whiten = whiten(strain,psd,dt)
        
        # # We need to suppress the high frequency noise (no signal!) with some bandpassing:
        # bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
        # normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
        
        # # final product here
        # strain_whitenbp = filtfilt(bb, ab, strain_whiten) / normalization

        return strain_whiten, time, psd, strain_timeseries, timeOfRecord, timeZero
        # return np.concatenate((time.reshape(-1,1), strain_whitenbp.reshape(-1,1)), axis=1)

    # Finds the spectrogram and plots it. Returns the figure and the spectrogram object
    if plot_spectrogram:
        timewindow = 0.10  # hardcoded to plot 5 seconds centered on merger - can change

        ## Using GWPY to make graph

        midpoint = (gpsStart + gpsEnd)/2
        window = (gpsEnd-gpsStart)*timewindow
        lowerbound = midpoint-window
        upperbound = midpoint+window
        
    
        og_spectrogram = strain_timeseries.q_transform(outseg=(lowerbound, upperbound), norm=True)
        time_bounds, freq_bounds = chauvenet_windowing(og_spectrogram, mad_thresh=68, buffer_scalers=[0.025, 0.05])
        if not (time_bounds or freq_bounds):
            fig = og_spectrogram.plot()
            ax = fig.gca()
            fig.colorbar(label="Energy")
            ax.set(xlim=(lowerbound, upperbound))
            ax.grid(False)
            ax.set_yscale('log')
            return fig, og_spectrogram
        hq = strain_timeseries.q_transform(outseg=(time_bounds[0], time_bounds[1]), norm=True)
        fig = hq.plot()
        # vmin=0, vmax=25
        ax = fig.gca()
        fig.colorbar(label="Energy")
        ax.set(xlim=(time_bounds[0], time_bounds[1]), ylim=(freq_bounds[0], freq_bounds[1]))
        
        ax.grid(False)
        ax.set_yscale('log')

        return fig, hq
