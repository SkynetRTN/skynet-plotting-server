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
    

    #----------------------------------------------------------------
    # Plot the Amplitude Spectral Density (ASD)
    #
    # The ASDs are the square root of the power spectral densities (PSDs), which are averages of the square of the fast fourier transforms (FFTs) of the data.
    #
    # They are an estimate of the "strain-equivalent noise" of the detectors versus frequency, which limit the ability of the detectors to identify GW signals.
    

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

    # create_spectrogram = 1
    if plot_spectrogram:
        deltat = 0.3  # hardcoded to plot 5 seconds centered on merger - can change

        ## Using GWPY to make graph
        strain_timeseries = TimeSeries(strain, times=time)

        midpoint = (gpsStart + gpsEnd)/2
        window = (gpsEnd-gpsStart)*0.05

        hq = strain_timeseries.q_transform(outseg=(midpoint-window, midpoint+window), norm=False)
        fig = hq.plot()
        # vmin=0, vmax=25
        ax = fig.gca()
        fig.colorbar(label="Energy")
        ax.set(xlim=(midpoint-window, midpoint+window))
        ax.grid(False)
        ax.set_yscale('log')


        return fig, hq


# figor, hq = get_data_from_file("L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5", plot_spectrogram=1)
# print(hq)
# # print(hq.name)
# # print(hq.x0)
# # print(hq.dx)
# # print(hq.dy)
# # print(hq.y0)
# # print(hq.yindex)
# # print(hq.xindex)
# figor.savefig("specplot_test.png")
# extracted_model = extract_model_from_spectrogram(1.175, 62.797, 1126259462.4, hq)
# fig = matplotlib.pyplot.figure(60)
# fig.gca().plot(extracted_model[:,0], extracted_model[:,1])
# fig.savefig("extraction_test.png")


