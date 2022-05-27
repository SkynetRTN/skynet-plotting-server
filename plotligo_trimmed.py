# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt
import h5py
import json

import matplotlib.mlab as mlab

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
def perform_whitening_on_file(file_name):

    # frequency band for bandpassing signal
    fband = [43.0, 400.0]
    fs = 4096
    
    
    #----------------------------------------------------------------
    # Load LIGO data from a single file.
    # FIRST, define the filenames fn_H1 and fn_L1, above.
    #----------------------------------------------------------------
    try:
        # read in data from H1 and L1, if available:
        strain, time, chan_dict, fs_fake = rl.loaddata(file_name)
    except Exception as e:
        print(e)
        print("Cannot find data files! - - " + file_name)
        print("Quitting.")
        quit()
    
    #----------------------------------------------------------------
    # Grab time info and print some stuff
    #----------------------------------------------------------------
    
    # the time sample interval (uniformly sampled!)
    dt = time[1] - time[0]
    
    # Let's look at the data and print out some stuff:
    # print('time: len, min, mean, max = ', \
    #     len(time), time.min(), time.mean(), time.max() )
    # print('strain: len, min, mean, max = ', \
    #     len(strain), strain.min(),strain.mean(),strain.max())
    
    
    #----------------------------------------------------------------
    # Plot the Amplitude Spectral Density (ASD)
    #
    # The ASDs are the square root of the power spectral densities (PSDs), which are averages of the square of the fast fourier transforms (FFTs) of the data.
    #
    # They are an estimate of the "strain-equivalent noise" of the detectors versus frequency, which limit the ability of the detectors to identify GW signals.
    
    make_psd = 1
    if make_psd:
        # number of sample for the fast fourier transform:
        NFFT = 4*fs
        Pxx, freqs = mlab.psd(strain, Fs = fs, NFFT = NFFT)
    
        # We will use interpolations of the ASDs computed above for whitening:
        psd = interp1d(freqs, Pxx)
    
    #----------------------------------------------------------------
    # Whitening bysuppressing the extra noise at low frequencies 
    # and a the spectral lines. Also, bandpass the high frequency noise.
    #----------------------------------------------------------------
    
    
    
    whiten_data = 1
    if whiten_data:
        # now whiten the data from H1 and L1, and the template (use H1 PSD):
        strain_whiten = whiten(strain,psd,dt)
        
        # We need to suppress the high frequency noise (no signal!) with some bandpassing:
        bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
        normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
        
        # final product here
        strain_whitenbp = filtfilt(bb, ab, strain_whiten) / normalization
        return time, strain_whitenbp

path = 'temp-grav-data'
file_name = 'H-H1_LOSC_4_V2-1126259446-32.hdf5'
print(perform_whitening_on_file(path + '/' + file_name))



# #----------------------------------------------------------------
# # Audio files
# #----------------------------------------------------------------
#
# # make wav (sound) files from the whitened data, +-2s around the event.
# from scipy.io import wavfile
#
# # function to keep the data within integer limits, and write to wavfile:
# def write_wavfile(filename,fs,data):
#     d = np.int16(data/np.max(np.abs(data)) * 32767 * 0.9)
#     wavfile.write(filename,int(fs), d)
#
# deltat_sound = 2. # seconds around the event
#
# # index into the strain time series for this time interval:
# indxd = np.where((time >= tevent-deltat_sound) & (time < tevent+deltat_sound))
#
#
# # # write the files:
# # write_wavfile(eventname+"_H1_whitenbp.wav",int(fs), strain_H1_whitenbp[indxd])
# # write_wavfile(eventname+"_L1_whitenbp.wav",int(fs), strain_L1_whitenbp[indxd])
#
#
# # function that shifts frequency of a band-passed signal
# def reqshift(data,fshift=100,sample_rate=4096):
#     """Frequency shift the signal by constant
#     """
#     x = np.fft.rfft(data)
#     T = len(data)/float(sample_rate)
#     df = 1.0/T
#     nbins = int(fshift/df)
#     # print T,df,nbins,x.real.shape
#     y = np.roll(x.real,nbins) + 1j*np.roll(x.imag,nbins)
#     y[0:nbins]=0.
#     z = np.fft.irfft(y)
#     return z
#
# # parameters for frequency shift
# # fs = 4096
# # fshift = 400.
# # speedup = 1.
# # fss = int(float(fs)*float(speedup))
# #
# # # shift frequency of the data
# # strain_H1_shifted = reqshift(strain_H1_whitenbp,fshift=fshift,sample_rate=fs)
# # strain_L1_shifted = reqshift(strain_L1_whitenbp,fshift=fshift,sample_rate=fs)
# #
# # # write the files:
# # write_wavfile(eventname+"_H1_shifted.wav",int(fs), strain_H1_shifted[indxd])
# # write_wavfile(eventname+"_L1_shifted.wav",int(fs), strain_L1_shifted[indxd])

