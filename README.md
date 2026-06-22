This suite of Matlab code performs brain-speech coherence (aka "entrainment") anlayses. 
 
Key functions include
 - CrossSpectralDensity.m
 - Apply2Dataset_CrossSpectralDensity.m
 - run_CSD_dataset.m

# running an analysis
- `manage_analysis_cross_spectral_density(eeg_data_filenames, speech_files_path, highpass_cutoff, lowpass_cutoff, options)`
- eeg_data_filenames:  a cell array of strings; each string is a filepath to one EEG data file.
 - speech_files_path:  a string, which provides the path the folder containing the .wav speech files
 - highpass_cutoff:  we will high pass filter the EEG and speech at this cutoff (e.g., 2)
 - lowpass_cutoff:  we will low pass filter the EEG and speech at this cutoff (e.g., 35). 
 - optional StartTimeOffset:  % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.
 - options.EpochDuration  % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
   ## how to use this function
    you'll create a list of filenames that contains the EEG data files
   for example:  eegfiles{1} = '~/EntrainmentData/S1_C1_eeg.mat'
   for example:  eegfiles{2} = '~/EntrainmentData/S2_C1_eeg.mat'
speech_files_path = '~/SpeechFilesEntrainment/';  % tells us where to find the speech .wav files. The names of the individual .wav files should match the audio_file field in the EEG data structures.
   you'll also need a folder that contains all speech files
   for example:  
Function usage:  
function [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, MSPC] = CrossSpectralDensity(sig_speech, sig_eeg, nfft) 
input arguments:
    sig_speech:  vector of speech data (should be same length and sample rate as eeg data--the speech should be downsampled). 
    sig_eeg:    vector of eeg data
    nfft:  the number of points used in the Fourier transform. we'll typically make this equal to the legnth of the signals. 
