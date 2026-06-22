# Brain Speech Coherence Functions
This suite of Matlab code performs brain-speech coherence (aka "entrainment") anlayses.  
 
Key functions include
 - run_CSD_dataset_list.m:
  - acts on a list of filenames, each containing an EEG dataset (one subject, one condition) and runs CSD analysis.  
 - Apply2Dataset_CrossSpectralDensity.m
 - CrossSpectralDensity.m

# running an analysis
- The script `TestScripts/runscript_MSC_example.m` provides a template for running MSC analyses on a group of datasets.
  - Users can create a script like this to manage CSD/MSC analyses. 
  - The main function of this script is to organize data for and run `run_CSD_dataset_list`
  - view the example script to understand how to run it.
  - `run_CSD_dataset_list(eeg_data_filenames, speech_files_path, analysis_ID_label, highpass_cutoff, lowpass_cutoff, options)`
    - eeg_data_filenames:  a cell array of strings; each string is a filepath to one EEG data file.
    - speech_files_path:  a string, which provides the path the folder containing the .wav speech files
    - highpass_cutoff:  lower bound of the pass band for filtering prior to CSD analyses (e.g., 2.6 Hz)
    - lowpass_cutoff:  upper bound of the pass band for filtering prior to CSD analysis (e.g., 35). 
    - optional StartTimeOffset:  % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.
    - optional EpochDuration     % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.
   ## how to write a run script like this.  
     - you'll create a list of filenames that contains the EEG data files
       -   for example:  eegfiles{1} = '~/EntrainmentData/S1_C1_eeg.mat'
       - for example:  eegfiles{2} = '~/EntrainmentData/S2_C1_eeg.mat'
     - speech_files_path = '~/MYDATA/SpeechFilesEntrainment/';
       - % tells us where to find the speech .wav files.
       - The names of the individual .wav files should match the audio_file field in the EEG data structures.
   

  [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, MSPC] = CrossSpectralDensity(sig_speech, sig_eeg, nfft) 
input arguments:
    sig_speech:  vector of speech data (should be same length and sample rate as eeg data--the speech should be downsampled). 
    sig_eeg:    vector of eeg data
    nfft:  the number of points used in the Fourier transform. we'll typically make this equal to the legnth of the signals. 
