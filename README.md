# Brain Speech Coherence Functions
This suite of Matlab code performs brain-speech coherence (aka "entrainment") anlayses.  
 
Key functions include
 - run_CSD_dataset_list.m:
  - acts on a list of filenames, each containing an EEG dataset (one subject, one condition) and runs CSD analysis.  
 - Apply2Dataset_CrossSpectralDensity.m
 - CrossSpectralDensity.m

# running an analysis
- See the script `TestScripts/runscript_MSC_example.m` for an example of how to run MSC analyses on a group of datasets.
  - The main function of this script is to organize data for and run `run_CSD_dataset_list`
    - `run_CSD_dataset_list(eeg_data_filenames, speech_files_path, analysis_ID_label, highpass_cutoff, lowpass_cutoff, options)`
      - eeg_data_filenames:  a cell array of strings; each string is a filepath to one EEG data file.
      - speech_files_path:  a string, which provides the path the folder containing the .wav speech files
      - highpass_cutoff:  lower bound of the pass band for filtering prior to CSD analyses (e.g., 2.6 Hz)
      - lowpass_cutoff:  upper bound of the pass band for filtering prior to CSD analysis (e.g., 35). 
      - optional StartTimeOffset:  % optional offset to the start time of the analysis epoch, in seconds, relative to the original EEG onset.
      - optional EpochDuration     % optional epoch duration, in seconds-- if NaN, use the full length of the EEG epochs as the analysis epoch duration.

  ## how to write a run script.  
     
     - you'll create a list of filenames that contains the EEG data files
 
```
eegfiles = {
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj01.mat', 
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj02.mat',
'~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj03.mat'}
```

  - alternative syntax:
```
eegfiles{1} = '~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj01.mat'`
eegfiles{2} = '~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj02.mat''
eegfiles{3} = '~/DATA/MyStudy_Data/Condition1epochs/Condition1_Subj03.mat''
```
     
  - You'll also need to specify the directory containing all speech files needed for the anlaysis.
     `speech_files_path = '~/MYDATA/SpeechFilesEntrainment/';`
       - The specified folder should contain all speech files (e.g., `audiosample21.wav`).
       - The names of the individual .wav files should be recorded in the field `EEG_struct.audio_file` of the EEG data structures.

  - The analysis code assumes that the critical data is in a matlab structure called `EEG_struct`, whose format should as exemplified below:

```
EEG_struct = 

  struct with fields:

            Data: [64×4000×14 single]
        Chanlocs: [1×64 struct]
              Fs: 1000
      Num_trials: 14
    Num_channels: 64
     Num_samples: 4000
         Subj_id: '01'
          TrigID: {1×14 cell}
       Condition: "MyCondition"
      audio_file: "AudioSample1"
    OnsetLatency: [0 4 16 24 32 36 52 64 68 80 88 92 100 104]
```

   - Saving the results of MSC analyses. The analysis function `runrun_CSD_dataset_list` will write return resuls as Matlab structures and will also write the results to files in a folder that includes the date and time of the analysis
     - For example, the results folder CSD_Analysis_Results_Condition1_matched_speech_2026_0620_1729_58 contains the results for the main analysis of the condition 'Condition1', with matched speech controls and was created by a run on 2026-06-20 at 17:29:58. 

# Understanding the analysis functions.

  [CSD_one_sided, PSD_speech_one_sided, PSD_eeg_one_sided, MSPC] = CrossSpectralDensity(sig_speech, sig_eeg, nfft) 
input arguments:
    sig_speech:  vector of speech data (should be same length and sample rate as eeg data--the speech should be downsampled). 
    sig_eeg:    vector of eeg data
    nfft:  the number of points used in the Fourier transform. we'll typically make this equal to the legnth of the signals. 
   
   # plotting results
   - see `TestScripts/plot_results_script.m` for an exampple of how to plot topographic maps of the results of MSC analyses
   - working now to make this script more transparent and user friendly
    
