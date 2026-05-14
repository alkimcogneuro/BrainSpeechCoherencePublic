function [ResultsCoherenceAnalysis] = Apply2Dataset_CrossSpectralDensity(eeg_struct, speech_rawdata, highpass_cutoff, lowpass_cutoff)
    % functionality to add:
    % think about how to flexibly adjust the width of the analysis epochs for speech and eeg
    % we may need to use sub-sections of the epochs that are provided. 

    % This function applies Cross Spectral Density analysis to one EEG dataset. 
    %
    % Arguments:
    %   eeg_struct:     structure containing EEG data and metadata
    %   speech_rawdata: structure containing raw speech data and metadata
    %   highpass_cutoff: high-pass filter cutoff frequency (Hz)
    %   lowpass_cutoff:  low-pass filter cutoff frequency (Hz)
    %   XX deleting this parameter.  we'll save the results in the function that calls this function.... results_filename:  filename to save the results
    %   we're not currently using this one:  randomize_onsets: boolean, if true adds random 3-10 sec shift to speech onset for control
    %
    % Returns:
    %   ResultsCoherenceAnalysis: structure containing phase coherence results
    
    %% highpass_cutoff = 35;   % default high-pass filter cutoff frequency in Hz (can be adjusted as needed)
    %% lowpass_cutoff = 2;     % default low-pass filter cutoff frequency in Hz (can be adjusted as needed)
    ResultsCoherenceAnalysis = struct();  % Initialize
    ResultsCoherenceAnalysis.analysis_type = 'cross_spectral_density';     
    ResultsCoherenceAnalysis.highpass_cutoff = highpass_cutoff;    
    ResultsCoherenceAnalysis.lowpass_cutoff = lowpass_cutoff;      
    ResultsCoherenceAnalysis.results_filename = results_filename;  
    ResultsCoherenceAnalysis.eeg_Fs = eeg_struct.Fs;                               
    ResultsCoherenceAnalysis.speech_Fs = speech_rawdata.Fs;                       
    ResultsCoherenceAnalysis.eeg_file_name = eeg_struct.eeg_file_name; 
    ResultsCoherenceAnalysis.label = eeg_struct.Subj_id;
    ResultsCoherenceAnalysis.Chanlocs = eeg_struct.Chanlocs;  

    fprintf('Running Cross Spectral Density analysis for subject %s\n', eeg_struct.Subj_id);
    fprintf('Filter parameters: HighPass = %d Hz, LowPass = %d Hz\n', highpass_cutoff, lowpass_cutoff);
    
    % Initialize results matrix: NumTrials x NumChannels
    cross_spectral_density_vals = cell(eeg_struct.Num_trials, eeg_struct.Num_channels);
    
    % Epoch duration in seconds
    speech_epoch_duration = size(eeg_struct.Data, 2) / eeg_struct.Fs;

    for eeg_trial_idx = 1:eeg_struct.Num_trials
        % Original EEG onset
        speech_onset_latency = eeg_struct.OnsetLatency(eeg_trial_idx);  
        speech_offset_latency = speech_onset_latency + speech_epoch_duration;

        % Convert seconds to sample indices
        speech_onset_idx = uint64(speech_onset_latency * speech_rawdata.Fs);
        speech_offset_idx = uint64(speech_offset_latency * speech_rawdata.Fs);

        % Ensure indices are within bounds
        speech_onset_idx = max(1, min(speech_onset_idx, length(speech_rawdata.Amplitudes)));
        speech_offset_idx = max(1, min(speech_offset_idx, length(speech_rawdata.Amplitudes)));

        % Extract speech epoch
        % Ren has already cut the speech into epochs that match the EEG epochs, 
        % so this should just be indexing into the correct segment of the speech data.

        speech_amplitudes = speech_rawdata.Amplitudes(speech_onset_idx:speech_offset_idx);

        % prepare the speech for coherence anaysis
        % Extract amplitude envelope, bandpass filter, downsample, save phase and magnitude information. 
        [Speech_Struct] = preprocess_speech_epoch(speech_amplitudes, speech_rawdata.Fs, ...
                                                 highpass_cutoff, lowpass_cutoff, eeg_struct.Num_samples);
        % Loop through EEG channels
        for ch_idx = 1:eeg_struct.Num_channels
            eeg_epoch = squeeze(eeg_struct.Data(ch_idx, :, eeg_trial_idx));
            if size(eeg_epoch, 2) ~= 1
                eeg_epoch = eeg_epoch';
            end
            % Band-pass filter EEG epoch
            eeg_epoch_bpf = band_pass_filt(eeg_epoch, eeg_struct.Fs, highpass_cutoff, lowpass_cutoff);
            % Save the cross spectral density value for this trial and channel
            cross_spectral_density_vals{eeg_trial_idx, ch_idx} = CrossSpectralDensity(Speech_Struct.envelope, eeg_epoch_bpf, Speech_Struct.fs, eeg_struct.Fs);
            %%        cross_spectral_density_vals(eeg_trial_idx, ch_idx) = CrossSpectralDensity(Speech_Struct.envelope, eeg_epoch_bpf, Speech_Struct.fs, eeg_struct.Fs);
        end
    end

    % Save trial-level CSD values in the main results structure.
    ResultsCoherenceAnalysis.cross_spectral_density_vals = cross_spectral_density_vals;
    % Average across trials per channel and save in the main results structure.
    

    %%% fix the calculation of by-channel means (across trials), so that it works for CSD 
    

%{
 
    ResultsCoherenceAnalysis.cross_spectral_density_chanmeans = mean(cross_spectral_density, 1);

    % Save maximum coherence and corresponding channel
    [ResultsCoherenceAnalysis.max_val, ResultsCoherenceAnalysis.max_chan]  = ...
        max(ResultsCoherenceAnalysis.cross_spectral_density_chanmeans);
 
 
%}

end
