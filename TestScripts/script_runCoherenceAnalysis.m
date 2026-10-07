% =====================================================================================
% Pipeline for Welch-style CSD/MSC analysis 
% =====================================================================================
% Purpose: run the new WelchCrossSpectralDensity.m-based pipeline 
%
% What to look for when you run this:
%   1) It completes without error for each subject (errors are caught and reported per
%      subject below, rather than stopping the whole run).
%   2) The per-subject "Welch analysis window" line printed by
%      PreprocessSpeechEEGForCoherence.m (during each SpeechEEGCoherence_ShuffleTest
%      call) shows WindowDurationSec=3, ~50% overlap, and a segment count > 1 -- if it's
%      reporting 1 segment/trial, the Welch change isn't actually doing anything for that
%      subject's trial length and something's off (e.g. WindowDurationSec not being passed
%      through, or that subject's trials are unexpectedly short).
%   3) The summary table below shows True_MSC values in [0, 1], no NaN/Inf, and channel
%      counts that make sense for your montage.
%   4) Compare the printed NumSegmentsPerTrial against what you'd expect from your trial
%      length -- e.g. for run RecommendWelchWindowLength(30, EEG_struct.Fs, 2,
%      WindowDurationSec=3, OverlapFraction=0.5) by hand for one subject's fs and see that
%      .NumSegments matches ShuffleTestResults.NumSegmentsPerTrial for that subject.
%
% If something looks wrong: PreprocessSpeechEEGForCoherence.m, ComputeCoherenceForPairing.m,
% SpeechEEGCoherence_ShuffleTest.m, and GroupLevel_MSC_Analysis.m were all edited on
% 2026-09-17 as part of this change, each with a same-day backup saved alongside it
% (<name>_backup_20260917.m) -- useful for a side-by-side diff if you need to track down
% where a discrepancy is coming from.
%
% Each subject's SpeechEEGCoherence_ShuffleTest results are saved
% to their OWN file as soon as that subject finishes (so a crash or a later subject's error
% never costs you the subjects already done):
%     <Path>/Coherence_bandMSC_control_analyses/coherence_control_analysis_S01.mat
% Each file holds: Coherence_Analysis_ShuffleTestResults (the results struct itself, same
% fields as before), Subj_id, eeg_data_file, and RunParameters (this run's parameters, so each
% file is self-describing). The "S01" tag is 'S' plus EEG_struct.Subj_id from that subject's
% own data file (a char array like '01') -- see subj_tag below. A subject that errors gets NO file (rather
% than a misleading empty one); the error is printed and logged.
% The same folder also gets coherence_control_analysis_parameters.txt, a human-readable record
% of all the main parameters, plus one outcome line per subject -- so a folder is self-describing
% about which analysis produced it, even after being renamed aside as a backup (see below).
% output_dir is set directly in the "USER SHOULD SET THESE PARAMETERS" block below, like every
% other analysis parameter. If output_dir already exists when this script starts (e.g. you forgot
% to update it from a previous run), the EXISTING folder is renamed aside by appending
% "_backup_<timestamp>" to its name, rather than being written into or overwritten (see
% open-items.md item 20) -- this run then always starts from a fresh, empty output_dir, so
% individual per-subject files are never at risk of collision within a run.
% To reload a previous run for plotting, use dev/LoadCoherenceResultsFromFolder.m (it does the loop
% below for you, across every per-subject file in a folder, and falls back to a cell array if the
% loaded structs don't all share identical fields):
%     ResultsArray = LoadCoherenceResultsFromFolder(output_dir);
%     PlotSubjectsAndGroupTopography(ResultsArray, 'True_MSC', OutputDir=fullfile(output_dir, 'plots'), ...
%         SubjectIDField='Subj_id', ColorbarLabel='True MSC');
% ('True_MSC' always exists; 'PValues'/'ZScores'/etc. only exist on files saved with NumShuffles >= 1.)
% To reload just ONE subject's file by hand instead:
%     S = load(fullfile(output_dir, 'coherence_control_analysis_S01.mat'));
%     S.Coherence_Analysis_ShuffleTestResults      % <- the results struct for S01
% =====================================================================================


% --- ANALYSIS PARAMETERS -----------------------------------------------------------
% 
% speech_files_path:   Specify the folder (full path) where the speech audio files (wav) for this analysis are located. 
%                       Each subject's EEG_struct.audio_file field should point to the correct file in that folder.
%                       The code assumes that all the speech audio files (wav) for this analysis are in one folder, 
%                       and that the EEG_struct.audio_file field for each subject points to the correct file in that folder. 
% 
% highpass_cutoff:      Prior to the coherence analysis, we'll apply a broadband filter to both the speech and EEG data. 
%                       Here we set the highpass_cutoff, in Hz.  Frequencies below this will be attenuated.
%                       This is important for removing slow drifts and low-frequency noise from the data.
% lowpass_cutoff:       low pass cutoff, in Hz, for broadband filter. Frequencies above this will be attenuated. 
%                       This is important for removing high-frequency noise from the data.
% 
% band_low:             after computing coherence at each of many frequency bins, 
%                       we will average the coherence values across a specific, narrower frequency band (e.g.,  3-7 Hz).
%                       here, we set the lower bound of that band, in Hz (e.g., 3).
% band_high:            upper bound of the frequency band for averaging coherence, in Hz (e.g., 7).
% 
% WindowDurationSec:    This parameter sets the duration (width) of those windows for Welch-style segmentation.
%                       The idea is that we'll divide each individual trial into smaller, overlapping windows of this width, 
%                       and compute coherence for each window.  Then we'll average the coherence values across 
%                       all windows within a trial to get a single coherence estimate for that trial.
%                       example width: 3 seconds (WindowDurationSec = 3).
%                       Note:  in addition to this within-trial averaging, we will also average across trials to get a final coherence estimate for the subject.
%                       The window duration should be long enough to capture the frequencies of interest,
%                       -- for example, if you are interested in 3-7 Hz, a window of at least 1 second is needed 
%                       to capture at least one full cycle of the lowest frequency (3 Hz).
% OverlapFraction:      fractional overlap of Welch windows (0.5 = 50% overlap, 0.25 = 25% overlap, etc.)
%                       In Welch-style segmentation, windows can overlap. This parameter sets the fraction of overlap 
%                       overlap should not be too high (e.g., > 0.5) or too low (e.g., < 0.25) to ensure a good balance between frequency resolution and statistical reliability.
%                       
% NumShuffles           Our shuffle (permutation) test shuffles the data a number of times to estimate 
%                       a null distribution of coherence values.
%                       This parameter sets the number of shuffles to perform. 
%                       A higher number of shuffles gives a more accurate estimate of the null distribution, but takes longer to compute.
%                       Setting NumShuffles = 0 skips the shuffle test, just compute the true MSC.
%                       Note:  when NumShuffles gets large, the analysis can take a long time, 
%                       so test with NumShuffles = 2 first to make sure everything is working, then increase to a larger number (e.g., 1000) for the final analysis.
% RandomSeed:           random seed for the shuffle test (NaN = don't seed, i.e. use a different random seed each run)
% AvgMethod:            This parameter describes how we will average the MSC across the targetted band [band_low - band_high].
%                       'power-weighted' is the default, which weights the MSC values by the power at each frequency bin before averaging.
%                       'simple' is another option, which simply averages the MSC values across the band without weighting by power.
% output_dir:           name the output folder where the per-subject .mat files and the parameters/run-log text file will be saved.
%                       If the selected folder already exists, it will be renamed aside as a backup (with a timestamp appended to its name) before this run starts, 
%                       so that the new run always starts with a fresh, empty output folder.
%                       A log file in each output folder records the parameters used for that run, 
%                       so you can always check what was actually used even if you forget to update this script.
%                       If output_dir is not set by the caller before running this script, it will default to '~/DATA/ACSEEG_Control30s/Coherence_bandMSC_control_analyses'.
%                       Note: i haven't checked whether this works on Windows, but it should.  The '~' is expanded via HOME so fopen/mkdir/save all agree on it.
% 
% --- Welch windowing parameters -----------------------------------------------------
% Currently the individual trials are ~30 sec.
% We will divide each longer trial into shorter windows (e.g., WindowDurationSec = 3).
% We'll analzye coherence in each trial with a Welch-style segmentation, then average the results across segments within each trial.
% 
% We'll use a ~50% overlap (comfortable margin above the
% frequency-resolution floor for this highpass_cutoff -- see RecommendWelchWindowLength.m).
% Deliberately explicit rather than left at the NaN/auto-floor default, which would pick
% a much shorter (1.5 sec) window instead -- see open-items.md item 12's "operational note".

% Small NumShuffles and a fixed RandomSeed so this first run is fast and reproducible
% while you're checking correctness. For a real analysis, bump NumShuffles back up
% (1000 is what the rest of this codebase uses) and leave RandomSeed unset (NaN).

% How MSC is averaged across the [band_low, band_high] band (passed to BandAverageMSC.m via
% SpeechEEGCoherence_ShuffleTest.m). Set explicitly here, rather than left to the function's
% default, so that the parameter record written below always states what was actually used.


% -------------------------------------------------------------------------------------------------------------------
% USER SHOULD SET THESE PARAMETERS BEFORE RUNNING THIS SCRIPT
% THE DEFAULT PARAMETERS HERE ARE FROM DEVELOPMENT, AND MAY NOT BE APPROPRIATE FOR YOUR ANALYSIS.
% FILE PATHS, IN PARTICULAR, MUST BE SET TO MATCH YOUR LOCAL FILE SYSTEM AND YOUR DATA ORGANIZATION.
% -------------------------------------------------------------------------------------------------------------------
speech_files_path = '~/DATA/ACSEEG_Control30s/Control30sParticipantAudioFiles'; % path to the speech audio files (wav) for this analysis
highpass_cutoff = 1.75;     % Hz, for the broadband filter applied to both speech and EEG before CSD estimation
lowpass_cutoff = 35;        % Hz, for the broadband filter applied to both speech and EEG before CSD estimation
band_low = 3;               % Hz, for the band over which MSC is averaged (passed to BandAverageMSC.m via SpeechEEGCoherence_ShuffleTest.m)
band_high = 7;              % Hz, for the band over which MSC is averaged (passed to BandAverageMSC.m via SpeechEEGCoherence_ShuffleTest.m)
WindowDurationSec = 3;      % window width in seconds for Welch segmentation of each trial's CSD/PSD estimate
OverlapFraction = 0.5;      % fractional overlap of Welch windows (0.5 = 50% overlap, 0.25 = 25% overlap, etc.)
NumShuffles = 2;            % number of shuffles in the shuffle test (0 = skip the shuffle test, just compute the true MSC).
RandomSeed = 1;             % random seed for the shuffle test (NaN = don't seed, i.e. use a different random seed each run)
AvgMethod = 'power-weighted';   % how to average the MSC across the [band_low, band_high] band (passed to BandAverageMSC.m via SpeechEEGCoherence_ShuffleTest.m)
output_dir = '~/DATA/ACSEEG_Control30s/Coherence_bandMSC_control_analyses';  % where to save the per-subject .mat files and the parameters/run-log text file

% ------------------------------------------------------------------------------------------------------------------------
% USER SHOULD SPECIFY THE EEG DATA FILES FOR THE SUBJECTS TO ANALYZE BELOW
% User should set eeg_data_filenames to the list of EEG_struct .mat files to analyze. 
% Each .mat file contains EEG data for one subject, and
% should contain a variable EEG_struct with fields Subj_id and audio_file, among others.
% The .audio_file field should point to the corresponding speech audio file in speech_files_path.
% 
% Note:  it is recommended to first test with NumShuffles = 2 to ensure everything is working correctly, then increase to a larger number (e.g., 1000) for the final analysis.
% Note:  Test this code on a short list of 1-2 subjects, before expanding to the full list once you are confident it is working correctly.
% Note:  The full analysis with NumShuffles = 1000 can take a long time (possibly >1 hour per subject).  
%           therefore, you might want to run this analysis for a few subjects at a time, 
%           and then combine the results later for group-level analysis.
% ------------------------------------------------------------------------------------------------------------------------
eeg_data_filenames = {
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct01.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct02.mat',
     '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct04.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct05.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct06.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct07.mat'}
%{
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct08.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct09.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct10.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct11.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct12.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct13.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct14.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct15.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct16.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct17.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct18.mat', 
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct19.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct20.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct21.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct22.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct23.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct24.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct25.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct26.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct27.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct28.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct29.mat',
    '~/DATA/ACSEEG_Control30s/EEG_structs/EEG_struct30.mat'};
 
%}


% -------------------------------------------------------------------------------------------------------------------

% --- OUTPUT ---------------------------------------------------------------------------
% One .mat file per subject, plus a parameters/run-log text file, all in a dedicated folder.
% ~ is expanded via HOME so fopen/mkdir/save all agree on it.
%
% output_dir itself is set above, in the "USER SHOULD SET THESE PARAMETERS" block, the same way
% as every other analysis parameter. If that folder already exists (see "Set up the output
% folder" below), the existing folder is renamed aside as a backup rather than written into --
% this run always gets a fresh, empty directory, so there's no need for a per-file overwrite
% check later on.
output_file_prefix = 'coherence_control_analysis_';            % + subj_tag + '.mat', e.g. coherence_control_analysis_S01.mat
params_log_filename = 'coherence_control_analysis_parameters.txt';


% =====================================================================================
% Run
% =====================================================================================
Num_test_subjects = length(eeg_data_filenames);

% --- Set up the output folder and the parameter record ----------------------------------
% If output_dir already exists (e.g. the caller forgot to change it from a previous run, or
% simply re-ran with the unchanged default), rename the EXISTING folder out of the way --
% appending "_backup_<timestamp>" to its name -- rather than writing into it or erroring. This
% is a safety net against silently losing or mixing in a previous run's results; the timestamp
% keeps each backup name unique even if this happens more than once in a row. Every folder
% (including each backup) gets its own coherence_control_analysis_parameters.txt below, so
% there's no ambiguity later about which analysis a given folder's contents came from. Once
% any rename is done, output_dir is (re)created fresh and EMPTY, so the per-subject saves in
% the loop below never need to check for an existing file -- everything in this output_dir is
% guaranteed to be from THIS run.
if isfolder(output_dir)
    backup_timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    [output_parent, output_name] = fileparts(output_dir);
    backup_dir = fullfile(output_parent, sprintf('%s_backup_%s', output_name, backup_timestamp));
    movefile(output_dir, backup_dir);
    fprintf('NOTE: %s already existed -- moved it to %s before starting this run.\n', output_dir, backup_dir);
end
mkdir(output_dir);
params_log_file = fullfile(output_dir, params_log_filename);

% Every main parameter of this run in one struct. It is (a) written to the human-readable
% log file just below and (b) saved inside every per-subject .mat file.
RunParameters = struct( ...
    'RunStarted',          char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'ScriptName',          mfilename, ...
    'MatlabRelease',       version('-release'), ...
    'speech_files_path',   speech_files_path, ...
    'eeg_data_filenames',  {eeg_data_filenames(:)'}, ...
    'highpass_cutoff',     highpass_cutoff, ...
    'lowpass_cutoff',      lowpass_cutoff, ...
    'band_low',            band_low, ...
    'band_high',           band_high, ...
    'AvgMethod',           AvgMethod, ...
    'WindowDurationSec',   WindowDurationSec, ...
    'OverlapFraction',     OverlapFraction, ...
    'NumShuffles',         NumShuffles, ...
    'RandomSeed',          RandomSeed, ...
    'output_dir',          output_dir);

% L is a cell array of lines to write to the parameters/run-log text file. 
% It is written in one go below, so that if MATLAB crashes or a later subject errors, 
% the log file still has a complete record of the run parameters and outcomes for all subjects that finished successfully.

L = {};
L{end+1} = repmat('=', 1, 90);
L{end+1} = sprintf('COHERENCE CONTROL ANALYSIS -- run started %s', RunParameters.RunStarted);
L{end+1} = repmat('=', 1, 90);
L{end+1} = sprintf('Script:                     %s.m  (MATLAB %s)', RunParameters.ScriptName, RunParameters.MatlabRelease);
L{end+1} =         'Analysis function:         SpeechEEGCoherence_ShuffleTest.m (Welch-segmented CSD -> band-averaged MSC, derangement shuffle test)';
L{end+1} = sprintf('Output folder:              %s', output_dir);
L{end+1} = sprintf('Per-subject files:          %s<subj_tag>.mat, e.g. %sS01.mat', output_file_prefix, output_file_prefix);
L{end+1} =         '  variables in each file:  Coherence_Analysis_ShuffleTestResults, Subj_id, eeg_data_file, RunParameters';
L{end+1} =         '';
L{end+1} =         'DATA';
L{end+1} = sprintf('  speech_files_path:        %s', speech_files_path);
L{end+1} = sprintf('  EEG data files (%d):', Num_test_subjects);
for k = 1:Num_test_subjects
    L{end+1} = sprintf('    %2d) %s', k, eeg_data_filenames{k}); %#ok<SAGROW>
end
L{end+1} =         'BROADBAND FILTER (applied before CSD estimation)';
L{end+1} = sprintf('  highpass_cutoff (Hz):     %g', highpass_cutoff);
L{end+1} = sprintf('  lowpass_cutoff (Hz):      %g', lowpass_cutoff);
L{end+1} =         'BAND FOR MSC AVERAGING';
L{end+1} = sprintf('  band_low / band_high (Hz):%g / %g', band_low, band_high);
L{end+1} = sprintf('  AvgMethod:                %s', AvgMethod);
L{end+1} =         'WELCH WINDOWING (within each trial)';
L{end+1} = sprintf('  WindowDurationSec:        %g', WindowDurationSec);
L{end+1} = sprintf('  OverlapFraction:          %g', OverlapFraction);
L{end+1} =         'SHUFFLE TEST';
L{end+1} = sprintf('  NumShuffles:              %d', NumShuffles);
L{end+1} = sprintf('  RandomSeed:               %g   (NaN = not seeded)', RandomSeed);
L{end+1} =         'NOT SET BY THIS SCRIPT (left at SpeechEEGCoherence_ShuffleTest.m defaults): StartTimeOffset, EpochDuration,';
L{end+1} =         '  MinFreqHz, MinCyclesPerWindow. The values each subject actually ended up with (sampling rate, nfft, segments';
L{end+1} =         '  per trial) are logged in that subject''s line below.';
L{end+1} =         '';
L{end+1} =         'SUBJECT OUTCOMES';
append_log(params_log_file, '%s\n', strjoin(L, newline));   % write the parameter record to the log file now, so it survives a crash or a later subject's error
clear L k

RunSummary = [];    % one small scalar-only struct per subject, for the summary table -- NOT the results themselves

for idx = 1:Num_test_subjects
    eeg_data_file = eeg_data_filenames{idx};
    fprintf('\n=== [%d/%d] %s ===\n', idx, Num_test_subjects, eeg_data_file);
    % Placeholder tag, used only in messages if this subject fails BEFORE its ID can be read
    % (e.g. the file won't load). It is replaced by the real tag, built from EEG_struct.Subj_id,
    % right after the load below; that real tag is what names the output file.
    [~, eeg_file_stem] = fileparts(eeg_data_file);
    subj_tag = regexprep(eeg_file_stem, '[^A-Za-z0-9]', '');
    subject_output_file = '';

    % Reset per-subject state so an errored subject can never save the PREVIOUS subject's results.
    Coherence_Analysis_ShuffleTestResults = [];
    Subj_id = '';
    thisSummary = struct('SubjTag', subj_tag, 'Subj_id', '', 'Status', 'ERRORED', 'SavedFile', '', ...
        'ElapsedSec', NaN, 'NumChannels', NaN, 'eeg_Fs', NaN, 'nfft', NaN, 'AvgMethod', '', ...
        'WindowDurationSec', NaN, 'OverlapFraction', NaN, 'NumSegmentsPerTrial', NaN, ...
        'MSC_min', NaN, 'MSC_max', NaN, 'MinP', NaN, 'NumSig', NaN, 'NumNonFinite', NaN, ...
        'MSCOutOfRange', false, 'Error', '');
    subject_timer = tic;

    try
        load(eeg_data_file);   % loads EEG_struct
        Subj_id = EEG_struct.Subj_id;

        % Output-file tag = 'S' + the subject ID stored in the data file itself
        % (EEG_struct.Subj_id, a char array such as '04' -> 'S04').
        if isnumeric(Subj_id)
            subj_id_str = sprintf('%02d', Subj_id);
        else
            subj_id_str = regexprep(char(Subj_id), '[^A-Za-z0-9]', '');
        end
        subj_tag = ['S' subj_id_str];
        thisSummary.Subj_id = subj_id_str;   % char version, safe to print with %s
        thisSummary.SubjTag = subj_tag;
        subject_output_file = fullfile(output_dir, [output_file_prefix subj_tag '.mat']);
        fprintf(' -- EEG data corresponds to AudioFile: %s\n', EEG_struct.audio_file);

        speech_full_filepath = fullfile(speech_files_path, EEG_struct.audio_file);
        [speech_amplitudes, fs_speech] = audioread(speech_full_filepath);
        Speech_RawData = struct;
        Speech_RawData.Amplitudes = speech_amplitudes;
        Speech_RawData.Fs = fs_speech;
        
        % this is the main call to the new Welch-style CSD/MSC pipeline, 
        % which does the shuffle test if NumShuffles >= 1
        Coherence_Analysis_ShuffleTestResults = SpeechEEGCoherence_ShuffleTest( ...
            EEG_struct, Speech_RawData, band_low, band_high, highpass_cutoff, lowpass_cutoff, ...
            WindowDurationSec=WindowDurationSec, OverlapFraction=OverlapFraction, ...
            NumShuffles=NumShuffles, RandomSeed=RandomSeed, AvgMethod=AvgMethod);

        % --- Save this subject's results to their own file --------------------------------
        % No overwrite check needed here: output_dir is guaranteed fresh and empty for this run
        % (see "Set up the output folder" above -- any pre-existing folder at this path was
        % already renamed aside as a backup before the loop started), so a collision within one
        % run would only happen if two subjects mapped to the same subj_tag, which would be a
        % real bug worth surfacing rather than something to silently paper over.
        save(subject_output_file, 'Coherence_Analysis_ShuffleTestResults', 'Subj_id', 'eeg_data_file', 'RunParameters');
        fprintf(' -- Saved %s\n', subject_output_file);
        
        % --- Scalars for the summary table (not the results themselves) ---------------------
        R = Coherence_Analysis_ShuffleTestResults;
        thisSummary.Status              = 'OK';
        thisSummary.SavedFile           = [output_file_prefix subj_tag '.mat'];
        thisSummary.NumChannels         = numel(R.True_MSC);
        thisSummary.eeg_Fs              = R.eeg_Fs;
        thisSummary.nfft                = R.nfft;
        thisSummary.AvgMethod           = R.AvgMethod;
        thisSummary.WindowDurationSec   = R.WindowDurationSec;
        thisSummary.OverlapFraction     = R.OverlapFraction;
        thisSummary.NumSegmentsPerTrial = R.NumSegmentsPerTrial;
        thisSummary.MSC_min             = min(R.True_MSC);
        thisSummary.MSC_max             = max(R.True_MSC);
        % Shuffle-test-derived fields (.PValues/.ZScores) only exist when NumShuffles >= 1 --
        % SpeechEEGCoherence_ShuffleTest.m returns early (no PValues/ZScores/Null_MSC/etc.) when
        % NumShuffles < 1, since no shuffle test was run. Leave thisSummary.MinP/.NumSig at their
        % NaN default (set above) in that case rather than reading a field that was never computed.
        if isfield(R, 'PValues')
            thisSummary.MinP             = min(R.PValues);
            thisSummary.NumSig           = sum(R.PValues < 0.05);
        end
        thisSummary.NumNonFinite        = sum(~isfinite(R.True_MSC));
        thisSummary.MSCOutOfRange       = any(R.True_MSC < 0 | R.True_MSC > 1);
        
    catch ME
        fprintf(2, ' !! ERROR for subject %s (%s): %s\n', subj_tag, thisSummary.Subj_id, ME.message);
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        thisSummary.Error = ME.message;
    end
    thisSummary.ElapsedSec = toc(subject_timer);
    RunSummary = [RunSummary, thisSummary]; %#ok<AGROW>
    
    % --- One line per subject in the parameter/run log (written immediately, so it survives a crash) ---
    if strcmp(thisSummary.Status, 'OK')
        append_log(params_log_file, ['  %s  Subj_id=%s  OK       saved %s | channels=%d, eeg_Fs=%g Hz, nfft=%d, ' ...
        'segments/trial=%d, AvgMethod=%s | %.0f s\n'], subj_tag, thisSummary.Subj_id, thisSummary.SavedFile, ...
        thisSummary.NumChannels, thisSummary.eeg_Fs, thisSummary.nfft, thisSummary.NumSegmentsPerTrial, ...
        thisSummary.AvgMethod, thisSummary.ElapsedSec);
    else
        append_log(params_log_file, '  %s  Subj_id=%s  ERRORED  (no file saved) %s | %.0f s\n', ...
        subj_tag, thisSummary.Subj_id, strrep(thisSummary.Error, newline, ' '), thisSummary.ElapsedSec);
    end
end
clear Coherence_Analysis_ShuffleTestResults R   % clear last subject's big struct for the next subject, so it doesn't get accidentally saved if the next subject errors

append_log(params_log_file, '\nRun finished %s -- %d/%d subjects saved.\n\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), sum(strcmp({RunSummary.Status}, 'OK')), Num_test_subjects);

% =====================================================================================
% Sanity-check summary -- one line per subject, so problems are easy to spot at a glance.
% (Built from RunSummary, the small per-subject scalar record above; the full results are in
% the per-subject .mat files in output_dir.)
% =====================================================================================
fprintf('\n===================================================================\n');
fprintf('SUMMARY (%d subjects, NumShuffles=%d, WindowDurationSec=%g, OverlapFraction=%g)\n', ...
    Num_test_subjects, NumShuffles, WindowDurationSec, OverlapFraction);
fprintf('===================================================================\n');
% The table has two shapes: with p-value/significance columns when the shuffle test actually ran
% (NumShuffles >= 1), and without them when it was skipped (NumShuffles == 0) -- there is nothing
% shuffle-test-derived to show in that case, so we don't print placeholder NaN columns for it.
if NumShuffles >= 1
    fprintf('%-10s %8s %8s %10s %10s %10s %10s %14s\n', ...
        'Subject', 'Window', 'Overlap', 'Segs/Trial', 'MSC_min', 'MSC_max', 'MinP', 'NumChanSig(p<.05)');
else
    fprintf('NumShuffles=0 -- shuffle test was skipped, so no p-value/significance columns below.\n');
    fprintf('%-10s %8s %8s %10s %10s %10s\n', ...
        'Subject', 'Window', 'Overlap', 'Segs/Trial', 'MSC_min', 'MSC_max');
end

for idx = 1:Num_test_subjects
    Ssum = RunSummary(idx);

    if ~strcmp(Ssum.Status, 'OK')
        fprintf('%-10s   ERRORED -- see message above (%s)\n', Ssum.SubjTag, Ssum.Error);
        continue;
    end

    if NumShuffles >= 1
        fprintf('%-10s %8.3f %7.0f%% %10d %10.4f %10.4f %10.4f %14d\n', ...
            Ssum.SubjTag, Ssum.WindowDurationSec, Ssum.OverlapFraction*100, Ssum.NumSegmentsPerTrial, ...
            Ssum.MSC_min, Ssum.MSC_max, Ssum.MinP, Ssum.NumSig);
    else
        fprintf('%-10s %8.3f %7.0f%% %10d %10.4f %10.4f\n', ...
            Ssum.SubjTag, Ssum.WindowDurationSec, Ssum.OverlapFraction*100, Ssum.NumSegmentsPerTrial, ...
            Ssum.MSC_min, Ssum.MSC_max);
    end

    if Ssum.NumNonFinite > 0
        fprintf('   !! %d/%d channels have non-finite True_MSC for %s -- investigate before trusting this subject.\n', ...
            Ssum.NumNonFinite, Ssum.NumChannels, Ssum.SubjTag);
    end
    if Ssum.MSCOutOfRange
        fprintf('   !! True_MSC out of [0,1] range for %s -- this should never happen; something is wrong.\n', Ssum.SubjTag);
    end
    if Ssum.NumSegmentsPerTrial <= 1
        fprintf(['   !! Only %d segment/trial for %s -- the Welch change isn''t doing anything for this ' ...
                 'subject''s trial length. Check trial duration vs. WindowDurationSec.\n'], Ssum.NumSegmentsPerTrial, Ssum.SubjTag);
    end
end
clear Ssum

fprintf('\nDone. Per-subject results are in:\n    %s\n', output_dir);
fprintf('Parameters for this run are recorded in:\n    %s\n', params_log_file);
fprintf('If this all looks sane, re-run with the full eeg_data_filenames list and your\n');
fprintf('normal NumShuffles (1000) for the real analysis. RunSummary (small per-subject scalars) is\n');
fprintf('left in the workspace; the full results live in the per-subject .mat files.\n');

% Example: reload this run's saved files and plot them (uncomment to run right after this script,
% or paste into a fresh session later -- output_dir just needs to point at the folder with the
% saved .mat files; LoadCoherenceResultsFromFolder.m does the per-file loading loop for you):
%
% ResultsArray = LoadCoherenceResultsFromFolder(output_dir);
% PlotSubjectsAndGroupTopography(ResultsArray, 'True_MSC', OutputDir=fullfile(output_dir, 'plots'), ...
%     SubjectIDField='Subj_id', ColorbarLabel='True MSC');
% % Only meaningful if this run used NumShuffles >= 1 (otherwise .PValues/.ZScores don't exist):
% % PlotSubjectsAndGroupTopography(ResultsArray, 'PValues', OutputDir=fullfile(output_dir, 'plots'), ...
% %     SubjectIDField='Subj_id', ColorbarLabel='p-value');


% =====================================================================================
% Local helper (script-local functions must come at the very end of the file)
% =====================================================================================
function append_log(logfile, fmt, varargin)
    % Append one formatted chunk of text to the parameters/run-log file. Opens and closes the
    % file on every call so nothing is lost if MATLAB crashes or a later subject errors.
    fid = fopen(logfile, 'a');
    if fid < 0
        warning('script_runTestWelchCoherence:LogWriteFailed', 'Could not open %s for appending.', logfile);
        return;
    end
    fprintf(fid, fmt, varargin{:});
    fclose(fid);
end
