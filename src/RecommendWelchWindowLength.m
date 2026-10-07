function WindowInfo = RecommendWelchWindowLength(trial_duration_sec, fs, min_freq_hz, options)
    % ================================================================================
    % Recommend / validate a Welch analysis-window length for a given trial duration
    % ================================================================================
    % Two uses:
    %   1) ADVISORY: call with just trial_duration_sec, fs, min_freq_hz (and optionally
    %      MinCyclesPerWindow) to get the FLOOR window length needed to resolve min_freq_hz
    %      reliably, plus how many segments that floor window would give you at <=50% overlap,
    %      for that trial length. Use this once per new dataset (e.g., before committing to a
    %      WindowDurationSec for a dataset whose trials are much shorter than 30 seconds) to see
    %      what you're actually getting, before wiring the value into
    %      PreprocessSpeechEEGForCoherence.m / ShuffleTest_SpeechEEGCoherence.m.
    %   2) VALIDATION: pass your chosen options.WindowDurationSec, and this function checks it
    %      against the floor (errors if it's shorter than the floor -- i.e., too short to resolve
    %      min_freq_hz), checks it fits within the trial (errors if longer than the trial), caps
    %      OverlapFraction at 0.5 (errors above that), and warns if the resulting segment count is
    %      low. This is what PreprocessSpeechEEGForCoherence.m calls internally, once per dataset
    %      (not per trial/shuffle), so a bad configuration fails fast with a clear message instead
    %      of quietly degrading -- or erroring deep inside a 1000-shuffle loop.
    %
    % Arguments:
    %   trial_duration_sec   : length of one trial's epoch, in seconds.
    %   fs                   : EEG sampling rate (Hz).
    %   min_freq_hz          : lowest frequency you need this analysis to resolve reliably
    %                          (e.g., your highpass_cutoff, or a specific band edge you care about).
    %   options.MinCyclesPerWindow (default 3) : minimum number of full cycles of min_freq_hz that
    %                          must fit inside one analysis window for that window to be considered
    %                          long enough to resolve it. 3 is a conservative, standard choice; 2 is
    %                          more aggressive (shorter windows, more segments, slightly less stable
    %                          per-segment estimate at the low-frequency edge).
    %   options.OverlapFraction (default 0.5) : must be in [0, 0.5] -- see WelchCrossSpectralDensity.m
    %                          header for why more than 50%% isn't accepted here.
    %   options.WindowDurationSec (default NaN = advisory mode) : if provided, this function
    %                          validates THIS window length instead of just reporting the floor.
    %   options.MinSegments (default 3) : segment counts below this trigger a warning (not an
    %                          error) -- the analysis can still run, but with little of the
    %                          averaging benefit Welch segmentation is meant to provide.
    %
    % Returns WindowInfo, a struct with fields:
    %   .FloorWindowDurationSec, .FloorWindowSamples : the minimum window length that resolves
    %                          min_freq_hz at MinCyclesPerWindow cycles, in seconds and samples.
    %   .WindowDurationSec, .WindowSamples : the window length actually used for the report below
    %                          -- equals options.WindowDurationSec if given, else the floor.
    %   .OverlapFraction, .StepSamples, .NumSegments : the overlap/step/segment-count that this
    %                          window length + trial_duration_sec would produce.
    %   .TrialDurationSec, .MinFreqHz, .MinCyclesPerWindow : echoed back for bookkeeping.
    % ================================================================================
    arguments
        trial_duration_sec (1,1) double {mustBePositive}
        fs (1,1) double {mustBePositive}
        min_freq_hz (1,1) double {mustBePositive}
        options.MinCyclesPerWindow (1,1) double {mustBePositive} = 3
        options.OverlapFraction (1,1) double {mustBeGreaterThanOrEqual(options.OverlapFraction, 0)} = 0.5
        options.WindowDurationSec (1,1) double {mustBeReal} = NaN
        options.MinSegments (1,1) double {mustBePositive, mustBeInteger} = 3
    end

    if options.OverlapFraction > 0.5
        error('RecommendWelchWindowLength:OverlapTooHigh', ...
            ['options.OverlapFraction = %.4f exceeds the 0.5 (50%%) cap. Segments overlapping more than ' ...
             '50%% share the majority of their samples with their neighbor and stop behaving like separate ' ...
             'observations in the average -- see WelchCrossSpectralDensity.m header.'], options.OverlapFraction);
    end

    FloorWindowDurationSec = options.MinCyclesPerWindow / min_freq_hz;
    FloorWindowSamples = ceil(FloorWindowDurationSec * fs);

    if FloorWindowDurationSec > trial_duration_sec
        error('RecommendWelchWindowLength:TrialTooShortForFrequency', ...
            ['Resolving %.3g Hz reliably (>= %g cycles per window) needs a window of at least %.3f sec, ' ...
             'but the trial itself is only %.3f sec long. This trial length cannot resolve %.3g Hz at all ' ...
             'with this method -- either accept a higher min_freq_hz (fewer cycles required), lower ' ...
             'options.MinCyclesPerWindow, or use a different trial/epoch length.'], ...
            min_freq_hz, options.MinCyclesPerWindow, FloorWindowDurationSec, trial_duration_sec, min_freq_hz);
    end

    if isnan(options.WindowDurationSec)
        % Advisory mode: report the floor itself as the recommended window length.
        WindowDurationSec = FloorWindowDurationSec;
    else
        WindowDurationSec = options.WindowDurationSec;
        if WindowDurationSec < FloorWindowDurationSec
            error('RecommendWelchWindowLength:WindowBelowFloor', ...
                ['options.WindowDurationSec = %.3f sec is shorter than the %.3f sec floor needed to resolve ' ...
                 '%.3g Hz at %g cycles/window. Either use a window >= %.3f sec, or relax min_freq_hz / ' ...
                 'MinCyclesPerWindow if you genuinely don''t need to resolve %.3g Hz.'], ...
                WindowDurationSec, FloorWindowDurationSec, min_freq_hz, options.MinCyclesPerWindow, ...
                FloorWindowDurationSec, min_freq_hz);
        end
        if WindowDurationSec > trial_duration_sec
            error('RecommendWelchWindowLength:WindowLongerThanTrial', ...
                'options.WindowDurationSec = %.3f sec is longer than the trial itself (%.3f sec).', ...
                WindowDurationSec, trial_duration_sec);
        end
    end

    WindowSamples = round(WindowDurationSec * fs);
    TrialSamples = round(trial_duration_sec * fs);
    StepSamples = max(round(WindowSamples * (1 - options.OverlapFraction)), 1);
    NumSegments = floor((TrialSamples - WindowSamples) / StepSamples) + 1;

    if NumSegments < options.MinSegments
        warning('RecommendWelchWindowLength:FewSegments', ...
            ['This window length (%.3f sec) gives only %d segment(s) from a %.3f sec trial at %.0f%% overlap ' ...
             '(minimum recommended: %d). With this few segments, the Welch estimate gets little of the ' ...
             'bias/variance benefit it''s meant to provide over a single whole-trial estimate. Consider a ' ...
             'shorter window if min_freq_hz/MinCyclesPerWindow allow it.'], ...
            WindowDurationSec, NumSegments, trial_duration_sec, options.OverlapFraction*100, options.MinSegments);
    end

    WindowInfo = struct();
    WindowInfo.FloorWindowDurationSec = FloorWindowDurationSec;
    WindowInfo.FloorWindowSamples = FloorWindowSamples;
    WindowInfo.WindowDurationSec = WindowDurationSec;
    WindowInfo.WindowSamples = WindowSamples;
    WindowInfo.OverlapFraction = options.OverlapFraction;
    WindowInfo.StepSamples = StepSamples;
    WindowInfo.NumSegments = NumSegments;
    WindowInfo.TrialDurationSec = trial_duration_sec;
    WindowInfo.MinFreqHz = min_freq_hz;
    WindowInfo.MinCyclesPerWindow = options.MinCyclesPerWindow;
end
