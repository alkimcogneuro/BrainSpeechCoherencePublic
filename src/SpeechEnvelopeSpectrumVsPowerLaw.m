function [Result, SpeechEnvelopes] = SpeechEnvelopeSpectrumVsPowerLaw(Speech_RawData, options)
    % ===============================================================================================================
    % SPEECH AMPLITUDE-ENVELOPE SPECTRUM vs. A 1/f (POWER-LAW) NULL
    % ===============================================================================================================
    % Question this answers: does the speech amplitude envelope contain MORE energy in a low-frequency band (classically
    % ~4-8 Hz, the syllable rate) than a smooth 1/f^beta spectrum would predict?
    %
    % This is a thin speech-specific front end for SpectrumVsPowerLaw.m (which does the spectral estimation, the
    % flank-only power-law fit, and the Monte-Carlo / bootstrap statistics -- read that header for the method and its
    % caveats). All this function does is produce the envelope epochs, using the SAME code path as the coherence pipeline
    % (PreprocessSpeechEEGForCoherence.m -> preprocess_speech_epoch.m -> band_pass_filt.m, hilbert, resample), so the
    % envelope whose spectrum you inspect is the same kind of envelope that goes into the brain-speech coherence.
    %
    % ONE DELIBERATE DIFFERENCE from the coherence pipeline: the pipeline band-pass filters the envelope at 2-35 Hz
    % (highpass_cutoff / lowpass_cutoff). That is fine for coherence (a zero-phase filter applied to both signals cancels
    % in the magnitude-squared coherence ratio) but it is WRONG for a spectral-shape analysis, because band_pass_filt's
    % effective response is |H(f)| = 1/(1+(fc/f)^8) per edge (filtfilt squares a 4th-order Butterworth, so the
    % nominal cutoff is the -6 dB point). With a 2 Hz high-pass the envelope is already down by 48 dB at 1 Hz and 21 dB
    % at 1.5 Hz, and with a 35 Hz low-pass by 2 dB at 30 Hz -- an artificial "spectral shape" right where a 1/f
    % fit needs its flanks. So here preprocess_speech_epoch is called with a much wider passband (defaults 0.5-64 Hz,
    % flat to within 0.1 dB from 1 to 35 Hz), and the fit is done well inside it.
    %
    % Arguments
    %   Speech_RawData : struct with .Amplitudes (raw waveform, samples x 1; multi-channel audio is averaged) and .Fs
    %                    (Hz) -- the same convention as in your run scripts.
    %
    % Options (name-value)
    %   EpochDurationSec (30)   Epoch length. The recording is cut into consecutive epochs of this length (the last
    %                           partial epoch is dropped) and each epoch is processed independently, exactly like the
    %                           coherence pipeline processes each trial. Inf = analyse the whole recording as ONE piece
    %                           (no epoch-edge filter/Hilbert transients, but no bootstrap-over-epochs either).
    %   OnsetLatency ([])       Optional epoch onsets, in seconds into the recording (e.g. EEG_struct.OnsetLatency),
    %                           to analyse exactly the speech epochs that entered the EEG analysis. Requires a finite
    %                           EpochDurationSec (which should equal the EEG trial duration).
    %   EnvelopeHighpassHz (0.5), EnvelopeLowpassHz (64)   Passed to preprocess_speech_epoch as its cutoffs (see above).
    %   TargetFs (256)          Sampling rate the envelope is resampled to (Hz). Use >= 2.5 x the top of FitRangeHz.
    %   SpectrumOptions (struct())   Name/value pairs forwarded to SpectrumVsPowerLaw, e.g.
    %                           struct('BandHz',[3 8],'NumSurrogates',500,'Exponent',1). See its header.
    %   MakePlot (true)         Call PlotSpectrumVsPowerLaw on the result.
    %   PlotTitle ('')          Title for the figure (default: auto-generated).
    %
    % Returns
    %   Result          the struct returned by SpectrumVsPowerLaw (plus .SpeechFs, .EpochDurationSec,
    %                   .NumEpochs, .EnvelopeHighpassHz, .EnvelopeLowpassHz, .EpochOnsetSec, .Figure)
    %   SpeechEnvelopes [Num_samples_per_epoch x Num_epochs] envelope epochs (at TargetFs) that were analysed.
    %                   NOTE: preprocess_speech_epoch's filter removes the envelope's mean (DC), which is irrelevant
    %                   here because the Welch estimate removes each segment's mean anyway.
    %
    % Sample-count bookkeeping: preprocess_speech_epoch takes the number of OUTPUT samples rather than a target rate and
    % calls resample(x, n_out, n_in), whose FIR length scales with max(n_out/g, n_in/g), g = gcd(n_out, n_in). That is
    % harmless when n_out/n_in reduces exactly to TargetFs/Fs (as for the whole-second epochs used in this repo), but
    % for an arbitrary-length input the ratio may not reduce and resample would try to build a filter with millions of
    % taps. So each epoch here is trimmed (by at most Fs/gcd(Fs,TargetFs) samples) to a length for which the ratio
    % reduces exactly to the rate ratio.
    % ===============================================================================================================
    arguments
        Speech_RawData (1,1) struct
        options.EpochDurationSec (1,1) double {mustBePositive} = 30
        options.OnsetLatency double = []
        options.EnvelopeHighpassHz (1,1) double {mustBePositive} = 0.5
        options.EnvelopeLowpassHz (1,1) double {mustBePositive} = 64
        options.TargetFs (1,1) double {mustBePositive} = 256
        options.SpectrumOptions (1,1) struct = struct()
        options.MakePlot (1,1) logical = true
        options.PlotTitle char = ''
    end

    x = Speech_RawData.Amplitudes;
    fs = Speech_RawData.Fs;
    if abs(fs - round(fs)) > 1e-9
        error('SpeechEnvelopeSpectrumVsPowerLaw:NonIntegerFs', 'Speech_RawData.Fs = %g is not an integer number of Hz.', fs);
    end
    fs = round(fs);
    if size(x, 2) > 1
        warning('SpeechEnvelopeSpectrumVsPowerLaw:MultiChannel', 'Audio has %d channels; averaging them.', size(x, 2));
        x = mean(x, 2);
    end
    x = x(:);
    % Same normalization as PreprocessSpeechEEGForCoherence.m (whole-recording z-score); only sets the PSD's units.
    x = (x - mean(x)) / std(x);
    num_speech_samples = numel(x);

    target_fs = options.TargetFs;
    if abs(target_fs - round(target_fs)) > 1e-9
        error('SpeechEnvelopeSpectrumVsPowerLaw:NonIntegerTargetFs', 'TargetFs must be an integer number of Hz.');
    end
    target_fs = round(target_fs);
    if options.EnvelopeLowpassHz >= target_fs/2
        error('SpeechEnvelopeSpectrumVsPowerLaw:LowpassAboveNyquist', ...
            'EnvelopeLowpassHz (%.1f) must be below the Nyquist frequency of TargetFs (%.1f Hz).', options.EnvelopeLowpassHz, target_fs/2);
    end

    % -------- epoch boundaries (in speech samples), trimmed so the resampling ratio reduces exactly ------------
    g = gcd(fs, target_fs);
    q0 = fs / g;            % input samples per "unit"
    p0 = target_fs / g;     % output samples per unit
    if isempty(options.OnsetLatency)
        if isinf(options.EpochDurationSec)
            onset_sec = 0;
            n_in = floor(num_speech_samples / q0) * q0;
        else
            n_epoch = floor(num_speech_samples / (options.EpochDurationSec * fs));
            if n_epoch < 1
                error('SpeechEnvelopeSpectrumVsPowerLaw:TooShort', 'Recording (%.2f s) is shorter than one epoch (%.2f s).', ...
                    num_speech_samples/fs, options.EpochDurationSec);
            end
            onset_sec = (0:n_epoch-1) * options.EpochDurationSec;
            n_in = floor(round(options.EpochDurationSec * fs) / q0) * q0;
        end
    else
        if isinf(options.EpochDurationSec)
            error('SpeechEnvelopeSpectrumVsPowerLaw:NeedFiniteEpoch', 'OnsetLatency requires a finite EpochDurationSec.');
        end
        onset_sec = options.OnsetLatency(:)';
        n_in = floor(round(options.EpochDurationSec * fs) / q0) * q0;
    end
    n_out = n_in / q0 * p0;
    if n_in < q0
        error('SpeechEnvelopeSpectrumVsPowerLaw:EpochTooShort', 'Epoch is shorter than %d samples; cannot be resampled to TargetFs.', q0);
    end
    Num_epochs = numel(onset_sec);
    fprintf(' -- Speech envelope spectrum: %d epoch(s) of %.3f s (%d samples @ %d Hz -> %d samples @ %d Hz); envelope band %.2f-%.1f Hz\n', ...
        Num_epochs, n_in/fs, n_in, fs, n_out, target_fs, options.EnvelopeHighpassHz, options.EnvelopeLowpassHz);

    SpeechEnvelopes = zeros(n_out, Num_epochs);
    for e = 1:Num_epochs
        i0 = round(onset_sec(e) * fs) + 1;         % 1-based index of the sample at time onset_sec(e)
        i1 = i0 + n_in - 1;
        if i0 < 1 || i1 > num_speech_samples
            error('SpeechEnvelopeSpectrumVsPowerLaw:EpochOutsideRecording', ...
                'Epoch %d (onset %.3f s) needs samples %d:%d but the recording has only %d.', e, onset_sec(e), i0, i1, num_speech_samples);
        end
        % Existing repo code: Hilbert envelope -> band_pass_filt -> resample to n_out samples.
        S = preprocess_speech_epoch(x(i0:i1), fs, options.EnvelopeHighpassHz, options.EnvelopeLowpassHz, n_out);
        if abs(S.fs - target_fs) > 1e-6 * target_fs
            error('SpeechEnvelopeSpectrumVsPowerLaw:UnexpectedRate', 'Envelope came back at %.6f Hz, expected %d Hz.', S.fs, target_fs);
        end
        SpeechEnvelopes(:, e) = S.envelope(:);
    end

    % -------- spectrum vs power law --------------------------------------------------------------------------
    sp_opts = options.SpectrumOptions;
    names = fieldnames(sp_opts);
    args = cell(1, 2 * numel(names));
    for k = 1:numel(names)
        args{2*k - 1} = names{k};
        args{2*k} = sp_opts.(names{k});
    end
    Result = SpectrumVsPowerLaw(SpeechEnvelopes, target_fs, args{:});
    Result.SpeechFs = fs;
    Result.EpochDurationSec = n_in / fs;
    Result.NumEpochs = Num_epochs;
    Result.EpochOnsetSec = onset_sec;
    Result.EnvelopeHighpassHz = options.EnvelopeHighpassHz;
    Result.EnvelopeLowpassHz = options.EnvelopeLowpassHz;
    Result.Figure = [];

    if options.MakePlot
        ttl = options.PlotTitle;
        Result.Figure = PlotSpectrumVsPowerLaw(Result, 'Title', ttl, 'YLabel', 'Envelope PSD (a.u.^2/Hz)');
    end
end
