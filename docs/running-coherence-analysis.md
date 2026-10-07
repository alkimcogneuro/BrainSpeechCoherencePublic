# Running the Speech–EEG Coherence Analysis

This page covers two things:

1. How to run `dev/TestScripts/script_runCoherenceAnalysis.m`, which computes brain–speech
   coherence (Welch-segmented cross-spectral density → band-averaged magnitude-squared
   coherence, with an optional permutation/shuffle significance test) for a list of subjects.
2. How to plot the results once a run has finished.

It assumes you already have EEG data loaded into the `EEG_struct` format this repo uses
(see `dev/ACSEEG_read_data.m` / `dev/load_data_acseeg.m` if you need to build one) and the
corresponding speech audio files as `.wav` files.

## 1. Before you run it

Open `dev/TestScripts/script_runCoherenceAnalysis.m` in MATLAB. Near the top, under the banner

```matlab
% USER SHOULD SET THESE PARAMETERS BEFORE RUNNING THIS SCRIPT
```

edit the following variables directly in the file. The defaults checked into the repo are
development values — they will very likely not be the right values for your data or your
machine's file paths.

| Parameter | What it controls |
|---|---|
| `speech_files_path` | Folder containing the speech `.wav` files. Each subject's `EEG_struct.audio_file` must name a file in this folder. |
| `highpass_cutoff`, `lowpass_cutoff` | Hz. The broadband filter applied to both EEG and speech before coherence is computed. |
| `band_low`, `band_high` | Hz. The narrower frequency band that the per-bin coherence is averaged over to produce one coherence value per channel (e.g. 3–7 Hz). |
| `WindowDurationSec`, `OverlapFraction` | Welch segmentation: each trial is split into overlapping windows of this duration (seconds) and this overlap fraction (0.5 = 50%), and coherence is averaged across those windows before being averaged across trials. |
| `NumShuffles` | Number of permutations for the shuffle/significance test. **Set this to a small number (e.g. 2) for your first test run** — a real run (`NumShuffles = 1000`) can take a long time per subject. Set to `0` to skip the shuffle test entirely and only compute the raw coherence (`True_MSC`); this is much faster and useful for a first look at the data. |
| `RandomSeed` | Seed for the shuffle test. `NaN` means "don't seed" (a different random seed every run). |
| `AvgMethod` | How per-bin coherence is averaged into one band value: `'power-weighted'` (default) or `'simple'`. |
| `output_dir` | Folder where this run's results are saved (see below). |

Below that, set `eeg_data_filenames` to the list of `EEG_struct` `.mat` files you want to
analyze — one file per subject. The repo ships this as a short list of a few subjects with a
longer list commented out underneath (wrapped in `%{ ... %}`); uncomment and edit that block,
or replace the whole list, to analyze a different set of subjects. **Test on 1–2 subjects with
`NumShuffles = 2` first** to confirm everything is configured correctly before scaling up to
the full subject list and `NumShuffles = 1000` — the full analysis can take well over an hour
per subject.

## 2. Running it

Make sure `dev/`, `dev/TestScripts/`, and `dev/PlottingFunctions/` are on your MATLAB path
(e.g. `addpath(genpath('path/to/BrainSpeechCoherence/dev'))`), then either open the script and
press **Run**, or from the Command Window:

```matlab
run('script_runCoherenceAnalysis.m')
```

For each subject, the script prints progress to the Command Window and either saves that
subject's results or prints an error and moves on to the next subject — one subject failing
(e.g. a missing data file) does not stop the rest of the run.

## 3. What gets written to `output_dir`

- **If `output_dir` does not already exist**, it is created.
- **If `output_dir` already exists** (for example, you forgot to change it since a previous
  run), the existing folder is renamed aside — a timestamp is appended to its name, e.g.
  `Coherence_bandMSC_control_analyses_backup_20261007_143533` — rather than being written into
  or overwritten. The script then creates a fresh, empty `output_dir` for this run. This means
  you can never silently lose a previous run's results by re-running with the same
  `output_dir`.
- Inside `output_dir`, each subject that completes successfully gets its own file,
  `coherence_control_analysis_S<id>.mat` (e.g. `coherence_control_analysis_S01.mat`). Each file
  contains:
  - `Coherence_Analysis_ShuffleTestResults` — that subject's full results struct.
  - `Subj_id`, `eeg_data_file` — which subject and input file this came from.
  - `RunParameters` — every parameter used for this run, so the file is self-describing even if
    you don't have the log below handy.
- `coherence_control_analysis_parameters.txt` — a human-readable log of every parameter for
  this run plus one outcome line per subject (saved / errored). This file is appended to, so if
  you look back at a folder you no longer remember the settings for, this is the place to check.

A subject that errors gets **no** `.mat` file (so you never end up with a misleading empty or
partial result) — the error is printed to the Command Window and recorded in the log.

## 4. Reading the summary table

After all subjects finish, the script prints a summary table — one row per subject — with
columns for window duration, overlap, segments per trial, and the min/max coherence across
channels (plus a minimum p-value and significant-channel count if `NumShuffles >= 1`; those
columns are omitted entirely when `NumShuffles == 0`, since there's nothing shuffle-test-derived
to show). A few things are flagged automatically if they look wrong:

- **`NumSegmentsPerTrial` of 1 or less** — the Welch segmentation isn't actually doing anything
  for that subject's trial length; check `WindowDurationSec` against how long that subject's
  trials actually are.
- **Non-finite values in `True_MSC`**, or **values outside `[0, 1]`** — coherence should always
  land in this range; either of these means something upstream is wrong and that subject's
  results shouldn't be trusted yet.

## 5. Plotting the results of a completed run

Once a run has finished, reload its per-subject files and plot them with
`dev/LoadCoherenceResultsFromFolder.m` and `dev/PlottingFunctions/PlotSubjectsAndGroupTopography.m`.
This works even in a brand-new MATLAB session — you just need to know (or look up, in that
run's `coherence_control_analysis_parameters.txt`) the `output_dir` that run used.

```matlab
% Reload every subject's saved results from one run:
ResultsArray = LoadCoherenceResultsFromFolder(output_dir);

% Plot the raw coherence (always present, whether or not the shuffle test ran):
PlotSubjectsAndGroupTopography(ResultsArray, 'True_MSC', ...
    OutputDir=fullfile(output_dir, 'plots'), ...
    SubjectIDField='Subj_id', ...
    ColorbarLabel='True MSC');
```

This produces one combined figure (every subject's topography plus the group-average
topography, all on one shared color scale) and saves an individual PNG per subject plus one
for the group average into `OutputDir`.

If that run used `NumShuffles >= 1`, you can plot the shuffle-test fields the same way:

```matlab
PlotSubjectsAndGroupTopography(ResultsArray, 'PValues', ...
    OutputDir=fullfile(output_dir, 'plots'), ...
    SubjectIDField='Subj_id', ...
    ColorbarLabel='p-value');
```

Trying to plot `'PValues'` or `'ZScores'` from a run that used `NumShuffles == 0` will raise a
clear "field does not exist" error — those fields are only computed when the shuffle test
actually ran.

To reload just one subject's file by hand instead of the whole folder:

```matlab
S = load(fullfile(output_dir, 'coherence_control_analysis_S01.mat'));
S.Coherence_Analysis_ShuffleTestResults   % that subject's results struct
```

### Other useful options for `PlotSubjectsAndGroupTopography`

- `OutputDir` is required; everything else is optional.
- `GroupLabel` / `FigureTitle` — titles for the group panel and the overall combined figure.
- `Clim` — `[min max]` color limits shared by every panel. Left empty, this is computed
  automatically from the data.
- `HighlightThreshold` — a scalar; values above it are rendered in a visually distinct
  "pops out" color scheme instead of a smooth gradient. Useful for a field like `ZScores`,
  e.g. `HighlightThreshold=2.0`.
- `Chanlocs` — only needed if your results structs don't already carry a `.Chanlocs` field.

See the header comments in `PlotSubjectsAndGroupTopography.m` for the full list of options.
