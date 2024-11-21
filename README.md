# Mysticete Utterance Processing and Parameter Extraction Tool (MUPPET)
## Work in progress  
Developed by Wilfried Beslin and Mike Adams  
Maritimes Team Whale
-------------------
Goal: A tool for use in cetacean research


-------------------

## Introduction
This tool is for use in calculating the spectral characteristics for cetacean vocalizations obtained using JASCO's PAMLAB annotation software.
Initial development was focused on Blue Whale audible calls. The tool currently exists as a helper script<sup>*</sup> with several underlying functions. The helper script (`MUPPET`) imports the output of JASCO's PAMLAB annotation files, extracts the needed inputs, and then matches the annotated calls to and imports the appropriate .wav files. These inputs are passed to the underlying functions to extract the data within the .wav files that corresponds to the annotated call and a sample of noise taken some time before the call. These clippings are bandpassed to the frequencies of interest using a Kaiser window-based FIR filter. After bandpass filtering, precise start and end times of the signals are determined based on a user-specified percentage of the energy within the signal window. Noise samples are isolated relative to the energy-based start time of the signal. Ideal duration of the noise window may either be set by the user to some common value (e.g., 10 secs), or be made equal to the signal duration. If a noise window happens to includes parts of other signals that were annotated in PAMLAB, then those parts will be eliminated from the noise window (thereby shortening the noise window). The clipped and bandpassed call and noise samples are then used to calculate both time and spectral parameters. 

<sup>*</sup>_The helper script is actually implemented as a function that can optionally accept input arguments for greater flexibility. This will be discussed further in the section "Running the Tool" below_.

## Set up

### Initial Start Up
1) Either clone or download the .zip and unzip the MUPPET repository to your local machine.
2) Add this new directory with all of its subfolders to your MATLAB path.
### Requirements

The tool requires inputs to calculate the parameters. These include:
  -  SNR_PARAMS.csv: a parameter file which contains the filtering and noise presets for each species' call type. This file has values for:
      - **Species** - The species of interest (e.g. Blue Whale)
      - **Call_Type** - The call type (e.g. Audible, Tonal, etc...) 
      - **Lower_Passband_Frequency** - Lower bound of the bandpass filter before which frequencies become attenuated (Hz). This value should correspond to the lowest frequency of interest.
      - **Upper_Passband_Frequency** - Upper bound of the bandpass filter before which frequencies become attenuated (Hz). This value should correspond to the highest frequency of interest.
      - **Stopband_Rolloff_Bandwidth** - The amount of frequency bandwidth that it takes for the passband frequencies to be attenuated to a level of 60 dB (Hz). The closer this value is to zero, the sharper the passband frequency cutoffs will be, but at the expense of increased filter order (and thus processing time/memory usage). Note that filter order also increases with sampling rate.
      - **Noise_Distance** - Value to determine how far before the ***signal*** the ***noise*** sample will be taken
      - **Ideal_Noise_Duration** - The target duration that the noise window should have (seconds). This parameter can also be left empty, in which case the target noise duration will be equal to the signal duration. Note that the actual noise duration may be shorter than ideal, if the noise window contains other signals that must be removed to avoid contaminated noise estimates.
      - **Signal_Energy_Percent** - The percentage of energy within a signal window that determines the start and stop times of a signal based on cumulative energy. For example, a value of 90 will mark the start time at 5% of the cumulative energy in the window, and the stop time at 95%.
- A directory containing the .csv annotation output generated by JASCO's PAMLAB. The inputs used from these files are:
  - **Filename** - used to identify the original .wav file containing the call
  - **Relative Start Time** - The manually selected start time in seconds of the call artefact relative to the start of the .wav file
- A directory containing all the .wav files for which PAMLAB annotation .csv files exist. *Note: This directory can contain additional .wav files, but **must** contain all .wav files for which PAMLAB annotations exist*   

## Running the Tool
### Basic Usage
To use the tool, open MATLAB and run the file _MUPPET.m_. The most basic way to run this file is by pressing the "Run" button in the MATLAB editor, or by typing `MUPPET` in the Command Window. This will load the parameters in the file _SNR_PARAMS.csv_, and you will be prompted to set the input and output file paths.

When running, the tool will also prompt the user to specify which species and call type to analyze. This dictates which row of the parameter file will be read. ***For each run of the tool, all annotations in a CSV file will be processed using the species and call type that were specified by the user at runtime using the PARAM file; the tool does not interpret species or call type information within the annotation files themselves.***

### Input Arguments
It is possible to pass certain input arguments when calling _MUPPET_ via the command window, and avoid having to set them manually or use defaults. The supported arguments are:

  - **PAMLAB_DATA_FOLDER** - Path to a folder with PAMlab output. If not specified, the tool will prompt the user to select the path manually.
  - **WAV_FILE_FOLDER** - Path to the folder containing WAV files for the dataset of interest. If not specified, the tool will prompt the user to select the path manually.
  - **OUTPUT_FOLDER_LOCATION** - Path where the tool's output folder will be created. If not specified, user will be prompted to select the path manually, or simply use the parent folder of the PAMlab output by cancelling the prompt.
  - **PARAMFILE** - Path to a CSV file of project parameters, as an alternative to the default _SNR_PARAMS.csv_ file. The specified file must still have the same format as _SNR_PARAMS.csv_.
  - **PLOT_TRACE_LINES** - true/false value that determines whether or not to save images of trace line plots for each annotation in the output folder. Default is false.

Input arguments are set by using MATLAB's Name-Value pair syntax. For example:
```matlab
MUPPET('PAMLAB_DATA_FOLDER',data_dir, 'WAV_FILE_FOLDER',wav_dir)
```
where in this case, `data_dir` and `wav_dir` are char string variables specifying the paths to the PAMLAB data folder and WAV file folder, respectively (e.g., `data_dir = 'C:\Users\Me\SNR_Tool_analysis\PAMLAB'`).

### Output Arguments
_MUPPET_ will always save its output as CSV files and, if requested, PNG files within a folder called _MUPPET_OUTPUT_INPUT_. However, it is also possible to pass the output into the MATLAB workspace by requesting it as output variables when running the tool. The syntax for this is:
```matlab
out = MUPPET
```

where `out` is the variable that will contain the output (it does not have to be called _out_ necessarily; give it any name you want). The output variable comes in the form of a MATLAB struct containing tables of PAMLAB annotations with the SNR data appended. Each field of the struct corresponds to one PAMLAB CSV file that was processed.


## Output Files
As the tool processes each PAMLAB annotation CSV file it finds, it will, by default, create two CSV files. Those files will be saved in the folder _MUPPET_OUTPUT_INPUT_, whose location is set by the user.

**SNR Output CSV contents:**
  - Copy of columns in input PAMlab file.
  - **SNR_Direct** - The "raw", uncorrected SNR value in dB re. noise power (see "SNR Calculation")
  - **SNR_Corrected** - The SNR value (in dB re. noise power) corrected for noise within the signal window (see "SNR Calculation")
  - **SNRCalc_SignalDuration** - Duration of the signal window used in the SNR calculations, as determined based on a user-specified percentage of the total energy within the signal annotation box, in seconds.
  - **SNRCalc_NoiseDuration** - Actual duration of the noise clip that was used to calculate SNR. This will always be equal to or less than the ideal noise duration, depending on whether the noise window contained other signals that needed to be removed or not, or if there were not enough samples in the time series to generate the ideal noise clip. ***Be wary of noise durations that are extremely small (e.g., < 1 sec). In such cases, the noise power estimate may not be an accurate representation of the noise power coinciding with the signal, which in turn means that the SNR estimates may be inaccurate.***

**CallParams Output CSV contents:**

| Base Parameter                | MUPPET Variable Name                 | Technical Definition (specific to the interpretation)                                                                                                                                                  | Obtained from                     |
| ----------------------------- | ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------- |
| Start time (s)                | StartTime_Waveform                   | Time (s) at which significant cumulative energy begins within a user-set interval containing a call, where "significance" is also determined by a user-set threshold (e.g., top 90% energy)            | Waveform (i.e., time domain)      |
| StartTime_Trace               | Time(s) at which a trace line begins | Spectrogram (trace line)                                                                                                                                                                               |
| End time (s)                  | EndTime_Waveform                     | Time (s) at which significant cumulative energy ends within a user-set interval containing a call, where "significance" is also determined by a user-set threshold (e.g., top 90% energy)              | Waveform (i.e., time domain)      |
| EndTime_Trace                 | Time(s) at which a trace line ends   | Spectrogram (trace line)                                                                                                                                                                               |
| Total duration (s)            | Duration_Waveform                    | Waveform end time (s) minus waveform start time (s)                                                                                                                                                    | Waveform (i.e., time domain)      |
|                               | Duration_Trace                       | Trace line end time (s) minus trace line start time (s)                                                                                                                                                | Spectrogram (trace line)          |
| Minimum frequency (Hz)        | MinFrequency_Trace                   | Lowest frequency value (Hz) that occurs in a trace line                                                                                                                                                | Spectrogram (trace line)          |
|                               | MinFrequency_PSD                     | Frequency (Hz) at which significant cumulative energy begins within the (truncated) power spectral density estimate, where "significance" is determined by a user-set threshold (e.g., top 90% energy) | Spectrum (i.e., frequency domain) |
| Maximum frequency (Hz)        | MaxFrequency_Trace                   | Highest frequency value (Hz) that occurs in a trace line                                                                                                                                               | Spectrogram (trace line)          |
|                               | MaxFrequency_PSD                     | Frequency (Hz) at which significant cumulative energy ends within the (truncated) power spectral density estimate, where "significance" is determined by a user-set threshold (e.g., top 90% energy)   | Spectrum (i.e., frequency domain) |
| Bandwidth (Hz)                | Bandwidth_Trace                      | Trace maximum frequency (Hz) minus trace minimum frequency (Hz)                                                                                                                                        | Spectrogram (trace line)          |
|                               | Bandwidth_PSD                        | PSD maximum frequency (Hz) minus PSD minimum frequency (Hz)                                                                                                                                            | Spectrum (i.e., frequency domain) |
| Time (s) of minimum frequency | TraceMinFreqTime                     | Time (s) corresponding to the point in the trace line that has the lowest frequency                                                                                                                    | Spectrogram (trace line)          |
| Time (s) of maximum frequency | TraceMaxFreqTime                     | Time (s) corresponding to the point in the trace line that has the highest frequency                                                                                                                   | Spectrogram (trace line)          |
| Peak frequency (Hz)           | PeakFrequency_Trace                  | Frequency (Hz) corresponding to the point in the trace line that has the highest magnitude                                                                                                             | Spectrogram (trace line)          |
|                               | PeakFrequency_PSD                    | Frequency (Hz) corresponding to the point in the power spectral density estimate that has the highest magnitude                                                                                        | Spectrum (i.e., frequency domain) |
| Center frequency (Hz)         | PSDCentreFrequency                   | Frequency (Hz) that divides the spectrum into two halves of equal energy                                                                                                                               | Spectrum (i.e., frequency domain) |
| Start frequency (Hz)          | TraceStartFrequency                  | Frequency (Hz) of the first point in the trace line                                                                                                                                                    | Spectrogram (trace line)          |
| End frequency (Hz)            | TraceEndFrequency                    | Frequency (Hz) of the last point in the trace line                                                                                                                                                     | Spectrogram (trace line)          |
| Call ID                       | Call_ID                              | Identification key to map calls across output results                                                                                                                                                  |                                   |

## Functions
MUPPET employs a toolbox of custom functions to calculate the parameters of mysticete calls. 

  - **isolateFilteredSNClip** - Isolate signal and associated noise samples from a larger audio time series vector, given a pre-determined signal location. Works by isolating a subsection of the audio vector and applyinging a bandpass filter to it before extracting signal and noise samples.
  - **calcEng** - Calculate the 90% energy start and stop position of a input signal.
  - **noDelayFilt** - Filter a signal using a digitalFilter object from the Signal Processing Toolbox, compensating for group delay introduced by the filter. This function only works if the delay is not frequency-dependent (usually the case with FIR filters). It will generally NOT work with IIR filters like the Butterworth filter.
  - **calculateSNR** - Calculate signal-to-noise ratio, given pre-isolated windows of signal and noise.
  - **computeSTFT** - Calculates the spectrogram of a given signal via the short-time Fourier transform (STFT). Also returns the Welch power spectral density (PSD) estimate, which is essentially just the average of STFT bins.
  - **processSpec** - Processes a signal spectrogram in a particular order requested by the user. Processing operations are "Smooth", "Log", and "Denoise". Any combination of the above may be specified, and some may be omitted (but at least one should be specified).
  - **smoothSpec** - Function for smoothing a spectrogram. Smoothing is done by convolving the spectrogram with a Gaussian kernel (based on Baumgartner and Mussoline). Note that special corrections must be applied to take away edge effects.
  - **getTraceLine** -  Find an ideal trace line through a call in a spectrogram by using a "shortest path" searching algorithm (Dijkstra's algorithm) followed by energy-based clipping.
  - **plotTraceLine** - Plots the trace line of a call against a spectrogram.
  - **extractCallParams** - Calculate several time, frequency, and time-frequency parameters for baleen whale calls.


## Parameter Calculation

## Tracelines

## SNR Calculation
The tool calculates signal-to-noise ratios in Decibels re. noise power, based on the following equations:
$$SNR = 10\log_{10}\left(\frac{P_{signal}}{P_{noise}}\right)$$
where $P_{signal}$ and $P_{noise}$ are the average power of signal and noise, respectively. The average power $P$ for a sampled time series $x$ can be calculated as:
$$P = \frac{1}{N}\displaystyle\sum_{i=1}^{N} x(i)^{2}$$
where $N$ is the total number of samples in $x$, and $i$ is the sample number.

### Direct vs. Corrected SNR
A true signal-to-noise ratio compares the average power of a _pure signal_ to that of noise. However, in virtually all marine mammal PAM analyses, calls are extracted from noisy time series and thus actually consist of _signal + noise mixtures_ rather than pure signals. To account for this, the tool returns a _Corrected SNR_ value in addition to the direct (uncorrected) value.

Given a noisy time series $x$ containing a signal of interest (i.e., a detected whale call), let $P_{xs}$ be the average power of the region in $x$ containing the signal, and $P_{xn}$ be the average power of a region in $x$ with no signal (i.e, noise). Thus, $P_{xs}$ is the power of a signal + noise mixture, and $P_{xn}$ is an estimate of noise power only. The direct SNR is calculated as:
$$SNR_{Direct} = 10\log_{10}\left(\frac{P_{xs}}{P_{xn}}\right)$$
However, since $P_{xs}$ corresponds to a signal + noise mixture, $SNR_{Direct}$ is not an accurate measure of the true SNR. To better approximate the true SNR, the corrected value is calculated by subtracting the estimated noise power from the power of the measured signal:
$$SNR_{Corrected} = 10\log_{10}\left(\frac{P_{xs} - P_{xn}}{P_{xn}}\right)$$
The table below shows how the expected SNR value differs between the direct vs. corrected calculations, based on how the energy of the true (pure) signal compares to that of the noise.

| True Signal vs. Noise Energy | Expected SNR Value<br>(Direct) | Expected SNR Value<br>(Corrected) |
| :--------------------------- | :----------------------------: | :-------------------------------: |
| *Signal Only; Noise Absent*  | +Infinity | +Infinity |
| *Signal > Noise*             | Positive  | Positive  |
| *Signal = Noise*             | Positive  | 0         |
| *Signal < Noise*             | Positive  | Negative  |
| *Signal Absent; Noise Only*  | 0         | -Infinity |

### Accounting for Inaccurate Noise Power Estimates
The average power of noise during a signal of interest can only be estimated. The SNR tool performs this estimate by taking a small clip of noise preceeding the signal. This assumes that the clip is representative of the noise coinciding with the signal. While this is a reasonable assumption in most cases, it may not always hold. There can be situations where the noise within the noise clip is contaminanted by other transient sounds, or is simply not representative of the noise level during the signal.

The SNR tool attempts to reduce poor noise power estimates caused by other signals by removing any annotated signals that may exist within the noise clip. However, the key here is that those signals must have been annotated by an analyst and are included in the PAMLAB CSV files. If there are transient sounds that were not annotated, then those will still be included in the noise clips and result in inflated estimates of noise power. Analysts should also be mindful of this signal removal feature when deciding where to set the start and stop bounds of their annotation boxes in PAMLAB.

Inflated noise power estimates can result in situations where the noise power appears greater than the power of signal + noise, which is of course not possible in reality and would produce anomalous SNR values (i.e., negative values for the direct SNR, and complex numbers for the corrected SNR). To prevent this from happening, the tool will limit noise power estimates such that any value greater than the signal + noise power will be capped at the signal + noise power level. Thus, any situation where the estimated noise power is excessively large will produce values of 0 for the direct SNR, and `-Inf` for the corrected SNR. 

If there are not enough samples available to extract a noise clip for a given signal (i.e., because the signal is too close to the beginning of the time series), then the tool will return `NaN` for both the direct and corrected SNR values for that signal.
