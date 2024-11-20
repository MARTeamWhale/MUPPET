function [xClipProcessed, j_targetSigEnergyPos, j_noisePos, j_targetSigBoxPos, tClipStart] = getProcessedSNClip(varargin)
%
% Extract a clip from a larger audio timeseries vector that contains a
% signal of interest and associated noise samples preceding the signal. The
% clip is bandpass-filtered and may be downsampled to a target sampling
% rate.
%
% This function can also accept the locations of other, non-target signals
% in the time series as an input argument. These other signals will be
% removed from the noise window, in which case the noise sample will
% consist of a truncated and possibly concatenated noise vector.
%
% NOTES:
% There are several types of relative indices at play here. To keep track
% of which index variables are using which reference values, the index
% variables are prefixed as follows:
% - i_ = relative to the full audio vector "x"
% - j_ = relative to the DOWNSAMPLED clip extracted from "x"
% - ks_ = relative to the signal annotation window (not the cumulative
%       energy limits) within the downsampled clip
%
% SYNTAX:
%   [xClipFilt, j_sigEnergyPos, j_noisePos, tClipStart] = isolateFilteredSNClip(x, fs, targetSigBoxPos, dFilter)
%   [__] = isolateFilteredSNClip(__, Name,Value)
%
% INPUT ARHUMENTS:
%   Required
%   .......................................................................
%   "x" - vector representing an audio time series (must be a column)
%   .......................................................................
%   "fsOriginal" - sampling rate of WAV file, in Hertz
%   .......................................................................
%   "fsResampled" - sampling rate if downsampling is required 
%   .......................................................................
%   "targetSigBoxPos" - 2-element vector representing the start and stop
%       times of a window containing the target signal, in seconds
%   .......................................................................
%   "dFilter" - Bandpass filter object, created using the "designfilt"
%       function from the Signal Processing Toolbox
%   .......................................................................
%
%   Optional (Name-Value Pairs)
%   .......................................................................
%   "EnergyPercent" - Percentage of the total energy within the signal
%       window that will determine the precise start and stop times of the
%       signal within that window. Default is 90%.
%   .......................................................................
%   "NoiseDistance" - Amount of separation that the noise window should
%       have before the energy-based start of the target signal, in
%       seconds. Default is 1 sec.
%   .......................................................................
%   "IdealNoiseSize" - The ideal duration of the noise sample, in seconds.
%       If not specified, then this parameter will be made equal to the
%       energy-based duration of the target signal. The actual noise
%       duration may be shorter than the ideal size if it contains other
%       signals that must be removed.
%   .......................................................................
%   "RemoveFromNoise" - N-by-2 matrix of start and stop times of non-target
%       signals to be removed from the noise window, should they overlap 
%       with it. N is the number of signals. Specified in seconds.
%   .......................................................................
%   "ClipBufferSize" - Amount of buffer samples to include at the beginning
%       and end of the bandpassed sub-clip, to ensure that the signal and
%       noise windows do not contain edge artefacts introduced by the
%       bandpass filter. Specified in seconds. The default is 0.032 secs.
%       NOTE: This parameter may be removed in the future and determined
%       automatically.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "xClipProcessed" - vector of bandpass-filtered and downsampled audio
%       subclip samples containing signal and noise
%   .......................................................................
%   "j_sigEnergyPos" - 2-element vector containing the start and stop
%       samples of the target signal based on cumulative energy threshold 
%       duration
%   .......................................................................
%   "j_noisePos" - N-by-2 matrix containing the start and stop samples
%       of signal-free sections within a noise window preceding the signal 
%   .......................................................................
%   "j_targetSigBoxPos" - 2-element vector containing the start and stop
%       samples of the target signal annotation box
%   .......................................................................
%   "tClipStart" - start time of the subclip, in seconds
%   .......................................................................
%
%
% Written by Wilfried Beslin and Mike Adams
% Last updated by Wilfried Beslin
% 2024-11-20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - [Wilfried] Things I might do:
%   -- [DONE] add noise range and/or duration as output arguments
%   -- [ABANDONED] remove buffer size as an argument and calculate it 
%       automatically
%           --- THIS IS NO LONGER RECOMMENDED, BECAUSE STFT PARAMETERS WILL
%           BE LIMITED BY BUFFER SIZE, AND THUS BUFFER SHOULD BE FLEXIBLE 
%           TO TAKE THIS INTO ACCOUNT
%   -- [DONE] shrink the noise window if there is not enough space in the 
%       time series to isolate a clip with the ideal noise size
%   -- extract noise samples that occur after the signal, either to
%       complement the preceding noise and give more samples to work with,
%       or as an alternative in case there is not enough space to isolate
%       a full noise window ahead of the target signal.
% - [Mike] Additional thing to consider
%   -- exporting 90% energy 
%   -- [DONE] add parameter to set cumulative energy thresehold

    
    % parse input
    p = inputParser();
    p.addRequired('x', @(val)validateattributes(val,{'numeric'},{'column'}))
    p.addRequired('fsOriginal', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
    p.addRequired('fsResampled', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
    p.addRequired('targetSigBoxPos', @(val)validateattributes(val,{'numeric'},{'numel',2}))
    p.addRequired('dFilter', @(val)validateattributes(val,{'digitalFilter'},{'scalar'}))
    p.addParameter('EnergyPercent', 90, @(val)validateattributes(val,{'numeric'},{'scalar','positive','<=',100}))
    p.addParameter('NoiseDistance', 1, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}))
    p.addParameter('IdealNoiseSize', [], @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
    p.addParameter('RemoveFromNoise', [], @(val)validateattributes(val,{'numeric'},{'ncols',2}))
    p.addParameter('ClipBufferSize', 0.032, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}))
    
    p.parse(varargin{:})
    %%% required input
    x = p.Results.x;
    fsOriginal = p.Results.fsOriginal;
    fsResampled = p.Results.fsResampled;
    targetSigBoxPos = p.Results.targetSigBoxPos;
    dFilter = p.Results.dFilter;
    %%% optional input
    enerPerc = p.Results.EnergyPercent;
    noiseDist = p.Results.NoiseDistance;
    idealNoiseSize = p.Results.IdealNoiseSize;
    otherSigBoxPos = p.Results.RemoveFromNoise;
    clipBufferSize = p.Results.ClipBufferSize;
    % end input parsing

    % get clip
    [xClip, tClipStart] = extractClip(x, fsOriginal, targetSigBoxPos, noiseDist, idealNoiseSize, clipBufferSize);

    % apply filter and downsample
    if ~isempty(xClip)
        xClipProcessed = filterDownsampleClip(xClip, fsOriginal, fsResampled, dFilter);
    else
        xClipProcessed = [];
    end

    % isolate final signal and noise windows within clip
    if ~isempty(xClipProcessed)
        [j_targetSigEnergyPos, j_noisePos, j_targetSigBoxPos] = isolateSN(xClipProcessed, fsResampled, tClipStart, targetSigBoxPos, enerPerc, noiseDist, idealNoiseSize, otherSigBoxPos);
    else
        %%% if it's not possible to generate the clip (i.e., because the
        %%% signal is too close to the beginning of the sequence or cannot
        %%% be downsampled), then return empties and NaNs
        j_targetSigEnergyPos = double.empty(0,2);
        j_noisePos = double.empty(0,2);
        j_targetSigBoxPos = double.empty(0,2);
        tClipStart = NaN;
    end
end


%% extractClip ------------------------------------------------------------
function [xClip, tClipStart] = extractClip(x, fs, targetSigBoxPos, noiseDist, idealNoiseSize, clipBufferSize)
% returns a subset of samples from a larger timeseries that contains a
% target signal plus some additional noise preceeding it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % convert all temporal arguments to samples
    secs2samples = @(t) round(t*fs);
    %%% signal box
    i_targetSigBoxPos = secs2samples(targetSigBoxPos) + 1;
    %%% noise distance and buffer
    numNoiseDistSamples = secs2samples(noiseDist);
    numClipBufferSamples = secs2samples(clipBufferSize);
    %%% noise samples - for this one, if user did not specify a duration,
    %%% then use the same number of samples as the target signal.
    useEqualSigNoiseSizes = isempty(idealNoiseSize) || isnan(idealNoiseSize);
    if useEqualSigNoiseSizes
        numIdealNoiseSamples = i_targetSigBoxPos(2) - i_targetSigBoxPos(1) + 1;
    else
        numIdealNoiseSamples = secs2samples(idealNoiseSize);
    end

    % generate shorter clip from audio vector
    %%% first, check if there are enough samples to generate the full clip,
    %%% making adjustments if needed and where possible. 
    i_clipStart = i_targetSigBoxPos(1) - numNoiseDistSamples - numIdealNoiseSamples - numClipBufferSamples;
    i_clipStop = i_targetSigBoxPos(2) + numClipBufferSamples;
    if i_clipStop > numel(x)
        %%% If i_clipStop is out of bounds, cancel completely
        goodClip = false;
    elseif i_clipStart < 1
        %%% If i_clipStart is out of bounds, try shrinking the noise window
        numAttemptedNoiseSamples = numIdealNoiseSamples + i_clipStart - 1;
        i_clipStart = 1;
        goodClip = numAttemptedNoiseSamples > 0; %** This would be a good place to use a minimum noise size parameter
    else
        %%% bounds of clip are within range
        %numAttemptedNoiseSamples = numIdealNoiseSamples;
        goodClip = true;
    end

    % set output (return empties if clip could not be extracted)
    if goodClip
        xClip = x(i_clipStart:i_clipStop);
        tClipStart = (i_clipStart-1)/fs;
    else
        xClip = [];
        tClipStart = [];
    end
end


%% filterDownsampleClip ---------------------------------------------------
function xClipProcessed = filterDownsampleClip(xClip, fsOriginal, fsResampled, dFilter)
% Applies a digital filter and downsamples an audio clip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MUPPET.noDelayFilt

    % start by making sure that the clip can be processed at all.
    %%% To allow downsampling, these rules must be followed:
    %%% - the resampling rate must be less than the original rate
    %%% - the downsampling factor must be an integer
    %%% - the filter's upper stopband frequency must be less than half the
    %%%   resampling rate
    do_downsampling = ~isempty(fsResampled) && fsResampled ~= fsOriginal;
    if do_downsampling
        dsFactor = fsOriginal/fsResampled;
        if dsFactor < 1
            dsErrMsg = sprintf('Upsampling from %d Hz to %d Hz is not supported.', fsOriginal, fsResampled);
        elseif dsFactor ~= round(dsFactor)
            dsErrMsg = sprintf('Cannot downsample evenly from %d Hz to %d Hz.', fsOriginal, fsResampled);
        elseif dFilter.StopbandFrequency2 > (fsResampled/2)
            dsErrMsg = sprintf('The filter has an upper stopband frequency of %d Hz; this is too high for downsampling to %d Hz.', dFilter.StopbandFrequency2, fsResampled);
        else
            dsErrMsg = '';
        end

        if ~isempty(dsErrMsg)
            warning(dsErrMsg)
            xClipProcessed = [];
            return
        end
    end

    % apply digital filter
    xClipFilt = noDelayFilt(dFilter, xClip);

    % downsample, if requested
    if do_downsampling
        xClipProcessed = downsample(xClipFilt, dsFactor);
    else
        xClipProcessed = xClipFilt;
    end
end


%% isolateSN --------------------------------------------------------------
function [j_targetSigEnergyPos, j_noisePos, j_targetSigBoxPos] = isolateSN(xClipProcessed, fs, tClipStart, targetSigBoxPos, enerPerc, noiseDist, idealNoiseSize, otherSigBoxPos)
% Returns sample indices (relative to clip) for energy-based signal limits
% and noise window(s).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MUPPET.calcEng

    % convert temporal arguments of signal positions to samples (noise
    % noise samples will be done later, once we have the signal energy
    % limits)
    secs2samples = @(t) round(t*fs);
    j_targetSigBoxPos = secs2samples(targetSigBoxPos - tClipStart) + 1;
    j_otherSigBoxPos = secs2samples(otherSigBoxPos - tClipStart) + 1;

    % get number of buffer samples
    numClipSamples = numel(xClipProcessed);
    numClipBufferSamples = numClipSamples - j_targetSigBoxPos(2);

    % isolate specific signal bounds based on a percentage of 
    % cumulative energy
    xSigInitial = xClipProcessed(j_targetSigBoxPos(1):j_targetSigBoxPos(2));
    [ks_sigEnergyStart, ks_sigEnergyStop] = calcEng(xSigInitial,enerPerc);

    % get signal position relative to clip
    j_targetSigEnergyPos = j_targetSigBoxPos(1) + [ks_sigEnergyStart,ks_sigEnergyStop] - 1;

    % isolate noise samples, this time based on signal energy bounds rather
    % than the annotation box
    %%% start by converting temporal positions to samples
    numNoiseDistSamples = secs2samples(noiseDist);
    useEqualSigNoiseSizes = isempty(idealNoiseSize) || isnan(idealNoiseSize);
    if useEqualSigNoiseSizes
        numIdealNoiseSamples = ks_sigEnergyStop - ks_sigEnergyStart + 1;
    else
        numIdealNoiseSamples = secs2samples(idealNoiseSize);
    end
    %%% get initial noise window, adjusting if needed
    j_noisePosInitial = j_targetSigEnergyPos(1) - numNoiseDistSamples - [numIdealNoiseSamples, 1];
    if j_noisePosInitial(1) < numClipBufferSamples
        bufferCreep = numClipBufferSamples - j_noisePosInitial(1) + 1;
        j_noisePosInitial(1) = j_noisePosInitial(1) + bufferCreep;
    end

    % get indices of clean (i.e., signal-free) noise sections.
    %%% This uses a loopless approach with logical indices and
    %%% cumulative sums.
    
    %%% start by finding where noise occurs
    clipIdcs = 1:numClipSamples;
    isInNoiseWin = clipIdcs >= j_noisePosInitial(1) & clipIdcs <= j_noisePosInitial(end);

    %%% find where non-target signals ("contaminants") occur.
    %%% Adjust the j-indices of the other signals if they are out of
    %%% bounds, or this won't work.
    j_otherSigBoxPosClipped = j_otherSigBoxPos(~all(j_otherSigBoxPos < 1, 2) & ~all(j_otherSigBoxPos > numClipSamples, 2), :);
    j_otherSigBoxPosClipped(j_otherSigBoxPosClipped < 1) = 1;
    j_otherSigBoxPosClipped(j_otherSigBoxPosClipped > numClipSamples) = numClipSamples;
    numContaminants = size(j_otherSigBoxPosClipped,1);
    contaminantPosAnalVec = cumsum(ismember(clipIdcs, j_otherSigBoxPosClipped+repmat([0,1],numContaminants,1))); % every odd number in this sequence indicates the presence of a signal
    isContaminated = rem(contaminantPosAnalVec, 2) ~= 0;

    %%% find sections of noise that are free of signal and create a matrix
    %%% of start/stop positions for each section
    isCleanNoise = isInNoiseWin & ~isContaminated;
    j_noisePos = [...
        find(isCleanNoise-[0,isCleanNoise(1:end-1)] > 0)',...
        find([isCleanNoise(2:end),0]-isCleanNoise < 0)'...
        ];
    
    %** DEBUG PLOT
    %{
    if size(j_otherSigBoxPosClipped, 1) > 0
        tClip = (((1:numel(xClipProcessed))-1)./fs)';
        figure;
        ax = axes();
        ax.NextPlot = 'add';
        plot(ax, tClip, xClipProcessed)
        plot(ax, tClip(j_targetSigBoxPos(1):j_targetSigBoxPos(2))', xClipProcessed(j_targetSigBoxPos(1):j_targetSigBoxPos(2)));
        plot(ax, tClip(j_targetSigEnergyPos(1):j_targetSigEnergyPos(2))', xClipProcessed(j_targetSigEnergyPos(1):j_targetSigEnergyPos(2)))
        plot(ax, tClip(isInNoiseWin)', xClipProcessed(isInNoiseWin));
        if sum(isContaminated) > 0
            stem(ax, tClip(isContaminated)', xClipProcessed(isContaminated));
        else
            stem(ax, NaN, NaN);
        end
        stem(ax, tClip(isCleanNoise)', xClipProcessed(isCleanNoise));
        legend(ax, {'Processed Clip', 'Signal Box', 'Energy Signal', 'Noise', 'Contaminating Signal', 'Signal-Free Noise'})
        grid(ax, 'on')
        box(ax, 'on')
        keyboard
    end
    %}
end