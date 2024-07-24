function [xClipFilt, j_targetSigEnergyPos, j_noisePos, tClipStart] = isolateFilteredSNClip(varargin)
%
% Isolate signal and associated noise samples from a larger audio time
% series vector, given a pre-determined signal location. Works by isolating
% a subsection of the audio vector and applyinging a bandpass filter to it 
% before extracting signal and noise samples.
%
% This function can also accept the locations of other, non-target signals
% in the time series as an input argument. These other signals will be
% removed from the noise window, in which case the noise sample will
% consist of a truncated and possibly concatenated noise vector.
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
%   "fs" - sampling rate, in Hertz
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
%   "xClipFilt" - vector of bandpass-filtered audio subclip samples
%       containing signal and noise
%   .......................................................................
%   "j_sigEnergyPos" - 2-element vector containing the start and stop
%       samples of the target signal based on cumulative energy threshold 
%       duration
%   .......................................................................
%   "j_noisePos" - N-by-2 matrix containing the start and stop samples
%       of signal-free sections within a noise window preceding the signal 
%   .......................................................................
%   "tClipStart" - start time of the subclip, in seconds
%   .......................................................................
%
%
% Written by Wilfried Beslin
% Last updated by Wilfried Beslin
% 2024-07-24
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

    
    import MUPPET.noDelayFilt
    import MUPPET.calcEng
    
    % parse input
    p = inputParser();
    p.addRequired('x', @(val)validateattributes(val,{'numeric'},{'column'}))
    p.addRequired('fs', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
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
    fs = p.Results.fs;
    targetSigBoxPos = p.Results.targetSigBoxPos;
    dFilter = p.Results.dFilter;
    %%% optional input
    enerPerc = p.Results.EnergyPercent;
    noiseDist = p.Results.NoiseDistance;
    idealNoiseSize = p.Results.IdealNoiseSize;
    otherSigBoxPos = p.Results.RemoveFromNoise;
    clipBufferSize = p.Results.ClipBufferSize;
    % end input parsing
    
    % convert all temporal arguments to samples
    %%% ...................................................................
    %%% NOTE: There will be several types of relative indices at play here,
    %%% so to keep track of which index variables are using which reference
    %%% values, I will prefix them as follows:
    %%% i_ = relative to the full audio vector "x"
    %%% j_ = relative to the clip that will be isolated from "x"
    %%% kn_ = relative to the noise sample isolated from the clip
    %%% ks_ = relative to the signal sample isolated from the clip 
    %%%     (from annotation box, not the cumulative energy limits)
    %%% ...................................................................
    secs2samples = @(t) round(t*fs);
    %%% signal boxes
    i_targetSigBoxPos = secs2samples(targetSigBoxPos);
    i_otherSigBoxPos = secs2samples(otherSigBoxPos);
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
    
    % generate shorter clip from audio vector (easier to filter)
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
        goodClip = numAttemptedNoiseSamples > 0;
    else
        %%% bounds of clip are within range
        numAttemptedNoiseSamples = numIdealNoiseSamples;
        goodClip = true;
    end
    
    if goodClip
        xClip = x(i_clipStart:i_clipStop);

        % get signal box start/stop samples relative to clip
        %%% for both target and non-target signals
        j_targetSigBoxPos = i_targetSigBoxPos - i_clipStart + 1;
        j_otherSigBoxPos = i_otherSigBoxPos - i_clipStart + 1;

        % apply digital filter
        xClipFilt = noDelayFilt(dFilter, xClip);

        % isolate specific signal bounds based on a percentage of 
        % cumulative energy
        xSigInitial = xClipFilt(j_targetSigBoxPos(1):j_targetSigBoxPos(2));
        [ks_sigEnergyStart, ks_sigEnergyStop] = calcEng(xSigInitial,enerPerc);
        %xSignal = xSigInitial(ks_sigEnergyStart:ks_sigEnergyStop);
        
        % update number of attempted noise samples if using equal signal
        % and noise sizes
        if useEqualSigNoiseSizes
            numAttemptedNoiseSamples = min([ks_sigEnergyStop - ks_sigEnergyStart + 1, numAttemptedNoiseSamples]);
        end
       
        % get signal and initial noise position relative to clip
        j_targetSigEnergyPos = j_targetSigBoxPos(1) + [ks_sigEnergyStart,ks_sigEnergyStop] - 1;
        j_noisePosInitial = j_targetSigEnergyPos(1) - numNoiseDistSamples - [numAttemptedNoiseSamples, 1];
        
        % get indices of clean (i.e., signal-free) noise sections.
        %%% This uses a loopless approach with logical indices and
        %%% cumulative sums.
        
        %%% start by finding where noise occurs
        numClipSamples = numel(xClipFilt);
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
        isContaminant = rem(contaminantPosAnalVec, 2) ~= 0;
        
        %%% find sections of noise that are free of signal and create a
        %%% matrix of start/stop positions for each section
        isCleanNoise = isInNoiseWin & ~isContaminant;
        j_noisePos = [...
            find(isCleanNoise-[0,isCleanNoise(1:end-1)] > 0)',...
            find([isCleanNoise(2:end),0]-isCleanNoise < 0)'...
            ];
        
        % get start time of clip
        tClipStart = i_clipStart/fs;
        
        %** old code for isolating noise
        %{
        % isolate noise
        xNoiseInitial = xClipFilt(j_noisePos(1):j_noisePos(2));
        
        % look for non-target signals that exist in the noise window
        kn_otherSigBoxPos = j_otherSigBoxPos - j_noisePos(1) + 1;
        otherSignalsInNoise = find(any(kn_otherSigBoxPos > 0 & kn_otherSigBoxPos <= numAttemptedNoiseSamples, 2));
        
        % remove parts of noise that contain other signals
        xNoise = xNoiseInitial;
        for ss = 1:numel(otherSignalsInNoise)
            idx_ss = otherSignalsInNoise(ss);
            kn_otherSignalStart_ss = max([1,kn_otherSigBoxPos(idx_ss,1)]);
            kn_otherSignalStop_ss = min([kn_otherSigBoxPos(idx_ss,2),numAttemptedNoiseSamples]);
            xNoise(kn_otherSignalStart_ss:kn_otherSignalStop_ss) = NaN;
        end
        xNoise = xNoise(~isnan(xNoise));
        %}
        
        %** DEBUG PLOT
        %{
        if size(j_otherSigBoxPosClipped, 1) > 0
            tClip = (((1:numel(xClip))-1)./fs)';
            figure;
            ax = axes();
            ax.NextPlot = 'add';
            %plot(ax, tClip, xClip)
            plot(ax, tClip, xClipFilt)
            plot(ax, tClip(j_targetSigBoxPos(1):j_targetSigBoxPos(2))', xClipFilt(j_targetSigBoxPos(1):j_targetSigBoxPos(2)));
            plot(ax, tClip(j_targetSigEnergyPos(1):j_targetSigEnergyPos(2))', xClipFilt(j_targetSigEnergyPos(1):j_targetSigEnergyPos(2)))
            plot(ax, tClip(isInNoiseWin)', xClipFilt(isInNoiseWin));
            if sum(isContaminant) > 0
                stem(ax, tClip(isContaminant)', xClipFilt(isContaminant));
            else
                stem(ax, NaN, NaN);
            end
            stem(ax, tClip(isCleanNoise)', xClipFilt(isCleanNoise));
            %legend(ax, {'Unfiltered Clip', 'Filtered Clip', 'Signal Box', 'Energy Signal', 'Noise', 'Intra-Noise Signal'})
            legend(ax, {'Full Filtered Clip', 'Signal Box', 'Energy Signal', 'Noise', 'Contaminating Signal', 'Signal-Free Noise'})
            grid(ax, 'on')
            box(ax, 'on')
            keyboard
        end
        %}
        

        %** DEBUG PLOT (OLD)
        %{
        %if numel(otherSignalsInNoise) > 0
            j_intraNoiseSigs = [];
            for ss = 1:numel(otherSignalsInNoise)
                idx_ss = otherSignalsInNoise(ss);
                j_intraNoiseSigs = [j_intraNoiseSigs, j_otherSigBoxPos(idx_ss,1):j_otherSigBoxPos(idx_ss,2)];
            end
            j_intraNoiseSigs(j_intraNoiseSigs < 1) = 1;
            j_intraNoiseSigs(j_intraNoiseSigs > numel(xClip)) = numel(xClip);
            %tWav = (((1:numel(x))-1)./fs)';
            tClip = (((1:numel(xClip))-1)./fs)';
            %tSig = (((1:numel(xSignal))-1)./fs)';
            %tNoise = (((1:numel(xNoiseInitial))-1)./fs)';
            figure;
            ax = axes();
            ax.NextPlot = 'add';
            %plot(ax, tClip, xClip)
            plot(ax, tClip, xClipFilt)
            plot(ax, tClip(j_targetSigBoxPos(1):j_targetSigBoxPos(2))', xClipFilt(j_targetSigBoxPos(1):j_targetSigBoxPos(2)));
            plot(ax, tClip(j_targetSigEnergyPos(1):j_targetSigEnergyPos(2))', xSignal)
            plot(ax, tClip(j_noisePos(1):j_noisePos(2))', xNoiseInitial);
            stem(ax, tClip(j_intraNoiseSigs)', xClipFilt(j_intraNoiseSigs));
            %legend(ax, {'Unfiltered Clip', 'Filtered Clip', 'Signal Box', 'Energy Signal', 'Noise', 'Intra-Noise Signal'})
            legend(ax, {'Full Filtered Clip', 'Signal Box', 'Energy Signal', 'Noise', 'Intra-Noise Signal'})
            grid(ax, 'on')
            box(ax, 'on')
            keyboard
        %end
        %}
        
    else
        %%% if it's not possible to generate the clip (i.e., because the
        %%% signal is too close to the beginning of the sequence), then
        %%% return empties and NaNs
        xClipFilt = [];
        j_targetSigEnergyPos = double.empty(0,2);
        j_noisePos = double.empty(0,2);
        tClipStart = NaN;
    end
end