function [xSignal, xNoise] = extractSN(x, fs, targetSigBoxPos, otherSigBoxPos, noiseDist, idealNoiseSize, clipBufferSize, dFilter, units)
% Isolate signal and associated noise samples from a larger audio time
% series vector, given a pre-determined signal location.
%
% This function can accept the locations of other, non-target signals in
% the time series as an input argument. These other signals will be removed
% from the noise window, in which case xNoise will consist of a truncated
% and possibly concatenated noise vector.
%
%
% Written by Wilfried Beslin
% Last updated by Wilfried Beslin
% 2024-05-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - [Wilfried] Things I might do:
%   -- add input parsing with inputParser
%   -- add noise range as an output argument (in samples)
%   -- ideal noise size could be an optional argument, where the default
%       makes it equal to the duration of the signal
% - [Mike] Additional thing to consider
%   -- exporting 90% energy 
%   -- add parameter to set cumulative energy thresehold

    
    import snr.noDelayFilt
    import snr.calcEng

    % define equation to get time series positional parameters in samples,
    % depending on the value of "units"
    switch units
        case 'seconds'
            samplesFromInput = @(a) round(a*fs);
        case 'samples'
            samplesFromInput = @(a) a;
        otherwise
            error('Invalid value for "units": must be either ''seconds'' or ''samples''.')
    end
    
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
    
    % get all positional and size parameters in samples
    i_targetSigBoxPos = samplesFromInput(targetSigBoxPos);
    i_otherSigBoxPos = samplesFromInput(otherSigBoxPos);
    numNoiseDistSamples = samplesFromInput(noiseDist);
    numIdealNoiseSamples = samplesFromInput(idealNoiseSize);
    numClipBufferSamples = samplesFromInput(clipBufferSize);
    
    %sigStartSample = samplesFromInput(sigStartTarget);
    %sigStopSample = samplesFromInput(sigStopTarget);
    %sigStartSamplesOther = samplesFromInput(sigStartOther);
    %sigStopSamplesOther = samplesFromInput(sigStopOther);
    % get SNR params
    %noiseDistSamples = samplesFromInput(noiseDist);
    %nNoiseSamplesIdeal = samplesFromInput(noiseIdealSize);
    %clipBufferSamples = samplesFromInput(clipBufferSize);
    
    % get size of signal annotation box in samples
    %numSigBoxSamples = sigStopSample - sigStartSample + 1;
    numTargetSigBoxSamples = i_targetSigBoxPos(2) - i_targetSigBoxPos(1) + 1;
    
    % generate shorter clip (easier to filter)
    %%% if it's not possible to generate the clip (i.e., because the signal
    %%% is too close to the beginning of the sequence), then return empties
    %clipStartSample = sigStartSample - noiseDistSamples - nSigSamples - clipBufferSamples;
    %clipStopSample = sigStopSample + clipBufferSamples;
    i_clipStart = i_targetSigBoxPos(1) - numNoiseDistSamples - numIdealNoiseSamples - numClipBufferSamples;
    i_clipStop = i_targetSigBoxPos(2) + numClipBufferSamples;
    %clipStartSample = sigStartSample - numNoiseDistSamples - numIdealNoiseSamples - numClipBufferSamples;
    %clipStopSample = sigStopSample + numClipBufferSamples;
    
    if i_clipStart > 0 && i_clipStop <= numel(x)
        xClip = x(i_clipStart:i_clipStop);

        % get signal box start/stop samples relative to clip
        %%% target signal
        %** Old complicated way
        j_targetSigBoxPos_OLD = numClipBufferSamples + numIdealNoiseSamples + numNoiseDistSamples + [1, numTargetSigBoxSamples];
        %sigStartSampleClip_OLD = numClipBufferSamples + numTargetSigBoxSamples + numNoiseDistSamples + 1;
        %sigStopSampleClip_OLD = sigStartSampleClip_OLD + numTargetSigBoxSamples - 1;
        %** new easy way [still need to verify that this works]
        j_targetSigBoxPos = i_targetSigBoxPos - i_clipStart + 1;
        %sigStartSampleClip = sigStartSample - clipStartSample + 1;
        %sigStopSampleClip = sigStopSample - clipStartSample + 1;
        %** DEBUG
        %fprintf('Old sig clip indices: [%d, %d]\n New sig clip indices: [%d, %d]\n', j_targetSigBoxPos_OLD(1), j_targetSigBoxPos_OLD(2), j_targetSigBoxPos(1), j_targetSigBoxPos(2))
        %keyboard
        %** END DEBUG
        %%% other signals
        j_otherSigBoxPos = i_otherSigBoxPos - i_clipStart + 1;
        %sigStartSamplesClipOther = sigStartSamplesOther - clipStartSample + 1;
        %sigStopSamplesClipOther = sigStopSamplesOther - clipStartSample + 1;

        % apply digital filter
        xClipFilt = noDelayFilt(dFilter, xClip);

        % isolate 90% energy signal
        xSigBoxFilt = xClipFilt(j_targetSigBoxPos(1):j_targetSigBoxPos(2));
        [ks_start90, ks_stop90] = calcEng(xSigBoxFilt,90);
        %xSigStart90 = clipStartSample + sigStartSampleClip + Start90;
        %xSigStop90 = clipStartSample + sigStartSampleClip + Stop90;
        xSignal = xSigBoxFilt(ks_start90:ks_stop90);
       
        % get relative noise position
        %nSig90Samples = length(xSignal);
        j_targetSigEnergyPos = j_targetSigBoxPos(1) + [ks_start90,ks_stop90] - 1;
        j_noisePos = j_targetSigEnergyPos(1) - numNoiseDistSamples - [numIdealNoiseSamples, 1];
        %noiseStartSampleClip = sigStartSampleClip - noiseDistSamples - nSig90Samples;
        %noiseStartSampleClip = sigStartSampleClip - numNoiseDistSamples - numIdealNoiseSamples;
        %noiseStopSampleClip = sigStartSampleClip - numNoiseDistSamples - 1;
        
        % isolate noise
        xNoiseInitial = xClipFilt(j_noisePos(1):j_noisePos(2));
        
        % look for non-target signals that exist in the noise window
        kn_otherSigBoxPos = j_otherSigBoxPos - j_noisePos(1) + 1;
        %otherSignalsInNoise = find(kn_otherSigBoxPos(:,1) > 0 & kn_otherSigBoxPos(:,2) <= numIdealNoiseSamples);
        otherSignalsInNoise = find(any(kn_otherSigBoxPos > 0 & kn_otherSigBoxPos <= numIdealNoiseSamples, 2));
        
        % remove parts of noise that contain other signals
        xNoise = xNoiseInitial;
        for ss = 1:numel(otherSignalsInNoise)
            idx_ss = otherSignalsInNoise(ss);
            kn_otherSignalStart_ss = max([1,kn_otherSigBoxPos(idx_ss,1)]);
            kn_otherSignalStop_ss = min([kn_otherSigBoxPos(idx_ss,2),numIdealNoiseSamples]);
            xNoise(kn_otherSignalStart_ss:kn_otherSignalStop_ss) = NaN;
        end
        xNoise = xNoise(~isnan(xNoise));
        
        % remove parts of noise that contain other signals
        %xNoise = xNoiseInitial;
        %numOtherSigs = numel(sigStartOther);
        %for ii = 1:numOtherSigs
        %    sigSamplesOther_ii = (sigStartSamplesClipOther(ii):sigStopSamplesClipOther(ii)) - noiseStartSampleClip + 1;
        %    xNoise(sigSamplesOther_ii) = NaN;
        %end
        %xNoise = xNoise(~isnan(xNoise));

        %** DEBUG PLOT
        %{
        if numel(otherSignalsInNoise) > 0
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
        end
        %}
        
        %** DEBUG PLOT (OLD)
        %{
        tClip = (((1:numel(xClip))-1)./fs)';
        figure;
        ax = axes();
        ax.NextPlot = 'add';
        plot(ax, tClip, xClip)
        plot(ax, tClip, xClipFilt)
        plot(ax, tClip(sigStartSampleClip:sigStopSampleClip), xSignal);
        plot(ax, tClip(noiseStartSampleClip:noiseStopSampleClip), xNoise);
        legend(ax, {'Unfiltered Clip', 'Filtered Clip', 'Signal', 'Noise'})
        grid(ax, 'on')
        box(ax, 'on')
        keyboard
        %}
    else
        xSignal = [];
        xNoise = [];
    end
end