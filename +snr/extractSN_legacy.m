function [xSignal, xNoise] = extractSN_legacy(x, fs, sigStart, sigStop, noiseDist, clipBufferSize, dFilter, units)
% Isolate signal and associated noise samples from a larger audio time
% series vector, given a pre-determined signal location.
%
% Written by Wilfried Beslin
% Last updated by Wilfried Beslin
% 2024-05-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - Things I might do:
%   -- combine sigStart and sigStop into a single variable (makes for less
%   documentation and fewer input arguments
%   -- add input parsing with inputParser
%   -- add noise range as an output argument (in samples)
% - Additional thing to consider
%   -- exporting 90% energy 
%   -- add parameter to set cumulative energy thresehold

    
    import snr.noDelayFilt
    import snr.calcEng

    % get sigStart, sigStop, and noiseDist as samples, based on "units"
    switch units
        case 'seconds'
            samplesFromInput = @(a) round(a*fs);
        case 'samples'
            samplesFromInput = @(a) a;
        otherwise
            error('Invalid value for "units": must be either ''seconds'' or ''samples''.')
    end
    sigStartSample = samplesFromInput(sigStart);
    sigStopSample = samplesFromInput(sigStop);
   % get SNR params
    noiseDistSamples = samplesFromInput(noiseDist);
    clipBufferSamples = samplesFromInput(clipBufferSize);
    
    % get unclipped signal and noise size
    nSigSamples = sigStopSample - sigStartSample + 1;
    
    % generate shorter clip (easier to filter)
    %%% if it's not possible to generate the clip (i.e., because the signal
    %%% is too close to the beginning of the sequence), then return empties
    clipStartSample = sigStartSample - noiseDistSamples - nSigSamples - clipBufferSamples;
    clipStopSample = sigStopSample + clipBufferSamples;
    
    if clipStartSample > 0 && clipStopSample <= numel(x)
        xClip = x(clipStartSample:clipStopSample);

        % get relative signal
        sigStartSampleClip = clipBufferSamples + nSigSamples + noiseDistSamples + 1;
        sigStopSampleClip = sigStartSampleClip + nSigSamples - 1;

        % apply digital filter
        xClipFilt = noDelayFilt(dFilter, xClip);

        % isolate 90% energy signal
        xSigClip = xClipFilt(sigStartSampleClip:sigStopSampleClip);
        [Start90, Stop90] = calcEng(xSigClip,90);
        xSigStart90 = clipStartSample + sigStartSampleClip + Start90;
        xSigStop90 = clipStartSample + sigStartSampleClip + Stop90;
        xSignal = xSigClip(Start90:Stop90);
       
        %get realtive noise
        nSig90Samples = length(xSignal);
        noiseStartSampleClip = sigStartSampleClip - noiseDistSamples - nSig90Samples;
        noiseStopSampleClip = sigStartSampleClip - noiseDistSamples - 1;
        
        % isolate noise
         xNoise = xClipFilt(noiseStartSampleClip:noiseStopSampleClip);

        %** DEBUG PLOT
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