function [xSignal, xNoise] = extractSN(x, fs, sigStart, sigStop, noiseDist, units)
% Isolate signal and associated noise samples from a larger audio time
% series vector, given a pre-determined signal location.
%
% Last updated by Wilfried Beslin
% 2024-02-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - Consider whether or not bandpass filtering should be implemented as an
% option within this function


    % get sigStart, sigStop, and noiseDist as samples, based on "units"
    switch units
        case 'time'
            samplesFromInput = @(a) round(a*fs);
        case 'samples'
            samplesFromInput = @(a) a;
        otherwise
            error('Invalid value for "units": must be either ''time'' or ''samples''.')
    end
    sigStartSample = samplesFromInput(sigStart);
    sigStopSample = samplesFromInput(sigStop);
    noiseDistSamples = samplesFromInput(noiseDist);
    
    % isolate signal
    xSignal = x(sigStartSample:sigStopSample);
    nSamples = numel(xSignal);
    
    % isolate noise
    noiseStartSample = sigStartSample - noiseDistSamples - nSamples;
    noiseStopSample = sigStartSample - noiseDistSamples - 1;
    xNoise = x(noiseStartSample:noiseStopSample);
end