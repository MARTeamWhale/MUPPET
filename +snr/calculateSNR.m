function snr_dB = calculateSNR(xSigWin, xNoiseWin, fs, varargin)
% Calculate signal-to-noise ratio, given pre-isolated windows of signal and
% noise.
%
% This function accepts matrices for the signal and noise windows, where 
% rows represent samples and columns represent channels. The number of
% samples may be different between signal and noise windows, but the number
% of channels must be identical.
%
% If xSigWin represents a signal + noise mixture (almost always the case),
% the optional parameter 'SubtractNoise' may be set to true to get a more
% accurate measure of signal-to-noise ratio.
%
%
% Last updated by Wilfried Beslin
% 2024-02-12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % parse input
    p = inputParser();
    p.addRequired('xSigWin', @(val)validateattributes(val,{'numeric'},{'2d'}))
    p.addRequired('xNoiseWin', @(val)validateattributes(val,{'numeric'},{'2d'}))
    p.addRequired('fs', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
    p.addParameter('SubtractNoise', false, @(val)validateattributes(val,{'logical'},{'scalar'}))
    p.parse(xSigWin,xNoiseWin,fs,varargin{:});
    subtractNoise = p.Results.SubtractNoise;
    
    assert(size(xSigWin,2) == size(xNoiseWin,2), 'Signal and noise windows must have the same number of channels!')
    
    % calculate the power of the signal and noise windows
    %%% Remember: "power" in DSP usually means the **average power** of a
    %%% signal over a period of time, i.e., 
    %%%     sum(E_t*dt)/T
    %%% where E_t is the energy at time t (measured as x^2)
    %%% Power itself is a rate: energy/s
    avepow = @(x,dt,dur) sum((x.^2).*dt,1)/dur; % equation for calculating average power
    dt = 1/fs;
    
    durSigWin = size(xSigWin,1)/fs;
    pSigWin = avepow(xSigWin, dt, durSigWin);
    
    durNoiseWin = size(xNoiseWin,1)/fs;
    pNoiseWin = avepow(xNoiseWin, dt, durNoiseWin);
    
    % calculate SNR
    if subtractNoise
        snr_linear = (pSigWin - pNoiseWin)./pNoiseWin;
    else
        snr_linear = pSigWin./pNoiseWin;
    end
    snr_dB = 10*log10(snr_linear);
end