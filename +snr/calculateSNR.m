function [snr_dB, snr_adjusted_dB] = calculateSNR(xSigWin, xNoiseWin, varargin)
% Calculate signal-to-noise ratio, given pre-isolated windows of signal and
% noise.
%
% This function accepts matrices for the signal and noise windows, where 
% rows represent samples and columns represent channels. The number of
% samples may be different between signal and noise windows, but the number
% of channels must be identical.
%
% The first output is the base SNR equation, that is:
%   SNR = powSig/powNoise
% where powSig and powNoise are the average power of the signal and noise
% windoes, respectively.
%
% The second output subtracts the estimated noise power from the numerator,
% that is:
%   SNR_Adjusted = (powSig - powNoise)/powNoise
% This representation of SNR is more true to the real definition of SNR
% when the "signal" window actually contains a signal + noise mixture,
% which is almost always the case.
%
%
% Last updated by Wilfried Beslin
% 2024-03-12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % parse input
    p = inputParser();
    p.addRequired('xSigWin', @(val)validateattributes(val,{'numeric'},{'2d'}))
    p.addRequired('xNoiseWin', @(val)validateattributes(val,{'numeric'},{'2d'}))
    %p.addRequired('fs', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}))
    %p.addParameter('SubtractNoise', false, @(val)validateattributes(val,{'logical'},{'scalar'}))
    %p.parse(xSigWin, xNoiseWin, fs, varargin{:});
    p.parse(xSigWin, xNoiseWin, varargin{:});
    %subtractNoise = p.Results.SubtractNoise;
    
    assert(size(xSigWin,2) == size(xNoiseWin,2), 'Signal and noise windows must have the same number of channels!')
    
    % calculate the power of the signal and noise windows
    avepow = @(x) sum((x.^2),1)./size(x,1); % equation for calculating RMS-based average power for each channel
    pSigWin = avepow(xSigWin);
    pNoiseWin = avepow(xNoiseWin);
    
    %%% Old time-based calculation
    %{
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
    %}
    
    % calculate SNR
    snr_dB = 10*log10(pSigWin./pNoiseWin);
    snr_adjusted_dB = 10*log10((pSigWin - pNoiseWin)./pNoiseWin);
    %{
    if subtractNoise
        snr_linear = (pSigWin - pNoiseWin)./pNoiseWin;
    else
        snr_linear = pSigWin./pNoiseWin;
    end
    snr_dB = 10*log10(snr_linear);
    %}
end