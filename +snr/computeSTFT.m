function [spectrogram_data, welchPSD_data] = computeSTFT(varargin)
%
% Calculates the spectrogram of a given signal via the short-time Fourier
% transform (STFT). Also returns the Welch power spectral density (PSD)
% estimate, which is essentially just the average of STFT bins.
%
% OUTPUT:
% spectrogram_data = struct with spectrogram data;
%   fields are "t", "f", and "psd"
%
% welchPSD_data = struct with Welch PSD estimate data;
%   fields are "f" and "psd"
%
% Written by Wilfried Beslin
% Last Updated 2024-05-21 by Wilfried Beslin
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % input parsing
    p = inputParser();
    
    p.addRequired('x', @(val)validateattributes(val,{'numeric'},{'vector'}))
    p.addRequired('fs', @(val)validateattributes(val,{'numeric'},{'positive','scalar'}))
    p.addRequired('sigPos', @(val)validateattributes(val,{'numeric'},{'positive','nondecreasing','numel',2}))
    p.addRequired('winSize', @(val)validateattributes(val,{'numeric'},{'positive','integer','scalar'}))
    p.addRequired('nOverlap', @(val)validateattributes(val,{'numeric'},{'positive','integer','scalar'}))
    p.addParameter('WindowType', 'hamming', @(val)validateattributes(val,{'char','function_handle'},{}))
    p.addParameter('NFFT', [], @(val)validateattributes(val,{'numeric'},{'positive','integer','scalar'}))
    
    p.parse(varargin{:})
    x = p.Results.x;
    fs = p.Results.fs;
    sigPos = p.Results.sigPos;
    winSize = p.Results.winSize;
    nOverlap = p.Results.nOverlap;
    winType = p.Results.WindowType;
    nfft = p.Results.NFFT;
    if isempty(nfft)
       %nfft = 2^nextpow2(winSize);
       nfft = length(x);
    else
        assert(nfft >= winSize, 'NFFT values smaller than the window size are not allowed')
    end
    
    % isolate the part of the time series to use in the STFT
    %%% Should start when just over half the window covers the signal, and 
    %%% end just before the window would contain more noise than signal
    i_stftStart = sigPos(1) - floor(winSize/2);
    stftStepSize = winSize - nOverlap;
    stftInitialBinStarts = i_stftStart:stftStepSize:sigPos(end);
    i_stftStop = stftInitialBinStarts(find(stftInitialBinStarts + ceil(winSize/2) - 1 <= sigPos(end), 1, 'last')) + winSize;
    
    %%% ensure that there are enough samples
    assert(i_stftStart > 0 && i_stftStop <= length(x), 'Sample indices for STFT are out of bounds')
    
    %%% extract clip for STFT
    xSTFT = x(i_stftStart:i_stftStop);
    
    % get windowing vector
    winVec = feval(winType, winSize);
    
    % calculate spectrogram and save to output
    [~, f_gram, t_gram, psd_gram] = spectrogram(xSTFT, winVec, nOverlap, nfft, fs);
    %%% The "spectrogram" function returns times corresponding to the
    %%% midpoint of the spectrogram bins, which is what we want. However,
    %%% since the input vector began a little before the signal, the bin
    %%% midpoint times should be adjusted to correct for this, such that
    %%% they are relative to the start of the signal. This correction
    %%% should make them start at about t=0.
    t_gram = t_gram - ((sigPos(1) - i_stftStart)/fs);
    spectrogram_data = struct(...
        'f', f_gram,...
        't', t_gram,...
        'psd', psd_gram);
    
    % calculate Welch PSD estimate and save to output
    [psd_welch, f_welch] = pwelch(xSTFT, winVec, nOverlap, nfft, fs);
    welchPSD_data = struct(...
        'f', f_welch,...
        'psd', psd_welch);
    
    % DEBUG PLOT
    %{
    fMax_plot = 300;
    fKeep_plot = spectrogram_data.f <= fMax_plot;
    figure;
    ax_gram = subplot(1, 2, 1);
    ax_welch = subplot(1, 2, 2);
    
    %%% plot spectrogram as waterfall plot with each strip corresponding to
    %%% a time bin
    [T, F] = meshgrid(spectrogram_data.t, spectrogram_data.f(fKeep_plot));
    waterfall(ax_gram, T', F', spectrogram_data.psd(fKeep_plot,:)');
    xlabel(ax_gram, 'Time [s]')
    ylabel(ax_gram, 'Frequency [Hz]')
    zlabel(ax_gram, 'PSD')
    title(ax_gram, 'Spectrogram');
    
    %%% plot Welch PSD
    plot(ax_welch, welchPSD_data.f(fKeep_plot), welchPSD_data.psd(fKeep_plot))
    xlabel(ax_welch, 'Frequency [Hz]')
    ylabel(ax_welch, 'PSD')
    title(ax_welch, 'Welch PSD Estimate')
    grid(ax_welch, 'on')
    %}
    
end