function snrVal = calculateSNR(xSignal, xNoise)
% Calculate Signal to Noise Ratio (SNR)
% 
%
% Last updated by Mike Adams
% 2024-02-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTES: Still need to decide which SNR calc to use and wether or not to
%remove mean noise from signal or not. 
% 
% Remove Mean Noise
RMN = 1; % change to 0 to skip
if RMN == 1
   xSignal = xSignal - mean(xNoise); 
end

% Temporary JASCO SNR calculation
snr1 = 20*log10(rssq(xSignal)/rssq(xNoise));

% Matlab function
snr2 = snr(xSignal,xNoise);

% Equation used in matlab function (same as JASCO...)
snr3 = mag2db(rssq(xSignal)/rssq(xNoise));

% Method adapted from BWD
%signal
nSignal = numel(xSignal);
rmsSignalLin = sqrt(sum(xSignal.^2)/nSignal);
rmsSignal = 20*log10(rmsSignalLin);
% noise
nNoise = numel(xNoise);  
rmsNoiseLin = sqrt(sum(xNoise.^2)/nNoise);
rmsNoise = 20*log10(rmsNoiseLin);
% SNR
snr4 = rmsSignal - rmsNoise;


snrVal = snr4; %temporary

end