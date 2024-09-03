function psdm_smooth = smoothSpec(psdm)
% function for smoothing a spectrogram. Smoothing is done by convolving the
% spectrogram with a Gaussian kernel (based on Baumgartner and Mussoline).
% Note that special corrections must be applied to take away edge effects.
%
% Written by Wilfried Beslin
% Last updated 2024-09-03
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create smoothing kernel
% (Baumgartner and Mussoline 2011, eq. 2)
spec_smooth_kernel = [1,2,1; 2,4,2; 1,2,1];

% duplicate the time and frequency edges of the spectrogram to prepare for
% edge effects
psdm_anal = psdm;
psdm_anal = [psdm_anal(:,1), psdm_anal, psdm_anal(:,end)];
psdm_anal = [psdm_anal(1,:); psdm_anal; psdm_anal(end,:)];

% apply smoothing kernel
psdm_anal_smooth = conv2(psdm_anal, spec_smooth_kernel);

% remove the edges introduced by both edge duplication and convolution 
% (2 levels)
psdm_smooth = psdm_anal_smooth(3:(end-2), 3:(end-2));