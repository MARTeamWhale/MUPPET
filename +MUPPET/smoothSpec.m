function psdm_smooth = smoothSpec(psdm)
% function for smoothing a spectrogram. Smoothing is done by convolving the
% spectrogram with a Gaussian kernel (based on Baumgartner and Mussoline).
% Note that special corrections must be applied to take away edge effects.
%
% Written by Wilfried Beslin
% Last updated 2024-09-04
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create smoothing kernel
% (Modified from Baumgartner and Mussoline 2011, eq. 2)
%%% Here the kernel matrix is divided by the sum of all elements - this is
%%% done because using the original matrix actually changes the scale of
%%% the magnitudes being smoothed, which is not good (this can impact
%%% spectrogram-based SNR calculations, for instance).
spec_smooth_kernel = [1,2,1; 2,4,2; 1,2,1];
spec_smooth_kernel = spec_smooth_kernel./sum(spec_smooth_kernel(:));

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

%** DEBUG
% create a plot to compare the magnitudes of the smoothed and unsmoothed
% spectrograms
%{
fig = figure();
ax1 = subplot(1,2,1);
surf(ax1, psdm, 'EdgeColor','none');
view(ax1,[0,0])

ax2 = subplot(1,2,2);
surf(ax2, psdm_smooth, 'EdgeColor','none')
view(ax2,[0,0])
keyboard
close(fig);
%}