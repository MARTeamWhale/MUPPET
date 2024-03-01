function xFilt = noDelayFilt(dFilter, x)
% Filter a signal using a digitalFilter object from the Signal Processing
% Toolbox, compensating for group delay introduced by the filter. This 
% function only works if the delay is not frequency-dependent (usually the
% case with FIR filters). It will generally NOT work with IIR filters like 
% the Butterworth filter.
%
% Last updated by Wilfried Beslin
% 2024-02-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% - consider supporting multichannel inputs (this would involve modifying
% how the 'zeros' function is used, as well as the indexing operation for
% subsetting xFilt)

    % get group delay
    gd = grpdelay(dFilter);
    delay = round(unique(gd));
    
    % make sure delay is not frequency-dependent
    assert(numel(delay) == 1, 'Delay is frequency-dependent, unable to compensate')
    
    % filter signal with zeros added
    xFilt = filter(dFilter, [x; zeros(delay,1)]);
    
    % remove first part
    xFilt = xFilt((delay+1):end);
end