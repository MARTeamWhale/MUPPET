function [psdm_out, psdm_noise] = processSpec(psdm_sig, psdmc_noise, ops)
% function for processing a signal spectrogram in a particular order
% requested by the user. Processing operations are "Smooth", "Log", and
% "Denoise". Any combination of the above may be specified, and some may be
% omitted (but at least one should be specified).
%
% The noise PSD matrix must be specified as a cell, since MUPPET may 
% produce multiple noise sources.
%
% Written by Wilfried Beslin
% Last updated 2024-09-10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = inputParser();
    p.addRequired('psdm_sig',@(val)validateattributes(val,{'numeric'},{'2d'}));
    p.addRequired('psdmc_noise',@(val)validateattributes(val,{'cell'},{}));
    p.addRequired('ops',@iscellstr);
    p.parse(psdm_sig, psdmc_noise, ops);
    
    % further validate operations
    valid_ops = {'log','smooth','denoise'};
    ops = lower(ops);
    assert(all(ismember(ops, valid_ops)), 'One or more operation is unrecognized')
    assert(numel(ops) == numel(unique(ops)), 'Cannot specify an operation more than once!')
    
    % initialize variables
    psdm_anal = psdm_sig;
    psdm_noise = [];
    psdmc_noise_anal = psdmc_noise;
    num_ops = numel(ops);
    num_noisespecs = numel(psdmc_noise);
    process_noise = ismember('denoise',ops);
    
    % loop through each operation
    for ii = 1:num_ops
        op_ii = ops{ii};
        switch op_ii
            case 'log'
                % convert spectrogram(s) to decibels
                psdm_anal = 10.*log10(psdm_anal);
                if process_noise
                    for jj = 1:num_noisespecs
                        try
                            psdmc_noise_anal{jj} = 10.*log10(psdmc_noise_anal{jj});
                        catch
                            psdmc_noise_anal{jj} = [];
                        end
                    end
                end
                
            case 'smooth'
                % smooth spectrogram(s)
                psdm_anal = MUPPET.smoothSpec(psdm_anal);
                if process_noise
                    for jj = 1:num_noisespecs
                        try
                            psdmc_noise_anal{jj} = MUPPET.smoothSpec(psdmc_noise_anal{jj});
                        catch
                            psdmc_noise_anal{jj} = [];
                        end
                    end
                end
                
            case 'denoise'
                % remove noise from signal spectrogram
                %%% concatenate noise spectrogram if possible
                bad_noise = cellfun('isempty',psdmc_noise);
                if ~all(bad_noise)
                    psdm_noise = horzcat(psdmc_noise{~bad_noise});
                else
                    warning('Insufficient noise spectrogram data available; will estimate noise from signal spectrogram instead')
                    psdm_noise = median(psdm_anal, 2);
                end
                
                %%% create a denoised version of the signal spectrogram
                %%% that shifts each frequency band such that 0 corresponds
                %%% to the average noise level in that band
                psdm_anal = psdm_anal - mean(psdm_noise, 2);
                
                %%% apply further corrections to account for time-varying
                %%% noise by equalizing the medians of each time bin
                psdm_anal = psdm_anal - median(psdm_anal, 1);
                
                %%% update noise processing flag
                process_noise = false;
                
            otherwise
                error('Unsupported operation "%s"',op_ii)
        end
    end
    
    % save output
    psdm_out = psdm_anal;
end