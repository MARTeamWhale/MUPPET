function call_params = extractCallParams(sigPosRel, fs, f_stft, psd, perEng, trace_line)
%
% Calculate several time, frequency, and time-frequency parameters for
% baleen whale calls.
%
% Written by Wilfried Beslin
% Last Updated by Wilfried Beslin 2024-07-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get handles of all local extraction functions
    param_fcns = localfunctions;
    num_call_params = numel(param_fcns);
    
    % initialize output
    call_params = struct();
    
    % loop through each parameter and calculate it
    for ii = 1:num_call_params
        param_fcn_ii = param_fcns{ii};
        param_name_ii = strrep(func2str(param_fcn_ii),'get_','');
        
        try
            param_val_ii = param_fcn_ii(sigPosRel, fs, f_stft, psd, perEng, trace_line, call_params);
            if isempty(param_val_ii)
                param_val_ii = NaN;
            end
        catch ME
            param_val_ii = NaN;
            % TEMPORARY
            %warning(ME.getReport)
            %keyboard 
            % END TEMPORARY
        end
        call_params.(param_name_ii) = param_val_ii;
    end
end


%% LOCAL FUNCTIONS ========================================================

% get_StartTime_Waveform --------------------------------------------------
function param = get_StartTime_Waveform(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the energy-based signal start time, in seconds relative to the
% start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = (sigPosRel(1) - 1)./fs;
end


% get_StartTime_Trace -----------------------------------------------------
function param = get_StartTime_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the signal start time based on the beginning time of the trace
% line, in seconds relative to the start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Remember that the trace times are relative to the energy-based signal
    % start time.
    param = trace_line.RelTime(1) + current_params.StartTime_Waveform;
end


% get_EndTime_Waveform ----------------------------------------------------
function param = get_EndTime_Waveform(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the energy-based signal end time, in seconds relative to the
% start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = (sigPosRel(2) - 1)./fs;
end


% get_StartTime_Trace -----------------------------------------------------
function param = get_EndTime_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the signal end time based on the ending time of the trace
% line, in seconds relative to the start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Remember that the trace times are relative to the energy-based signal
    % start time.
    param = trace_line.RelTime(end) + current_params.StartTime_Waveform;
end


% get_Duration_Waveform ---------------------------------------------------
function param = get_Duration_Waveform(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the energy-based duration of a call, in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = current_params.EndTime_Waveform - current_params.StartTime_Waveform;
end


% get_Duration_Trace ------------------------------------------------------
function param = get_Duration_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the duration of a call based on the trace line, in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = current_params.EndTime_Trace - current_params.StartTime_Trace;
end


% get_MinFrequency_Trace --------------------------------------------------
function param = get_MinFrequency_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the call's minimum significant frequency based on the trace line,
% in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = min(trace_line.Freq);
end


% get_MinFrequency_PSD ----------------------------------------------------
function param = get_MinFrequency_PSD(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the call's minimum significant frequency as measured from the
% full PSD, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [i_min, ~] = MUPPET.calcEng(sqrt(psd), perEng);
    param = f_stft(i_min);
end


% get_MaxFrequency_Trace --------------------------------------------------
function param = get_MaxFrequency_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the call's maximum significant frequency based on the trace line,
% in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = max(trace_line.Freq);
end


% get_MaxFrequency_PSD ----------------------------------------------------
function param = get_MaxFrequency_PSD(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the call's maximum significant frequency as measured from the
% full PSD, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, i_max] = MUPPET.calcEng(sqrt(psd), perEng);
    param = f_stft(i_max);
end


% get_Bandwidth_Trace -----------------------------------------------------
function param = Bandwidth_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the difference between the minimum and maximum significant
% frequencies as measured from the trace line, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = current_params.MaxFrequency_Trace - current_params.MinFrequency_Trace;
end


% get_Bandwidth_PSD -------------------------------------------------------
function param = Bandwidth_PSD(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the difference between the minimum and maximum significant
% frequencies as measured from the PSD, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = current_params.MaxFrequency_PSD - current_params.MinFrequency_PSD;
end


% get_TraceMinFreqTime ----------------------------------------------------
function param = get_TraceMinFreqTime(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the time at which the minimum frequency in the trace line occurs,
% in seconds relative to the start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Because of frequency quantization, there is a significant chance
% that there may be multiple instances of the same minimum frequency. If
% that happens, here I will take the mean time of the longest consecutive
% sequence of min frequencies. Since time is also highly quantized, if
% there happen to be multiple min frequency sequences with the same number
% of samples, then return NaN.

    % identify all instances where the max frequency occurs
    f_min = current_params.MinFrequency_Trace;
    is_min = trace_line.Freq == f_min;
    
    % find where sequences of the max start and end
    min_seq_starts = find(diff([0;is_min]) > 0);
    min_seq_stops = find(diff([is_min;0]) < 0);
    
    % get the length of each sequence and find the longest one
    seq_durs = min_seq_stops - min_seq_starts + 1;
    longest_seq = max(seq_durs);
    i_longest_seq = find(seq_durs == longest_seq);
    
    % get the mean time of the longest sequence, if there is only one
    if i_longest_seq == 1
        i_mean_min = mean([min_seq_starts(i_longest_seq),min_seq_stops(i_longest_seq)]);
        param = interp1((1:numel(trace_line.RelTime))', trace_line.RelTime, i_mean_min) + current_params.StartTime_Waveform;
    else
        param = NaN;
    end
end


% get_TraceMaxFreqTime ----------------------------------------------------
function param = get_TraceMaxFreqTime(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the time at which the maximum frequency in the trace line occurs,
% in seconds relative to the start of the annotation box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Because of frequency quantization, there is a significant chance
% that there may be multiple instances of the same maximum frequency. If
% that happens, here I will take the mean time of the longest consecutive
% sequence of max frequencies. Since time is also highly quantized, if
% there happen to be multiple max frequency sequences with the same number
% of samples, then return NaN.

    % identify all instances where the max frequency occurs
    f_max = current_params.MaxFrequency_Trace;
    is_max = trace_line.Freq == f_max;
    
    % find where sequences of the max start and end
    max_seq_starts = find(diff([0;is_max]) > 0);
    max_seq_stops = find(diff([is_max;0]) < 0);
    
    % get the length of each sequence and find the longest one
    seq_durs = max_seq_stops - max_seq_starts + 1;
    longest_seq = max(seq_durs);
    i_longest_seq = find(seq_durs == longest_seq);
    
    % get the mean time of the longest sequence, if there is only one
    if i_longest_seq == 1
        i_mean_max = mean([max_seq_starts(i_longest_seq),max_seq_stops(i_longest_seq)]);
        param = interp1((1:numel(trace_line.RelTime))', trace_line.RelTime, i_mean_max) + current_params.StartTime_Waveform;
    else
        param = NaN;
    end
end


% get_PeakFrequency_Trace -------------------------------------------------
function param = get_PeakFrequency_Trace(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the time at which the frequency with greatest energy in the trace
% line occurs, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, i_max] = max(trace_line.Power_dB);
    param = trace_line.Freq(i_max);
end


% get_PeakFrequency_PSD ---------------------------------------------------
function param = get_PeakFrequency_PSD(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the time at which the frequency with greatest energy in the
% spectrum occurs, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, i_max] = max(psd);
    param = f_stft(i_max);
end


% get_PSDCentreFrequency --------------------------------------------------
function param = get_PSDCentreFrequency(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the frequency which divides the spectrum into two halves of equal
% energy, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This calculation is based on Equation 12.14 in:
%   Au and Hastings (2008), "Principles of Marine Bioacoustics"
%
% Here I am limiting the calculation to be between the min and max PSD
% frequencies calculated earlier.

    % remove frequencies that are insignificant
    is_sig_freq = (f_stft >= current_params.MinFrequency_PSD) & (f_stft <= current_params.MaxFrequency_PSD);
    f_stft_sig = f_stft(is_sig_freq);
    psd_sig = psd(is_sig_freq);
    
    % calculate centroid frequency
    df = mean(unique(diff(f_stft)));
    param = (sum(f_stft_sig.*psd_sig).*df)./(sum(psd_sig).*df);
end


% get_TraceStartFrequency -------------------------------------------------
function param = get_TraceStartFrequency(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the frequency at which the trace line begins, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = trace_line.Freq(1);
end 


% get_TraceEndFrequency -------------------------------------------------
function param = get_TraceEndFrequency(sigPosRel, fs, f_stft, psd, perEng, trace_line, current_params)
% Returns the frequency at which the trace line ends, in Hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param = trace_line.Freq(end);
end 