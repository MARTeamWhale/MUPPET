function PARAMS = importInputParams(paramFilePath)
%
% Returns a struct of parameters read from a text file of parameters for
% MUPPET.
%
% Written by Wilfried Beslin
% Last updated by Wilfried Beslin
% 2024-11-22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % load the file as a single char vector
    paramFileText = fileread(paramFilePath);
    %paramFileText = splitlines(fileread(paramFilePath));

    % initialize parameter struct
    PARAMS = struct();

    %%% Bandpass filter parameters
    PARAMS.Lower_Passband_Frequency = readParam(paramFileText, 'Lower_Passband_Frequency', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.Upper_Passband_Frequency = readParam(paramFileText, 'Upper_Passband_Frequency', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.Stopband_Rolloff_Bandwidth = readParam(paramFileText, 'Stopband_Rolloff_Bandwidth', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});

    %%% Signal and noise isolation
    PARAMS.Noise_Distance = readParam(paramFileText, 'Noise_Distance', {@(var)validateattributes(var,{'numeric'},{'scalar','nonnegative'})});
    PARAMS.Ideal_Noise_Duration  = readParam(paramFileText, 'Ideal_Noise_Duration', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})}); 
    PARAMS.Signal_Energy_Percent  = readParam(paramFileText, 'Signal_Energy_Percent', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','<=',100})});

    %%% STFT parameters
    PARAMS.STFT_Win_Size  = readParam(paramFileText, 'STFT_Win_Size', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.STFT_Step_Size  = readParam(paramFileText, 'STFT_Step_Size', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.NFFT_8kHz  = readParam(paramFileText, 'NFFT_8kHz', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','integer'})});
    PARAMS.Smooth_Spec  = readParam(paramFileText, 'Smooth_Spec', {@(var)validateattributes(var,{'logical'},{'scalar'})});

    %%% Resampling
    PARAMS.Downsampled_Sampling_Rate = readParam(paramFileText, 'Downsampled_Sampling_Rate', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','integer'})});

    %%% Trace line calculation
    PARAMS.Trace_Penalty_Coefficients = readParam(paramFileText, 'Trace_Penalty_Coefficients', {@(var)validateattributes(var,{'numeric'},{'numel',2})});
    PARAMS.Trace_Energy_Percent = readParam(paramFileText, 'Trace_Energy_Percent', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','<=',100})});
    PARAMS.Trace_Threshold_Type = readParam(paramFileText, 'Trace_Threshold_Type', {@(var)validateattributes(var,{'char'},{'vector'})});
    PARAMS.Trace_Threshold_Val = readParam(paramFileText, 'Trace_Threshold_Val', {@(var)validateattributes(var,{'numeric'},{'scalar'})});

    %%% Trace line plot parameters
    PARAMS.Spec_Plot_Colour_Map = readParam(paramFileText, 'Spec_Plot_Colour_Map', {@(var)validateattributes(var,{'char'},{'vector'})});
    PARAMS.Log_Spec_Plot_Colours = readParam(paramFileText, 'Log_Spec_Plot_Colours', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.Trace_Plot_Line_Colour  = readParam(paramFileText, 'Trace_Plot_Line_Colour', {@(var)validateattributes(var,{'char'},{'scalar'}),@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1})});
    PARAMS.Trace_Plot_Line_Width  = readParam(paramFileText, 'Trace_Plot_Line_Width', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.Trace_Plot_Marker_Type  = readParam(paramFileText, 'Trace_Plot_Marker_Type', {@(var)validateattributes(var,{'char'},{'scalar'})});
    PARAMS.Trace_Plot_Marker_Size  = readParam(paramFileText, 'Trace_Plot_Marker_Size', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});

    %%% Advanced parameters
    PARAMS.Cap_Noise = readParam(paramFileText, 'Cap_Noise', {@(var)validateattributes(var,{'logical'},{'scalar'})});
end


%% readParam --------------------------------------------------------------
% This function was copied from BAIT
function var = readParam(filestr, varname, validation_fcns)
% Read a parameter value from a parameter file.
% The validation_fcns argument is a cell array of function handles for
% ensuring the extracted variable is valid. It should be populated with
% functions that throw errors when conditions are not met, such as 'assert'
% and 'validateattributes'. The variable is considered valid if it passes
% at least one condition.
%
%   Written by Wilfried Beslin
%   Last updated 2023-12-06 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define standard regex for reading parameters
    param_prefix_regex = 'VAR[ \t]*=[ \t]*';  % regular expression part to isolate parameter name and the equal sign that comes after it
    param_std_value_regex = '[\w.-]+';  % regular expression part to isolate standard (scalar) parameter values
    param_array_value_regex = '\[[\d.-, ]*\]';  % regular expression part to isolate multivalue parameters enclosed in braces
    param_regex = ['(?<=',param_prefix_regex,')((',param_std_value_regex,')|(',param_array_value_regex,'))'];

    % read variable string
    varstr_raw = char(regexp(filestr, strrep(param_regex,'VAR',varname), 'match'));
    varstr_num = regexpi(varstr_raw, '[\d.-]+|NaN', 'match');
    varstr_bool = regexpi(varstr_raw, 'true|false', 'match');

    % check what kind of variable it is
    if ~isempty(varstr_num)
        % numeric inputs
        var = str2double(varstr_num);
    elseif ~isempty(varstr_bool)
        % logical inputs
        var = strcmpi(varstr_bool, 'true');
    elseif ~isempty(varstr_raw) && ~strcmp(varstr_raw,'[]')
        % char input
        var = varstr_raw;
    else
        % empty
        var = [];
    end

    % validate input
    valid_input = true;
    for ii = 1:numel(validation_fcns)
        try
            feval(validation_fcns{ii}, var)
            valid_input = true;
            break
        catch
            valid_input = false;
        end
    end
    if ~valid_input
        error('The value of parameter "%s" is invalid or could not be read', varname)
    end
end