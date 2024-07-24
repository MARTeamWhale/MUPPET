function varargout = Baleen_SNR_Tool(varargin)
%
% Baleen_SNR_Tool.m
%
% Process Pamlab output for use in SNR tool.
%
% SYNTAX:
%   Baleen_SNR_Tool
%   Baleen_SNR_Tool(Name,Value)
%   out = Baleen_SNR_Tool(__)
%
% OPTIONAL INPUT ARGUMENTS (Name-Value Pairs):
%   .......................................................................
%   "PAMLAB_DATA_FOLDER" - Path to folder with PAMlab output. If not
%       specified, user will be prompted to select the path manually.
%   .......................................................................
%   "WAV_FILE_FOLDER" - Path to folder with WAV files. If not specified,
%       user will be prompted to select the path manually.
%   .......................................................................
%   "OUTPUT_FOLDER_LOCATION" - Path where the tool's output folder will be
%       created. If not specified, user will be prompted to select the path
%       manually, or simply use the parent folder of the PAMlab output by
%       cancelling the prompt.
%   .......................................................................
%   "PARAMFILE" - Path to CSV file of SNR parameters. If not specified, the
%       tool will try to load the default parameters from the file
%       SNR_PARAMS.csv.
%   .......................................................................
%
% OPTIONAL OUTPUT ARGUMENTS:
%   .......................................................................
%   "out" = struct containing the output annotation tables that were saved
%       to CSV files. Mostly useful for debugging or otherwise working with
%       the results directly in the MATLAB workspace. Only returned if
%       explicitly requested.
%   .......................................................................
%
%
% Written by Mike Adams
% Last updated by Wilfried Beslin
% 2024-07-24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTE: https://www.mathworks.com/help/matlab/ref/listdlg.html

    % check if output arguments were requested
    nargoutchk(0,1) % allows the function to return exactly one or zero output parameters
    returnOutputVar = nargout > 0;

    % process optional input parameters (I/O paths)
    ip = inputParser();
    ip.addParameter('PAMLAB_DATA_FOLDER', '', @isfolder)
    ip.addParameter('WAV_FILE_FOLDER', '', @isfolder)
    ip.addParameter('OUTPUT_FOLDER_LOCATION', '', @isfolder)
    ip.addParameter('PARAMFILE', '', @isfile)
    ip.parse(varargin{:})
    PATH2INPUT = ip.Results.PAMLAB_DATA_FOLDER;
    PATH2DATA = ip.Results.WAV_FILE_FOLDER;
    PATH2OUTPUTDIRECTORY = ip.Results.OUTPUT_FOLDER_LOCATION;
    PARAMFILE = ip.Results.PARAMFILE;
    
    % prompt user to specify I/O folder paths interactively if they were
    % not entered via the command line
    pathNotSpecified = @(pth) isnumeric(pth) && pth == 0;
    %%% PAMLAB data
    if isempty(PATH2INPUT)
        PATH2INPUT = uigetdir('','SELECT FOLDER WITH PAMLAB OUTPUT');
        if pathNotSpecified(PATH2INPUT)
            disp('No PAMLAB data folder specified; cancelling')
            return
        end
    end
    %%% WAV files
    if isempty(PATH2DATA)
        PATH2DATA = uigetdir('','SELECT FOLDER WITH WAV FILES');
        if pathNotSpecified(PATH2DATA)
            disp('No WAV file folder specified; cancelling')
            return
        end
    end
    %%% Output folder location
    if isempty(PATH2OUTPUTDIRECTORY)
        PATH2OUTPUTDIRECTORY = uigetdir('','SELECT DIRECTORY TO CREATE OUTPUT FOLDER');
        
        % if output directory location was not specified, set it to the
        % parent directory of the PAMLAB output
        if pathNotSpecified(PATH2OUTPUTDIRECTORY)
            [PATH2OUTPUTDIRECTORY, ~, ~] = fileparts(PATH2INPUT);
            disp(['No location specified for output folder; will use: ',PATH2OUTPUTDIRECTORY])
        end
        %%% OLD CODE
        %{
        if PATH2OUTPUTDIRECTORY == 0
            temp1 = split(PATH2INPUT,'\');
            temp2 = temp1(1:end-1)';
            PATH2OUTPUTDIRECTORY = fullfile(strjoin(temp2,'\'));
            clear(temp1,temp2)
        end
        %}
    end

    %Get list of Pamlab Output csv and path to wav files
    PAMLAB_ANNOTATIONS = dir(fullfile(PATH2INPUT, '\*.csv'));
    WAVFILES = struct2table(dir(fullfile(PATH2DATA, '**\*.wav')));
    numAnnotations = length(PAMLAB_ANNOTATIONS);
    
    % create output folder
    PATH2OUTPUT = fullfile(PATH2OUTPUTDIRECTORY,'SNR_OUTPUT');
    if ~exist(PATH2OUTPUT, 'dir')
       mkdir(PATH2OUTPUT)
    end

    % read parameter file
    %%% Look for default file if one was not specified in the command line
    if isempty(PARAMFILE)
        toolScriptPath = mfilename('fullpath');
        [tooldir, ~, ~] = fileparts(toolScriptPath);
        PARAMFILE = fullfile(tooldir, 'SNR_PARAMS.csv');
        disp('Using default parameter file')
    end
    
    SNR_PARAMS = readtable(PARAMFILE);
    paramnames = {'Species', 'Call_Type', 'Lower_Passband_Frequency', 'Upper_Passband_Frequency', 'Stopband_Rolloff_Bandwidth', 'Noise_Distance', 'Ideal_Noise_Duration', 'Signal_Energy_Percent'};
    assert(all(ismember(SNR_PARAMS.Properties.VariableNames, paramnames)), 'Parameter file does not include all expected parameters')
    specieslist = {SNR_PARAMS.Species};
    specieslist = unique(horzcat(specieslist{:}),'stable');
    sp_idx = listdlg('PromptString','Select a species:',...
                      'SelectionMode','single',...
                      'ListString',specieslist);
    species = specieslist(sp_idx,1);
    calltypelist = SNR_PARAMS.Call_Type(strcmp(SNR_PARAMS.Species,string(species{1,1})),:);
    calltypelist = {calltypelist};
    calltypelist = unique(horzcat(calltypelist{:}),'stable');
    %if length(calltypelist)== 1
    %    calltypelist = string(calltypelist);
    %end
    ct_idx = listdlg('PromptString','Select a call type:',...
                      'SelectionMode','single',...
                      'ListString',calltypelist);
    calltype = calltypelist(ct_idx,1);
    SNR_PARAMS_filtered = SNR_PARAMS(strcmp(SNR_PARAMS.Species,string(species{1,1})) & strcmp(SNR_PARAMS.Call_Type,string(calltype{1,1})),:);

    %%% extract or derive variables from PARAMS table
    LowerStopbandFreq = SNR_PARAMS_filtered.Lower_Passband_Frequency - SNR_PARAMS_filtered.Stopband_Rolloff_Bandwidth;
    LowerPassbandFreq = SNR_PARAMS_filtered.Lower_Passband_Frequency;
    UpperPassbandFreq = SNR_PARAMS_filtered.Upper_Passband_Frequency;
    UpperStopbandFreq = SNR_PARAMS_filtered.Upper_Passband_Frequency + SNR_PARAMS_filtered.Stopband_Rolloff_Bandwidth;
    NoiseDistance = SNR_PARAMS_filtered.Noise_Distance; 
    NoiseSize = SNR_PARAMS_filtered.Ideal_Noise_Duration;
    EnergyPercent = SNR_PARAMS_filtered.Signal_Energy_Percent;

    %%% create empty variable to store bandpass filter object
    bandpass_filter = [];
    
    %%% initialize output variable if output was requested
    if returnOutputVar
        out = struct();
    end

    %%% process each artefact file
    for p = 1:numAnnotations %read in in Pamlab csv (Loop) Possibly redundant...
        file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
        PLA = readtable(file);
        PLA.SNR_Direct = NaN(height(PLA),1); %create location to save SNR
        PLA.SNR_Corrected = NaN(height(PLA),1); %create location to save SNR with noise power subtracted from numerator
        PLA.SNRCalc_SignalDuration = NaN(height(PLA),1); %create location to save the energy-based signal duration used to calculate SNR
        PLA.SNRCalc_NoiseDuration = NaN(height(PLA),1); %create location to save the noise duration used to calculate SNR
        x = [];
        FileName =[];

        %%% initialize waitbar
        num_annotations = height(PLA);
        waitmsg = 'Processing PAMlab annotations...';
        waitfig = waitbar(0, waitmsg);
        tic

        %%% process each annotation
        for w = 1:num_annotations %Start rows loop
            %%% update waitbar
            t_elapsed = toc;
            t_rem = t_elapsed.*((num_annotations-(w-1))/(w-1));
            waitbar(w/num_annotations, waitfig, sprintf('%s\nEstimated time remaining: %s', waitmsg, duration(0,0,t_rem)))

            %%% get wav file and read it in, if not already done
            temp = PLA.filename(w);
            if isempty(x)||~strcmp(temp, FileName) %check if first time running
                FileName = temp;
                for i = 1:length(WAVFILES.name)
                    if contains(WAVFILES.name(i), FileName)
                       PATH2WAV = char(fullfile(WAVFILES.folder(i),WAVFILES.name(i)));
                       continue
                    end
                end
                if ~exist('PATH2WAV','var')
                    disp("File not found in directory") 
                    return
                end
                [x,Fs] = audioread(PATH2WAV);
                [M,q] = size(x); %get size length of audio
                dt = 1/Fs;      %time between samples in seconds
                t = dt*(0:M-1)';%get time index in seconds
                xt = [x t];       
                %%% create bandpass filter object if it doesn't exist already
                if isempty(bandpass_filter) || Fs ~= bandpass_filter.SampleRate
                    bandpass_filter = designfilt(...
                        'bandpassfir',...
                        'StopbandFrequency1', LowerStopbandFreq,...
                        'PassbandFrequency1', LowerPassbandFreq,...
                        'PassbandFrequency2', UpperPassbandFreq,...
                        'StopbandFrequency2', UpperStopbandFreq,...
                        'StopbandAttenuation1', 60,...
                        'StopbandAttenuation2', 60,...
                        'PassbandRipple', 1,...
                        'DesignMethod', 'kaiserwin',...
                        'SampleRate', Fs...
                        );
                end
                
                %%% define spectrogram parameters (TEMPORARILY HARDCODED)
                spectype = 'custom'; % CHANGE AS NEEDED
                switch spectype
                    case 'MERIDIAN'
                        %%%%%% based on MERIDIAN settings for NARW;
                        %%%%%% tends to be very hi-res for blue calls, and
                        %%%%%% slow to run/render
                        stftWinSize = 2^nextpow2(round(0.256*Fs));
                        stftOverlap = stftWinSize - round(0.032*Fs);
                        stftN = stftWinSize;
                    case 'longcalls'
                        %%%%%% based on PAMlab settings for long calls;
                        %%%%%% very coarse temporal resolution
                        stftWinSize = round(2*Fs);
                        stftOverlap = stftWinSize - round(0.5*Fs);
                        stftN_candidates = 2.^(nextpow2(Fs/0.4) + [-1,0]);
                        stftN = interp1(stftN_candidates, stftN_candidates, Fs/0.4, 'nearest');
                    case 'custom'
                        %%%%%% finer resolution than longcalls
                        stftWinSize = 2^nextpow2(round(Fs));
                        stftOverlap = stftWinSize - round(0.1*Fs);
                        stftN = stftWinSize;
                end
            end

            %%% Get Start of annotation and End of annotation
            PLA_Start = PLA.annotation_relative_start_time_sec(w);
            PLA_Stop = PLA.annotation_relative_end_time_sec(w);

            % Start90 = PLA_StartTime90 + RelativeStartTime;
            % End90 = PLA_StopTime90 + RelativeStartTime;

            %%% identify other annotations in same audio file and get their
            %%% start and stop times
            others_in_wav = strcmp(PLA.filename, FileName);
            others_in_wav(w) = false;
            PLA_Start_other = PLA.annotation_relative_start_time_sec(others_in_wav);
            PLA_Stop_other = PLA.annotation_relative_end_time_sec(others_in_wav);

            %%% extract bandpass-filtered signal and noise samples.
            %%% NOTE: the number of buffer samples at the ends of the clip
            %%% must be large enough to accommodate STFT windows for
            %%% generation of spectrograms and Welch spectra; thus, buffer
            %%% size is dependent on STFT parameters.
            [xClip, sigPos, noisePos, tClipStart] = snr.isolateFilteredSNClip(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter,...
                'NoiseDistance', NoiseDistance,...
                'IdealNoiseSize', NoiseSize,...
                'RemoveFromNoise', [PLA_Start_other,PLA_Stop_other],...
                'EnergyPercent', EnergyPercent,...
                'ClipBufferSize', (stftWinSize*0.75)/Fs ... % the factor of 75% is probably overkill, but will stick with this for now
                );
            %** old method
            %{
            [xSignal, xNoise] = snr.extractSN(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter,...
                'NoiseDistance', NoiseDistance,...
                'IdealNoiseSize', NoiseSize,...
                'RemoveFromNoise', [PLA_Start_other,PLA_Stop_other],...
                'EnergyPercent', EnergyPercent);
            %}
            
            %%% calculate SNR and other features
            %%% (leave NaN if not possible because signal is too close to 
            %%% endpoints)
            %if ~isempty(xSignal)
            if ~isempty(xClip)
                %%% extract signal and noise from clip.
                %%% Separated noise sections will be concatenated.
                xSignal = xClip(sigPos(1):sigPos(2));
                xNoise = [];
                for ii = 1:size(noisePos,1)
                    xNoise = [xNoise; xClip(noisePos(ii,1):noisePos(ii,2))];
                end
                
                %[PLA.SNR(w), PLA.SNR_Adjusted(w)] = snr.calculateSNR(xSignal, xNoise);
                [PLA.SNR_Direct(w), PLA.SNR_Corrected(w)] = snr.calculateSNR(xSignal, xNoise, 'CapNoise',true);
                PLA.SNRCalc_SignalDuration(w) = numel(xSignal)/Fs;
                PLA.SNRCalc_NoiseDuration(w) = numel(xNoise)/Fs;
                
                %%% get spectrogram and full PSD estimate of signal
                %%% (truncated to passband frequencies)
                [t_stft, f_stft, psdm, psd] = snr.computeSTFT(xClip, Fs, sigPos, stftWinSize, stftOverlap, 'NFFT',stftN, 'FRange',[LowerStopbandFreq,UpperStopbandFreq]);
                
                %%% smooth the spectrogram
                spec_smooth_kernel = [1,2,1; 2,4,2; 1,2,1];
                psdm_smooth = conv2(psdm, spec_smooth_kernel);
                psdm_smooth = psdm_smooth(2:(end-1), 2:(end-1));
                
                % get subset of spectrogram that includes the user-defined
                % frequency bounds
                f_min_ann = PLA.annotation_fmin_hz(w);
                f_max_ann = PLA.annotation_fmax_hz(w);
                is_f_in_range = f_stft >= f_min_ann & f_stft <= f_max_ann;
                psdm_trace = psdm_smooth(is_f_in_range,:);
                f_stft_trace = f_stft(is_f_in_range);
                
                %%% get noise threshold for pitch tracing
                psdm_anal = 10*log10(psdm_smooth);
                psdm_anal = psdm_anal - max(psdm_anal(:));
                th_psdm = prctile(psdm_anal(:),95);
                
                %%% find the best trace line
                %%% (hard-code trace parameters for now)
                %%% (also trying out different parameters)
                trace_penalty_coeff = 0.01; %0.008;
                trace_penalty_exp = 3; %2;

                t_trace = struct();
                f_trace = struct();
                try
                    % "final" set
                    [t_trace.unconstrained, f_trace.unconstrained] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp);
                    [t_trace.constrained, f_trace.constrained] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.333, 'MaxFreqGap',10);
                    
                    % Testing set
                    %[t_trace.unbound, f_trace.unbound] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp);
                    %[t_trace.bound, f_trace.bound] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm);
                    %[t_trace.gapcut_t033, f_trace.gapcut_t033] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.333);
                    %[t_trace.gapcut_t033_f10, f_trace.gapcut_t033_f10] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.333, 'MaxFreqGap',10);
                    %[t_trace.gapcut_t025_f10, f_trace.gapcut_t025_f10] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.25, 'MaxFreqGap',10);
                    %[t_trace.gapcut_t025_f5, f_trace.gapcut_t025_f5] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.25, 'MaxFreqGap',5);
                    %[t_trace.penalized, f_trace.penalized] = snr.getTraceLine(t_stft, f_stft_trace, psdm_trace, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',0.25, 'MaxFreqGap',5, 'GapPenalty',true);
                    
                    bad_trace = false;
                catch ME
                    warning('Failed to find a trace line for call No. %d:\n%s', w, ME.message)
                    bad_trace = true;
                end
                
                %** DEBUG
                %{
                testtrace_sizes = structfun(@numel, t_trace);
                if ~bad_trace && numel(unique(testtrace_sizes(1:end))) > 1
                    fig = figure(1);
                    ax = gca();
                    cla(ax);

                    patch_tshift = -mean(unique(diff(t_stft)))/2;
                    line_fshift = mean(unique(diff(f_stft)))/2;

                    surf(t_stft-mean(diff(t_stft))/2, f_stft, psdm_anal, 'EdgeColor','none')

                    axis(ax, 'xy');
                    view(ax, 0, 90);
                    ylim(ax,[LowerPassbandFreq,UpperPassbandFreq])
                    xlabel(ax, 'Time [s]')
                    ylabel(ax, 'Frequency [Hz]')
                    title(sprintf('Pathfinding Trace Line (Cost = %g\\it\\Deltaf\\rm^{%g} + 1)\nNo. %d', trace_penalty_coeff, trace_penalty_exp, w))

                    caxis(ax, [-60, -3])
                    set(ax, 'ColorScale', 'log')

                    hold on
                    
                    trace_fields = fieldnames(t_trace);
                    num_traces = numel(trace_fields);
                    lin = cell(1,num_traces);
                    %marktypes = {'s', 'o', 'x', '^', 'v', '*'};
                    marktypes = {'wx', 'ro'};
                    for ii = 1:num_traces
                        field_ii = trace_fields{ii};
                        lin_ii = plot3(t_trace.(field_ii)', f_trace.(field_ii)+line_fshift, ones(size(f_trace.(field_ii))), [marktypes{ii},'-'], 'LineWidth',1.5, 'DisplayName',field_ii);
                        lin{ii} = lin_ii;
                    end
                    plot3(t_stft([1,end]), [f_min_ann,f_min_ann]+line_fshift, [1,1], 'w:')
                    plot3(t_stft([1,end]), [f_max_ann,f_max_ann]+line_fshift, [1,1], 'w:')
                    
                    %%% draw a partially transparent plane corresponding to
                    %%% the trace clipping threshold
                    %{
                    patch(...
                        t_stft([1,1,end,end,1]) + patch_tshift,...
                        f_stft([1,end,end,1,1]),...
                        repelem(th_psdm,1,5),...
                        'r',...
                        'FaceAlpha', 0.25,...
                        'EdgeColor', 'none')
                    %}
                    
                    legend(ax, [lin{:}], 'Interpreter','none')
                    
                    keyboard
                end
                %}
                
            end

            %%%
            %pass: raw wav,Start90, End90,[frequency band],buffer size,and noiseDistance to BP_clip.m
            %output: bandpassed wav clip + buffer
            %BP_clip = snr.BP_clip(x,Start90,End90,Freq_band,NoiseDistance,BP_buffer_samples);
            %%%
            %%%
            %pass: start90, stop90, noiseDistance,and bandpassed wav clip + buffer to extractSN.m
            %output: signal clip and noise clip
            %Temporary function test:
            %[xSignal, xNoise] = snr.extractSN(x,Fs,Start90,End90,NoiseDistance,Units);
            %[xSignal, xNoise] = snr.extractSN(x, fs, sigStart, sigStop, noiseDist, units)
            %
            %%%
            %%%
            %pass: signal clip and noise clip to calculateSNR.m
            %output: SNR

            %%%
        end %call loop
        close(waitfig)
        
        % create output file
        temp_name = split(PAMLAB_ANNOTATIONS(p).name,'.');
        temp_filename = [char(temp_name(1)) '_SNR.csv'];
        final_filename = generateUniqueName(PATH2OUTPUT,temp_filename);
        PATH2OUTPUT_FILENAME = fullfile(PATH2OUTPUT,final_filename);
        writetable(PLA,PATH2OUTPUT_FILENAME);
        
        % save table to output variable if needed
        if returnOutputVar
            outFieldName = temp_name{1};
            outFieldNameAlt = sprintf('Annotations_%d',p);
            try
                out.(outFieldName) = PLA;
            catch
                warning('"%s" is not a valid field name for the output variable; it will be replaced with "%s".', outFieldName, outFieldNameAlt)
                out.(outFieldNameAlt) = PLA;
            end
        end
        
    end % end PAMLAB annotations loop
    
    % return the output variable if requested
    if returnOutputVar
        varargout{1} = out;
    end
    
end % end main function

            
%% generateUniqueName
function newFileName = generateUniqueName(dirPath, fileName)
%
% Checks if a file in a directory exists, and proposes a unique alternative
% filename if it does. Use this when you want to avoid overwriting files or
% folders that may already exist. The new proposed names consist of the
% original plus a unique (incremental) integer appended at the end,
% starting with 2.
% 
% If the input file name does not yet exist, this function will return the
% original name.
%
%
% Written by Wilfried Beslin
% Last updated 2024-01-10, using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some possible options here:
% - Allow zero padding
% - Choose whether or not to force an integer for names that are NOT taken
% - Choose which number to start with
% - specify a delimiter string (default is underscore)
% - Maybe: support letters rather than numbers

    % get full path
    filePath = fullfile(dirPath,fileName);

    % check if file already exists
    if logical(exist(filePath,'file')) % also includes folders
        
        % if so, generate unique name
        [~,fileNameBase,fileExt] = fileparts(fileName);
        fileNum = 2;
        haveUnique = false;
        while ~haveUnique
            newFileName = [fileNameBase,'_',num2str(fileNum),fileExt];
            newFilePath = fullfile(dirPath,newFileName);
            if logical(exist(newFilePath,'file'))
                fileNum = fileNum + 1;
            else
                haveUnique = true;
            end
        end
        
    else
        % if file doesn't exist, then keep old name
        newFileName = fileName;
    end

end
