function varargout = MUPPET(varargin)
%
% MUPPET.m
%
% Process Pamlab output for use with MUPPET. Runs MUPPET code to extract
% call parameters from Pamlab annotations.
%
% SYNTAX:
%   MUPPET
%   MUPPET(Name,Value)
%   out = MUPPET(__)
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
%   "PARAMFILE" - Path to text file of parameters. If not specified, user
%       will be prompted to select a file manually.
%   .......................................................................
%   "PLOT_TRACE_LINES" - true/false value that determines whether or not to
%       save images of trace line plots for each annotation in the output
%       folder. Default is false.
%   .......................................................................
%
% OPTIONAL OUTPUT ARGUMENTS:
% These arguments are only returned when requested in the command window.
% MUPPET already saves all output as files, so these arguments are mostly
% useful for debugging or otherwise working with the results directly in
% the MATLAB workspace. Output arguments all take the form of structs,
% where each field corresponds to one PAMLab annotation file.
%   .......................................................................
%   "out1" = data from the first MUPPET output file, i.e., expanded
%   versions of the PAMLab annotation tables that include SNR information
%   at the end.
%   .......................................................................
%   "out2" - data from the second MUPPET output file, i.e., tables of call
%   parameters.
%   .......................................................................
%   "out3" - data from the third MUPPET output file, i.e., tables of
%   concatenated trace line points.
%   .......................................................................
%
% Written by Mike Adams and Wilfried Beslin
% Last updated by Wilfried Beslin
% 2024-12-16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTE: https://www.mathworks.com/help/matlab/ref/listdlg.html

    % check if output arguments were requested
    nargoutchk(0,3) % allows the function to return up to three optional output parameters

    % process optional input parameters (I/O paths)
    ip = inputParser();
    ip.addParameter('PAMLAB_DATA_FOLDER', '', @isfolder)
    ip.addParameter('WAV_FILE_FOLDER', '', @isfolder)
    ip.addParameter('OUTPUT_FOLDER_LOCATION', '', @isfolder)
    ip.addParameter('PARAMFILE', '', @isfile)
    ip.addParameter('PLOT_TRACE_LINES', false, @(v)validateattributes(v,{'logical'},{'scalar'}))
    ip.addParameter('EXPORT_SPECTROGRAMS', false, @(v)validateattributes(v,{'logical'},{'scalar'})) % for debugging
    ip.parse(varargin{:})
    PATH2INPUT = ip.Results.PAMLAB_DATA_FOLDER;
    PATH2DATA = ip.Results.WAV_FILE_FOLDER;
    PATH2OUTPUTDIRECTORY = ip.Results.OUTPUT_FOLDER_LOCATION;
    PARAMFILE = ip.Results.PARAMFILE;
    plot_trace = ip.Results.PLOT_TRACE_LINES;
    export_specs = ip.Results.EXPORT_SPECTROGRAMS;
    
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
    %PATH2OUTPUT = fullfile(PATH2OUTPUTDIRECTORY,'SNR_OUTPUT');
    [~,dataDirName] = fileparts(PATH2INPUT);
    PATH2OUTPUT = fullfile(PATH2OUTPUTDIRECTORY, generateUniqueName(PATH2OUTPUTDIRECTORY,['MUPPET_OUTPUT_',dataDirName]));
    if ~exist(PATH2OUTPUT, 'dir')
       mkdir(PATH2OUTPUT)
    end

    % read parameter file
    %%% prompt user for file if one was not specified in the command line
    if isempty(PARAMFILE)
        [paramfile_name, paramfile_dir] = uigetfile('*.txt', 'SELECT PARAMETER FILE');
        PARAMFILE = fullfile(paramfile_dir, paramfile_name);
        %toolScriptPath = mfilename('fullpath');
        %[tooldir, ~, ~] = fileparts(toolScriptPath);
        %PARAMFILE = fullfile(tooldir, 'SNR_PARAMS.csv');
        %disp('Using default parameter file')
    end
    PARAMS = MUPPET.importInputParams(PARAMFILE);
    
    %%% OLD CSV PARAMFILE PARSING
    %{
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
    %}

    % extract or derive parameters from the imported PARAMS struct, to
    % avoid having to access the struct multiple times
    %%% Bandpass filter settings
    LowerStopbandFreq = PARAMS.Lower_Passband_Frequency - PARAMS.Stopband_Rolloff_Bandwidth;
    LowerPassbandFreq = PARAMS.Lower_Passband_Frequency;
    UpperPassbandFreq = PARAMS.Upper_Passband_Frequency;
    UpperStopbandFreq = PARAMS.Upper_Passband_Frequency + PARAMS.Stopband_Rolloff_Bandwidth;
    %%% Signal and noise isolation
    NoiseDistance = PARAMS.Noise_Distance;
    NoiseSize = PARAMS.Ideal_Noise_Duration;
    EnergyPercent = PARAMS.Signal_Energy_Percent;
    %%% Resampling
    FsResampled = PARAMS.Downsampled_Sampling_Rate;
    %%% STFT parameters
    stftWinSize = PARAMS.STFT_Win_Size;
    stftOverlap = stftWinSize - PARAMS.STFT_Step_Size;
    stftN = PARAMS.NFFT;
    smoothSpec = PARAMS.Smooth_Spec;
    %%% Trace line calculation
    trace_penalty_sigma = PARAMS.Trace_Penalty_Sigma;
    trace_penalty_at_sigma = PARAMS.Trace_Penalty_At_Sigma;
    trace_energy_percent = PARAMS.Trace_Energy_Percent;
    trace_th_type = validatestring(lower(PARAMS.Trace_Threshold_Type), {'fixed','percentile','xnoisesd'});
    trace_th_val = PARAMS.Trace_Threshold_Val;
    %%% Trace line plot parameters
    colmap = lower(PARAMS.Spec_Plot_Colour_Map);
    log_specplot_cols = PARAMS.Log_Spec_Plot_Colours;
    trace_plot_data = struct(...
        'Color', {PARAMS.Trace_Plot_Line_Colour},...
        'Marker', {PARAMS.Trace_Plot_Marker_Type},...
        'MarkerSize', {PARAMS.Trace_Plot_Marker_Size},...
        'LineWidth', {PARAMS.Trace_Plot_Line_Width}...
        );
    %%% Advanced
    cap_noise = PARAMS.Cap_Noise;

    % create empty variable to store bandpass filter object
    bandpass_filter = [];
    
    % initialize output variables if output was requested
    varargout = cell(1,nargout);
    for outvarnum = 1:nargout
        varargout{outvarnum} = struct();
    end
    
    % set trace line table variable names
    trace_table_vars = {'RelTime','Freq','RelPow_dB','Call_ID'};

    % process each artefact file
    %** Still need to apply error-checking for unsupported CSV files
    for p = 1:numAnnotations %read in in Pamlab csv (Loop) Possibly redundant...
        file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
        PLA = readtable(file);
        PLA.Call_ID = [1:height(PLA)]';
        PLA.SNR_Direct = NaN(height(PLA),1); %create location to save SNR
        PLA.SNR_Corrected = NaN(height(PLA),1); %create location to save SNR with noise power subtracted from numerator
        PLA.SNRCalc_SignalDuration = NaN(height(PLA),1); %create location to save the energy-based signal duration used to calculate SNR
        PLA.SNRCalc_NoiseDuration = NaN(height(PLA),1); %create location to save the noise duration used to calculate SNR
        x = [];
        FileName =[];
        
        outfile_refname = erase(PAMLAB_ANNOTATIONS(p).name, '.csv');
        
        %%% initialize empty table to store features
        call_params = table();
        
        %%% initialize empty table to store trace lines
        trace_lines = table();
        
        %%% set up plot folder, if relevant
        if plot_trace
            PATH2OUTPUT_TRACEPLOTS = fullfile(PATH2OUTPUT, [outfile_refname,'_TracePlots']);
            mkdir(PATH2OUTPUT_TRACEPLOTS)
        end
        
        %%% set up spectrogram export folder if relevant
        if export_specs
            PATH2OUTPUT_SPECS = fullfile(PATH2OUTPUT, [outfile_refname,'_Spectrograms']);
            mkdir(PATH2OUTPUT_SPECS)
        end

        %%% initialize waitbar
        num_annotations = height(PLA);
        waitmsg = 'Processing PAMlab annotations...';
        waitfig = waitbar(0, waitmsg);
        tic

        %%% process each annotation
        for w = 1:num_annotations %Start rows loop
            if PLA.Species(w) == "NN"
                continue
            end
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
                
                %%% define spectrogram and downsampling parameters 
                %%% (TEMPORARILY HARDCODED)
                %{
                FsResampled = 8000;
                spectype = 'custom'; % CHANGE AS NEEDED
                switch spectype
                    case 'MERIDIAN'
                        %%%%%% based on MERIDIAN settings for NARW;
                        %%%%%% tends to be very hi-res for blue calls, and
                        %%%%%% slow to run/render
                        stftWinSize = 2^nextpow2(round(0.256*FsResampled));
                        stftOverlap = stftWinSize - round(0.032*FsResampled);
                        stftN = stftWinSize;
                    case 'longcalls'
                        %%%%%% based on PAMlab settings for long calls;
                        %%%%%% very coarse temporal resolution
                        stftWinSize = round(2*FsResampled);
                        stftOverlap = stftWinSize - round(0.5*FsResampled);
                        stftN_candidates = 2.^(nextpow2(FsResampled/0.4) + [-1,0]);
                        stftN = interp1(stftN_candidates, stftN_candidates, FsResampled/0.4, 'nearest');
                    case 'custom'
                        %%%%%% finer resolution than longcalls
                        stftWinSize = 2^nextpow2(round(FsResampled));
                        stftOverlap = stftWinSize - round(0.1*FsResampled);
                        stftN = stftWinSize;
                end
                %}
            end

            %%% Get Start of annotation and End of annotation
            PLA_Start = PLA.LeftTime_sec_(w);
            PLA_Stop = PLA.RightTime_sec_(w);

            % Start90 = PLA_StartTime90 + RelativeStartTime;
            % End90 = PLA_StopTime90 + RelativeStartTime;

            %%% identify other annotations in same audio file and get their
            %%% start and stop times
            others_in_wav = strcmp(PLA.filename, FileName);
            others_in_wav(w) = false;
            PLA_Start_other = PLA.LeftTime_sec_(others_in_wav);
            PLA_Stop_other = PLA.RightTime_sec_(others_in_wav);

            %%% extract bandpass-filtered signal and noise samples.
            %%% NOTE: the number of buffer samples at the ends of the clip
            %%% must be large enough to accommodate STFT windows for
            %%% generation of spectrograms and Welch spectra; thus, buffer
            %%% size is dependent on STFT parameters.
            [xClip, sigPos, noisePos, annotPos, tClipStart] = MUPPET.getProcessedSNClip(x, Fs, FsResampled, [PLA_Start,PLA_Stop], bandpass_filter,...
                'NoiseDistance', NoiseDistance,...
                'IdealNoiseSize', NoiseSize,...
                'RemoveFromNoise', [PLA_Start_other,PLA_Stop_other],...
                'EnergyPercent', EnergyPercent,...
                'ClipBufferSize', (stftWinSize*0.75)/FsResampled ... % the factor of 75% is probably overkill, but will stick with this for now
                );
            
            %%% calculate SNR and other features
            %%% (leave NaN if not possible because signal is too close to 
            %%% endpoints)
            %if ~isempty(xSignal)
            if ~isempty(xClip)
                %%% extract samples of signal and noise from clip.
                %%% Separated noise sections will be concatenated.
                xSignal = xClip(sigPos(1):sigPos(2));
                xNoise = [];
                for ii = 1:size(noisePos,1)
                    xNoise = [xNoise; xClip(noisePos(ii,1):noisePos(ii,2))];
                end
                %xAnnot = xClip(annotPos(1):annotPos(2));
                
                [PLA.SNR_Direct(w), PLA.SNR_Corrected(w)] = MUPPET.calculateSNR(xSignal, xNoise, 'CapNoise',cap_noise);
                PLA.SNRCalc_SignalDuration(w) = numel(xSignal)/FsResampled;
                PLA.SNRCalc_NoiseDuration(w) = numel(xNoise)/FsResampled;
                
                %%% get spectrogram and full PSD estimate of signal
                %%% (truncated to passband frequencies)
                [t_stft, f_stft, psdm, psd] = MUPPET.computeSTFT(xClip, FsResampled, annotPos, stftWinSize, stftOverlap, 'NFFT',stftN, 'FRange',[LowerStopbandFreq,UpperStopbandFreq]);
                
                %%% do the same for noise
                %%% Here I am calculating different spectrograms for every
                %%% noise piece seperately, rather than piecing the noise
                %%% together. Why? Because piecing disjointed noise
                %%% sections together may affect the frequency spectrum
                %%% (think of sudden spikes causing broadband impulses).
                %%% NOTE: for now, there is no special processing to fix
                %%% cases where there are not enough samples at the ends to
                %%% generate a noise spectrum - those cases are simply
                %%% removed. But it might be good to try shrinking the
                %%% noise windows in the future if this happens often.
                num_noiseparts = size(noisePos,1);
                psdmc_noise = cell(1, num_noiseparts);
                for ii = 1:num_noiseparts
                    try
                        [~, ~, psdm_noise_ii, ~] = MUPPET.computeSTFT(xClip, FsResampled, [noisePos(ii,1), noisePos(ii,2)], stftWinSize, stftOverlap, 'NFFT',stftN, 'FRange',[LowerStopbandFreq,UpperStopbandFreq]);
                        psdmc_noise{ii} = psdm_noise_ii;
                    catch
                        warning('Annotation %d: Failed to create noise spectrogram for noise part %d of %d', w, ii, num_noiseparts)
                    end
                end
                
                %%% set spectrogram processing order 
                if smoothSpec
                    spec_ops = {'log','smooth','denoise1'}; % this is the order used in LFDCS
                else
                    spec_ops = {'log','denoise1'};
                end
                
                %%% process signal and noise spectrograms
                [psdm_anal, psdm_noise_anal] = MUPPET.processSpec(psdm, psdmc_noise, spec_ops);
                psdm_noise_anal_demeaned = psdm_noise_anal - mean(psdm_noise_anal, 2);
                
                % get subset of signal spectrogram that includes the user-defined
                % frequency bounds
                f_min_ann = PLA.BottomFreq_Hz_(w);
                f_max_ann = PLA.TopFreq_Hz_(w);
                is_f_in_annot_range = f_stft >= f_min_ann & f_stft <= f_max_ann;
                psdm_anal_annwin = psdm_anal(is_f_in_annot_range,:);
                f_stft_annwin = f_stft(is_f_in_annot_range);
                
                %%% export raw (unprocessed) spectrograms, if requested
                %** This is only meant to be used for debugging
                if export_specs
                    PATH2OUTPUT_SPEC_FILE = fullfile(PATH2OUTPUT_SPECS, [outfile_refname,'_spectrogram_',num2str(w),'.mat']);
                    save(PATH2OUTPUT_SPEC_FILE, 'f_stft', 't_stft', 'psdm', 'psdmc_noise', 'is_f_in_annot_range');
                end
                
                % set power-over-noise threshold for pitch tracing
                switch trace_th_type
                    case 'fixed'
                        %%% set threshold according to a fixed dB value
                        trace_power_th = trace_th_val;
                    case 'percentile'
                        %%% set threshold based on a certain percentile
                        %%% of all values in the signal spectrogram
                        trace_power_th = prctile(psdm_anal_annwin(:), trace_th_val);
                    case 'xnoisesd'
                        %%% set threshold based on a certain number of
                        %%% standard deviations above the noise spectrogram
                        trace_power_th = trace_th_val*std(psdm_noise_anal_demeaned(:));
                    otherwise
                        error(sprintf('Unrecognized trace threshold type "%s"', trace_th_type))
                end
                
                %%% find the best trace line
                try
                    [t_trace, f_trace] = MUPPET.getTraceLine(t_stft, f_stft_annwin, psdm_anal_annwin, 'PenaltySigma',trace_penalty_sigma, 'PenaltyAtSigma',trace_penalty_at_sigma, 'PowerThreshold',trace_power_th, 'EnergyPercent',trace_energy_percent);
                    
                     %** parameters for alternate trace lines
                     alt_trace_args = {};
                     %** CHANGE AS NEEDED *******************************
                     %alt_trace_args = [alt_trace_args, {{'PenaltySigma',trace_penalty_sigma, 'PenaltyAtSigma',trace_penalty_at_sigma, 'PowerThreshold',3.5*std(psdm_noise_anal_demeaned(:)), 'EnergyPercent',trace_energy_percent}}];
                     %alt_trace_args = [alt_trace_args, {{'PenaltySigma',trace_penalty_sigma, 'PenaltyAtSigma',trace_penalty_at_sigma, 'PowerThreshold',4*std(psdm_noise_anal_demeaned(:)), 'EnergyPercent',trace_energy_percent}}];
                     %alt_trace_args = [alt_trace_args, {{'PenaltySigma',trace_penalty_sigma, 'PenaltyAtSigma',trace_penalty_at_sigma, 'PowerThreshold',trace_power_th, 'EnergyPercent',90}}];
                     %** END CHANGE *************************************

                     % compile all trace line types
                     %%% (data will only be extracted from the first)
                     t_trace_all = {t_trace};
                     f_trace_all = {f_trace};
                     for ii = 1:numel(alt_trace_args)
                         [t_trace_ii, f_trace_ii] = MUPPET.getTraceLine(t_stft, f_stft_annwin, psdm_anal_annwin, alt_trace_args{ii}{:});
                         t_trace_all = [t_trace_all, {t_trace_ii}];
                         f_trace_all = [f_trace_all, {f_trace_ii}];
                     end
                    
                    %%% override plot data if there are multiple trace
                    %%% lines
                    if numel(t_trace_all) > 1
                        %%% CHANGE AS NEEDED
                        trace_plot_data = struct(...
                            'Color', {'r','c'},...
                            'Marker', {'x','o'},...
                            'MarkerSize', {4,4},...
                            'LineWidth', {0.5,0.5},...
                            'DisplayName', {'E100','E90'}...
                            );
                    end
                    
                    %%% get relative power values of the trace line
                    %%% (using the processed spectrogram)
                    %psdm_tracedata = 10.*log10(psdm(is_f_in_annot_range,:));
                    [~, i_f_trace] = ismember(f_trace, f_stft_annwin);
                    [~, i_t_trace] = ismember(t_trace, t_stft);
                    logp_trace = psdm_anal_annwin(sub2ind(size(psdm_anal_annwin), i_f_trace, i_t_trace'));
                    
                    %%% store trace data in table
                    trace_line_w = table(t_trace', f_trace,  logp_trace, repelem(w,numel(t_trace),1), 'VariableNames',trace_table_vars);
                    trace_lines = [trace_lines; trace_line_w];
                    
                catch ME
                    %%% issue warning if unable to build trace line
                    warning('Failed to find a trace line for call No. %d:\n%s', w, ME.message)
                    trace_line_w = [];
                end
                
                %%% plot trace line if specified
                if plot_trace && ~isempty(trace_line_w)
                    if log_specplot_cols
                        % process LOG colour scale
                        plot_zshift = trace_power_th*1.5 + 1;
                        psdm_plot = psdm_anal - plot_zshift;
                    
                        caxis_val = trace_power_th.*[-3,1.5] - plot_zshift;
                        log_specplot_cols = true;
                    else
                        % process LINEAR colour scale
                        psdm_plot = psdm_anal;
                    
                        caxis_val = trace_power_th.*[-2,2];
                        log_specplot_cols = false;
                    end
                    
                    % setup figure and axes
                    fig = gcf();
                    fig.Visible = 'off';
                    clf(fig);
                    ax = axes();

                    % make plot
                    MUPPET.plotTraceLine(ax, t_stft, f_stft, psdm_plot, t_trace_all, f_trace_all, [f_min_ann,f_max_ann], 'LineData',trace_plot_data, 'CAxis',caxis_val, 'LogCols',log_specplot_cols);
                    ylim(ax,[LowerPassbandFreq,UpperPassbandFreq])
                    xlim(ax,t_stft([1,end]))

                    % add lines representing fullband %-energy limits
                    tSigStart_re_annot = (sigPos(1) - annotPos(1))./FsResampled;
                    tSigStop_re_annot = (sigPos(2) - annotPos(1))./FsResampled;
                    fPlotMin = ax.YLim(1);
                    fPlotMax = ax.YLim(2);
                    plot3(ax, [tSigStart_re_annot,tSigStart_re_annot,NaN,tSigStop_re_annot,tSigStop_re_annot], [fPlotMin,fPlotMax,NaN,fPlotMin,fPlotMax], repelem(ax.ZLim(2),5), 'w:')

                    % update colormap
                    try
                        colormap(ax, colmap);
                    catch ME
                        warning(sprintf('Invalid colourmap "%s"; will use "parula" instead',colmap))
                        colormap(ax, 'parula');
                    end

                    % set plot title
                    datasetname_printf = strrep(outfile_refname, '_', '\_');
                    wavname_printf = strrep(PLA.filename{w}, '_', '\_');
                    plot_title = sprintf('\\bf %s - Call ID = %d\n\\rm %s, t = %.2f s', datasetname_printf, w, wavname_printf, PLA_Start);
                    title(ax, plot_title)

                    % save the plot
                    PATH2OUTPUT_TRACEPLOT_FILE = fullfile(PATH2OUTPUT_TRACEPLOTS, ['Trace_',num2str(w),'.png']);
                    figRes = 150;
                    figDims = [1440, 810];
                    fig.PaperUnits = 'inches';
                    fig.PaperPosition = [0, 0, figDims/figRes];
                    fig.PaperPositionMode = 'manual';
                    print(fig, PATH2OUTPUT_TRACEPLOT_FILE, '-dpng', ['-r',num2str(figRes)])
                else
                    fig = [];
                end
                
                %%% extract call parameters and store to table
                sigPosRel = sigPos - annotPos(1) + 1; % energy-based signal position relative to start of annotation
                call_params_w = MUPPET.extractCallParams(sigPosRel, FsResampled, f_stft, psd, EnergyPercent, trace_line_w);
            else
                % cannot calculate parameters for unsuitable calls, so just
                % pass NaNs to parameter extraction routine
                call_params_w = MUPPET.extractCallParams(NaN, NaN, NaN, NaN, NaN, NaN);
            end
            
            %%% add call parameters to running table
            call_params_w.Call_ID = w; %keep track of call ID
            call_params = [call_params; struct2table(call_params_w)];
            
        end %call loop
        if ~isempty(fig)
            close(fig);
        end
        close(waitfig)
        
        % Create output files
        
        %%% SNR-expanded annotations output file
        %remove additional noise annotations from SNR output
        PLA_calls = PLA(~ismember(PLA.Species,'NN'),:);
        out1_filename = [outfile_refname, '_SNR.csv'];
        PATH2OUTPUT_SNR_FILE = fullfile(PATH2OUTPUT, out1_filename);
        writetable(PLA_calls, PATH2OUTPUT_SNR_FILE);
        %%% OLD CODE
        %temp_name = split(PAMLAB_ANNOTATIONS(p).name,'.');
        %temp_filename = [char(temp_name(1)) '_SNR.csv'];
        %final_filename = generateUniqueName(PATH2OUTPUT,temp_filename);
        %PATH2OUTPUT_FILENAME = fullfile(PATH2OUTPUT,final_filename);
        %writetable(PLA,PATH2OUTPUT_FILENAME);
        
        %%% call parameter output file
        out2_filename = [outfile_refname, '_CallParams.csv'];
        PATH2OUTPUT_CALLPARAMS_FILE = fullfile(PATH2OUTPUT, out2_filename);
        writetable(call_params, PATH2OUTPUT_CALLPARAMS_FILE);
        
        %%% trace line output table
        out3_filename = [outfile_refname, '_TraceLines.csv'];
        PATH2OUTPUT_TRACE_FILE = fullfile(PATH2OUTPUT, out3_filename);
        writetable(trace_lines, PATH2OUTPUT_TRACE_FILE);

        %%% paramfile
        [~,paramfile_name,paramfile_ext] = fileparts(PARAMFILE);
        PATH2OUTPUT_PARAMFILE = fullfile(PATH2OUTPUT, [paramfile_name,'.',paramfile_ext]);
        copyfile(PARAMFILE, PATH2OUTPUT_PARAMFILE);
        
        % save tables to output variables if needed
        if nargout > 0
            outFieldName_primary = outfile_refname;
            outFieldName_secondary = sprintf('Annotations_%d',p);
            outFieldName = outFieldName_primary;
            for outvarnum = 1:nargout
                switch outvarnum
                    case 1
                        outvar = PLA_calls;
                    case 2
                        outvar = call_params;
                    case 3
                        outvar = trace_lines;
                    otherwise
                        warning('Ignoring extra output arguments')
                        break
                end
                try
                    varargout{outvarnum}.(outFieldName) = outvar;
                catch
                    warning('"%s" is not a valid field name for the output variable(s); it will be replaced with "%s".', outFieldName, outFieldName_secondary)
                    outFieldName = outFieldName_secondary;
                    varargout{outvarnum}.(outFieldName) = outvar;
                end
            end
        end
        
    end % end PAMLAB annotations loop
    
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
