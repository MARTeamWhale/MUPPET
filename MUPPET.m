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
%   "PARAMFILE" - Path to CSV file of SNR parameters. If not specified, the
%       tool will try to load the default parameters from the file
%       SNR_PARAMS.csv.
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
% Written by Mike Adams
% Last updated by Wilfried Beslin
% 2024-09-04
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
    ip.parse(varargin{:})
    PATH2INPUT = ip.Results.PAMLAB_DATA_FOLDER;
    PATH2DATA = ip.Results.WAV_FILE_FOLDER;
    PATH2OUTPUTDIRECTORY = ip.Results.OUTPUT_FOLDER_LOCATION;
    PARAMFILE = ip.Results.PARAMFILE;
    plot_trace = ip.Results.PLOT_TRACE_LINES;
    
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
    
    %%% initialize output variables if output was requested
    varargout = cell(1,nargout);
    for outvarnum = 1:nargout
        varargout{outvarnum} = struct();
    end
    
    %%% set trace line table variable names
    trace_table_vars = {'RelTime','Freq','Power_dB','CallNumber'};

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
            [xClip, sigPos, noisePos, annotPos, tClipStart] = MUPPET.isolateFilteredSNClip(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter,...
                'NoiseDistance', NoiseDistance,...
                'IdealNoiseSize', NoiseSize,...
                'RemoveFromNoise', [PLA_Start_other,PLA_Stop_other],...
                'EnergyPercent', EnergyPercent,...
                'ClipBufferSize', (stftWinSize*0.75)/Fs ... % the factor of 75% is probably overkill, but will stick with this for now
                );
            %** old method
            %{
            [xSignal, xNoise] = MUPPET.extractSN(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter,...
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
                %%% extract samples of signal and noise from clip.
                %%% Separated noise sections will be concatenated.
                xSignal = xClip(sigPos(1):sigPos(2));
                xNoise = [];
                for ii = 1:size(noisePos,1)
                    xNoise = [xNoise; xClip(noisePos(ii,1):noisePos(ii,2))];
                end
                %xAnnot = xClip(annotPos(1):annotPos(2));
                
                %[PLA.SNR(w), PLA.SNR_Adjusted(w)] = MUPPET.calculateSNR(xSignal, xNoise);
                [PLA.SNR_Direct(w), PLA.SNR_Corrected(w)] = MUPPET.calculateSNR(xSignal, xNoise, 'CapNoise',true);
                PLA.SNRCalc_SignalDuration(w) = numel(xSignal)/Fs;
                PLA.SNRCalc_NoiseDuration(w) = numel(xNoise)/Fs;
                
                %%% get spectrogram and full PSD estimate of signal
                %%% (truncated to passband frequencies)
                [t_stft, f_stft, psdm, psd] = MUPPET.computeSTFT(xClip, Fs, sigPos, stftWinSize, stftOverlap, 'NFFT',stftN, 'FRange',[LowerStopbandFreq,UpperStopbandFreq]);
                dt_stft = median(diff(t_stft));
                
                %%% convert spectrogram to log scale **(NEW)**
                logpsdm = 10.*log10(psdm);
                
                %%% smooth the spectrogram
                logpsdm_smooth = MUPPET.smoothSpec(logpsdm);
                % OLD METHOD
                %spec_smooth_kernel = [1,2,1; 2,4,2; 1,2,1];
                %logpsdm_smooth = conv2(logpsdm, spec_smooth_kernel);
                %logpsdm_smooth = logpsdm_smooth(2:(end-1), 2:(end-1));
                
                %** NEW: get noise spectrogram
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
                logpsdm_noise_smooth_cell = cell(1, num_noiseparts);
                for ii = 1:num_noiseparts
                    try
                        [~, ~, psdm_noise_ii, ~] = MUPPET.computeSTFT(xClip, Fs, [noisePos(ii,1), noisePos(ii,2)], stftWinSize, stftOverlap, 'NFFT',stftN, 'FRange',[LowerStopbandFreq,UpperStopbandFreq]);
                        logpsdm_noise_ii = 10.*log10(psdm_noise_ii);
                        logpsdm_noise_smooth_ii = MUPPET.smoothSpec(logpsdm_noise_ii);
                        logpsdm_noise_smooth_cell{ii} = logpsdm_noise_smooth_ii;
                    catch
                        warning('Annotation %d: Failed to create noise spectrogram for noise part %d of %d', w, ii, num_noiseparts)
                    end
                end
                % remove any noise sections that could not be processed
                logpsdm_noise_smooth_cell = logpsdm_noise_smooth_cell(~cellfun('isempty',logpsdm_noise_smooth_cell));
                if isempty(logpsdm_noise_smooth_cell)
                    % Handle cases where there is absolutely no noise data
                    % available
                    %%% For now, this is done by getting a noise estimate
                    %%% straight from the signal spectrogram instead. But
                    %%% this is not ideal because it can result in a few
                    %%% calls using a different method than all the others.
                    %%% Still, there are not many options. The only other
                    %%% real solution I can think of is to reject the calls
                    %%% outright, which seems excessive. I think the best
                    %%% improvement that can be done would be to reduce the
                    %%% number of cases which run into this issue as much
                    %%% as possible, e.g. by getting better noise windows
                    %%% from the rear if there aren't enough samples in
                    %%% front.
                    warning('No spectral noise data available for this annotation; will estimate noise from the signal spectrogram instead.')
                    logpsdm_noise_smooth = median(logpsdm_smooth,2);
                else
                    % concatenate the noise spectrograms
                    logpsdm_noise_smooth = horzcat(logpsdm_noise_smooth_cell{:});
                end
                
                % create a denoised version of the signal spectrogram that
                % shifts each frequency band such that 0 dB corresponds to
                % the average noise level in that band
                logpsdm_denoised = logpsdm_smooth - mean(logpsdm_noise_smooth,2);
                
                % get subset of signal spectrogram that includes the user-defined
                % frequency bounds
                f_min_ann = PLA.annotation_fmin_hz(w);
                f_max_ann = PLA.annotation_fmax_hz(w);
                is_f_in_annot_range = f_stft >= f_min_ann & f_stft <= f_max_ann;
                logpsdm_denoised_annwin = logpsdm_denoised(is_f_in_annot_range,:);
                f_stft_annwin = f_stft(is_f_in_annot_range);
                
                % set noise threshold for pitch tracing
                spec_snr_th = 10;
                %** OLD METHOD
                %{
                %%% (hard-coded to 95th percentile for now)
                psdm_anal = 10*log10(logpsdm_smooth);
                %%% PROBLEM: the threshold is based on a shifted version of
                %%% the full passband spectrogram, but the spectrogram used
                %%% to calculate trace lines is truncated to the annotation
                %%% box and shifted relative to that. Try fixing this.
                %psdm_anal = psdm_anal - max(psdm_anal(:));
                psdm_log_offset = max(psdm_anal(is_f_in_annot_range,:),[],'all');
                psdm_anal = psdm_anal - psdm_log_offset;
                th_psdm = prctile(psdm_anal(:),95);
                %}
                
                %%% set other pitch tracing parameters
                %%% (hard-coded for now)
                trace_penalty_coeff = 0.01; %0.008;
                trace_penalty_exp = 3; %2;
                trace_max_t_gap = Inf; %0.333;
                trace_max_f_gap = Inf; %10;
                
                %%% find the best trace line
                try
                    [t_trace, f_trace] = MUPPET.getTraceLine(t_stft, f_stft_annwin, logpsdm_denoised_annwin, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',spec_snr_th, 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap);
                    
                    %** TESTING                    
                    %[t_trace_nolog, f_trace_nolog] = MUPPET.getTraceLine(t_stft, f_stft_trace, psdm_tracecalc, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',10.^(th_psdm./10), 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap, 'LogWeights',false);
                    %[t_trace_average, f_trace_average] = MUPPET.getTraceLine(t_stft, f_stft_trace, psdm_tracecalc, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap, 'Method','dijkstra_averaged');
                    %[t_trace_reg1, f_trace_reg1] = MUPPET.getTraceLine_regression(t_stft, f_stft_trace, psdm_tracecalc, round(1./dt_stft), 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'Order',1);
                    %[t_trace_reg2, f_trace_reg2] = MUPPET.getTraceLine_regression(t_stft, f_stft_trace, psdm_tracecalc, round(1./dt_stft), 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_psdm, 'Order',2);
                    %[t_trace_dn, f_trace_dn] = MUPPET.getTraceLine(t_stft, f_stft_trace, psdm_denoised_tracecalc, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',th_denoised_psdm, 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap);
                    %[t_trace_relth1, f_trace_relth1] = MUPPET.getTraceLine(t_stft, f_stft_trace, psdm_tracecalc, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',3, 'ThresholdType','timebin', 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap);
                    %[t_trace_relth, f_trace_relth] = MUPPET.getTraceLine(t_stft, f_stft_trace, psdm_tracecalc, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',10, 'ThresholdType','timebin', 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap);
                    [t_trace_th3, f_trace_th3] = MUPPET.getTraceLine(t_stft, f_stft_annwin, logpsdm_denoised_annwin, 'PenaltyCoefficient',trace_penalty_coeff, 'PenaltyExponent',trace_penalty_exp, 'ClippingThreshold',3, 'MaxTimeGap',trace_max_t_gap, 'MaxFreqGap',trace_max_f_gap);
                    
                    t_trace_all = {t_trace_th3, t_trace};
                    f_trace_all = {f_trace_th3, f_trace};
                    
                    trace_plot_data = struct(...
                        'Color', {'w', 'r'},...
                        'Marker', {'x', 'o'},...
                        'MarkerSize', {4, 4},...
                        'LineWidth', {1, 1,},...
                        'DisplayName', {'Th 3', 'Th 10'}...
                        );
                    %** END TEST
                    %%% get the original (logged but unsmoothed) power values of the 
                    %%% trace line 
                    logpsdm_tracedata = logpsdm(is_f_in_annot_range,:);
                    [~, i_f_trace] = ismember(f_trace, f_stft_annwin);
                    [~, i_t_trace] = ismember(t_trace, t_stft);
                    logp_trace = logpsdm_tracedata(sub2ind(size(logpsdm_tracedata), i_f_trace, i_t_trace'));
                    %%% OLD
                    %{
                    psdm_tracedata = 10.*log10(psdm(is_f_in_annot_range,:));
                    [~, i_f_trace] = ismember(f_trace, f_stft_annwin);
                    [~, i_t_trace] = ismember(t_trace, t_stft);
                    p_trace = psdm_tracedata(sub2ind(size(psdm_tracedata), i_f_trace, i_t_trace'));
                    %}
                    
                    %%% store trace data in table
                    trace_line_w = table(t_trace', f_trace,  logp_trace, repelem(w,numel(t_trace),1), 'VariableNames',trace_table_vars);
                    %trace_line_w = table(t_trace', f_trace,  p_trace, repelem(w,numel(t_trace),1), 'VariableNames',trace_table_vars);
                    trace_lines = [trace_lines; trace_line_w];
                    
                catch ME
                    %%% issue warning if unable to build trace line
                    warning('Failed to find a trace line for call No. %d:\n%s', w, ME.message)
                    trace_line_w = [];
                end
                
                %%% plot trace line if specified
                if plot_trace && ~isempty(trace_line_w)
                    logpsdm_plot = logpsdm_denoised - max(logpsdm_denoised);
                    %logpsdm_plot = logpsdm_smooth - max(logpsdm_smooth);
                    
                    fig = gcf();
                    fig.Visible = 'off';
                    clf(fig);
                    ax = axes();

                    %MUPPET.plotTraceLine(ax, t_stft, f_stft, psdm_smooth, t_trace, f_trace, [f_min_ann,f_max_ann]);
                    MUPPET.plotTraceLine(ax, t_stft, f_stft, logpsdm_plot, t_trace_all, f_trace_all, [f_min_ann,f_max_ann], 'LineData',trace_plot_data);
                    ylim(ax,[LowerPassbandFreq,UpperPassbandFreq])
                    xlim(t_stft([1,end]))

                    %title(sprintf('%s - No. %d\n\\rmCost = %g\\it\\Deltaf\\rm^{%g} + 1;   Th=%.2fdB,  MaxTGap=%.2fs,  MaxFGap=%gHz', strrep(outfile_refname,'_','\_'), w, trace_penalty_coeff, trace_penalty_exp, th_psdm, trace_max_t_gap, trace_max_f_gap))
                    title(sprintf('%s - No. %d', strrep(outfile_refname,'_','\_'), w))

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
                sigPosRel = sigPos - annotPos(1); % energy-based signal position relative to start of annotation
                call_params_w = MUPPET.extractCallParams(sigPosRel, Fs, f_stft, psd, EnergyPercent, trace_line_w);
            else
                % cannot calculate parameters for unsuitable calls, so just
                % pass NaNs to parameter extraction routine
                call_params_w = MUPPET.extractCallParams(NaN, NaN, NaN, NaN, NaN, NaN);
            end
            
            %%% add call parameters to running table
            call_params = [call_params; struct2table(call_params_w)];
            
        end %call loop
        if ~isempty(fig)
            close(fig);
        end
        close(waitfig)
        
        % Create output files
        
        %%% SNR-expanded annotations output file
        out1_filename = [outfile_refname, '_SNR.csv'];
        PATH2OUTPUT_SNR_FILE = fullfile(PATH2OUTPUT, out1_filename);
        writetable(PLA, PATH2OUTPUT_SNR_FILE);
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
        
        % save tables to output variables if needed
        if nargout > 0
            outFieldName_primary = outfile_refname;
            outFieldName_secondary = sprintf('Annotations_%d',p);
            outFieldName = outFieldName_primary;
            for outvarnum = 1:nargout
                switch outvarnum
                    case 1
                        outvar = PLA;
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
