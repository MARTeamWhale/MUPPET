function Baleen_SNR_Tool(varargin)
%
% Baleen_SNR_Tool.m
%
% Process Pamlab output for use in SNR tool.
%
%
% Written by Mike Adams
% Last updated by Wilfried Beslin
% 2024-05-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTE: https://www.mathworks.com/help/matlab/ref/listdlg.html

    % process optional input parameters (I/O paths)
    p = inputParser();
    p.addParameter('PAMLAB_DATA_FOLDER', '', @isfolder)
    p.addParameter('WAV_FILE_FOLDER', '', @isfolder)
    p.addParameter('OUTPUT_FOLDER_LOCATION', '', @isfolder)
    p.parse(varargin{:})
    PATH2INPUT = p.Results.PAMLAB_DATA_FOLDER;
    PATH2DATA = p.Results.WAV_FILE_FOLDER;
    PATH2OUTPUTDIRECTORY = p.Results.OUTPUT_FOLDER_LOCATION;
    
    % prompt user to specify I/O folder paths interactively if they were
    % not entered via the command line
    pathNotSpecified = @(p) isnumeric(p) && p == 0;
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
    
    % create output folder
    PATH2OUTPUT = fullfile(PATH2OUTPUTDIRECTORY,'SNR_OUTPUT');
    if ~exist(PATH2OUTPUT, 'dir')
       mkdir(PATH2OUTPUT)
    end

    
    %%% read in call type params: 
                             % Species x
                             % Call type x
                             %[frequency band] x
                             % noiseDistance x
                             % buffer size (maybe) x
    SNR_PARAMS = readtable("SNR_PARAMS.csv"); %Change PARAM file here if required
    specieslist = {SNR_PARAMS.Species};
    specieslist = unique(horzcat(specieslist{:}),'stable');
    sp_idx = listdlg('PromptString','Select a species:',...
                      'SelectionMode','single',...
                      'ListString',specieslist);
    species = specieslist(sp_idx,1);
    calltypelist = SNR_PARAMS.CallType(strcmp(SNR_PARAMS.Species,string(species{1,1})),:);
    calltypelist = {calltypelist};
    calltypelist = unique(horzcat(calltypelist{:}),'stable');
    %if length(calltypelist)== 1
    %    calltypelist = string(calltypelist);
    %end
    ct_idx = listdlg('PromptString','Select a call type:',...
                      'SelectionMode','single',...
                      'ListString',calltypelist);
    calltype = calltypelist(ct_idx,1);
    SNR_PARAMS_filtered = SNR_PARAMS(strcmp(SNR_PARAMS.Species,string(species{1,1})) & strcmp(SNR_PARAMS.CallType,string(calltype{1,1})),:);

    %%% extract variables from PARAMS table
    %Freq_band = [SNR_PARAMS_filtered.LowerFrequency SNR_PARAMS_filtered.UpperFrequency];
    NoiseDistance = SNR_PARAMS_filtered.NoiseDistance; 
    NoiseSize = 10; % hardcoded to 10 sec for now
    BP_buffer = SNR_PARAMS_filtered.BP_Buffer;
    Units = string(SNR_PARAMS_filtered.Units);


    %%% create empty variable to store bandpass filter object
    bandpass_filter = [];

    %%% process each artefact file
    for p = 1:length(PAMLAB_ANNOTATIONS)%read in in Pamlab csv (Loop) Possibly redundant...
        file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
        PLA = readtable(file);
        PLA.SNR = NaN(height(PLA),1); %create location to save SNR
        PLA.SNR_Adjusted = NaN(height(PLA),1); %create location to save SNR with noise power subtracted from numerator
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
                        'StopbandFrequency1', SNR_PARAMS_filtered.LowerStopbandFrequency,...
                        'PassbandFrequency1', SNR_PARAMS_filtered.LowerPassbandFrequency,...
                        'PassbandFrequency2', SNR_PARAMS_filtered.UpperPassbandFrequency,...
                        'StopbandFrequency2', SNR_PARAMS_filtered.UpperStopbandFrequency,...
                        'StopbandAttenuation1', 60,...
                        'StopbandAttenuation2', 60,...
                        'PassbandRipple', 1,...
                        'DesignMethod', 'kaiserwin',...
                        'SampleRate', Fs...
                        );
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

            %%% extract bandpass-filtered signal and noise samples
            %%{
            [xSignal, xNoise] = snr.extractSN(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter,...
                'NoiseDistance', NoiseDistance,...
                'IdealNoiseSize', NoiseSize,...
                'RemoveFromNoise', [PLA_Start_other,PLA_Stop_other],...
                'ClipBufferSize', BP_buffer);
            %}
            %** Testing defaults
            %[xSignal, xNoise] = snr.extractSN(x, Fs, [PLA_Start,PLA_Stop], bandpass_filter);
            %** legacy code
            %[xSignal, xNoise] = snr.extractSN_legacy(x, Fs, PLA_Start, PLA_Stop, NoiseDistance, BP_buffer, bandpass_filter, Units);

            %%% calculate SNR 
            %%% (leave NaN if not possible because signal is too close to 
            %%% endpoints)
            %** Consider adding to the output the duration of the signal and
            %** noise estimates, as well as a comment if the signal is too
            %** close to the beginning or the noise is too powerful.
            %** Or, instead of a comment, rig calculateSNR such that the noise
            %** power can never be greater than the signal power (make them
            %** equal if that happens) - this will return zeros and -Infs
            %** instead of negatives and complex numbers.
            if ~isempty(xSignal)
                [PLA.SNR(w), PLA.SNR_Adjusted(w)] = snr.calculateSNR(xSignal, xNoise);
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
            temp_name = split(PAMLAB_ANNOTATIONS(p).name,'.');
            temp_filename = [char(temp_name(1)) '_SNR.csv'];
            final_filename = generateUniqueName(PATH2OUTPUT,temp_filename);
            PATH2OUTPUT_FILENAME = fullfile(PATH2OUTPUT,final_filename);
            writetable(PLA,PATH2OUTPUT_FILENAME);
    end % end PAMLAB annotations loop
    
end
            
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
          

    
