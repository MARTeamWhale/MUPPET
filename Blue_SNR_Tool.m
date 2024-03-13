% Blue_SNR_Tool.m
%
% Process Pamlab output for use in SNR tool.
%
%
% Written by Mike Adams
% Last updated by Wilfried Beslin
% 2024-03-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTE: https://www.mathworks.com/help/matlab/ref/listdlg.html

clear
close all

%Get list of Pamlab Output csv and path to wav files
PATH2INPUT = uigetdir('','SELECT FOLDER WITH PAMLAB OUTPUT');
PAMLAB_ANNOTATIONS = dir(fullfile(PATH2INPUT, '\*.csv'));
PATH2DATA = uigetdir('','SELECT FOLDER WITH WAV FILES');
WAVFILES = struct2table(dir(fullfile(PATH2DATA, '**\*.wav')));
PATH2OUTPUTDIRECTORY = uigetdir('','SELECT DIRECTORY TO CREATE OUTPUT FOLDER');

if PATH2OUTPUTDIRECTORY == 0
    temp1 = split(PATH2INPUT,'\');
    temp2 = temp1(1:end-1)';
    PATH2OUTPUTDIRECTORY = fullfile(strjoin(temp2,'\'));
    clear(temp1,temp2)
end

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
BP_buffer = SNR_PARAMS_filtered.BP_Buffer;
Units = string(SNR_PARAMS_filtered.Units);

%%% create empty variable to store bandpass filter object
bandpass_filter = [];

%%% process each artefact file
for p = 1:length(PAMLAB_ANNOTATIONS)%read in in Pamlab csv (Loop) Possibly redundant...
    file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
    %opts = detectImportOptions(file, 'NumHeaderLines',2, 'Delimiter',',');
    %opts = detectImportOptions(file);
    %opts.VariableNamesLine = 3;
    %opts.Delimiter = ",";
    %PLA = readtable(file,opts);
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
    
    %%% get wav file and read it in 
    for w = 1:num_annotations %Start rows loop
        %%% update waitbar
        t_elapsed = toc;
        t_rem = t_elapsed.*((num_annotations-(w-1))/(w-1));
        waitbar(w/num_annotations, waitfig, sprintf('%s\nEstimated time remaining: %s', waitmsg, duration(0,0,t_rem)))
        
        %%% process 
        
        temp = split(PLA.filename(w),'.');
        temp(end) = {'wav'};    
        if isempty(x)||~strcmp(strjoin(temp, '.'), FileName) %check if first time running
           FileName = strjoin(temp, '.');
           for i = 1:length(WAVFILES.name)
                if contains(WAVFILES.name(i), FileName)
                   PATH2WAV = char(fullfile(WAVFILES.folder(i),WAVFILES.name(i)));
                   continue
                end
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
  
       %%% Get Start90 and End90 RelativeStartTime
       %%% Transform Start90 and End90 with RelativeStartTime

       RelativeStartTime = PLA.RelativeStartTime(w);
       if ~isa(RelativeStartTime,'double')
             RelativeStartTime = str2double(PLA.RelativeStartTime(w));
        end

        PLA_StartTime90 = PLA.StartTime90(w);
        if ~isa(PLA_StartTime90,'double')
            PLA_StartTime90 = str2double(PLA.StartTime90(w));
        end

        PLA_StopTime90 = PLA.StopTime90(w);
        if ~isa(PLA_StopTime90,'double')
            PLA_StopTime90 = str2double(PLA.StopTime90(w));
        end

        Start90 = PLA_StartTime90 + RelativeStartTime;
        End90 = PLA_StopTime90 + RelativeStartTime;
        
        %%% extract bandpass-filtered signal and noise samples
        [xSignal, xNoise] = snr.extractSN(x, Fs, Start90, End90, NoiseDistance, BP_buffer, bandpass_filter, Units);
        
        %%% calculate SNR 
        %%% (leave NaN if not possible because signal is too close to 
        %%% endpoints)
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
          

    
