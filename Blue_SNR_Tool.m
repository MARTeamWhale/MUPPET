% Blue_SNR_Tool.m
%
% Process Pamlab output for use in SNR tool.
%
%
% Last updated by Mike Adams
% 2024-02-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTE: https://www.mathworks.com/help/matlab/ref/listdlg.html

clear
close all

%Get list of Pamlab Output csv and path to wav files
PATH2INPUT = uigetdir('','SELECT FOLDER WITH PAMLAB OUTPUT');
PAMLAB_ANNOTATIONS = dir(fullfile(PATH2INPUT, '**\*.csv'));
PATH2DATA = uigetdir('','SELECT FOLDER WITH WAV FILES');

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
Freq_band = [SNR_PARAMS_filtered.LowerFrequency SNR_PARAMS_filtered.UpperFrequency];
NoiseDistance = SNR_PARAMS_filtered.NoiseDistance; 
BP_buffer = SNR_PARAMS_filtered.BP_Buffer;
Units = string(SNR_PARAMS_filtered.Units);
%%%
%%%
for p = 1:length(PAMLAB_ANNOTATIONS)%read in in Pamlab csv (Loop)
    file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
    opts = detectImportOptions(file, 'NumHeaderLines',2, 'Delimiter',',');
    %opts = detectImportOptions(file);
    %opts.VariableNamesLine = 3;
    %opts.Delimiter = ",";
    
    PLA = readtable(file,opts);
    PLA.SNR = NaN(height(PLA),1); %create location to save SNR
    
    %%% get wav file and read it in
    temp = split(PAMLAB_ANNOTATIONS(p).name,'.');
    temp(end) = {'wav'};
    FileName = strjoin(temp, '.');
    [x,Fs] = audioread(fullfile(PATH2DATA,FileName));
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    xt = [x t];
    
    %%%
    
    %%%  Loop through blue whale calls
    for i = 1:height(PLA)
    %%% Get Start90 and End90 RelativeStartTime
    %%% Transform Start90 and End90 with RelativeStartTime
    
    RelativeStartTime = PLA.RelativeStartTime(i);
    if ~isa(RelativeStartTime,'double')
         RelativeStartTime = str2double(PLA.RelativeStartTime(i));
    end
    
    PLA_StartTime90 = PLA.StartTime90(i);
    if ~isa(PLA_StartTime90,'double')
        PLA_StartTime90 = str2double(PLA.StartTime90(i));
    end
    
    PLA_StopTime90 = PLA.StopTime90(i);
    if ~isa(PLA_StopTime90,'double')
        PLA_StopTime90 = str2double(PLA.StopTime90(i));
    end
    
    Start90 = PLA_StartTime90 + RelativeStartTime;
    End90 = PLA_StopTime90 + RelativeStartTime;
    
    
    %%%
    %pass: raw wav,Start90, End90,[frequency band],buffer size,and noiseDistance to BP_clip.m
    %output: bandpassed wav clip + buffer
    %BP_clip = snr.BP_clip(x,Start90,End90,Freq_band,NoiseDistance,BP_buffer_samples);
    %%%
    %%%
    %pass: start90, stop90, noiseDistance,and bandpassed wav clip + buffer to extractSN.m
    %output: signal clip and noise clip
    %Temporary function test:
    [xSignal, xNoise] = snr.extractSN(x,Fs,Start90,End90,NoiseDistance,Units);
    %[xSignal, xNoise] = snr.extractSN(x, fs, sigStart, sigStop, noiseDist, units)
    %
    %%%
    %%%
    %pass: signal clip and noise clip to calculateSNR.m
    %output: SNR
    %PLA.SNR(i)  = snr.calculateSNR(xSignal, xNoise,Fs);
    %%%
    end           
end % end PAMLAB annotations loop
            
%OUTPUT: filename RelativeStartTime Start90 End90 SNR
          

    
