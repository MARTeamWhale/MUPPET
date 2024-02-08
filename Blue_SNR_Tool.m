%Blue_SNR_Tool.m
%
%Process Pamlab output for use in SNR tool.

clear
close all

%Get list of Pamlab Output csv
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
specieslist = table(SNR_PARAMS.Species);
specieslist.Properties.VariableNames(1) = "Species";
specieslist.index = (1:height(specieslist)).';
disp(specieslist)
s = 0;
while s == 0
    sp_idx = input("Choose species by entering index number:");
    if sp_idx > height(specieslist)
        disp("Invalid species index")
    elseif sp_idx < 1
        disp("Invalid species index")
    elseif isinteger(sp_idx)
        disp("Invalid species index")
    else
        s = 1;
    end
end
species = specieslist(sp_idx,1);
calltypelist = SNR_PARAMS.CallType(strcmp(SNR_PARAMS.Species,string(species{1,1})),:);
calltypelist = table(calltypelist);
calltypelist.Properties.VariableNames(1) = "Call Type";
calltypelist.index = (1:height(calltypelist)).';
disp(calltypelist)
c = 0;
while c == 0
    ct_idx = input("Choose call type by entering index number:");
    if ct_idx > height(calltypelist)
        disp("Invalid call type index")
    elseif ct_idx < 1
        disp("Invalid call type index")
    elseif isinteger(ct_idx)
        disp("Invalid call type index")
    else
        c = 1;
    end
end
calltype = calltypelist(ct_idx,1);
SNR_PARAMS_filtered = SNR_PARAMS(strcmp(SNR_PARAMS.Species,string(species{1,1})) & strcmp(SNR_PARAMS.CallType,string(calltype{1,1})),:);
Freq_band = [SNR_PARAMS_filtered.LowerFrequency SNR_PARAMS_filtered.UpperFrequency];
NoiseDistance = SNR_PARAMS_filtered.NoiseDistance;
BP_buffer_samples = SNR_PARAMS_filtered.BufferSamples;
%%%
%%%
for p = 1:length(PAMLAB_ANNOTATIONS)%read in in Pamlab csv (Loop)
    file = fullfile(PAMLAB_ANNOTATIONS(p).folder,PAMLAB_ANNOTATIONS(p).name);
    opts = detectImportOptions(file);
    opts.VariableNamesLine = 3;
    opts.Delimiter = ",";
    PLA = readtable(file,opts);
    
    %%% get wav file and read it in
    temp = split(PAMLAB_ANNOTATIONS(p).name,'.');
    FileStartTime = readDateTime(PAMLAB_ANNOTATIONS(p).name);
    temp(end) = {'wav'};
    FileName = strjoin(temp, '.');
    x = audioread(fullfile(PATH2DATA,FileName));
    %%%
   
    %%% Get Start90 and End90 RelativeStartTime
    Start90 = str2double(PLA.StartTime90);
    End90 = str2double(PLA.StopTime90);
    RelativeStartTime = str2doublePLA.RelativeStartTime);
    %%%
    %%% Transform Start90 and End90 with RelativeStartTime 
    

        %Loop through blue whale calls
            %pass: raw wav,Start90, End90,[frequency band],buffer size,and noiseDistance to BP_clip.m
                %output: bandpassed wav clip + buffer 
            %pass: start90, stop90, noiseDistance,and bandpassed wav clip + buffer to extractSN.m
                %output: signal clip and noise clip
            %pass: signal clip and noise clip to calculateSNR.m
                %output: SNR
                
end % end PAMLAB annotations loop
            
%OUTPUT: filename RelativeStartTime Start90 End90 SNR
          

    
