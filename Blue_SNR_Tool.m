%Blue_SNR_Tool.m
%
%Process Pamlab output for use in SNR tool.

clear
close all

%Get list of Pamlab Output csv
PATH2INPUT = uigetdir('','SELECT FOLDER WITH PAMLAB OUTPUT');
PAMLAB_ANNOTATIONS = dir(fullfile(PATH2INPUT, '**\*.csv'));

%read in call type params: 
                         % Species 
                         % Call type
                         %[frequency band]
                         % noiseDistance
                         % buffer size (maybe)
SNR_PARAMS = readtable("SNR_PARAMS.csv");
fig = uifigure;
species = SNR_PARAMS.Species;
uidropdown(fig,"Species",species)
%
%read in in Pamlab csv (Loop)
    % get wav file and read it in
    % Get Start90 and End90 RelativeStartTime
    % Transform Start90 and End90 with RelativeStartTime 
    
        %Loop through blue whale calls
            %pass: raw wav,Start90, End90,[frequency band],buffer size,and noiseDistance to BP_clip.m
                %output: bandpassed wav clip + buffer 
            %pass: start90, stop90, noiseDistance,and bandpassed wav clip + buffer to extractSN.m
                %output: signal clip and noise clip
            %pass: signal clip and noise clip to calculateSNR.m
                %output: SNR
            
%OUTPUT: filename RelativeStartTime Start90 End90 SNR
           

        
    
