function [Start, Stop] = calcEng(x,perEng)
% Calculate the 90% energy start and stop position of a input signal
% x = data
% perEng =  percentage energy
% Last updated by Mike Adams
% 2024-03-21
temp = ((100 - perEng)/2)/100;
upperbound = 1-temp;
lowerbound = temp;
x2 = x.^2; %square signal to approx. energy

cum_FreqPing = cumsum(x2); % 
maxFreqPing = max(cum_FreqPing);
minFreqPing = min(cum_FreqPing);
max95Ping = maxFreqPing*upperbound;
max05Ping = maxFreqPing*lowerbound;
            
%find closest value to 95%
V = max95Ping;
C = cum_FreqPing;
A = repmat(C,[1 length(V)]);
[~,closestIndex] = min(abs(A-V'));
FreqPing95 = C(closestIndex); 
            
%find closest value to 5%
V = max05Ping;
C = cum_FreqPing;
A = repmat(C,[1 length(V)]);
[~,closestIndex] = min(abs(A-V'));
FreqPing05 = C(closestIndex); 
%%%
            
Stop90 = find(cum_FreqPing == FreqPing95);
Start90 = find(cum_FreqPing == FreqPing05);

           