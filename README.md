# SNR_Tool_Development
Work in progress
Developing a SNR tool for use in cetacean research

## Working outline for SNR functions
### function "snr.extractSN"
[xSignal, xNoise] = snr.extractSN(x, fs, sigStart, sigStop, noiseDist, units)
* x = data vector
* fs = sampling rate
* sigStart = signal start time or sample
* sigStop = signal stop time or sample
* noiseDist = distance from signal from which to sample noise, in time or samples
* units = string specifying if start, stop, and distance inputs represent time or samples

### function "snr.calculateSNR"
[snrVal] = calculateSNR(xSignal, xNoise)
* xSignal = signal samples
* xNoise = noise samples
