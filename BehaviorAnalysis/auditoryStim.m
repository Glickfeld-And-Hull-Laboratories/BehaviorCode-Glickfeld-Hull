function [ output_args ] = auditoryStim(maxFreqHz, stepsPerOctave, nLevels, durationS)
startup
%%AUDITORYSTIM: generates an octave of wave files exactly stepsPerOctave/8
%%log units away from each other. Calls sineWAV.m
%MaxFreqHz = max frequency (specified in Hz)
%stepsPerOctave = whole range of one octave
%nLevels = number of stimuli to be created
%durationS = duration of stimulus (in seconds)

minFreqHz = log2(maxFreqHz)-stepsPerOctave;
tFreqSpace = logspace(log2(maxFreqHz), minFreqHz, 8);
freqSpace=tFreqSpace(1:nLevels);
for i=1:length(freqSpace);
    freq = freqSpace(i);
    fileName = strcat('sound-level', mat2str(i), '.wav')
    sineWAV(freq, durationS, fileName)
end