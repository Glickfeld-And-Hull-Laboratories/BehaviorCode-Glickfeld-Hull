function [file]=sineWAV(freq,durationS, fileName)
%Generates a sine wave .wav file for a given duration and frequency.
%AM 14064
audioFolder = '~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/wavs';
if nargin<3
    fileName = strcat('sound-',round(mat2str(freq)), 'Hz-', mat2str(durationS), 's.wav');
end

Fs = 44800;
fullName = fullfile(audioFolder, fileName);
t = [0: 1/Fs : durationS-1/Fs]; %1 second, length 44800
f1 = sin(2*pi*freq*t);
sound(f1,Fs)

%write to disk
wavwrite(f1,Fs, fullName)
end