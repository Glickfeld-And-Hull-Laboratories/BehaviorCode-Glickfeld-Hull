function [file]=sineWAV(freq,durationS, fileName, Fs)
%Generates a sine wave .wav file for a given duration and frequency.
%AM 140522
audioFolder = '~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode';
if nargin<3
    fileName = strcat('sound-',mat2str(freq), 'Hz-', mat2str(durationS), 's');
if nargin<4,
    Fs = 44800;
end
fullName = fullfile(audioFolder, fileName);
t = [0:durationS/Fs:1-1/Fs]; %1 second, length 44800
f1 = sin(2*pi*freq*t);
sound(f1,Fs)
f2 = sin(2*pi*(2*freq)*t);
sound(f2,Fs)
%write to disk
wavwrite(f2,Fs, fullName)

end

