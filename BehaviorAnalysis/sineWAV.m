function [file]=sineWAV(freq,duration)
%Generates a sine wave .wav file for a given duration and frequency.
%AM 140522
audioFolder = '~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode';
fileName = strcat('sound-',mat2str(freq), '-', mat2str(duration));

rfreq=2*pi*freq;
d=[linspace(0,duration,8192*duration)];
signal=cos(rfreq*d);
soundsc(signal);

end

