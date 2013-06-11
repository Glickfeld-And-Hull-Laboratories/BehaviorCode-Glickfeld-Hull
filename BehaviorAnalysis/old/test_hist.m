%%
addpath('/Library/Monkeyworks/Matlab');
cd '~/MonkeyWorks/trunk/MatlabToolbox';
%fname = 'TestDataFiles/100112-hold-test.mwk';
fname = 'TestDataFiles/100108-hold-b5.mwk';

%% read file
cos = getCodecs(fname);
codec = cos.codec;
unqCodes = unique([codec.code]);

%% t
tic;
ev = getEvents(fname, codec_tag2code(codec, 'trialStart'));
toc

%% get all user events
codeNames = {codec.tagname};
userIx = cellfun(@(x)~isempty(x),regexp(codeNames, '^[^#].*', 'start'));
userCodes = codeNames(userIx);
tic
ev = getEvents(fname, codec_tag2code(codec, userCodes));
toc

    
    
%% get all events
%getCodes = setdiff(unqCodes, codec_tag2code(codec, '#codec')); %drop reserved code
tic
ev = getEvents(fname, unqCodes);
codes = [ev.event_code];
time_us = [ev.time_us];
code_data = {ev.data};
toc

%% plot trial times
trialStartTimes = time_us(codec_tag2code(codec, 'trialStart') == codes);
trialStartIx = codec_tag2code(codec, 'trialStart') == codes;

%% plot a test of codes
cRange = [1:1000]; 
figH = figure;
subplot(1,2,1);
plot(time_us(cRange), 'x');
subplot(1,2,2);
plot(diff(time_us(cRange), 'x');
for iC = cRange
    tagname = codec(codec_code2idx(codec, codes(iC))).tagname;
    if iC == 1, lastTime = 0; else lastTime = time_us(iC-1); end
    fprintf(1, '* %-25s - abs %8d, diff %8d\n', ...
            tagname, time_us(iC), time_us(iC)-lastTime);
    disp(code_data{iC}); 
end


%% extract all announceMessages
ev = getEvents(fname, codec_tag2code(codec, '#announceMessage'));
nEv = length(ev);
for iE = 1:nEv
    if iE == 1, lastTime = 0; else lastTime = ev(iE-1).time_us; end    
    fprintf(1, '%10d diff %10d: %s\n', ...
            ev(iE).time_us, ev(iE).time_us-lastTime, ev(iE).data.message);
end

%% get 

