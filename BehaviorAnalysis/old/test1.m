%%
addpath('/Library/MonkeyWorks/Matlab');
cd '~/MonkeyWorks/trunk/MatlabToolbox';
%fname = 'TestDataFiles/100112-hold-test.mwk';
%fname = '/Library/MonkeyWorks/DataFiles/test-100114.mwk'; % broken
fname = '/Library/MonkeyWorks/DataFiles/100115-test.mwk';

[tp tn te] = fileparts(fname);
if exist(fname, 'dir') == 7
    disp('removing index');
    unix('rm -rf /tmp/t');
    unix(sprintf('/bin/mv %s /tmp/t', fname));
    unix(sprintf('mv /tmp/t/%s /Library/MonkeyWorks/DataFiles/%s', [tn te], ...
                 [tn te]));
    if ~exist(fname, 'file')
        error('!');
    end
end

% read file
cos = getCodecs(fname);
%codec = cos.codec;
%unqCodes = unique([codec.code]);

%% t

%% t
tic;
ev = getEvents(fname, 1:53);
toc

