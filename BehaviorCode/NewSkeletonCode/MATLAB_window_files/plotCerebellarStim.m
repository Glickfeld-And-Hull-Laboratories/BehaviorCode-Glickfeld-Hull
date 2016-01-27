function input = plotCerebellarStim(data_struct, input)

% essential consts
figNum = 1;
name = 'CerebellarStim';
smoothType = 'moving';

% set up environment
tHostname = lower(hostname);
tUsername = getenv('USER');
homeDir = fullfile('/Users',tUsername);

% defaults
consts.centralDataPath = fullfile(homeDir, 'Documents/MWorks');
consts.dataPath = fullfile(homeDir, 'Documents/MWorks/Data');
consts.behavPdfPath = fullfile(homeDir, 'Documents/MWorks/BehavOutputPdfs');

%smoothType = 'lowess';

%nTrial = input.trialSinceReset;

wheelIx = input.quadratureValues;
emptyIndex = cellfun(@isempty,wheelIx);       %# Find indices of empty cells
wheelIx(emptyIndex) = {int64(0)};                    %# Fill empty cells with 0
wheelIx = cell2mat(wheelIx);
wheelIx = diff(wheelIx);

wheelTime = input.quadratureTimesUs;
emptyIndex = cellfun(@isempty,wheelTime);       %# Find indices of empty cells
wheelTime(emptyIndex) = {int64(1)};                    %# Fill empty cells with 0
wheelTime = cell2mat(wheelTime);
wheelTime_diff = diff(wheelTime)/1000;     % convert time to ms
time_zero = find(wheelTime_diff==0);       % find non zero time difference
wheelIx(time_zero) = [];
wheelTime_diff(time_zero) = [];

time_one = find(wheelTime==1);
wheelTime(time_one(2:end))=[];


wheelSpeed  = double(wheelIx)./double(wheelTime_diff);
avgWheelN   = 10;

if length(wheelSpeed) >= 10
    WheelTimePoint  = wheelTime;
    WheelTimePoint  = WheelTimePoint(1:avgWheelN:end)/1e6; 
    meanWheelSpeed  = arrayfun(@(i) mean(wheelSpeed(i:i+avgWheelN-1)),1:avgWheelN:length(wheelSpeed)-avgWheelN+1);
    plot(WheelTimePoint(1:length(meanWheelSpeed)),smooth(WheelTimePoint(1:length(meanWheelSpeed)),meanWheelSpeed,smoothType));
    %plot(WheelTimePoint(1:length(meanWheelSpeed)),meanWheelSpeed);
    %plot(meanWheelSpeed);
% lH = plot(double(wheelIx), nTrial);
% set(lH, 'Color', 'r', ...
%         'LineWidth', 3);
% lH2 = plot(double(wheelIx), 100);
% set(lH2, 'Color', 'k', ...
%         'LineWidth', 2);
        
    title('Mean Running rate every 10 movements registered')
    ylabel('Running Rate (pulses/ms)');
    xlabel('Timestamp (s)')
% else
%     plot(wheelSpeed);
%     title('Mean Running rate for every 1 movement registered')
%     ylabel('Running Rate');
%     %xlabel('Second')
end

%% Add a save button
if ~exist(consts.behavPdfPath, 'file')
  mkdir(consts.behavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  consts.behavPdfPath, ...
                  datestr(now, 'yymmdd'), input.subjectNum);
epParams = { figNum, outName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF figure : %s', outName), ...
               'Units', 'pixels', ...
               'Position', [5 5 650 20], ...
               'Callback', { @saveButtonCb, epParams });
           
% subfunciton

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfig_print(epParamsIn{:});
