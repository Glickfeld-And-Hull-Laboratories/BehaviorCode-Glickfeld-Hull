function input = plotCerebellarStim(data_struct, input)

% essential consts
figNum = 1;
name = 'CerebellarStim';

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
%mylogicalarray = logical(cell2mat(mycellarray));  %# Convert the cell array
wheelIx = diff(cell2mat(wheelIx));

wheelTime = input.quadratureTimesUs;
emptyIndex = cellfun(@isempty,wheelTime);       %# Find indices of empty cells
wheelTime(emptyIndex) = {int64(0)};                    %# Fill empty cells with 0
%mylogicalarray = logical(cell2mat(mycellarray));  %# Convert the cell array
wheelTime = diff(cell2mat(wheelTime))/1000;     % convert time to ms

wheelSpeed  = double(wheelIx)./double(wheelTime);
avgWheelN   = 100;

if length(wheelSpeed) >= 100
    meanWheelSpeed  = arrayfun(@(i) mean(wheelSpeed(i:i+avgWheelN-1)),1:avgWheelN:length(wheelSpeed)-avgWheelN+1);
    plot(meanWheelSpeed);
% lH = plot(double(wheelIx), nTrial);
% set(lH, 'Color', 'r', ...
%         'LineWidth', 3);
% lH2 = plot(double(wheelIx), 100);
% set(lH2, 'Color', 'k', ...
%         'LineWidth', 2);
        
    title('Mean Running rate every 1 second')
    ylabel('Running Rate');
    xlabel('Second')
else
    plot(wheelSpeed);
    title('Mean Running rate for about 10 milisecond')
    ylabel('Running Rate');
    xlabel('Second')
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
bH = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF figure : %s', outName), ...
               'Units', 'pixels', ...
               'Position', [5 5 650 20], ...
               'Callback', { @saveButtonCb, epParams });
           
% subfunciton

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfig_print(epParamsIn{:});
