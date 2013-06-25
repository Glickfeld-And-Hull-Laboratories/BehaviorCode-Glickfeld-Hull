function consts = reward_delay_constants

tHostname = lower(hostname);
tUsername = getenv('USER');
homeDir = fullfile('/Users',tUsername);

% defaults
consts.centralDataPath = fullfile(homeDir, 'Desktop/CentralData');
consts.dataPath = fullfile(homeDir, 'Documents/MWorks/Data');
consts.behavPdfPath = fullfile(homeDir, 'Documents/MWorks/BehavOutputPdfs');
