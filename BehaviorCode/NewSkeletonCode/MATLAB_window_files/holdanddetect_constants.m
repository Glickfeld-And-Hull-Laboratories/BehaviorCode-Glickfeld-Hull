function cs = holdanddetect_constants

tHostname = lower(hostname);
tUsername = getenv('USER');

homeDir = fullfile('/Users',tUsername);

% defaults
cs.centralDataPath = fullfile(homeDir, 'Desktop/CentralData');
cs.dataPath = fullfile(homeDir, 'Documents/MWorks/Data');
cs.behavPdfPath = fullfile(homeDir, 'Documents/MWorks/BehavOutputPdfs');
