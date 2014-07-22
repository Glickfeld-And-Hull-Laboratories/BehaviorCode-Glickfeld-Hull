function consts = lego_constants

tHostname = lower(hostname);
tUsername = getenv('USER');
homeDir = fullfile('/Users',tUsername);

% defaults
consts.centralDataPath = fullfile(homeDir, 'Documents/MWorks');
consts.dataPath = fullfile(homeDir, 'Documents/MWorks/Data');
consts.behavPdfPath = fullfile(homeDir, 'Documents/MWorks/BehavOutputPdfs');
consts.dbPath = fullfile(homeDir, 'Dropbox/LegoPDFs');
