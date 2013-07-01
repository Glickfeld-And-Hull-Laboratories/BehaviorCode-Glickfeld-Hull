function cs = exptConstants

cs.trialOutcomes = { ...
    'success', ...
    'early', ...
    'late', ...
    'norelease'};

%% directories

tHostname = lower(hostname);

cs.userName = getenv('USER');
cs.homeDir = fullfile('/Users', cs.userName);

% defaults
cs.centralDataPath = '';
cs.dataPath = fullfile(cs.homeDir, 'Documents/MWorks/Data');
cs.behavPdfPath = fullfile(cs.homeDir, 'Documents/MWorks/BehavOutputPdfs');

switch tHostname
  case {'maunsellmouse1', 'maunsellmouse2', 'maunsellmouse3', ...
        'maunsellmouse4' }
    cs.centralDataPath = '/Users/holdanddetect/MWDataCentral/Data';
  case 'maunsellmousetest'
    % use defaults, pass
  case 'mambo' % MH laptop
    cs.centralDataPath = '/Users/histed/data/mus-behavior/SnowLeopard/Data';
  case 'andrews-macbook-pro-7'
    cs.centralDataPath = '/Users/andrewmckinney/Desktop/Data';    
end


