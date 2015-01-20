function cs = holdanddetect_constants

tHostname = lower(hostname);
tUsername = getenv('USER');

homeDir = fullfile('/Users',tUsername);

% defaults
cs.centralDataPath = '';
cs.dataPath = fullfile(homeDir, 'Documents/MWorks/Data');
cs.behavPdfPath = fullfile(homeDir, 'Documents/MWorks/BehavOutputPdfs');

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


