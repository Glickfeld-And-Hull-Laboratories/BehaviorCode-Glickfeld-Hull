function cs = holdanddetect_constants

tHostname = lower(hostname);

switch tHostname
  case {'maunsellmouse1', 'maunsellmouse2', 'maunsellmouse3', ...
        'maunsellmouse4' }
    cs.centralDataPath = '/Users/holdanddetect/MWDataCentral/Data';
    cs.dataPath = '/Users/holdanddetect/Documents/MWorks/Data';
    cs.behavPdfPath = '/Users/holdanddetect/Documents/MWorks/BehavOutputPdfs';
 case 'maunsellmousetest'
    cs.centralDataPath = '';
    cs.dataPath = '/Users/lindsey/Documents/MWorks/Data';
    cs.behavPdfPath = '/Users/lindsey/Documents/MWorks/BehavOutputPdfs';
  case 'mambo' % MH laptop
    cs.centralDataPath = '/Users/histed/data/mus-behavior/SnowLeopard/Data';
  case 'andrews-macbook-pro-7'
    cs.centralDataPath = '/Users/andrewmckinney/Desktop/Data';    
end

