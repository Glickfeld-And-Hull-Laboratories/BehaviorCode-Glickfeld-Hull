function map = cmap_unique13(m)
%CMAP_POSNEG_RDBU (ca-utils): a colormap to represent pos/neg values
%   black is the middle value
%
%   MH 080501: support rescaling zero point / asymmetric map
%
%$id: cmap_posneg_rdbu.m 200 2006-04-01 00:24:05Z histed $ 

if nargin < 1 || isempty(m), m = 13; end
if m ~= 13, error('Only 13 colors in this colormap'); end

uniqueColors = { ...
    '007D34', %Vivid Green
    'F6768E', %Strong Purplish Pink
    '00538A', %Strong Blue
    '93AA00', %Vivid Yellowish Green
    '53377A', %Strong Violet
    'FF8E00', %Vivid Orange Yellow
    'B32851', %Strong Purplish Red
    'F4C800', %Vivid Greenish Yellow
    '7F180D', %Strong Reddish Brown
    'FF7A5C', %Strong Yellowish Pink
    '593315', %Deep Yellowish Brown
    'F13A13', %Vivid Reddish Orange
    '232C16', %Dark Olive Green
    };

nColors = length(uniqueColors); 
map = repmat(NaN, [nColors 3]);
for iC = 1:nColors
    tC = uniqueColors{iC};
    
    for iN=1:3
        stN = (iN-1)*2+1;
        tStr = tC(stN:stN+1);
        tNum = hex2dec(tStr);
        map(iC, iN) = tNum/255;
    end
end




