function lock_handle_properties(objH, propList, functionId)
%LOCK_HANDLE_PROPERTIES (mh): install a listener to prevent prop changes
%
%    function lock_handle_properties(objH, propList, functionId)
%
%    After changing a graphics property statically, so that changes in other
%    properties are not reflected, this function will cause an error on
%    changing those properties.
%
% histed 120514

%% reg listener
axHH = handle(objH);  % get handle obj
for iC = 1:length(propList)
    tCN = propList{iC};
    hProp = findprop(axHH,tCN);
    hListener = handle.listener(axHH, hProp, 'PropertyPostSet', @errorCallback);
    setappdata(objH, [mfilename 'Listener'], hListener);
end
setappdata(objH, [ 'lock_handle_propertiesId'], functionId);

%% subf
function errorCallback(h, eventData)
disp(h)
disp(eventData)
axH = eventData.affectedObject;
functionId = getappdata(axH, ['lock_handle_propertiesId']);
tStr = sprintf('Function %s prevents changes to property %s; set this property BEFORE calling %s', ...
    functionId, h.Name, functionId);
error(tStr);

