function cam_struct = camera_save(axisHandle)
%CAMERA_SAVE (ps-utils): save position of camera so it can be loaded later
%
%   CAM_STRUCT = CAMERA_SAVE(AXISHANDLE)
%
%   See also CAMERA_APPLY
%
%  MH - http://github.com/histed/tools-mh

if nargin < 1, axisHandle = gca; end

cs.CameraPosition = get(axisHandle, 'CameraPosition');
cs.CameraTarget = get(axisHandle, 'CameraTarget');
cs.CameraUpVector = get(axisHandle, 'CameraUpVector');
cs.CameraViewAngle = get(axisHandle, 'CameraViewAngle');

cam_struct = cs;