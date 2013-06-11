function camera_apply(cam_struct, axisHandle)
%CAMERA_APPLY (ps-utils): restore position of camera previously saved
%
%   CAMERA_APPLY(CAM_STRUCT, AXISHANDLE)
%
%  MH - http://github.com/histed/tools-mh

if nargin < 2, axisHandle = gca; end

cs = cam_struct;

set(axisHandle, 'CameraPosition', cs.CameraPosition, ...
                'CameraTarget', cs.CameraTarget, ...
                'CameraUpVector', cs.CameraUpVector, ...
                'CameraViewAngle', cs.CameraViewAngle);
