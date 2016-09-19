function [ initGeoParam, initAffParam] = initState(gtCenters,gtCorners,objParam)
% This fucntion find intial state of tracker
% -------------------------------------------------------------------------
%--- input---
% gtCorners = 
% objSize   =
% 
%--- output---
% initGeoParam = [x_center, y_center, scale, rotation, aspect_ratio, skew]; 
% initAffParam = initial 6 affine parameters
% 
% -- other----
% initP = [x_center, y_center, width, height, rotation]; 
%          initial position of object:
%
% Copyright (C) Ali Zarezade.  All rights reserved.
% -------------------------------------------------------------------------
objSize = objParam.size;
initP = paramCornersAndCenter2Loc(gtCorners(1,:), gtCenters(1,:));
initGeoParam = paramLoc2Geom(initP, objSize(1));
initAffParam = paramGeom2Aff(initGeoParam);
end

