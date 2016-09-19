function L = paramCornersAndCenter2Loc(corners, center)
   
   width = norm(corners(1:2) - corners(3:4));
   height = norm(corners(3:4) - corners(5:6));
   
   tmp = (corners(3:4)+corners(5:6))/2 - center;
   theta = angle(tmp(1) + tmp(2)*i);
   
   L = [center(1), center(2), width, height, theta];
end