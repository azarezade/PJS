function ftrs=generateharrftrs(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to generate set of harr-like features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftrnum = 250;

for i = 1:ftrnum
    ftrs{1,i}=generateftr(param);
end
