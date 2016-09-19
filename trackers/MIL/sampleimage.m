function samples = sampleimage(params,state,flag,image)
%   ---------------------------------> x
%   |
%   |
%   |
%   |
%   |
%   |
%   |
%   \/
%   y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to sample negative and positive instances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag == 0  %negative samples in the first frame
    inrad = params.srchwinsz*2;
    outrad = round(params.initpostrainrad*1.5);
    maxnum = params.negnumtrain;
end
if flag ==1  %positive samples in the first frame
    inrad = params.initpostrainrad;
    outrad = 0;
    maxnum = params.posmaxtrain;
end
if flag ==2  %negative samples for updating classifiers in run-time stage
    inrad = round(params.srchwinsz*1.5);
    outrad = params.posradtrain+5;
    maxnum = params.negnumtrain;
end
if flag ==3   %positive samples for updating classifiers in run-time stage
    inrad = params.posradtrain;
    outrad = 0;
    maxnum = params.posmaxtrain; 
end
if flag ==4  %sampling candidates for object estimation
    inrad = params.srchwinsz;
    outrad = 0;
    maxnum = params.posmaxtrain; 
end
x = state(1,1);
y = state(1,2);
w = state(1,3);
h = state(1,4);

[row,col]=size(image);

rowsz = row - h;
colsz = col - w;

inradsq = inrad*inrad;
outradsq = outrad*outrad;
	
minrow = max(1,y-inrad+1);
maxrow = min(rowsz,y+inrad+1);
mincol = max(1,x-inrad+1);
maxcol = min(colsz,x+inrad+1);

samplesize = (maxrow-minrow+1)*(maxcol-mincol+1);
prob = maxnum/samplesize;

i = 1;
for r = minrow:maxrow
		for c = mincol:maxcol
			dist = (y-r)*(y-r) + (x-c)*(x-c);
			if( rand(1)<prob && dist < inradsq && dist >= outradsq )
      
				samples(i,1) = c;  %x
				samples(i,2) = r;  %y      
				samples(i,3) = w;  %width
				samples(i,4) = h;  %height
				i=i+1;
	
            end
        end
end