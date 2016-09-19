function mat2LatexTable(data,horiz,vertic,grayLevel,numFormat,roundDigit,filename,best)
% data = [ 99    68     8    15    99;
%          26    74    22    82     7;
%          65    45    91    53    44];
% horiz = {'OAB','MIL','IVT','PBGM','PBGS'};
% vertic = {'david','faceocc','trellis'};
% filename = '/home/ali/Documents/Thesis/reports/Paper/result tables/test.tex';
% grayLevel = [0.6 0.8]

%% Initalize
data=roundn(data,roundDigit);

firstGrayLevel = grayLevel(1);
secondGrayLevel = grayLevel(2);

%% Create Rows
rows = cell(length(horiz)+3,1);
rows{1,1} = '\hline';

horizNames=[];
for i=1:length(horiz)
    horizNames = [horizNames, ' & ', horiz{i}];
end
horizNames = [horizNames, ' \\'];
rows{2,1} = horizNames;

rows{3,1} = '\hline \hline';

for i=1:size(data,1)
    currentRow = vertic{i};
    tmpRowData = data(i,:);
    
    if isequal(best,'max')
        valueFirst = max(tmpRowData);
        indxFirst = find(tmpRowData==valueFirst);
        tmpRowData(indxFirst) = -inf;
        valueSecond = max(tmpRowData);
        indxSecond = find(tmpRowData==valueSecond);
    elseif isequal(best,'min')
        valueFirst = min(tmpRowData);
        indxFirst = find(tmpRowData==valueFirst);
        tmpRowData(indxFirst) = inf;
        valueSecond = min(tmpRowData);
        indxSecond = find(tmpRowData==valueSecond);        
    else
        disp('Error:  Set best')
    end
    
    for j=1:size(data,2)
        switch j
            case num2cell(indxFirst)
                currentRow = [currentRow, ' & ', '\cellcolor[gray]','{',num2str(firstGrayLevel),'}', sprintf(numFormat,data(i,j)) ];    
            case num2cell(indxSecond)
                currentRow = [currentRow, ' & ', '\cellcolor[gray]','{',num2str(secondGrayLevel),'}', sprintf(numFormat,data(i,j)) ];    
            otherwise
                currentRow = [currentRow, ' & ', sprintf(numFormat,data(i,j)) ];    
        end
    end
    if i==size(data,1)-1
        rows{i+3,1} = [currentRow ' \\  \hline \hline'];
    else
        rows{i+3,1} = [currentRow ' \\  \hline'];
    end
end

%% Creat Column Type
col = ['|c||',repmat('c|',[1,length(horiz)])];

%% Print to File
fid = fopen(filename,'w'); 
fprintf(fid,'%s\n',['\begin{tabular}{',col,'}']);
for i = 1:length(rows)
    fprintf(fid,'%s\n',rows{i}); 
end
fprintf(fid,'%s\n','\end{tabular}');
fclose(fid);
