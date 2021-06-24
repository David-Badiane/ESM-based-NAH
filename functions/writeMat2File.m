function [dataTable] = writeMat2File(data,dstFileName, name, numVars, singleTitles)
% WRITEMAT2FILE this function saves data into a .csv or .txt file
%   INPUTS
%   data         (2Darray) = data to write ;
%   dstFileName  (string) = filename, can be .txt or .csv;
%   name         (cell array) = containing the variables names;
%   numVars      (2Darray) = num variable names -- length(names);
%   singleTitles (boolean) = if name contains all the names of the file - true 
%                            if you want to numerate name for all the cols
%                            of the file, false 
%   ex1 - name = {'x' 'y'}, numVars = 2, singleTitles = false    --> variableNames = {'x1' 'y1' ... 'xnCols' 'ynCols'}
%   ex2 - name = {'x' 'y'}, numVars = 2, singleTitles = true --> variableNames = {'x' 'y'}                           
%   OUTPUT
%   dataTable    (table) = written table in the file ;

columns = length(data(1,:));
names = cell(1,columns);
nIndexes = ceil(columns /numVars);
indexes = 1:nIndexes;

select = 1;
for ii = 1:columns
    names{ii} = [name{1+mod(ii-1, numVars)}, int2str(indexes(select))]; 
    
    if mod(ii,numVars) == 0
        select = select + 1;
    end
end

if singleTitles == true 
    names = cell(size(name));
        for ii = 1:length(name)
            names{ii} = name{ii};
        end
end

dataTable = array2table(data, 'VariableNames', names);
writetable(dataTable, dstFileName);
end

