function [dataTable] = writeMat2File(data,dstFileName, name, numVars, singleTitles)
% WRITEMAT2FILE this function save data into csv


%   INPUTS
%   data         (2Darray) = ;
%   dstFileName  (2Darray) = ;
%   name         (2Darray) = ;
%   numVars      (2Darray) = ;
%   singleTitles (2Darray) = ;

%   OUTPUT
%   dataTable    (2Darray) = ;

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

