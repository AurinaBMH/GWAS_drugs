function writeOut(proteinNames,filter,fileName)
% Writes a list of names (proteinNames) with a given filter applied (filter)
%-------------------------------------------------------------------------------

fid = fopen(fileName,'w');

theNames = proteinNames(filter);
for i = 1:length(theNames)
    fprintf(fid,'%s\n',theNames{i});
end

fclose(fid);

fprintf(1,'Wrote %u protein names to %s!\n',length(theNames),fileName);


end
