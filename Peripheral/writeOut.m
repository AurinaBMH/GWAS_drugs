function writeOut(proteinNames,filter,fileName)


fid = fopen(fileName,'w');

theNames = proteinNames(filter);
for i = 1:length(theNames)
    fprintf(fid,'%s\n',theNames{i});
end

fclose(fid);

fprintf(1,'Wrote %u protien names to %s!\n',length(theNames),fileName);


end
