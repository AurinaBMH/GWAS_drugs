function LDRelateTable = LDImport()
% Import SNP-SNP LD relationships
%-------------------------------------------------------------------------------

fid = fopen('2_2_LD_r2.csv','r');
C = textscan(fid,'%s%f%s','Delimiter',',','HeaderLines',1);
fclose(fid);
SNP_id_1 = C{1};
SNP_id_2 = C{3};
LD = C{2};
LDRelateTable = table(SNP_id_1,SNP_id_2,LD);

end
