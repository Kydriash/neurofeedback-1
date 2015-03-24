function [data] = ReadEEGData(fname)
f = fopen(fname,'r','n','US-ASCII');
% string_size = fread(f,1,'int');
% cols = fread(f,string_size,'char');
data_sz = fread(f,2,'int');
data = fread(f,data_sz','double'); 
fclose(f);


end

 