function f0 =  SaveData( fname,data)

f = fopen(fname,'w');
fwrite(f,size(data),'int'); 
fwrite(f,data,'float'); 
fclose(f);
f0 = f;
%f = fopen('cm_136_5.bin','r'), sz = fread(f,2,'int'); A = fread(f,sz,'float'); fclose(f);

%f = fopen('cm_136_5.bin','r'), sz = fread(f,2,'int'); A = fread(f,sz,'float'); fclose(f);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

