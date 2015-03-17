function [A,sz] =  ReadData( fname)

f = fopen(fname,'r');
sz = fread(f,2,'int'); 
A = fread(f,sz','float'); 
fclose(f);
f0 = f;

end

