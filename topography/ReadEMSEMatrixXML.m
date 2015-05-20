function M = ReadEMSEMatrixXML(XMLFileName)

a = xmlread(XMLFileName);
a.getElementsByTagName('EMSE_Labeled_Matrix').item(0);
lm = a.getElementsByTagName('EMSE_Labeled_Matrix').item(0);
em = lm.getElementsByTagName('EMSE_Mtx').item(0);
c = char(em.getFirstChild.getData);
k=1;
i=1;
while( i < length(c) )
    nm = [];
    while((double(c(i))~=10)) 
        nm = [nm c(i)];
        i = i+1;
    end;
    
    if(~isempty(nm))
        v{k} = str2num(nm);
        k = k+1;
    end;
    i = i+1;
end;
nrow = v{1}(1);        
ncol = v{1}(2);
M = zeros(nrow,ncol);
for i=1:nrow
    M(i,:) = v{i+1};
end;
