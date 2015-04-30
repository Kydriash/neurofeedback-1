function [protocols, durations, channels] = GetDataProperties(pathname,fnames)
protocols = {};
durations = [];
channels = {};

if ischar(fnames)
    [a, b, c]= fileparts(strcat(pathname,fnames)); %#ok<ASGLU>
    if verLessThan('matlab','8.1')
        s = regexp(b,' ','split');
    else
        s = strsplit(b);
    end
    protocols{end+1} = s{2};
    durations(end+1) = str2double(s{end});
else
    for f = fnames
        %files{end+1} = strcat(pathname,f{1});
        
        [a, b, c]= fileparts(strcat(pathname,f{1})); %#ok<ASGLU>
        if verLessThan('matlab','8.1')
            s = regexp(b,' ','split');
        else
            s = strsplit(b);
        end
        protocols{end+1} = s{2};
        durations(end+1) = str2double(s{end});
    end
end

all_fs = dir(pathname);
for n = 1:length(all_fs)
    if strcmp(all_fs(n).name,'ssd_exp_info.hdr')
        sh = fopen(strcat(pathname,'\',all_fs(n).name),'r');
        if verLessThan('matlab','8.1')
            channels = regexp(fscanf(sh,'%c'),',','split');
        else
            channels = strsplit(fscanf(sh,'%c'),',');
        end
        fclose(sh);
    end
end

end