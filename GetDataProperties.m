function [protocols, protocols_show_as,durations, channels,settings_file] = GetDataProperties(pathname,fnames)
protocols = {};
protocols_show_as = {};
durations = [];
channels = {};

% if nargin < 2
%     fnames = 
% end
if ischar(fnames)
    [a, b, c]= fileparts(strcat(pathname,fnames)); %#ok<ASGLU>
    if verLessThan('matlab','8.1')
        s = regexp(b,' ','split');
    else
        s = strsplit(b);
    end
    protocols{end+1} = s{2};
    try
    protocols_show_as{end+1} = s{3};
    end
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
        if length(s) > 1
        protocols{end+1} = s{2};
        try %#ok<TRYNC>
        protocols_show_as{end+1} = s{3};
        end
        durations(end+1) = str2double(s{end});
        end
    end
end

all_fs = dir(pathname);
for n = 1:length(all_fs)
    if strcmp(all_fs(n).name,'ssd_exp_info.hdr') || strcmp(all_fs(n).name,'csp_exp_info.hdr')
        sh = fopen(strcat(pathname,'\',all_fs(n).name),'r');
        ch_str = fscanf(sh,'%c');
        if verLessThan('matlab','8.1')
            channels = regexp(ch_str,',','split');
        else
            channels = strsplit(ch_str,',');
        end
        fclose(sh);
    elseif strcmp(all_fs(n).name,'Exp_design.xml')
        settings_file = strcat(pathname,'\',all_fs(n).name);
    end
end

end