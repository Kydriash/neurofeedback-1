function [filenames, protocols, durations, header, chs] = GetDataLength(dirname) %#ok<INUSD>
%dirname = 'D:\neurofeedback\results\Gorin_Alexey\2015-04-24\15-14-48';
files = dir(dirname);
filenames = {};
protocols = {};
durations = []; %seconds
header = {};


for f = 1:length(files)
    [folder, fname, ext] = fileparts(strcat(dirname,'\',files(f).name)); %#ok<ASGLU>
    if strcmp(ext,'.bin')
        if verLessThan('matlab','8.1')
            s = regexp(fname,' ','split');
        else
            s = strsplit(fname);
        end
        filenames{end+1} = strcat(dirname,'\',files(f).name); %#ok<AGROW>
        protocols{end+1} = s{2}; %#ok<AGROW>
        durations(end+1) = str2double(s{end}); %#ok<AGROW>
    elseif strcmp(fname,'exp_info')
        h = fopen(strcat(dirname,'\',files(f).name),'r');
        header = fscanf(h,'%c');
        fclose(h);
    elseif strcmp(fname,'ssd_exp_info')
        sh = fopen(strcat(dirname,'\',files(f).name),'r');
        if verLessThan('matlab','8.1')
            chs = regexp(fscanf(sh,'%c'),',','split');
        else
        chs = strsplit(fscanf(sh,'%c'),',');
        end
        fclose(sh);
    end
end


end