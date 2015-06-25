function eeg = load_exp()
%select folder
folder_name = uigetdir('','Select directory with experiment');
if isempty(folder_name)
    return
end
listing = dir(folder_name);

%load exp settings
[protocols, protocols_show_as, durations, channels,settings_file] = GetDataProperties(folder_name,fnames);

nfs = NeurofeedbackSession;
nfs.LoadFromFile(settings_file);
%load data

%create experiment
eeg = EEGLSL;

%push data


end