close all;
clear();
clear classes();

warning('on'); %#ok<WNON>
warning('off','backtrace');
eeg = EEGLSL;
eeg.RunInterface('name','NVX136_Data');
%eeg.Connect('type', 'Data');
%eeg.Connect('type', 'EEG');
%eeg.Connect('name','Mitsar');
%eeg.InitTimer; 