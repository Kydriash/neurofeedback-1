close all;
clear();
clear classes();

warning('on'); %#ok<WNON>
eeg = EEGLSL;
eeg.Connect('type', 'Data');
%eeg.Connect('type', 'EEG');
%eeg.Connect('name','Mitsar');
%eeg.InitTimer; 