close all;
clear();
clear classes();

eeg = EEGLSL;
eeg.Connect('type', 'Data');
%eeg.Connect('type', 'EEG');
%eeg.Connect('name','Mitsar');
%eeg.InitTimer; 