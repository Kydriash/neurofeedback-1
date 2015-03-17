close all;
clear();
clear classes();

eeg = EEGLSL;
eeg.Connect('type', 'EEG');
%eeg.Connect('name','Mitsar');
%eeg.InitTimer; 