close all;
clear();
clear classes();

eeg = EEGLSL;
eeg.Connect('type', 'Data');
%eeg.Connect('name','Mitsar');
%eeg.InitTimer; 