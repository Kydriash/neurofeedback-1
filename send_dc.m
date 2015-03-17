lsllib = lsl_loadlib();
eeg_info = lsl_streaminfo(lsllib,'Mitsar', 'EEG',136, 500,'cf_float32','test');
outlet = lsl_outlet(eeg_info);
DC=randn(136,1);

while true
    outlet.push_sample(randn(136,1)+ DC);
    pause(0.002);
end