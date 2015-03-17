lsllib = lsl_loadlib();
eeg_info = lsl_streaminfo(lsllib,'Mitsar', 'EEG',136, 512,'cf_float32','test');
outlet = lsl_outlet(eeg_info);


while true
    outlet.push_sample(randn(136,1));
    pause(0.01);
end
