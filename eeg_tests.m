function eeg_tests()
%% derived signal tests
DS = DerivedSignal;

%prepare channels
hardware_channels = {};
used_by_signal_channels = {};
%prepare dummy struct signal
signal = struct();
signal.channels = used_by_signal_channels ;
signal.sSignalName = 'test_signal';
signal.filters = {};
%create derived signal
DS = DerivedSignal(1,signal,100,4000,hardware_channels,200);


end