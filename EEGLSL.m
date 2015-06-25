classdef EEGLSL < handle
    % >> self = EEGLSL;
    % >> self.sampling_frequency = 500;
    % >> self.channel_labels = {'a', 'b', 'c', 'd', 'e' ,'f'};
    % >> self.exp_data_length = 10000;
    % >> self.current_window_size = 10;
    % >> self.derived_signals{1} = self.CreateNewDS('raw');
    % >> self.derived_signals{2} = self.CreateNewDS('ds', ones(1,6),[10 15]);
    % >> size(self.derived_signals{1}.collect_buff.raw,2)
    %
    %ans =
    %
    % 6
    %
    % >> self.derived_signals{2}.Apply(ones(6,10),1);
    % >> self.derived_signals{2}.ring_buff.raw(self.derived_signals{2}.ring_buff.fst:self.derived_signals{2}.ring_buff.lst)
    %
    %ans =
    %
    %    0.0001    0.0006    0.0021    0.0051    0.0101    0.0171    0.0262    0.0372    0.0496    0.0629
    %
    % >> self.UpdateFeedbackSignal;
    % >> self.feedback_manager.feedback_vector(1)
    %
    %ans =
    %
    % 0.0211
    %
    %
    % >> self.derived_signals{2}.Apply(zeros(6,2),1);
    %
    % >> self.derived_signals{2}.ring_buff.raw(self.derived_signals{2}.ring_buff.fst:self.derived_signals{2}.ring_buff.lst)
    %
    %ans =
    %
    %    0.0001    0.0006    0.0021    0.0051    0.0101    0.0171    0.0262    0.0372    0.0496    0.0629
    %
    % >> self.derived_signals{2}.Apply([ones(6,4) zeros(6,2)],1);
    % >> self.derived_signals{2}.ring_buff.raw(self.derived_signals{2}.ring_buff.lst-9:self.derived_signals{2}.ring_buff.lst)
    %
    %ans =
    %
    % 0.0262    0.0372    0.0496    0.0629 0.0763    0.0890    0.0998    0.1080    0.1124    0.1117
    %
    % >> self.UpdateFeedbackSignal;
    % >> self.feedback_manager.feedback_vector(1)
    %
    %ans =
    %
    % 0.0773
    %
    % >> self.feedback_protocols{1} = self.CreateNewRP('Baseline',16);
    % >> self.current_protocol = 1;
    % >> self.Update_Statistics;
    % >> self.feedback_manager.average(1)
    %
    %ans =
    %
    % 0.0505
    % >> self.feedback_manager.standard_deviation(1)
    %
    %ans =
    %
    % 0.0436
    %
    % >> self.derived_signals{3} = self.CreateNewDS('composite', repmat(ones(1,6),3,1),[10 15]);
    % >> self.derived_signals{3}.signal_type = 'Composite';
    % >> self.derived_signals{3}.Apply(ones(6,16),1);
    % >> self.derived_signals{3}.ring_buff.raw(self.derived_signals{2}.ring_buff.fst:self.derived_signals{2}.ring_buff.lst)
    %
    % ans =
    %
    % 0.0002    0.0011    0.0037    0.0089    0.0174    0.0296    0.0454    0.0644    0.0860    0.1090    0.1322    0.1541    0.1729    0.1870    0.1948    0.1946
    %
    
    properties
        %%%%%%%%%%%%%%%%%%%%%% plot & user interface related members %%%%%%%%%%
        %%figures
        fig_interface %settings
        raw_and_ds_figure %plots
        fig_feedback %fb bar
        %%settings
        plot_length %sec
        plot_size %samples
        plot_refresh_rate %sec
        fb_refresh_rate %sec
        %%plots and axes
        %%raw eeg subplot
        raw_subplot
        raw_plot
        raw_shift
        r_ytick_labels
        raw_ydata_scale
        raw_line %text
        raw_plot_min
        raw_plot_max
        raw_plot_shift
        raw_fit_plot
        %%derived_signals subplot
        ds_subplot
        ds_shift
        ds_plot
        ds_ytick_labels
        ds_ydata_scale
        ds_line
        ds_plot_min
        ds_plot_max
        ds_plot_shift
        ds_fit_plot
        %%feedback subplot
        fbplot_handle %bar
        feedback_axis_handle
        fb_stub %text
        show_fb %bool
        fb_type %string
        
        %%uicontrol
        connect_button
        disconnect_button
        curr_protocol_text
        sn_to_fb_dropmenu
        sn_to_fb_string
        y_limit %fb plot
        status_text
        status_field
        log_text
        add_notes_window
        add_notes_field
        write_notes %pushbutton
        raw_scale_slider
        ds_scale_slider
        montage_fname
        montage_fname_text
        %%data
        channel_labels
        %%%%%%%%%%%% LSL and data input related objects %%%%%%
        streams
        inlet
        data_receive_rate
        path_text
        settings_file
        settings_file_text
        nd %temporarily stores new data
        %%%%%%%%%%%%%%%%%%%%% Data info %%%%%%%%%%%%%%%%%%%%%%
        sampling_frequency
        channel_count
        max_ch_count
        %%%%%%%%%%%%%%%%%%%%% Timers and callbacks %%%%%%%%%%%
        timer_new_data
        timer_disp
        timer_fb
        
        timer_new_data_function
        timer_disp_function
        timer_fb_function
        %%%%%%%%%%%%%%%%%%%% Status members
        connected
        recording
        finished
        raw_ylabels_fixed
        ds_ylabels_fixed
        yscales_fixed
        raw_yscale_fixed
        ds_yscale_fixed
        fb_statistics_set
        %%%%%% Signal processing and feedback related members%
        composite_montage %used_ch by used_ch matrix
        allxall% all ch by all ch matrix
        signals %parameters
        derived_signals
        feedback_protocols
        feedback_manager
        current_protocol
        next_protocol
        protocol_sequence
        protocol_types
        samples_acquired
        signal_to_feedback
        used_ch
        last_to_fb
        to_proceed
        count
        to_fb
        %%%%%% subject data management related
        subject_record
        path
        exp_data_length
        %other
        tstop
        sizes
        window %counts fb windows
        default_window_size
        buffer_length
        ssd %if an ssd protocol exists
        program_path
        subjects_dropmenu
        edit_protocols_button
        protocol_duration_text
        exp_design
        fb_manager_set
        protocol_indices
        fb_sigmas
        from_file
        fnames
        files_pathname
        looped
        run_protocols
        settings
        bad_channels
        raw_data_indices
        paused
        current_window_size
        init_band
        csp_settings
        bluetooth_connections
    end
    
    methods
        
        function self = EEGLSL(self) %#ok<INUSD>
            
            
            
            
            self.plot_length = 4;
            self.sampling_frequency = -1;
            self.streams = {};
            self.plot_refresh_rate = 0.2;
            self.data_receive_rate = 0.005;
            self.fb_refresh_rate = 0.1;
            self.show_fb = 1;
            
            self.max_ch_count = -1; % -1 to get all the channels
            self.connected = false;
            self.fig_feedback = figure('Visible', 'off');
            self.channel_labels = {};
            
            self.current_protocol = 0;
            self.feedback_protocols = {};
            self.exp_data_length = 0;
            self.samples_acquired = 0;
            self.y_limit = [-1 7];
            self.subject_record = SubjectRecord;
            [self.program_path, ~, ~] = fileparts(which(mfilename));
            self.path = strcat(self.program_path,'\results');
            self.signal_to_feedback = 2;
            
            self.settings_file_text = 'LeftVsRightMu.nss.xml';
            self.settings_file =  'settings\LeftVsRightMu.nss.xml';
            self.recording = 0;
            self.next_protocol = 1;
            self.finished = 0;
            self.raw_ylabels_fixed = 0;
            self.ds_ylabels_fixed = 0;
            self.yscales_fixed = 0;
            self.fb_statistics_set = 0;
            self.raw_ydata_scale = 1000;
            self.ds_ydata_scale = 1000;
            self.nd = [];
            
            
            self.samples_acquired = 0;
            %self.montage_fname = 'C:\Users\user1\AppData\Local\MCS\NeoRec\nvx136.nvx136.monopolar-Pz';
            self.montage_fname = 'D:\neurofeedback\settings\nvx136.nvx136.monopolar-Pz.xml';
            self.montage_fname_text = 'nvx136.nvx136.monopolar-Pz';
            self.raw_shift = 1;
            self.ds_shift = 1;
            %self.sizes = [0]; %#ok<NBRAK>
            self.window = 0;
            self.fb_type = 'Color';
            self.default_window_size = 0;
            self.buffer_length = 0;
            self.ssd = 0;
            self.raw_yscale_fixed = 0;
            self.fb_manager_set = 0;
            self.protocol_indices = 0;
            self.fb_sigmas = 8;
            self.from_file = 0;
            
            self.settings = struct();
            self.settings.subject = 'Null';
            self.settings.montage_fname = 'D:\neurofeedback\settings\nvx136.nvx136.monopolar-Pz.xml';
            self.settings.settings_file =  'settings\LeftVsRightMu.nss.xml';
            self.bad_channels = {};
            self.raw_data_indices = [];
            self.paused = 0;
            self.init_band = [7 9];
            
            
        end
        function UpdateFeedbackSignal(self)
            if length(self.derived_signals) > 1
                if ~self.fb_manager_set
                    self.feedback_manager.average = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.standard_deviation = ones(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                    self.fb_manager_set=1;
                end
                
                %                 if self.current_protocol >0 && self.current_protocol <= length(self.feedback_protocols)
                %                     if self.current_protocol <= length(self.feedback_manager.window_size)
                %                         n = self.feedback_manager.window_size(self.current_protocol);
                %                     else
                %                         n = self.default_window_size;
                %                     end
                %                 elseif self.default_window_size == 0 %zero protocol at the beginning
                %                     n = self.feedback_manager.window_size(1);
                %                 else
                %                     n = self.default_window_size;
                %                 end
                
                for s = 2:length(self.derived_signals)
                    dat = self.derived_signals{s}.ring_buff.raw(self.derived_signals{s}.ring_buff.lst-self.current_window_size+1:self.derived_signals{s}.ring_buff.lst);
                    avg  = self.feedback_manager.average(s-1);
                    sdev = self.feedback_manager.standard_deviation(s-1);
                    val = sum(abs(dat))/self.current_window_size;
                    self.feedback_manager.feedback_vector(s-1)  = (val-avg)/sdev;
                end
                
                if self.recording
                    fb = zeros(self.current_window_size,6);
                    self.window = self.window + 1;
                    fb(:,1) = self.signal_to_feedback-1;
                    fb(:,2) = self.feedback_manager.feedback_vector(self.signal_to_feedback-1);
                    fb(:,3) = self.feedback_manager.average(self.signal_to_feedback-1);
                    fb(:,4) = self.feedback_manager.standard_deviation(self.signal_to_feedback-1);
                    try
                        if self.current_protocol <= length(self.feedback_manager.window_size)
                            fb(:,5) = self.feedback_manager.window_size(self.current_protocol);
                        else
                            fb(:,5) = self.default_window_size;
                        end
                    catch
                        'An error occured while accessing self.feedback_manager.window_size, function UpdateFeedbackSignal' %#ok<NOPRT>
                    end
                    fb(:,6) = self.window;
                    
                    self.feedback_manager.feedback_records.append(fb);
                end
            end
        end
        function Update_Statistics(self)
            global selected
            if(self.current_protocol>0 && self.current_protocol <= length(self.feedback_protocols))
                if strfind(lower(self.feedback_protocols{self.current_protocol}.protocol_name),'ssd')
                    self.inlet.close_stream();
                    stop(self.timer_new_data);
                    stop(self.timer_disp);
                    stop(self.timer_fb);
                    
                    if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                        mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                    end
                    pth = (strcat(self.path,'\',self.subject_record.subject_name)); %#ok<NASGU>
                    
                    
                    peaks_found = 0; %#ok<NASGU>
                    ds = 0;
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if N > self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                        self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1;
                    else
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                    end
                    %save x_raw2.mat x_raw;
                    if length(self.derived_signals) > 1
                        self.derived_signals(2:end) = [];
                    end
                    %%while ~peaks_found
                    ds = ds + 1;
                    try
                        
                        
                        %find the peak
                        
                        i = 1;
                        dF = self.feedback_protocols{self.current_protocol}.band;
                        %dF= 1;
                        clear C SSD V Fr G B A Rng chs;
                        %                    x_raw = x_raw+randn(size(x_raw))*0.01*mean(abs(x_raw(:)));
                        for f0 = 2.5:0.5:25
                            pp(1,:) = [f0-3*dF,f0-dF];
                            if(pp(1,1)<=0)
                                continue;
                            end;
                            pp(2,:) = [f0-dF,f0+dF];
                            pp(3,:) = [f0+dF,f0+3*dF];
                            for f = 1:3
                                [z, p, k] = cheby1(3,1,pp(f,:)/(0.5*self.sampling_frequency),'bandpass');
                                [b(f,:),a(f,:)] = zp2tf(z,p,k);
                                x = filtfilt(b(f,:),a(f,:),x_raw)';
                                C{f} = x*x'/size(x,2);
                                C{f} = C{f}+ 0.1*trace(C{f})/size(C{f},1)*eye(size(C{f}));
                            end;
                            try
                                %[v e] = eig(C{2},0.5*(C{1}+C{3}),'chol');
                                [u, s, v] = svd((C{1}+C{3})^(-1)*C{2}); %#ok<ASGLU>
                            catch
                                13 %#ok<NOPRT>
                            end
                            % [mxv, mxi] = max(diag(e));
                            
                            mxi = 1;
                            mxv = s(1,1);
                            SSD(i) = mxv;
                            V(:,i) = u(:,mxi);
                            G(:,i) = u(mxi,:);
                            Fr(i) = f0; %frequencies
                            Rng(i,:) = pp(2,:); % band
                            i = i+1;
                        end;
                    catch
                        16 %#ok<NOPRT>
                    end
                    %%draw the picture
                    hh = figure;
                    plot(Fr,SSD); xlabel('frequency, Hz');
                    annotation('textbox', [0.2,0.8,0.1,0.1],'String', {'Two left clicks and Enter to select range.','Backspace to delete last point','Right click to finish'});
                    % F = []; %#ok<NASGU>
                    F = getpts(hh);
                    
                    
                    if length(F) == 1
                        peaks_found = 1; %#ok<NASGU>
                        close(hh);
                        %break;
                    end
                    close(hh);
                    
                    left_point = min(F);
                    right_point = max(F);
                    
                    ind = (find(Fr>=left_point & Fr<=right_point));
                    [~, ind_max] = max(SSD(ind));
                    middle_point = ind(ind_max);
                    disp(strcat(num2str(Fr(middle_point)),' Hz'));
                    channel_mapping(ds) = figure;
                    stem(G(:,middle_point));
                    set(get(channel_mapping,'Children'),'XTick',(1:1:length(self.used_ch)));
                    set(get(channel_mapping,'Children'),'XTickLabel',self.used_ch(:,1));
                    %find the largest SSD for the middle point
                    w_ssd  = V(:,middle_point);
                    %                         w =1-diag(w_ssd);
                    %                         if w == ones(size(w))
                    %                             break
                    %                         end
                    %x_raw = x_raw * w;
                    hh1 = figure; %#ok<NASGU>
                    StandChannels = load('StandardEEGChannelLocations.mat');
                    rearranged_map = rearrange_channels(G(:,middle_point),self.used_ch, StandChannels.channels);
                    topoplot(rearranged_map, StandChannels.channels, 'electrodes', 'labelpoint', 'chaninfo', StandChannels.chaninfo);
                    
                    %b_ssd = B(middle_point,:);
                    %a_ssd = A(middle_point,:);
                    for i=1:length(w_ssd)
                        %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                        chan_w{i,1} = self.derived_signals{1}.channels{i};
                        chan_w{i,2} = w_ssd(i);
                    end;
                    
                    
                    %%write to file
                    try %write
                        %convert cell chan_w to structure
                        
                        chs(size(chan_w,1)) = struct();
                        for ch = 1:size(chan_w,1)
                            chs(ch).channel_name = chan_w{ch,1};
                            chs(ch).coefficient = chan_w{ch,2};
                        end
                        chan_structure = struct();
                        chan_structure.channels.channel = chs;
                        x = struct2xml(chan_structure);
                        sn = {''};
                        while isempty(sn{1})
                            sn = inputdlg('Enter derived signal name','Derived signal name',1,{strcat('DS_',num2str(ds))});
                        end
                        %                         dummy_signal = struct();
                        %                         dummy_signal.sSignalName = sn;
                        %                         dummy_signal.channels = self.derived_signals{1}.channels;
                        %                         dummy_signal.filters = cell(0,0);
                        if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                            mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                        end
                        pth = (strcat(self.path,'\',self.subject_record.subject_name));
                        full_name = [pth '\SSD_' sn{1} '.xml'];
                        if exist(full_name,'file')
                            choice = questdlg('The spatial filter for this person already exists. Rewrite the filter?','Rewrite?','Yes','No, use the old one','No, use the old one');
                            switch choice
                                case 'Yes'
                                    f = fopen(full_name,'w');
                                    fwrite(f,x);
                                    fclose(f);
                                    spatial_filter = chs;
                                case 'No, use the old one'
                                    
                                    s = xml2struct(full_name);
                                    channels_coeff = cell(length(s.channels.channel),2);
                                    for i = 1:length(s.channels.channel)
                                        channels_coeff{i,1} = strtrim(s.channels.channel{i}.channel_name.Text);
                                        channels_coeff{i,2} = str2double(strtrim(s.channels.channel{i}.coefficient.Text));
                                    end
                                    spatial_filter = channels_coeff;
                            end
                        else
                            spatial_filter = chs;
                            f = fopen(full_name,'w');
                            fwrite(f,x);
                            fclose(f);
                        end
                    catch
                        'Error while writing to file, function UpdateStatistics' %#ok<NOPRT>
                    end
                    try
                        
                        
                        %                         NewDS= DerivedSignal(1,dummy_signal, self.sampling_frequency, self.exp_data_length ,self.channel_labels,self.plot_length);
                        %                         NewDS.signal_name = sn{1};
                        %                         NewDS.ring_buff = circVBuf(self.plot_size,1,0);
                        %                         NewDS.collect_buff = circVBuf(self.exp_data_length,1,0);
                        %                         NewDS.UpdateSpatialFilter(spatial_filter,self.derived_signals{1},self.bad_channels);
                        %                         NewDS.UpdateTemporalFilter(Rng(middle_point,:));
                        %                         NewDS.channels_file = full_name;
                        self.derived_signals{end+1} = CreateNewDS(self,sn{1},spatial_filter, Rng(middle_point,:),full_name);
                    catch
                        'Error while creating a new derived signal, function UpdateStatistics' %#ok<NOPRT>
                    end
                    % end
                    % end
                    if exist('channel_mapping','var')
                        close(channel_mapping);
                    end
                    
                    %%calculate feedback avg,std
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if N > self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                        self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1;
                    else
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                    end
                    j = 1;
                    x = zeros(size(x_raw,1),length(self.derived_signals{1}.spatial_filter));
                    for i = 1:length(self.derived_signals{1}.spatial_filter)
                        if self.derived_signals{1}.spatial_filter(i)
                            x(:,i) = x_raw(:,j);
                            j = j+1;
                        else
                            x(:,i) = 0;
                        end
                    end
                    
                    %%filter
                    temp_derived_signal = CreateNewDS(self,'Temp',spatial_filter, Rng(middle_point,:),full_name);
                    %                     dummy_signal = struct();
                    %                     dummy_signal.sSignalName = 'Temp';
                    %                     channels = cell(length(self.derived_signals{1}.channels),2);
                    %                     for i = 1:length(channels)
                    %                         channels{i,1} = self.derived_signals{1}.channels{i};
                    %                         channels{i,2} = 1;
                    %                     end
                    %                     dummy_signal.filters = cell(0,0);
                    %                     dummy_signal.channels = channels;
                    %                     temp_derived_signal = DerivedSignal(1,dummy_signal, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                    %                     temp_derived_signal.ring_buff = circVBuf(self.plot_size,1,0);
                    %                     temp_derived_signal.collect_buff = circVBuf(self.exp_data_length,1,0);
                    %                     temp_derived_signal.UpdateSpatialFilter(spatial_filter,self.derived_signals{1},self.bad_channels);% = w_ssd(self.derived_signals{1}.channels_indices)';
                    %                     temp_derived_signal.UpdateTemporalFilter(Rng(middle_point,:));
                    temp_derived_signal.Apply(x',1);
                    values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                    %calcuate feedback stats
                    values = abs(values);
                    self.feedback_manager.average(1) = mean(values);
                    self.feedback_manager.standard_deviation(1) = std(values);
                    
                    self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                    self.fb_manager_set = 1;
                elseif strfind(lower(self.feedback_protocols{self.current_protocol}.protocol_name), 'csp')
                    
                    self.inlet.close_stream();
                    stop(self.timer_new_data);
                    stop(self.timer_disp);
                    stop(self.timer_fb);
                    
                    if length(self.derived_signals) > 1
                        self.derived_signals(2:end) = [];
                    end
                    
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if N > self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                        self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1;
                    else
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                    end
                    %%determine Number of Components (vectors to tell one
                    %pattern from another)
                    if ~isempty(self.feedback_protocols{self.current_protocol}.n_comp)
                        Ncomp = self.feedback_protocols{self.current_protocol}.n_comp;
                    else
                        Ncomp = 2;
                    end
                    if ~isempty(self.feedback_protocols{self.current_protocol}.init_band)
                        init_band = self.feedback_protocols{self.current_protocol}.init_band;
                    else
                        init_band = self.init_band;
                    end
                    for ib = 1:4
                        band = init_band +ib-1 ;
                        [z, p, k] = cheby1(3,1,band/(0.5*self.sampling_frequency),'bandpass');
                        [b,a] = zp2tf(z,p,k);
                        x = filtfilt(b,a,x_raw)';
                        C10 = x(:,1:fix(end/2))* x(:,1:fix(end/2))'/fix(size(x,2)/2);
                        C20 = x(:,fix(end/2)+1:end)* x(:,fix(end/2)+1:end)'/fix(size(x,2)/2);
                        
                        nchan = size(C10,1);
                        
                        %%regularize covariances
                        Lambda = 0.1;%%%%%%%%%%%%%%%%%%%%%%
                        C1 = C10 + Lambda * trace(C10) * eye(nchan) / nchan;
                        C2 = C20 + Lambda * trace(C20) * eye(nchan) / nchan;
                        %%do generalized eigenvalue decomp
                        [V, d] = eig(C1,C2); %#ok<ASGLU>
                        iV = inv(V);
                        M12{ib} = V(:,[1:Ncomp, end-Ncomp+1:end])';
                        G12{ib} = iV([1:Ncomp, end-Ncomp+1:end],:);
                        
                    end
                    %%present the heads
                    hh1 = figure;
                    StandChannels = load('StandardEEGChannelLocations.mat');
                    chan_labels = self.used_ch(:,1)';
                    Nbands = length(G12);
                    PlotIndex = 1;
                    selected = {};
                    for ib = 1:Nbands
                        rearranged_map = rearrange_channels(G12{ib}',chan_labels, StandChannels.channels);
                        Nmaps = size(rearranged_map,2);
                        for tpm=1:Nmaps
                            sp(PlotIndex) = subplot(Nbands,Nmaps,PlotIndex);
                            
                            %                            topoplot(rearranged_map(:,tpm), StandChannels.channels, 'electrodes', 'labelpoint', 'chaninfo', StandChannels.chaninfo);
                            topoplot(rearranged_map(:,tpm), StandChannels.channels,  'chaninfo', StandChannels.chaninfo);
                            hold on;
                            sibs = get(sp(PlotIndex),'Children');
                            for k = 1:length(sibs)
                                set(sp(PlotIndex).Children(k), 'ButtonDownFcn', @(src,event)toggleplot(src,event));
                            end
                            title(num2str(PlotIndex));
                            %add legend
                            PlotIndex = PlotIndex+1;
                        end;
                    end
                    okay_btn = uicontrol('Parent',hh1,'Style', 'pushbutton', 'String', 'OK', 'Callback', 'uiresume', 'Position', [230 10 100 35],'Tag','SelectHeadsBtn','enable','off');  %#ok<NASGU>
                    
                    
                    uiwait();
                    
                    
                    %                     mi = '';
                    %                     while isempty(mi)
                    %                         mi = inputdlg;
                    %                     end
                    %                     mi = strsplit(mi{1});
                    %
                    %                     MapIndex =  str2double(mi);
                    MapIndex = cellfun(@str2num,selected);
                    BandNumber = fix((MapIndex-1)/(2*Ncomp))+1;
                    CompNumber = mod(MapIndex-1,(2*Ncomp))+1;
                    close(hh1);
                    w_ssd = zeros(length(MapIndex),length(M12{1}));
                    for bn = 1:length(MapIndex)
                        w_ssd(bn,:) = M12{BandNumber(bn)}(CompNumber(bn),:);
                    end
                    %%if the band is the same for all inputs, a choice is
                    %offered whether to make one combined signal or
                    %several; If the bands are different, no such
                    %choice is offered
                    %determine whether the band is the same
                    if isempty(nonzeros(BandNumber - max(BandNumber))) %if it is empty, then all the bands are the same
                        button = questdlg(['One combined DS or ', num2str(length(MapIndex)), ' different DS?'],'Choose the number of DS', 'One', num2str(length(MapIndex)),'One');
                        if strcmp(button, 'One')
                            %make one DS with 2-D sp_filter
                            band = init_band + BandNumber(1)-1;
                            chan_w = cell(length(self.derived_signals{1}.channels),size(w_ssd,1)+1);
                            w_ssd = w_ssd';
                            % 2-d spatial_filter
                            for i=1:length(w_ssd)
                                %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                                chan_w{i,1} = self.derived_signals{1}.channels{i};
                                for ind = 2:size(w_ssd,2)+1
                                    chan_w{i,ind} = w_ssd(i,ind-1);
                                end
                                
                            end;
                            clear chs
                            chs(size(chan_w,1)) = struct();
                            for ch = 1:size(chan_w,1)
                                chs(ch).channel_name = chan_w{ch,1};
                                coeff = '';
                                for ind = 2:size(w_ssd,2)+1
                                    coeff = [coeff ',' num2str(chan_w{ch,ind})];
                                end
                                
                                chs(ch).coefficient = coeff(2:end);
                            end
                            
                            
                            
                            sn = {''};
                            while isempty(sn{1})
                                sn = inputdlg('Enter derived signal name','Derived signal name',1,{strcat('CombinedDS_',num2str(band(1)),'-',num2str(band(2)))});
                            end
                            
                            %%write spatial_filter to file
                            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                                mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                            end
                            pth = (strcat(self.path,'\',self.subject_record.subject_name));
                            full_name = [pth '\CSP_' sn{1} '.xml'];
                            
                            if exist(full_name,'file')
                                choice = questdlg('The spatial filter for this person already exists. Rewrite the filter?','Rewrite?','Yes','No, use the old one','No, use the old one');
                                switch choice
                                    case 'Yes'
                                        chan_structure = struct();
                                        chan_structure.channels.channel = chs;
                                        x = struct2xml(chan_structure);
                                        f = fopen(full_name,'w');
                                        fwrite(f,x);
                                        fclose(f);
                                        spatial_filter = chan_w;
                                    case 'No, use the old one'
                                        s = xml2struct(full_name);
                                        channels_coeff = cell(length(s.channels.channel),2);
                                        for i = 1:length(s.channels.channel)
                                            channels_coeff{i,1} = strtrim(s.channels.channel{i}.channel_name.Text);
                                            coeffs = str2num(strtrim(s.channels.channel{i}.coefficient.Text)); %#ok<ST2NM>
                                            for j = 2:length(coeffs)+1
                                                channels_coeff{i,j} = coeffs(j-1);
                                            end
                                        end
                                        %convert struct to cell!!
                                        spatial_filter = channels_coeff;
                                end
                            else
                                spatial_filter = chs;
                                f = fopen(full_name,'w');
                                fwrite(f,chs);
                                fclose(f);
                            end
                            
                            %%create DS
                            %
                            self.derived_signals{end+1} = CreateNewDS(self,sn{1},spatial_filter, band,full_name);
                            self.derived_signals{end}.signal_type = 'combined';
                            
                            %%set fb
                            
                            %
                            temp_derived_signal = CreateNewDS(self,'Temp',spatial_filter, band);
                            %
                            temp_derived_signal.signal_type = 'combined';
                            %
                            N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                            if N > self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1
                                x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                                self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1;
                            else
                                x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                            end
                            temp_derived_signal.Apply(x_raw',1);
                            values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                            %calcuate feedback stats
                            values = abs(values);
                            self.feedback_manager.average(1) = mean(values);
                            self.feedback_manager.standard_deviation(1) = std(values);
                            
                            
                            self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                            self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                            self.fb_manager_set = 1;
                            
                        end
                    end
                    
                    if length(self.derived_signals)< 2 %if no signals were added so far
                        for ind = 1:length(MapIndex)
                            band = init_band + BandNumber(ind)-1;
                            for i=1:length(w_ssd)
                                %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                                chan_w{i,1} = self.derived_signals{1}.channels{i};
                                chan_w{i,2} = w_ssd(ind,i);
                            end;
                            
                            %write to file
                            try %write
                                %convert cell chan_w to structure
                                clear chs
                                chs(size(chan_w,1)) = struct();
                                for ch = 1:size(chan_w,1)
                                    chs(ch).channel_name = chan_w{ch,1};
                                    chs(ch).coefficient = chan_w{ch,2};
                                end
                                chan_structure = struct();
                                chan_structure.channels.channel = chs;
                                x = struct2xml(chan_structure);
                                
                                
                                sn = {''};
                                if CompNumber(ind) <= Nmaps/2
                                    n = CompNumber(ind);
                                else
                                    n = CompNumber(ind)-Nmaps-1;
                                end
                                while isempty(sn{1})
                                    sn = inputdlg('Enter derived signal name','Derived signal name',1,{strcat('DS_', num2str(band(1)),'-',num2str(band(2)),'_',num2str(n))});
                                end
                                
                                if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                                    mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                                end
                                pth = (strcat(self.path,'\',self.subject_record.subject_name));
                                full_name = [pth '\CSP_' sn{1} '.xml'];
                                if exist(full_name,'file')
                                    choice = questdlg('The spatial filter for this person already exists. Rewrite the filter?','Rewrite?','Yes','No, use the old one','No, use the old one');
                                    switch choice
                                        case 'Yes'
                                            f = fopen(full_name,'w');
                                            fwrite(f,x);
                                            fclose(f);
                                            spatial_filter = chan_w;
                                        case 'No, use the old one'
                                            
                                            s = xml2struct(full_name);
                                            channels_coeff = cell(length(s.channels.channel),2);
                                            for i = 1:length(s.channels.channel)
                                                channels_coeff{i,1} = strtrim(s.channels.channel{i}.channel_name.Text);
                                                channels_coeff{i,2} = str2double(strtrim(s.channels.channel{i}.coefficient.Text));
                                            end
                                            spatial_filter = channels_coeff;
                                    end
                                else
                                    spatial_filter = chan_w;
                                    f = fopen(full_name,'w');
                                    fwrite(f,x);
                                    fclose(f);
                                end
                                
                                
                            catch
                                'Error while writing to file, function UpdateStatistics' %#ok<NOPRT>
                            end
                            
                            
                            try
                               
                                self.derived_signals{end+1} = CreateNewDS(self,sn{1},spatial_filter, band,full_name);
                                
                                %%temp derived signal to calculate stats
                                temp_derived_signal = CreateNewDS(self,sn{1},spatial_filter, band);
                               
                                %%calculating stats
                                %j = 1;
                                x = zeros(size(x_raw,1),length(self.derived_signals{1}.spatial_filter));
                                for i = 1:length(self.derived_signals{1}.spatial_filter)
                                    if self.derived_signals{1}.spatial_filter(i)
                                        x(:,i) = x_raw(:,i);
                                       % j = j+1;
                                    else
                                        x(:,i) = 0;
                                    end
                                end
                                temp_derived_signal.Apply(x',1);
                                values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                                self.feedback_manager.average(ind) = mean(values);
                                self.feedback_manager.standard_deviation(ind) = std(values);
                                
                                
                            catch
                                'Error while creating a new derived signal, function UpdateStatistics' %#ok<NOPRT>
                            end
                        end
                    end
                    
                    
                    
                    self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                    self.fb_manager_set = 1;
                    
                else %not csp, not ssd, mb simple baseline
                    % fetches the data from all derived_signals except raw and
                    % calculates their statistics
                    % though why do it if we need stats of feedback values, not of the
                    % 'raw' derived signals
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if(N>0) && length(self.derived_signals) > 1
                        for s = 2:length(self.derived_signals)
                            if self.derived_signals{s}.collect_buff.lst - N+1 < self.derived_signals{s}.collect_buff.fst
                                values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.fst:self.derived_signals{s}.collect_buff.lst,:);
                            else
                                values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.lst - N+1:self.derived_signals{s}.collect_buff.lst,:);
                            end
                            self.feedback_manager.average(s-1) = mean(values);
                            self.feedback_manager.standard_deviation(s-1) = std(values);
                        end
                        
                        self.fb_statistics_set = 1;
                        self.SetRawYTicks;
                        self.SetDSYTicks;
                        self.yscales_fixed = 1;
                        self.raw_yscale_fixed = 1;
                        self.ds_yscale_fixed = 1;
                    end;
                    
                end
                
                
            end
            
            
        end
        function Receive(self,timer_obj, event) %#ok<INUSD>
            
            [sample, timestamp] = self.inlet.pull_chunk(); %#ok<ASGLU>
            self.nd = [self.nd sample];
            sz = size(self.nd,2);
            
            if (sz >= self.current_window_size)
                for ds = 1:length(self.derived_signals)
                    self.derived_signals{ds}.Apply(self.nd(:,1:self.current_window_size),self.recording);
                end
                self.nd =self.nd(:,self.current_window_size+1:end);
                try
                    self.UpdateFeedbackSignal;
                catch
                    'Error while updating Feedback' %#ok<NOPRT>
                end
                
                
                self.samples_acquired = self.samples_acquired+self.current_window_size;
                if(self.current_protocol>0 && self.current_protocol <= length(self.feedback_protocols))
                    self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.feedback_protocols{self.current_protocol}.actual_protocol_size +self.current_window_size;
                    
                    if self.feedback_protocols{self.current_protocol}.actual_protocol_size + self.current_window_size >= self.feedback_protocols{self.current_protocol}.protocol_size
                        self.protocol_indices(self.current_protocol+1,:) = self.derived_signals{1}.collect_buff.lst -self.derived_signals{1}.collect_buff.fst +1;
                        try
                            if self.current_protocol == self.ssd
                                self.Prepare_CSP();
                            elseif self.feedback_protocols{self.current_protocol}.to_update_statistics
                                self.Update_Statistics();
                            end
                            temp_log_str = get(self.log_text,'String');
                            temp_log_str{end+1} = self.feedback_protocols{self.current_protocol}.show_as;
                            set(self.log_text,'String', temp_log_str);
                            
                        catch
                            'Error while updating statistics or setting log_string, function Receive' %#ok<NOPRT>
                        end
                        
                        try
                            if self.feedback_protocols{self.current_protocol}.stop_after
                                
                                set(self.connect_button, 'String', 'Start recording');
                                set(self.connect_button, 'Callback',@self.StartRecording);%%%%%
                                self.StopRecording();
                            else
                                
                                self.current_protocol = self.next_protocol;
                                self.next_protocol = self.next_protocol + 1;
                                
                                if self.current_protocol > length(self.feedback_protocols)
                                    if self.looped
                                        self.current_protocol = 1;
                                        self.next_protocol = 2;
                                        for pr = 1:length(self.feedback_protocols)
                                            self.feedback_protocols{pr}.actual_protocol_size = 0;
                                        end
                                    else
                                        self.StopRecording();
                                    end
                                end
                            end
                        catch
                            'Error while StopRecording, function Receive' %#ok<NOPRT>
                        end
                        if strcmp(self.timer_new_data.Running,'off')
                            self.inlet = lsl_inlet(self.streams{1});
                            self.InitTimer();
                        end
                    end;
                end
            end;
            
        end
        function Connect(self,predicate, value)
            if self.from_file && any([~isempty(self.fnames),~isempty(self.files_pathname)])
                RecordedStudyToLSL(self.fnames,self.files_pathname,self.looped);
            elseif self.from_file
                warning('No files selected')
                return;
            end
            lsllib = lsl_loadlib();
            disp('Connecting...')
            fl = 1;
            while fl<5
                self.streams = lsl_resolve_byprop(lsllib,predicate, value);%found
                if ~isempty(self.streams)
                    break;
                end
                disp(strcat('Streams were not found yet, attempt ', num2str(fl), '/', num2str(5)) )
                pause(0.5);
                fl = fl+1;
            end
            
            if fl == 5 && isempty(self.streams)
                %                      disp('Please make sure that the hardware is plugged and software running.'
                %                      'Or try to change pair "predicate/value" and re-run the script')
                %                 end
                in = questdlg('Select the source ','Choose the source', 'Hardware', 'File','File');
                switch in
                    case 'Hardware'
                        disp('Please make sure that the hardware is plugged and software running./nOr try to change pair "predicate/value" and re-run the script')
                        self.RunInterface(predicate,value);
                    case 'File'
                        self.from_file = 1;
                        
                        self.RunInterface(predicate,value);
                        %                         while isempty(self.streams)
                        %                             self.streams = lsl_resolve_byprop(lsllib,predicate, value);%found
                        %                         end
                end
                
            end
            disp('Connected')
            if length(self.streams) > 1
                warning('The pair predicate/value matches more than one channel. The results might be inconsistent. You might want to restart MATLAB');
                
            end
            self.sampling_frequency = self.streams{1}.nominal_srate();
            
            self.inlet = lsl_inlet(self.streams{1});
            if exist(strcat(self.program_path,'\','channels.txt'), 'file')
                delete(strcat(self.program_path,'\','channels.txt'))
            end
            
            disp('Trying to read the channels... ');
            cd(self.program_path);
            if strcmp(self.streams{1}.name,'File')
                winopen('channels_shcut.lnk') %%channels_shcut contain 'type' and 'Data' in 'object' field
                % change it to obtain channels from another source
                while ~exist('channels.txt','file')
                    pause(0.01)
                end
                self.channel_labels = read_channel_file('channels.txt');
                
            else
                command = ['resolve_channels8.exe' ' ' predicate ' ' value];
                status = system(command);
                if ~status
                    self.channel_labels = read_channel_file('channels.txt');
                    
                end
            end
            
            
            
            
            
            
            
            
            %channels = derive_channel_labels(self.streams{1});
            
            
            self.plot_size = self.plot_length * self.sampling_frequency;
            
            %%set durations and window size based on sampling frequency
            for pr = 1:length(self.feedback_protocols)
                self.feedback_protocols{pr}.Recalculate(self.sampling_frequency);
            end
            if ~self.from_file
                for i = 1:length(self.feedback_protocols)
                    if any([strfind(lower(self.feedback_protocols{i}.protocol_name),'ssd'), strfind(lower(self.feedback_protocols{i}.protocol_name),'csp')]);
                        self.ssd = i;
                    end
                    
                    
                    try %#ok<TRYNC>
                        if self.feedback_protocols{i}.window_size
                            self.feedback_manager.window_size(i) = self.feedback_protocols{i}.window_size;
                        end
                    end
                end
            end
            for j = 1:length(self.feedback_protocols)
                self.exp_data_length = self.exp_data_length + self.feedback_protocols{j}.protocol_size;
            end
            self.exp_data_length = fix(self.exp_data_length * 1.1); %just in case
            for i = 1:length(self.feedback_manager.window_size)
                if self.feedback_manager.window_size(i) <=5
                    warning('Given that the sampling frequency is %d ,the window length of protocol %s is less than 6 samples. Set the window size at least %d ms', self.sampling_frequency,self.feedback_protocols{i}.protocol_name, (5000/self.sampling_frequency+1))
                end
                if self.feedback_manager.window_size(i)/self.sampling_frequency <= self.data_receive_rate
                    warning('The window size of protocol %s is too small. Increase the window size or decrease data receive rate', self.feedback_protocols{i}.protocol_name)
                end
            end
            %set ds
            if self.from_file
                %                 dummy_signal = struct();
                %                 dummy_signal.sSignalName = 'Raw';
                %                 channels = cell(size(self.channel_labels,2));
                %                 for ch = 1:length(channels)
                %                     channels{ch,1} = self.channel_labels{ch};
                %                     channels{ch,2} = 1;
                %                 end
                %                 dummy_signal.channels = channels;
                %                 dummy_signal.filters = [];
                
                if self.ssd
                    self.derived_signals{1} = self.CreateNewDS('Raw',ones(length(self.channel_labels),1));
                else
                     for i = 1: length(self.signals)
                    
                    self.derived_signals{i} = DerivedSignal(1,self.signals{i}, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                    
                    self.derived_signals{i}.UpdateSpatialFilter(self.signals{i}.channels,self.derived_signals{1},self.bad_channels);
                    if ~isempty(self.signals{i}.filters)
                    self.derived_signals{i}.UpdateTemporalFilter(size(self.signals{i}.channels,2)-1,self.signals{i}.filters.range,self.signals{i}.filters.order,self.signals{i}.filters.mode);
                    end
                     end
                end
                %self.derived_signals{1} = DerivedSignal(1,dummy_signal, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                %self.derived_signals{1}.UpdateSpatialFilter(ones(length(self.channel_labels),1),self.channel_labels);
            else
                if self.ssd
                    self.derived_signals = cell(1,1);
                else
                    self.derived_signals = cell(1,length(self.signals));
                end
                for i = 1: length(self.derived_signals)
                    
                    self.derived_signals{i} = DerivedSignal(1,self.signals{i}, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                    
                    self.derived_signals{i}.UpdateSpatialFilter(self.signals{i}.channels,self.derived_signals{1},self.bad_channels);
                    
                end
                
            end
            for ds = 1:length(self.derived_signals)
                if strcmpi(self.derived_signals{ds}.signal_name, 'raw')
                    raw = self.derived_signals{ds};
                    self.used_ch = raw.channels;
                end
            end
            self.current_window_size = self.feedback_manager.window_size(1);
            %self.RunInterface;
            if self.from_file
                self.StartRecording();
            end
            %%bluetooth
            try %#ok<TRYNC>
                b  = Bluetooth('NXT');
                self.bluetooth_connection = fopen(b);
                self.fb_refresh_rate = 0.001;
            end
            %run the timers
            if ishandle(self.fig_interface)
                tic
                self.InitTimer();
            end
        end
        function InitTimer(self)
            if strcmp(self.timer_new_data.Running,'off')
                start(self.timer_new_data);
            end
            if strcmp(self.timer_disp.Running,'off')
                start(self.timer_disp);
            end
            if strcmp(self.timer_fb.Running,'off')
                start(self.timer_fb);
            end
            
            set(self.connect_button, 'String', 'Start recording');
            set(self.connect_button, 'Callback',@self.StartRecording);
        end
        function RunInterface(self,predicate,value)
            %read subjects file
            if exist(strcat(self.program_path,'\subjects.txt'),'file')
                subjects_file = fopen(strcat(self.program_path,'\subjects.txt'));
                subjects = {};
                
                subjects{end+1} = fgetl(subjects_file);
                while ischar(subjects{end})
                    subjects{end+1} = fgetl(subjects_file);
                end
                subjects = subjects(1:end-1);
            else
                subjects = {};
            end
            subjects = [{'Null','Add a new subject'} sort(subjects)];
            if verLessThan('matlab','8.4.0')
                self.fig_interface = figure('CloseRequestFcn',@self.DoNothing);
            else
                self.fig_interface = figure;
                self.fig_interface.CloseRequestFcn = @self.DoNothing;
            end
            
            prr_text = uicontrol('Parent',self.fig_interface, 'Style', 'text', 'String', 'Plot refresh rate, s', 'Position',[20 250 120 30],'HorizontalAlignment','left'); %#ok<NASGU>
            prr = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', num2str(self.plot_refresh_rate), 'Position', [125 260 50 20]);
            drr_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Data receive rate, s', 'Position', [20 210 100 30],'HorizontalAlignment','left'); %#ok<NASGU>
            drr = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', num2str(self.data_receive_rate), 'Position', [125 220 50 20]);
            frr_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Feedback refresh rate, s', 'Position', [20 190 100 30],'HorizontalAlignment','left'); %#ok<NASGU>
            frr = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', num2str(self.fb_refresh_rate), 'Position', [125 190 50 20]);
            self.path_text =uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.path,'Position', [120 125 200 35],'HorizontalAlignment','left');
            path_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select path', 'Callback', @self.SetWorkpath, 'Position', [20 135 100 35]); %#ok<NASGU>
            self.settings_file_text =uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.settings_file,'Position', [120 90 200 35],'HorizontalAlignment','left');
            settings_file_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select exp.design', 'Callback', @self.SetDesignFile, 'Position', [20 100 100 35]); %#ok<NASGU>
            set_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Run the experiment', 'Position', [100 20 200 40],'Callback','uiresume','Tag','set_button'); %#ok<NASGU>
            %montage_file_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select exp. montage', 'Callback', @self.SetMontageFile, 'Position', [20 60 100 35]); %#ok<NASGU>
            %self.montage_fname_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.montage_fname_text,'Position', [120 60 200 35],'HorizontalAlignment','left');
            show_feedback = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Show feedback to subject','Position', [20 290 135 20],'HorizontalAlignment','left'); %#ok<NASGU>
            show_fb_check = uicontrol('Parent', self.fig_interface, 'Style', 'checkbox' ,'Position', [160 295 20 20],'HorizontalAlignment','left','Value',1);
            self.subjects_dropmenu = uicontrol('Parent', self.fig_interface,'Style','popupmenu','Position',[170 320 100 20],'String',subjects,'Callback',@self.SetSubject);
            sn_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Choose/Enter subject name', 'Position',[20 315 140 20],'HorizontalAlignment','left'); %#ok<NASGU>
            subj_folder_button = uicontrol('Parent', self.fig_interface,'Style','pushbutton','Position',[285 320 150 20],'String','Or select subject folder','Callback',@self.SetSubjectFolder); %#ok<NASGU>
            from_file_chb =  uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 273 200 20],'Tag','from_file_chb','String','From file');
            %from_file_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[285 270 50 20],'String','From file'); %#ok<NASGU>
            loop_replay_chb = uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 253 200 20],'Tag','From_file_loop','String','Loop the recording?');
            %loop_replay_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[220 250 120 20],'String','Loop the recording?'); %#ok<NASGU>
            use_protocols_chb = uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 233 200 20],'Value',1,'Tag','Protocols from files','String','Use protocols from files?');
            %use_protocols_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[200 230 140 20],'String','Use protocols from files?'); %#ok<NASGU>
            uiwait();
            if ishandle(self.fig_interface) %if the window is not closed
                if verLessThan('matlab','8.4.0')
                    self.plot_refresh_rate = str2double(get(prr,'String'));
                    self.data_receive_rate = str2double(get(drr,'String'));
                    self.fb_refresh_rate = str2double(get(frr,'String'));
                    
                    set(self.fig_interface,'Visible', 'off');
                    self.show_fb = get(show_fb_check, 'Value');
                    self.from_file = get(from_file_chb,'Value');
                    if self.from_file
                        
                        self.looped = get(loop_replay_chb,'Value');
                        self.run_protocols = get(use_protocols_chb,'Value');
                    end
                    
                else
                    self.plot_refresh_rate = str2double(prr.String);
                    self.data_receive_rate = str2double(drr.String);
                    self.fb_refresh_rate = str2double(frr.String);
                    set(self.fig_interface,'Visible','off');
                    self.show_fb = get(show_fb_check,'Value');
                    self.from_file = get(from_file_chb,'Value');
                    if self.from_file
                        self.looped = loop_replay_chb.Value;
                        self.run_protocols = use_protocols_chb.Value;
                    end
                end
                subjects = [subjects {self.subject_record.subject_name}];
                subjects = sort(unique(subjects));
                subjects_file = fopen(strcat(self.program_path,'\subjects.txt'),'wt');
                for s = 1:length(subjects)
                    if ~(strcmp(subjects{s},'Add a new subject')|| strcmp(subjects{s},'Null'))
                        fprintf(subjects_file,'%s\n',subjects{s});
                    end
                end
                fclose(subjects_file);
                self.timer_new_data_function = @self.Receive;
                self.timer_new_data = timer('Name','receive_data','TimerFcn', self.timer_new_data_function,'ExecutionMode','fixedRate',...
                    'Period',self.data_receive_rate);
                self.timer_disp_function = @self.PlotEEGData;
                self.timer_disp = timer('Name','plot_data','TimerFcn', self.timer_disp_function,'ExecutionMode','fixedRate',...
                    'Period', self.plot_refresh_rate);
                self.timer_fb_function = @self.RefreshFB;
                self.timer_fb = timer('Name','refresh_feedback','TimerFcn',self.timer_fb_function,'ExecutionMode','fixedRate',...
                    'Period',self.fb_refresh_rate);
                self.feedback_manager = FeedbackManager;
                if self.from_file
                    predicate = 'name';
                    value = 'File';
                    self.looped = loop_replay_chb.Value;
                    self.run_protocols = use_protocols_chb.Value;
                    %self.channel_labels = get_channel_labels(self.inlet);
                    
                    [self.fnames, self.files_pathname, filterindex] = uigetfile('.bin','Select files to play','MultiSelect','on');
                    if any([~isempty(self.fnames), ~isempty(self.files_pathname), filterindex] )
                        %subject_folder = self.streams{1}.source_id;
                        [protocols, protocols_show_as, durations, channels,settings_file] = GetDataProperties(self.files_pathname,self.fnames);
                        self.channel_labels = channels;
                        if self.run_protocols
                            self.settings_file = settings_file;
                            self.protocol_sequence = protocols;
                            self.feedback_protocols = [];
                            for pr = 1:length(protocols)
                                
                                self.feedback_protocols{pr} = RealtimeProtocol;
                                self.feedback_protocols{pr}.protocol_name = protocols{pr};
                                self.feedback_protocols{pr}.protocol_duration = durations(pr);
                                if ~isempty(protocols_show_as)
                                    self.feedback_protocols{pr}.show_as = protocols_show_as{pr};
                                else
                                    self.feedback_protocols{pr}.show_as = protocols{pr};
                                end
                                if strfind(lower(protocols{pr}),'ssd')
                                    self.ssd = pr;
                                    self.feedback_protocols{pr}.band = 1;
                                elseif strfind(lower(protocols{pr}),'csp')
                                    self.ssd = pr;
                                    self.feedback_protocols{pr}.band = 1;
                                elseif strcmpi(protocols{pr},'baseline')
                                    self.feedback_protocols{pr}.to_update_statistics = 1;
                                elseif strfind(lower(protocols{pr}),'feedback')
                                    self.feedback_protocols{pr}.fb_type = protocols{pr};
                                end
                            end
                            nfs = NeurofeedbackSession;
                            nfs.LoadFromFile(self.settings_file);
                            self.protocol_types = nfs.protocol_types;
                            %self.feedback_protocols = nfs.feedback_protocols;
                            self.csp_settings = nfs.csp_settings;
                            self.signals = nfs.derived_signals;
                            
                        end
                    end
                    
                else
                    %self.channel_labels = read_montage_file(self.montage_fname);
                    nfs = NeurofeedbackSession;
                    nfs.LoadFromFile(self.settings_file);
                    self.protocol_types = nfs.protocol_types;
                    self.feedback_protocols = nfs.feedback_protocols;
                    self.signals = nfs.derived_signals;
                    self.protocol_sequence = nfs.protocol_sequence;
                    self.csp_settings = nfs.csp_settings;
                    self.feedback_manager.window_size = zeros(length(self.feedback_protocols),1);
                end
                self.protocol_indices = zeros(length(self.feedback_protocols)+1,2); %catch actual data length since 'act_protocol_size' can lie
                
                figure(self.fig_feedback);
                set(self.fig_feedback, 'OuterPosition', [-1280 0 1280 1024]);
                self.feedback_axis_handle = axes;
                self.fb_stub = uicontrol('Parent', self.fig_feedback, 'String', 'Baseline acquisition', 'Style', 'text', 'ForegroundColor',[0 1 0],'Position', [200 500 900 250], 'FontSize', 75, 'BackgroundColor',[1 1 1], 'FontName', 'Courier New', 'Visible', 'off' );
                self.fbplot_handle = bar(self.feedback_axis_handle,[0 1 0],'FaceColor',[1 1 1]);
                self.Connect(predicate,value);
            end
        end
        function PlotEEGData(self,timer_obj, event) %#ok<INUSD>
            
            if ~self.connected
                self.raw_and_ds_figure = figure('Tag','raw_and_ds_figure'); %add Tag
                set(self.raw_and_ds_figure,'ResizeFcn',@self.FitFigure);
                self.connect_button =  uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [10 10 150 20], ...
                    'String', 'Start recording','Tag','connect_button');
                self.disconnect_button = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [420 10 130 20], ...
                    'String', 'Disconnect', 'Callback', @self.Disconnect,'Tag','disconnect_button');
                self.log_text = uicontrol('Parent', self.raw_and_ds_figure  ,'Style', 'Text','String', {'Log'}, 'Position', [0 300 50 100],'Tag','log_text');
                self.status_text = uicontrol('Parent', self.raw_and_ds_figure,'Style', 'text', 'String', 'Status: ', 'Position', [0 210 200 20],'HorizontalAlignment','left','Tag','status_text');
                self.curr_protocol_text = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'text','String', 'Current protocol: ', 'Position', [0 40  190 100],'Tag','curr_protocol_text');
                %self.edit_protocols_button = uicontrol('Parent',self.raw_and_ds_figure,'Style','pushbutton','Callback',@self.EditProtocols,'Tag','edit_protocols_button','String','Edit protocols');
                select_bad_channels_button = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton', ...
                    'String', 'Select bad channels', 'Callback', @self.SelectBadChannels,'Tag','select_bad_channels_button'); %#ok<NASGU>
                %                 bad_channels_text = uicontrol('Parent', self.raw_and_ds_figure,'Style', 'text', 'String', '',...
                %                     'HorizontalAlignment','left','Tag','bad_channels_text');
                
                self.raw_subplot = subplot(2,1,1);
                set(self.raw_subplot,'YLim', [0, self.raw_shift*(length(self.used_ch)+1)]);
                self.raw_scale_slider = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'slider', 'String','Raw scale', 'Value', 0, 'Position', [520 300 10 100], 'Max', 24, 'Min',-24,'SliderStep',[1 1],'Callback',@self.SetYScale,'Tag','raw_slider');
                
                self.raw_data_indices = 1:length(self.used_ch);
                r_temp = zeros(length(self.raw_data_indices),fix(self.plot_size));
                self.raw_plot = plot(r_temp', 'Parent', self.raw_subplot);
                self.raw_line = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [480 320 100 25],'Tag', 'raw_line');
                
                if ~self.raw_ylabels_fixed
                    self.r_ytick_labels = {' '};
                    for i = 1:length(self.used_ch)
                        self.r_ytick_labels{end+1} = self.used_ch{i};
                    end
                    self.r_ytick_labels{end+1} = ' ';
                    for i = 1:length(self.used_ch)
                        set(self.raw_plot(i),'DisplayName', self.used_ch{i});
                    end
                    
                    
                    self.raw_ylabels_fixed = 1;
                end
                
                
                
                self.ds_subplot = subplot(2,1,2);
                set(self.ds_subplot,'YLim', [0 self.raw_shift*length(self.derived_signals)]);
                
                
                self.FitFigure;
                self.connected = 1;
                
            elseif self.connected
                set(self.raw_subplot,'YLim', [0, self.raw_shift*(length(self.raw_data_indices)+1)]);
                set(self.ds_subplot,'YLim', [0 self.raw_shift*(length(self.derived_signals))]);
                r_sp = get(self.raw_subplot);
                ds_sp = get(self.ds_subplot);
                
                try
                    
                    %plot filtered data
                    if length(self.derived_signals) > 1
                        if ~self.ds_ylabels_fixed && length(self.derived_signals)>1
                            
                            ds_temp = zeros(length(self.derived_signals)-1,fix(self.plot_size));
                            self.ds_plot = plot(ds_temp', 'Parent', self.ds_subplot);
                            self.ds_line = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [480 120 100 25],'Tag','ds_line');
                            self.ds_scale_slider= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'slider', 'String','DS scale', 'Value', 0, 'Position', [520 100 10 100], 'Max', 24, 'Min',-24,'SliderStep',[1 1],'Callback',@self.SetYScale,'Tag','ds_slider');
                            
                            self.ds_ytick_labels = {' '};
                            for i = 2:length(self.derived_signals)
                                self.ds_ytick_labels{end+1} = self.derived_signals{i}.signal_name;
                            end
                            self.ds_ytick_labels{end+1} = ' ';
                            for i = 2:length(self.derived_signals)
                                set(self.ds_plot(i-1),'DisplayName', self.derived_signals{i}.signal_name);
                            end
                            self.ds_ylabels_fixed = 1;
                        end
                        %set dropmenu with a choice of derived signals
                        if ~max(size(findobj('Tag','sn_to_fb_dropmenu')))
                            self.sn_to_fb_string = '';
                            for i = 2:length(self.derived_signals)
                                if i == 2
                                    self.sn_to_fb_string = self.derived_signals{i}.signal_name;
                                else
                                    self.sn_to_fb_string = strcat(self.sn_to_fb_string,'|',self.derived_signals{i}.signal_name);
                                end
                            end
                            self.sn_to_fb_dropmenu = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'popupmenu', 'String', self.sn_to_fb_string, 'Position',[300 10 100 20], 'Callback', @self.SelectSignalToFeedback,'Tag','sn_to_fb_dropmenu');
                        end
                        
                        %plot the ds data
                        ds_first_to_show = self.derived_signals{self.signal_to_feedback}.ring_buff.lst-self.plot_size;
                        ds_last_to_show = self.derived_signals{self.signal_to_feedback}.ring_buff.lst;
                        if ds_first_to_show < ds_last_to_show
                            for i = 2:length(self.derived_signals)
                                %pulled = self.derived_signals{i}.ring_buff.raw(self.derived_signals{self.signal_to_feedback}.ring_buff.fst:self.derived_signals{self.signal_to_feedback}.ring_buff.lst,:);
                                pulled = self.derived_signals{i}.ring_buff.raw(ds_first_to_show+1:ds_last_to_show,:);
                                pulled = pulled';
                                set(self.ds_plot(i-1), 'YData',pulled(1,:)*self.ds_ydata_scale+self.ds_shift*(i-1));
                            end
                        end
                    end
                    
                    %plot raw data
                    raw_first_to_show = self.derived_signals{1}.ring_buff.lst-self.plot_size;
                    raw_last_to_show = self.derived_signals{1}.ring_buff.lst;
                    if raw_last_to_show > raw_first_to_show
                        raw_data = self.derived_signals{1}.ring_buff.raw(self.derived_signals{1}.ring_buff.lst-self.plot_size+1:self.derived_signals{1}.ring_buff.lst,:);
                        raw_data = raw_data';
                        raw_data = raw_data(self.raw_data_indices,:);
                        for i = 1:size(raw_data,1)
                            set(self.raw_plot(i),'YData', raw_data(i:i,:)*self.raw_ydata_scale+self.raw_shift*i);
                        end
                    end
                    if ~self.raw_yscale_fixed
                        self.SetRawYTicks;
                        self.raw_yscale_fixed = 1;
                    end
                    
                    %set x limits of the plots
                    xlim(self.ds_subplot, [0 self.plot_size]);
                    xlim(self.raw_subplot, [0 self.plot_size]);
                    set(self.raw_subplot, 'XTick', [0:self.sampling_frequency:self.plot_size]); %#ok<NBRAK>
                    set(self.ds_subplot, 'XTick', [0:self.sampling_frequency:self.plot_size]); %#ok<NBRAK>
                    if self.samples_acquired > self.plot_size
                        set(self.ds_subplot, 'XTickLabel', [self.samples_acquired - ds_sp.XTick(end):ds_sp.XTick(2):self.samples_acquired]); %#ok<NBRAK>
                        set(self.raw_subplot, 'XTickLabel', [self.samples_acquired - r_sp.XTick(end):r_sp.XTick(2):self.samples_acquired]); %#ok<NBRAK>
                    end
                catch
                    'Error while plotting, function PlotEEGData' %#ok<NOPRT>
                end
            end
            %set recording status
            try
                if(self.current_protocol> 0 && self.current_protocol<=length(self.feedback_protocols)) %non-zero protocol
                    if verLessThan('matlab','8.4.0')
                        set(self.curr_protocol_text, 'String', [strcat('Current protocol: ',self.feedback_protocols{self.current_protocol}.show_as),...
                            strcat('Samples acquired', num2str(self.feedback_protocols{self.current_protocol}.actual_protocol_size),'/', num2str(self.feedback_protocols{self.current_protocol}.protocol_size)),...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            ]);
                    else
                        self.curr_protocol_text.String = [strcat('Current protocol: ',self.feedback_protocols{self.current_protocol}.show_as),...
                            strcat('Samples acquired', num2str(self.feedback_protocols{self.current_protocol}.actual_protocol_size),'/', num2str(self.feedback_protocols{self.current_protocol}.protocol_size)),...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            ];
                    end
                else %zero protocol
                    if verLessThan('matlab','8.4.0')
                        set(self.curr_protocol_text, 'String', ['Current protocol: idle, ',...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            ]);
                    else
                        self.curr_protocol_text.String = ['Current protocol: idle, ',...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate)), ' s'];
                    end
                end
            catch
                'Error while setting the status string, function PlotEEGData' %#ok<NOPRT>
            end
            
            self.SetRecordingStatus;
        end
        function RefreshFB(self,timer_obj,event) %#ok<INUSD>
            %feedback
            if length(self.derived_signals) > 1
                self.fig_feedback;
                set(self.feedback_axis_handle,'Visible','off'); %axes
                if self.current_protocol && self.current_protocol <= length(self.feedback_protocols)
                    if isprop(self.feedback_protocols{self.current_protocol},'fb_type') && ~isempty(self.feedback_protocols{self.current_protocol}.fb_type)
                        try
                            self.fb_type = self.feedback_protocols{self.current_protocol}.fb_type;
                        catch
                            'Error while getting fb_type, function RefreshFB' %#ok<NOPRT>
                        end
                    end
                end
                try
                    if (self.current_protocol> 0 && self.current_protocol<=length(self.feedback_protocols))
                        set(self.fb_stub,'String',self.feedback_protocols{self.current_protocol}.string_to_show);
                        if isempty(get(self.fb_stub,'String')) && self.show_fb %feedback
                            set(self.fb_stub, 'Visible', 'off'); %string
                            if strfind(lower(self.fb_type),'bar')
                                set(self.fbplot_handle,'Visible','on'); %feedback if bar
                                %self.y_limit = [-1 7];
                                xlim(self.feedback_axis_handle, [1 3]);
                                ylim(self.feedback_axis_handle,self.y_limit);
                                set(self.fbplot_handle,'FaceColor',[1 0 0]);
                                set(self.fbplot_handle,'EdgeColor',[0 0 0]);
                                if ~isnan(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))
                                    set(self.fbplot_handle,'YData',[0 self.feedback_manager.feedback_vector(self.signal_to_feedback-1) 0]);
                                end
                                
                            elseif strfind(lower(self.fb_type),'color')
                                set(self.fbplot_handle,'Visible', 'off');
                                if ~isnan(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))
                                    set(self.fig_feedback,'Color',[1 1-(1/16+1/16*self.feedback_manager.feedback_vector(self.signal_to_feedback-1)) 1-(1/16+1/16*self.feedback_manager.feedback_vector(self.signal_to_feedback-1))]);
                                end
                            elseif strfind(lower(self.fb_type),'mock')
                                set(self.fbplot_handle,'Visible','off');
                                if ~isnan(self.feedback_manager.average(self.signal_to_feedback-1) && ~isnan(self.feedback_manager.standard_deviation(self.signal_to_feedback-1)))
                                    mock_feedback = random('Normal',self.feedback_manager.average + 3*self.feedback_manager.standard_deviation,self.feedback_manager.standard_deviation,1,1);
                                    set(self.fig_feedback,'Color',[1 1-(1/16+1/16*mock_feedback) 1-(1/16+1/16*mock_feedback)]);
                                end
                            end
                        elseif strcmp(self.feedback_protocols{self.current_protocol}.protocol_name,'Rest')
                            set(self.fig_feedback,'Color',[0.94 0.94 0.94]);
                            set(self.fb_stub, 'Visible', 'off'); %string
                            set(self.fbplot_handle,'Visible','off');
                            
                        else %strings must be visible
                            set(self.fig_feedback,'Color',[0.94 0.94 0.94]);
                            set(self.fbplot_handle,'Visible','off');
                            set(self.fb_stub,'Visible', 'on'); %string
                            set(self.fbplot_handle,'FaceColor',[1 1 1]);
                            set(self.fbplot_handle,'EdgeColor','none');
                        end
                    elseif self.fb_statistics_set && self.show_fb %zero protocol after baseline recorded
                        set(self.fb_stub, 'Visible', 'off'); %string
                        if strfind(lower(self.fb_type),'bar')
                            set(self.fbplot_handle,'Visible','on'); %feedback if bar
                            %self.y_limit = [-1 7];
                            xlim(self.feedback_axis_handle, [1 3]);
                            ylim(self.feedback_axis_handle,self.y_limit);
                            set(self.fbplot_handle,'FaceColor',[1 0 0]);
                            set(self.fbplot_handle,'EdgeColor',[0 0 0]);
                            if ~isnan(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))
                                set(self.fbplot_handle,'YData',[0 self.feedback_manager.feedback_vector(self.signal_to_feedback-1) 0]);
                            end
                        elseif strfind(lower(self.fb_type),'color')
                            set(self.fbplot_handle,'Visible','off'); %feedback if bar
                            if ~isnan(self.feedback_manager.feedback_vector(self.signal_to_feedback-1)) && self.fb_sigmas*self.feedback_manager.feedback_vector(self.signal_to_feedback-1) > 0
                                set(self.fig_feedback,'Color',[1 1-(1/self.fb_sigmas+1/self.fb_sigmas*self.feedback_manager.feedback_vector(self.signal_to_feedback-1)) 1-(1/self.fb_sigmas+1/self.fb_sigmas*self.feedback_manager.feedback_vector(self.signal_to_feedback-1))]);
                            end
                        else
                            set(self.fig_feedback,'Color',[0.94 0.94 0.94]);
                            set(self.fb_stub, 'Visible', 'off'); %string
                            set(self.fbplot_handle,'Visible','off'); %fb if bar
                        end
                    else %zero protocol before baseline recorded
                        set(self.fig_feedback,'Color',[0.94 0.94 0.94]);
                        set(self.fbplot_handle,'Visible','off');
                        %set(self.fbplot_handle,'FaceColor',[1 1 1],'EdgeColor','none');
                        set(self.fb_stub,'Visible','off'); %string
                    end
                    
                catch
                    self.feedback_protocols{self.current_protocol}.protocol_name
                    'Error while setting feedback' %#ok<NOPRT>
                    1/self.fb_sigmas+1/self.fb_sigmas*self.feedback_manager.feedback_vector(self.signal_to_feedback-1) %#ok<NOPRT>
                end
            elseif self.current_protocol
                if strfind(lower(self.feedback_protocols{self.current_protocol}.protocol_name),'ssd')
                    set(self.fb_stub,'String',self.feedback_protocols{self.current_protocol}.string_to_show);
                    set(self.fb_stub,'Visible','on');
                    set(self.feedback_axis_handle,'Visible','off');
                    set(self.fbplot_handle,'Visible','off');
                end
            else%zero protocol before baseline recorded
                set(self.feedback_axis_handle,'Visible','off');
                set(self.fbplot_handle,'Visible','off');
                set(self.fb_stub,'Visible','off'); %string
            end
            %bluetooth
            self.TransmitToBluetooth();
        end
        function StartRecording(self,obj,event) %#ok<INUSD>
            if self.from_file && ~self.run_protocols
                self.current_protocol = self.next_protocol;
                self.next_protocol = self.next_protocol + 1;
            else
                self.current_protocol = self.next_protocol;
                self.next_protocol = self.next_protocol + 1;
            end
            if self.current_protocol > 0 && self.current_protocol <= length(self.feedback_protocols)
                if self.current_protocol <= length(self.feedback_manager.window_size)
                    self.current_window_size = self.feedback_manager.window_size(self.current_protocol);
                    self.default_window_size =  self.current_window_size;
                else
                    self.current_window_size = self.default_window_size;
                end
            elseif self.default_window_size == 0 %zero protocol at the beginning
                self.current_window_size = self.feedback_manager.window_size(1);
            else  %zero protocol after some data was recorded
                self.current_window_size = self.default_window_size;
            end
            self.paused = 0;
            self.recording = 1;
            self.InitTimer();
            set(self.connect_button, 'String', 'Stop recording');
            set(self.connect_button, 'Callback', @self.StopRecording);
        end
        function StopRecording(self, obj,event) %#ok<INUSD>
            self.recording = 0;
            if self.current_protocol >= length(self.feedback_protocols)
                if self.current_protocol == length(self.feedback_protocols) && self.feedback_protocols{self.current_protocol}.stop_after
                    self.current_protocol = 0;
                    self.tstop = toc;
                elseif self.current_protocol > length(self.feedback_protocols)
                    self.current_protocol = 0;
                end
                self.finished = 1;
                set(self.connect_button,'enable', 'off');
                temp_log_text = get(self.log_text,'String');
                temp_log_text{end+1} = 'Finished';
                set(self.log_text,'String',temp_log_text);
                set(self.disconnect_button,'String','Disconnect and write');
                set(self.connect_button, 'String', 'Recording finished');
                set(self.connect_button,'Callback','');
                %self.PlotFB;
                self.PlotErrorBar;
            end
            if  ~self.finished
                if self.feedback_protocols{self.current_protocol}.actual_protocol_size*1.1 < self.feedback_protocols{self.current_protocol}.protocol_size
                    self.next_protocol = self.current_protocol;
                    set(self.connect_button, 'String', 'Continue recording');
                    self.paused = 1;
                else
                    self.next_protocol = self.current_protocol + 1;
                    set(self.connect_button, 'String', 'Start recording');
                end
                self.current_protocol = 0;
                
                set(self.connect_button, 'Callback', @self.StartRecording);
            end
        end
        function Disconnect(self, obj,event) %#ok<INUSD>
            stop(self.timer_new_data);
            stop(self.timer_disp);
            set(self.fb_stub,'Visible','off');
            set(self.status_text, 'String', 'Status: disconnected');
            set(self.connect_button, 'String', 'Resume recording');
            set(self.connect_button, 'Callback',{@self.Connect});
            self.inlet.close_stream();
            self.subject_record.time_stop = datestr(now,13);
            if self.finished
                self.WriteToFile;
            end
        end
        function PlotErrorBar(self)
            
            %protocols:
            protocols = self.ssd+1:length(self.feedback_protocols);
            if self.ssd
                const_shift = self.protocol_indices(self.ssd+1,2);
            else
                const_shift = 0;
            end
%             if ~any([strfind(lower(self.feedback_protocols{1}.protocol_name),'ssd'),strfind(lower(self.feedback_protocols{1}.protocol_name),'csp')])
%                 pr = length(self.feedback_protocols);
%                 pr_shift = 0;
%                 const_shift = 0;
%             else
%                 pr = length(self.feedback_protocols)-1;
%                 pr_shift = 1;
%                 const_shift = self.protocol_indices(2,2);
%             end

            %signals
            
            averages = zeros(length(protocols),length(self.derived_signals)-1);
            deviations =  zeros(length(protocols),length(self.derived_signals)-1);
            names = cell(length(protocols),1);
            plot_size = length(protocols);
            for i = protocols
                idx11 = self.protocol_indices(i,1);
                idx12= self.protocol_indices(i+1,1);
                idx21 = self.protocol_indices(i,2);
                idx22 = self.protocol_indices(i+1,2);
                if any([strfind(lower(self.feedback_protocols{i}.protocol_name),'ssd'),strfind(lower(self.feedback_protocols{i}.protocol_name),'csp')])
                    if idx22-idx21 ~= idx12-idx11
                        warning('Something went wrong... Function PlotErrorBar')
                    end
                    plot_size = length(self.feedback_protocols)-1;
                else
                    names{i-self.ssd} = self.feedback_protocols{i}.show_as;
                    for ds = 2:length(self.derived_signals)
                        dat = self.derived_signals{ds}.collect_buff.raw(self.derived_signals{ds}.collect_buff.fst+idx21-const_shift: self.derived_signals{ds}.collect_buff.fst+idx22-1-const_shift);
                        fb = self.Recalculate(dat,self.feedback_protocols{i}.window_size, self.feedback_manager.average(ds-1), self.feedback_manager.standard_deviation(ds-1));
                        averages(i-self.ssd,ds-1) = mean(fb);
                        deviations(i-self.ssd,ds-1) = std(fb);
                    end
                    
                    %averages(i-self.ssd) = mean(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst+idx21-const_shift:self.feedback_manager.feedback_records.fst+idx22-1-const_shift,2));
                    %deviations(i-self.ssd) = std(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst+idx21-const_shift:self.feedback_manager.feedback_records.fst+idx22-1-const_shift,2));
                    
                end
            end
            %write to file
            curr_date = datestr(date,29);
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start));
            end
            filename = strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start,'\','stats.txt');
            f = fopen(filename,'w');
            fprintf(f, self.derived_signals{self.signal_to_feedback}.signal_name);
            fprintf(f,'\n');
            fprintf(f,'Protocol average stddev\n');
            for i = 1:length(names)
                fprintf(f, [names{i} ' ' num2str(averages(i)) ' ' num2str(deviations(i)) '\n']);
            end
            fclose(f);
            %show plot
            names = [names' ' '];
            legend_str = {};
            
            f = figure; %#ok<NASGU>
            for ds = 1:length(self.derived_signals)-1
                e = errorbar(averages(:,ds), deviations(:,ds)); %#ok<NASGU>
                hold on;
                legend_str{end+1} = self.derived_signals{ds+1}.signal_name;
            end
            set(gca,'XTick',1:plot_size);
            set(gca,'XTickLabel',names);
            xlabel('Protocols');
            ylabel('Normalized values of calculated feedback (Mean +/- std)');
            legend(legend_str);
            
        end
        function WriteToFile(self)
            curr_date = datestr(date,29);
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start));
            end
            cd(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start));
            string = '';
            for c = 1:length(self.used_ch)
                if c == 1
                    string = self.used_ch{c,1};
                else
                    string = strcat(string, ',', self.used_ch{c,1});
                end
            end
            
            
            for i = 1:length(self.feedback_protocols)
                idx11 = self.protocol_indices(i,1);%create two separate sets of indices in case the first fb_protocol
                idx12= self.protocol_indices(i+1,1);%was ssd or csp and fb during them was not recorded
                idx21 = self.protocol_indices(i,2);
                idx22 = self.protocol_indices(i+1,2);
                
                if any([strfind(lower(self.feedback_protocols{i}.protocol_name),'ssd'),strfind(lower(self.feedback_protocols{i}.protocol_name),'csp')])
                    %correct second indices
                    start = self.protocol_indices(2,2);
                    for j = 2:length(self.feedback_protocols)+1
                        
                        self.protocol_indices(j,2) = self.protocol_indices(j,2) - start;
                    end
                    
                    raw_data_matrix = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx11:self.derived_signals{1}.collect_buff.fst+idx12-1,:);
                    whole_data = raw_data_matrix;
                    if strfind(lower(self.feedback_protocols{i}.protocol_name),'ssd')
                        inf_file = fopen('ssd_exp_info.hdr','w');
                    elseif strfind(lower(self.feedback_protocols{i}.protocol_name),'csp')
                        inf_file = fopen('csp_exp_info.hdr','w');
                    end
                    fprintf(inf_file,string); %basically writes channels' names only
                    fclose(inf_file);
                    
                else
                    if idx22-idx21 ~= idx12-idx11
                        warning('Something went wrong... Function WriteToFile')
                    end
                    raw_data_matrix = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx11:self.derived_signals{1}.collect_buff.fst+idx12-1,:);
                    fb_matrix = self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst+idx21:self.feedback_manager.feedback_records.fst+idx22-1, :);
                    data_matrix = zeros(size(fb_matrix,1), length(self.derived_signals)-1);
                    for j = 2:length(self.derived_signals)
                        data_matrix(:,j-1) = self.derived_signals{j}.collect_buff.raw(self.derived_signals{j}.collect_buff.fst+idx21:self.derived_signals{j}.collect_buff.fst+idx22-1, :);
                    end
                    whole_data = [raw_data_matrix data_matrix fb_matrix];
                    
                end
                
                filename = [num2str(i) ' ' self.feedback_protocols{i}.protocol_name ' ' self.feedback_protocols{i}.show_as ' ' num2str(self.feedback_protocols{i}.actual_protocol_size/self.sampling_frequency) '.bin'];
                %write data
                f = fopen(filename,'w');
                fwrite(f,size(whole_data),'int');
                fwrite(f,whole_data, 'double');
                fclose(f);
                
                
            end
            %write header
            for j = 2:length(self.derived_signals)
                string = strcat(string,',',self.derived_signals{j}.signal_name);
            end
            string = strcat(string,',','Feedbacked signal', ',','Fb values',',','Average',',','Stddev',',','Window size',',','Window num');
            inf_file = fopen('exp_info.hdr','w');
            fprintf(inf_file,string);
            fclose(inf_file);
            %prepare the struct with exp.design and write it
            self.ExpDesignToXML;
            %write notes
            self.AddNotes;
        end
        function AddNotes(self)
            if ~isempty(self.bad_channels)
                notes_string = ['The channels ' strjoin(self.bad_channels) ' were excluded from analysis.'];
            else
                notes_string = '';
            end
            
            self.add_notes_window = figure;
            self.add_notes_field = uicontrol('Parent', self.add_notes_window, 'Style', 'edit', 'Position',...
                [ 10 30 300 200],'String', notes_string); %there's no such thing as VerticalAlignment in uicontrols
            self.write_notes = uicontrol('Parent', self.add_notes_window, 'Style', 'pushbutton', 'Position',[ 150 10 100 20], 'Callback', 'uiresume', 'String', 'Save notes');
            uiwait;
            if verLessThan('matlab','8.4.0')
                notes = get(self.add_notes_field, 'String');
            else
                notes = self.add_notes_field.String;
            end
            f = fopen('notes.txt','w');
            fwrite(f, notes,'char');
            fclose(f);
            close(self.add_notes_window);
            close(self.raw_and_ds_figure);
        end
        function SetRecordingStatus(self)
            if verLessThan('matlab','8.4.0')
                
                if self.from_file && isempty(self.feedback_protocols)
                    set(self.status_text,'String', 'Playing from file');
                elseif self.current_protocol == 0 || self.current_protocol > length(self.feedback_protocols)
                    set(self.status_text, 'String','Status: receiving');
                else
                    set(self.status_text,'String',strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.show_as, ' : ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration)));
                end
                if self.paused
                    set(self.status_text,'String', ['Protocol ' self.feedback_protocols{self.next_protocol}.show_as ' paused.'...
                        num2str(round(self.feedback_protocols{self.next_protocol}.actual_protocol_size/self.sampling_frequency)) '/'...
                        num2str(self.feedback_protocols{self.next_protocol}.protocol_duration)]);
                end
            else
                
                if self.from_file && isempty(self.feedback_protocols)
                    self.status_text.String = 'Playing from file';
                elseif self.current_protocol == 0 || self.current_protocol > length(self.feedback_protocols)
                    self.status_text.String = 'Status: receiving';
                else
                    self.status_text.String = strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.show_as, ': ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration));
                end
                if self.paused
                    set(self.status_text,'String',['Protocol ' self.feedback_protocols{self.next_protocol}.show_as ' paused.'...
                        num2str(round(self.feedback_protocols{self.next_protocol}.actual_protocol_size/self.sampling_frequency)) '/'...
                        num2str(self.feedback_protocols{self.next_protocol}.protocol_duration)]);
                end
            end
        end
        function SetYScale(self,obj,event) %#ok<INUSD>
            if strcmp(get(obj,'String'),'DS scale')
                self.ds_ydata_scale =  2^get(obj,'Value');
                self.SetDSYTicks;
            elseif strcmp(get(obj,'String'),'Raw scale')
                self.raw_ydata_scale = 2^get(obj,'Value');
                self.SetRawYTicks;
            end
        end
        function SetWorkpath(self,~,~)
            p = uigetdir;
            if p
                self.path = uigetdir;
            end
            if verLessThan('matlab','8.4.0')
                set(self.path_text,'String',self.path);
            else
                self.path_text.String = self.path;
            end
        end
        function SetDesignFile(self,~,~)
            [fname, fpath, fspec] = uigetfile('*.*');
            if ~isempty(nonzeros([fname fpath fspec]))
                self.settings_file = strcat(fpath,fname);
                if verLessThan('matlab','8.4.0')
                    set(self.settings_file_text, 'String',self.settings_file);
                else
                    self.settings_file_text.String = self.settings_file;
                end
            end
        end
        function SetMontageFile(self,~,~)
            [fname, fpath, fspec] = uigetfile('*.*');
            if ~isempty(nonzeros([fname fpath fspec]))
                self.montage_fname = strcat(fpath,fname);
                if verLessThan('matlab','8.4.0')
                    set(self.montage_fname_text, 'String',fname);
                else
                    self.montage_fname_text.String = self.montage_fname;
                end
            end
        end
        function SelectSignalToFeedback(self,obj,event) %#ok<INUSD>
            if verLessThan('matlab','8.4.0')
                self.signal_to_feedback = get(obj,'Value')+1;
            else
                self.signal_to_feedback = obj.Value+1;
            end
        end
        function FitFigure(self,obj, event) %#ok<INUSD>
            f = findobj('Tag','raw_and_ds_figure');
            fp = get(f,'Position');
            cb = findobj('Tag','connect_button');
            db = findobj('Tag','disconnect_button');
            dss = findobj('Tag','ds_slider');
            rss = findobj('Tag','raw_slider');
            dm = findobj('Tag','sn_to_fb_dropmenu');
            lt = findobj('Tag','log_text');
            st = findobj('Tag','status_text');
            cpt = findobj('Tag','curr_protocol_text');
            rl = findobj('Tag', 'raw_line');
            dsl = findobj('Tag','ds_line');
            sb = findobj('Tag','settings_button');
            sbch = findobj('Tag','select_bad_channels_button');
            epb = findobj('Tag','edit_protocols_button');
            try %#ok<TRYNC>
                ok = findobj('Tag', 'add_bad_channel_button');
                fin = findobj('Tag','finish adding bad channels button');
                bcht = findobj('Tag','bad_channels_text');
                cht = findobj('Tag','ch_text');
            end
            
            set(db,'Position',[0.85*fp(3), 0.02*fp(4), 0.12*fp(3), 0.04*fp(4)]);
            set(cb,'Position',[0.03*fp(3), 0.02*fp(4), 0.12*fp(3), 0.04*fp(4)]);
            set(dss,'Position',[0.93*fp(3),0.12*fp(4) , 0.02*fp(3), 0.3*fp(4)]);
            set(rss,'Position',[0.93*fp(3),0.60*fp(4) , 0.02*fp(3), 0.3*fp(4)]);
            set(lt,'Position',[0 0.6*fp(4) 0.1*fp(3), 0.4*fp(4)]);
            set(st,'Position', [0 0.49*fp(4), 0.3*fp(3), 0.05*fp(4)]);
            set(cpt, 'Position',[0.2*fp(3) 0.49*fp(4), 0.7*fp(3), 0.05*fp(4)]);
            set(dm,'Position', [0.45*fp(3), 0.015*fp(4),0.12*fp(3),0.04*fp(4)]);
            set(rl,'Position', [0.8 * fp(3), 0.62 *fp(4), 0.05*fp(3), 0.02*fp(4)]);
            set(dsl,'Position', [0.8 * fp(3), 0.15 *fp(4), 0.05*fp(3), 0.02*fp(4)]);
            set(epb,'Position', [0.7*fp(3), 0.945*fp(4), 0.2*fp(3), 0.05*fp(4)]);
            set(sb,'Position', [0.1*fp(3), 0.94*fp(4),0.1*fp(3), 0.05*fp(4)]);
            set(ok,'Position', [0.65*fp(3), 0.94*fp(4), 0.1*fp(3),0.05*fp(4)]);
            set(fin,'Position', [0.75*fp(3), 0.94*fp(4), 0.1*fp(3),0.05*fp(4)]);
            set(bcht, 'Position',[fp(3)*0.2,fp(4)*0.7,fp(3)*0.05,fp(4)*0.2]);
            set(cht,'Position', [fp(3)*0.35, fp(4)*0.94, fp(3)*0.3, fp(4)*0.05]);
            set(sbch,'Position',[fp(3)*0.14, fp(4)*0.945, fp(3)*0.2, fp(4)*0.05]);
            self.SetRawYTicks;
            self.SetDSYTicks;
        end
        function SetRawYTicks(self)
            try %#ok<TRYNC>
                r_sp = get(self.raw_subplot);
                r_yticks = [r_sp.YLim(1):(r_sp.YLim(2)-r_sp.YLim(1))/(length(self.raw_data_indices)+1):r_sp.YLim(2)]; %#ok<NBRAK>
                set(self.raw_subplot, 'YTick', r_yticks);
                set(self.raw_subplot, 'YTickLabel', self.r_ytick_labels);
                set(self.raw_line,'String',num2str((r_sp.YLim(2)-r_sp.YLim(1))/(length(self.raw_data_indices)+1)/self.raw_ydata_scale));
            end
            %self.FitFigure;
        end
        function SetDSYTicks(self)
            try %#ok<TRYNC>
                ds_sp = get(self.ds_subplot);
                ds_yticks = [ds_sp.YLim(1):self.ds_shift:ds_sp.YLim(2)]; %#ok<NBRAK>
                set(self.ds_subplot, 'YTick', ds_yticks);
                set(self.ds_subplot, 'YTickLabel', self.ds_ytick_labels);
                set(self.ds_line,'String',num2str((ds_sp.YLim(2)-ds_sp.YLim(1))/(length(self.derived_signals))/self.ds_ydata_scale));
            end
        end
        function SetSubject(self,obj,event) %#ok<INUSD>
            if strcmp(obj.String{obj.Value},'Add a new subject')
                p = get(obj,'Position');
                p(2) = p(2)-2;
                p(4) = p(4)+2;
                subj_text = uicontrol('Parent',self.fig_interface,'Style','edit','String', '','Position',p);
                subj_tip = uicontrol('Parent',self.fig_interface,'Style','text','String', 'Press Enter when finished', 'Position', [p(1) + p(3), p(2), 150, 20]);
                waitfor(subj_text,'String');
                self.subject_record.subject_name = strtrim(subj_text.String); %remove leading and trailing spaces
                set(obj,'String',[{subj_text.String} obj.String']);
                delete(subj_text);
                delete(subj_tip);
                set(obj,'Value', 1);
            else
                self.subject_record.subject_name = obj.String{obj.Value};
            end
        end
        function SetSubjectFolder(self,obj,event) %#ok<INUSD>
            subj_directory = uigetdir;
            if subj_directory
                [~, b, ~ ] = fileparts(subj_directory);
                if b
                    self.subject_record.subject_name = b;
                    set(self.subjects_dropmenu,'String',[{b} get(self.subjects_dropmenu,'String')']);
                    set(self.subjects_dropmenu,'Value',1);
                end
            end
        end
        function EditProtocols(self,obj,event) %#ok<INUSD>
            if min(length(findobj('Tag', 'EditProtocolFigure'))) %if it already exists, bring it to front
                uistack(findobj('Tag', 'EditProtocolFigure'));
            end
            protocol_figure = figure('Tag','EditProtocolFigure');
            delta_y = protocol_figure.Position(4)/(length(self.feedback_protocols)+3);
            max_height = protocol_figure.Position(4);
            existing_prs = cell(length(self.feedback_protocols),1);
            for p = 1:length(self.feedback_protocols)
                bgr = 0.94-[0.1 0.1 0.1] * mod(p-1,2);
                existing_prs{p} = self.feedback_protocols{p}.protocol_name;
                protocol_count{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.01,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(p),'HorizontalAlignment','left','Tag','Protocol count','BackgroundColor',bgr); %#ok<NASGU>
                protocol_name{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.04,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.25, protocol_figure.Position(4)*0.04],'String', existing_prs{p},'HorizontalAlignment','left','Tag','Protocol name text','BackgroundColor',bgr);
                
                if p < self.next_protocol %already recorded; duration cannot be changed
                    self.protocol_duration_text{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.3, max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(self.feedback_protocols{p}.protocol_duration),'Tag','Protocol duration text','HorizontalAlignment','right');
                else
                    self.protocol_duration_text{p} = uicontrol('Parent',protocol_figure,'Style','edit','Position', [protocol_figure.Position(3)*0.3, max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(self.feedback_protocols{p}.protocol_duration),'HorizontalAlignment','right','Tag', 'Protocol duration text');
                    
                end
                ms{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.35,max_height-protocol_figure.Position(4)*0.05*p,protocol_figure.Position(3)*0.02 , protocol_figure.Position(4)*0.04],'String', 's','HorizontalAlignment','left','BackgroundColor',bgr,'Tag','s text'); %#ok<NASGU>
            end
            
            
            prs_types = {};
            for pt = 1:length(self.protocol_types)
                prs_types = [prs_types {self.protocol_types{pt}.sProtocolName}];
            end
            pr_dpmenu = {};
            for p = self.next_protocol-1:length(existing_prs)
                if p < 1
                    continue;
                end
                
                pr_dpmenu = [pr_dpmenu {[num2str(p) ' ' protocol_name{p}.String]}];
            end
            if self.next_protocol == 1
                
                ins_pr_dpmenu = [{'0'} pr_dpmenu];
                del_pr_dpmenu = pr_dpmenu;
            else
                ins_pr_dpmenu = pr_dpmenu;
                del_pr_dpmenu = pr_dpmenu(2:end);
            end
            
            
            add_protocol_text = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.4,protocol_figure.Position(4)*0.92,protocol_figure.Position(3)*0.20 , protocol_figure.Position(4)*0.05],'String', 'Add a protocol','HorizontalAlignment','right'); %#ok<NASGU>
            insert_protocol_text = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.4,protocol_figure.Position(4)*0.85,protocol_figure.Position(3)*0.20 , protocol_figure.Position(4)*0.05],'String', 'Insert after','HorizontalAlignment','right'); %#ok<NASGU>
            
            add_protocol_dropmenu = uicontrol('Parent',protocol_figure,'Style','popupmenu','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.925,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', prs_types,'Tag','Add protocol dropmenu'); %#ok<NASGU>
            insert_protocol_dropmenu = uicontrol('Parent',protocol_figure,'Style','popupmenu','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.855,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', ins_pr_dpmenu,'Tag', 'Insert protocol dropmenu'); %#ok<NASGU>
            add_protocol_pushbutton =uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.78, protocol_figure.Position(3)*0.2, protocol_figure.Position(4)*0.06],'String', 'Add','Tag','add_protocol_button','Callback',@self.AddProtocol); %#ok<NASGU>
            
            delete_protocol_text = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.4,protocol_figure.Position(4)*0.6,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', 'Delete a protocol','HorizontalAlignment','right'); %#ok<NASGU>
            delete_protocol_pushbutton = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.53, protocol_figure.Position(3)*0.2, protocol_figure.Position(4)*0.06],'String', 'Delete','Tag','delete_protocol_button','Callback',@self.DeleteProtocol);
            if isempty(del_pr_dpmenu)
                
                set(delete_protocol_pushbutton,'enable','off');
                delete_protocol_dropmenu = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.6,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', 'No protocols to delete','Tag','Delete protocol dropmenu'); %#ok<NASGU>
            else
                delete_protocol_dropmenu = uicontrol('Parent',protocol_figure,'Style','popupmenu','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.61,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', del_pr_dpmenu,'Tag','Delete protocol dropmenu'); %#ok<NASGU>
                set(delete_protocol_pushbutton,'enable','on');
            end
            
            
            okay_button = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.75,delta_y/2, protocol_figure.Position(3)*0.09,protocol_figure.Position(4)*0.12],'String', 'Apply','Tag','okay_button','Callback',@self.ChangeProtocols); %#ok<NASGU>
            cancel_button = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.85,delta_y/2, protocol_figure.Position(3)*0.09, protocol_figure.Position(4)*0.12],'String', 'Cancel','Tag','cancel_button','Callback',@self.DoNothing); %#ok<NASGU>
            
        end
        
        function AddProtocol(self,obj,event) %#ok<INUSD>
            
            add_obj = findobj('Tag','Add protocol dropmenu');
            insert_obj = findobj('Tag','Insert protocol dropmenu');
            protocol_to_add = add_obj.String(add_obj.Value);
            
            %add the protocol
            for p = 1:length(self.protocol_types)
                if strcmp(protocol_to_add, self.protocol_types{p}.sProtocolName)
                    new_protocol = RealtimeProtocol(1,self.protocol_types{p},self.sampling_frequency);
                    break;
                end
            end
            
            
            %update the figure
            protocols_names_obj = findobj('Tag', 'Protocol name text');
            protocols_names = {};
            for pn = length(protocols_names_obj):-1:1
                protocols_names = [protocols_names {protocols_names_obj(pn).String}];
            end
            
            protocols_durations_obj = findobj('Tag', 'Protocol duration text');
            protocols_durations = {};
            for pd = length(protocols_durations_obj):-1:1
                protocols_durations = [protocols_durations {protocols_durations_obj(pd).String}];
            end
            %!!!!
            pr_idx = strsplit(insert_obj.String{insert_obj.Value});
            idx = str2num(pr_idx{1}); %#ok<ST2NM>
            protocols_names = [protocols_names(end:-1:idx+1) new_protocol.protocol_name protocols_names(idx:-1:1)];
            protocols_names = protocols_names(end:-1:1);
            protocols_durations = [protocols_durations(end:-1:idx+1) num2str(new_protocol.protocol_duration) protocols_durations(idx:-1:1)];
            protocols_durations = protocols_durations(end:-1:1);
            self.UpdateEditProtocolsFigure(protocols_names,protocols_durations);
        end
        function DeleteProtocol(self,obj,event) %#ok<INUSD>
            delete_obj = findobj('Tag','Delete protocol dropmenu');
            
            protocols_names_obj = findobj('Tag', 'Protocol name text');
            protocols_names = {};
            for pn = length(protocols_names_obj):-1:1
                protocols_names = [protocols_names {protocols_names_obj(pn).String}];
            end
            
            protocols_durations_obj = findobj('Tag', 'Protocol duration text');
            protocols_durations = {};
            for pd = length(protocols_durations_obj):-1:1
                protocols_durations = [protocols_durations {protocols_durations_obj(pd).String}];
            end
            
            protocols_names = [protocols_names(1:delete_obj.Value+self.next_protocol-2) protocols_names(delete_obj.Value+self.next_protocol:end)];
            protocols_durations = [protocols_durations(1:delete_obj.Value+self.next_protocol-2) protocols_durations(delete_obj.Value+self.next_protocol:end)];
            self.UpdateEditProtocolsFigure(protocols_names,protocols_durations);
        end
        function UpdateEditProtocolsFigure(self,protocols_names,protocols_durations)
            protocol_figure = findobj('Tag', 'EditProtocolFigure');
            old_names = findobj('Tag', 'Protocol name text');
            old_durations = findobj('Tag', 'Protocol duration text');
            old_count = findobj('Tag', 'Protocol count');
            old_ms = findobj('Tag','s text');
            
            delete (old_names);
            delete (old_durations);
            delete (old_count);
            delete(old_ms);
            
            max_height = protocol_figure.Position(4) - 0.05;
            %mind the numbers
            
            self.protocol_duration_text = {};
            for p = 1:length(protocols_names)
                bgr = 0.94-[0.1 0.1 0.1] * mod(p-1,2);
                protocol_count{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.01,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(p),'HorizontalAlignment','left','Tag','Protocol count','BackgroundColor',bgr); %#ok<NASGU>
                protocol_name{end+1} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.04,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.25, protocol_figure.Position(4)*0.04],'String', protocols_names{p},'HorizontalAlignment','left','Tag','Protocol name text','BackgroundColor',bgr); %#ok<NASGU>
                if p < self.next_protocol %already recorded; duration cannot be changed
                    self.protocol_duration_text{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.3, max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(protocols_durations{p}),'Tag','Protocol duration text','HorizontalAlignment','right');
                else
                    self.protocol_duration_text{p} = uicontrol('Parent',protocol_figure,'Style','edit','Position', [protocol_figure.Position(3)*0.3, max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(protocols_durations{p}),'HorizontalAlignment','right','Tag', 'Protocol duration text');
                end
                ms{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.35,max_height-protocol_figure.Position(4)*0.05*p,protocol_figure.Position(3)*0.02 , protocol_figure.Position(4)*0.04],'String', 's','HorizontalAlignment','left','BackgroundColor',bgr,'Tag','s text'); %#ok<NASGU>
                
            end
            %update dropmenus
            insert_protocol_dropmenu = findobj('Tag','Insert protocol dropmenu');
            delete_protocol_dropmenu = findobj('Tag','Delete protocol dropmenu');
            new_protocols_names = cell(1,length(protocols_names));
            for i = self.next_protocol-1:length(protocols_names)
                if i < 1
                    continue;
                end
                pn = protocols_names(i);
                new_protocols_names(i) = {[num2str(i) ' ' pn{1}]};
            end
            if self.next_protocol == 1
                insert_protocol_dropmenu.String = [{'0'} new_protocols_names];
            else
                insert_protocol_dropmenu.String = new_protocols_names; %%%%delete prt dropmenu
            end
            delete_protocol_dropmenu.String = new_protocols_names;
            
            
            
        end
        function ChangeProtocols(self,obj,event) %#ok<INUSD>
            %get the protocols
            self.feedback_protocols(self.next_protocol:end) = [];
            protocols_names_obj = findobj('Tag', 'Protocol name text');
            protocols_durations_obj = findobj('Tag', 'Protocol duration text');
            for j = length(protocols_names_obj)-self.next_protocol+1:-1:1
                for i = 1:length(self.protocol_types)
                    if strcmp(protocols_names_obj(j).String,self.protocol_types{i}.sProtocolName)
                        rtp = RealtimeProtocol(1,self.protocol_types{i});
                        rtp.protocol_duration = str2double(protocols_durations_obj(j).String);  %%%%%duration is taken from the figure
                        rtp.Recalculate(self.sampling_frequency);
                        self.feedback_protocols{end+1} = rtp;
                        break;
                    end
                end
            end
            %check if we reserved enough space
            data_length = 0;
            for p = 1:length(self.feedback_protocols)
                data_length = data_length + self.feedback_protocols{p}.protocol_size;
            end
            if fix(data_length*1.1)> self.exp_data_length
                
                self.exp_data_length = fix(data_length*1.1);
                for ds = 1:length(self.derived_signals)
                    new_circbuff = circVBuf(self.exp_data_length, size(self.derived_signals{ds}.collect_buff.raw,2),0);
                    new_circbuff.append(self.derived_signals{ds}.collect_buff.raw(self.derived_signals{ds}.collect_buff.fst:self.derived_signals{ds}.collect_buff.lst,:));
                    self.derived_signals{ds}.collect_buff = new_circbuff;
                end
                new_fb_records = circVBuf(self.exp_data_length,6,0);
                if ~isempty(self.feedback_manager.feedback_records)
                    new_fb_records.append(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst:self.feedback_manager.feedback_records.lst,:));
                    self.feedback_manager.feedback_records = new_fb_records;
                end
            end
            
            %delete the figure
            delete(obj.Parent);
        end
        function ExpDesignToXML(self)
            curr_date = datestr(date,29);
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date));
            end
            if ~isdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start))
                mkdir(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start));
            end
            
            %prepare xml
            self.exp_design = struct();
            vSignals = struct();
            for ds = 1:length(self.derived_signals)
                if ds > 1
                    vSignals.DerivedSignal(ds).fAverage = self.feedback_manager.average(ds-1);
                    vSignals.DerivedSignal(ds).fStdDev = self.feedback_manager.standard_deviation(ds-1);
                end
                vSignals.DerivedSignal(ds).sSignalName = self.derived_signals{ds}.signal_name;
                
                try %#ok<TRYNC>
                    vSignals.DerivedSignal(ds).fBandpassLowHz = self.derived_signals{ds}.temporal_filter{1}.range(1);
                    vSignals.DerivedSignal(ds).fBandpassHighHz = self.derived_signals{ds}.temporal_filter{1}.range(2);
                end
                
                try %#ok<TRYNC>
                    vSignals.DerivedSignal(ds).sType = self.derived_signals{ds}.signal_type;
                end
                %prepare spatial filter matrix
                spatial_filter_matrix_struct = struct();
                for ch = 1:length(self.derived_signals{ds}.channels)
                    s = '';
                    for i = 2:size(self.derived_signals{ds}.channels,2)
                        s = [s, ',', num2str(self.derived_signals{ds}.channels{ch,i})];
                    end
                    s = s(2:end);
                    spatial_filter_matrix_struct.channels.(self.derived_signals{ds}.channels{ch,1}) =  s;
                end
                sp_filter_matrix = struct2xml(spatial_filter_matrix_struct);
                spf_filename = strcat(vSignals.DerivedSignal(ds).sSignalName,'.xml');
                vSignals.DerivedSignal(ds).SpatialFilterMatrix = spf_filename;
                %write the file
                spffile = fopen(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start,'\',spf_filename),'w');
                fwrite(spffile,sp_filter_matrix);
                fclose(spffile);
                %vSignals.DerivedSignal(ds) = DS;
            end
            
            vProtocols = struct();
            for p = 1:length(self.protocol_types)
                %copy all the fields
                fields = fieldnames(self.protocol_types{p});
                for j = 1: numel(fields)
                    %try %#ok<TRYNC>
                    vProtocols.FeedbackProtocol(p).(fields{j}) = self.protocol_types{p}.(fields{j});
                    %end
                end
            end
            vPSequence = struct();
            for fp = 1:length(self.feedback_protocols)
                vPSequence.s{fp} = self.feedback_protocols{fp}.protocol_name;
            end
            self.exp_design.NeurofeedbackSignalSpecs.vSignals = vSignals;
            self.exp_design.NeurofeedbackSignalSpecs.vProtocols = vProtocols;
            self.exp_design.NeurofeedbackSignalSpecs.vPSequence = vPSequence;
            % necessary for formatting xml structures
            a = struct2xml(self.exp_design);
            
            design_file = fopen(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start,'\Exp_design.xml'),'w');
            fwrite(design_file,a);
            fclose(design_file);
            disp(['Design file successfully written to the location ' strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start,'\Exp_design.xml')]);
        end
        function DoNothing(self,obj,event) %#ok<INUSL,INUSD>
            %do nothing and destroy the window
            if strcmp(get(obj,'Tag'),'cancel_button')
                delete(obj.Parent);
            else
                delete(obj);
            end
        end
        function SelectBadChannels(self,obj,event) %#ok<INUSD>
            global finished;
            global ok;
            f = findobj('Tag','raw_and_ds_figure');
            epb = findobj('Tag','edit_protocols_button');
            set(epb,'Visible','off');
            bad_channels_text =  uicontrol('Parent',f,'Style','text','String', self.bad_channels,'HorizontalAlignment','right','Tag','bad_channels_text');
            ok_button = uicontrol('Parent',f,'style','pushbutton', 'String', 'Add','Tag','add_bad_channel_button','Callback','global ok; ok = 1;');
            finished_button= uicontrol('Parent',f,'style','pushbutton', 'String', 'Finish','Tag','finish adding bad channels button','Callback','global ok; global finished; ok = 1;finished = 1;');
            text =  uicontrol('Parent',f,'Style','text','String', 'Select a channel and press Add','HorizontalAlignment','right','Tag','ch_text','HorizontalAlignment','center');
            datacursormode('off');
            r_sp = get(self.raw_subplot);
            self.FitFigure;
            finished = 0;
            while ~finished
                ok = 0;
                dcm_obj = datacursormode(f);
                set(dcm_obj,'UpdateFcn',@cursor_callback);
                %if 'OK'
                while ~(~isempty(getCursorInfo(dcm_obj)) && ok)
                    datacursormode('on');
                    pause(1);
                    if ok && isempty(getCursorInfo(dcm_obj))
                        ok = 0;
                    end
                    if finished
                        break; %#ok<UNRCH>
                    end
                    
                end
                if ~finished
                    datacursormode('off');
                    cursor_info = getCursorInfo(dcm_obj);
                    self.bad_channels{end+1} = cursor_info.Target.DisplayName;
                    self.bad_channels = unique(self.bad_channels);
                    set(bad_channels_text,'String',self.bad_channels);
                    dcm_obj.removeAllDataCursors();
                end
            end
            delete(ok_button);
            delete(finished_button);
            delete(findobj('Tag','bad_channels_text'));
            delete(text);
            datacursormode('off');
            set(epb,'Visible','on');
            %update derived signals
            
            for b_ch = 1:length(self.bad_channels)
                for child = 1:length(r_sp.Children)
                    try
                        if strcmp(self.bad_channels{b_ch},r_sp.Children(length(r_sp.Children)-child+1).DisplayName)
                            self.raw_data_indices = self.raw_data_indices(self.raw_data_indices ~= child);
                            
                        end
                    catch
                        length(r_sp.Children)
                        child %#ok<NOPRT>
                    end
                end
            end
            for ds = 1:length(self.derived_signals)
                self.derived_signals{ds}.ZeroOutBadChannels(self.bad_channels);
            end
            %update  YTickLabels
            for b_ch = 1:length(self.bad_channels)
                self.r_ytick_labels(strcmp(self.bad_channels{b_ch},self.r_ytick_labels)) = [];
            end
            %update shift
            r_temp = zeros(length(self.raw_data_indices),fix(self.plot_size));
            self.raw_plot = plot(r_temp', 'Parent', self.raw_subplot);
            for i = 2:length(self.r_ytick_labels)-1
                set(self.raw_plot(i-1),'DisplayName',self.r_ytick_labels{i});
            end
            %and set them
            set(self.raw_subplot,'YLim',[0 self.raw_shift*length(self.r_ytick_labels)]);
            self.FitFigure;
            
        end
        function [baseline_data, data_sets, data_names] = Prepare_CSP(self)
            %stop the timers
            self.inlet.close_stream();
            stop(self.timer_new_data);
            stop(self.timer_disp);
            stop(self.timer_fb);
            %%list the previous protocols and allow to checkbox those
            %%that are supposed to be calculated
            csp_figure = figure('Tag','CSP choice');
            max_height = csp_figure.Position(4)-csp_figure.Position(4)*0.1;
            select_baseline_protocol = uicontrol('Parent', csp_figure,'Style','text','Position',[csp_figure.Position(3)*0.12,csp_figure.Position(4)*0.85, csp_figure.Position(3)*0.2, csp_figure.Position(4)*0.1],...
                'Tag','select_baseline_protocol_text','String','Select Baseline','FontSize',11,'HorizontalAlignment','center'); %#ok<NASGU>
            select_other_protocols = uicontrol('Parent', csp_figure,'Style','text','Position',[csp_figure.Position(3)*0.62,csp_figure.Position(4)*0.85, csp_figure.Position(3)*0.2, csp_figure.Position(4)*0.1],...
                'Tag','select_baseline_protocol_text','String','Select protocols to CSP','FontSize',11,'HorizontalAlignment','center'); %#ok<NASGU>
            
            protocol_rbgroup = uibuttongroup('Parent',csp_figure,'Position',[0.08 0.1 .3 0.76]);
            for pr = 1:self.next_protocol-1
                bgr = ones(1,3) * 0.94; %default background colo
                %bgr = 0.94-[0.1 0.1 0.1] * mod(pr-1,2); %background color
                protocol_rb{pr} = uicontrol(protocol_rbgroup,'Style','radiobutton','String',self.feedback_protocols{pr}.show_as,'Position',[csp_figure.Position(3)*0.05,max_height-csp_figure.Position(4)*0.15-csp_figure.Position(4)*0.09*pr, csp_figure.Position(3)*0.2, csp_figure.Position(4)*0.05],'HandleVisibility','off');
                protocol_chb{pr} = uicontrol('Parent',csp_figure,'Style','checkbox','Position',[csp_figure.Position(3)*0.45,max_height-csp_figure.Position(4)*0.04-csp_figure.Position(4)*0.09*pr, csp_figure.Position(3)*0.3, csp_figure.Position(4)*0.05],'Tag','protocols_chb','BackgroundColor',bgr,'Callback',@self.CheckIfSelected,'String',self.feedback_protocols{pr}.show_as);
                %edit_name{pr} = uicontrol('Parent',csp_figure,'Style','edit','Position', [csp_figure.Position(3)*0.65,max_height-csp_figure.Position(4)*0.04-csp_figure.Position(4)*0.09*pr, csp_figure.Position(3)*0.3, csp_figure.Position(4)*0.05],'String', self.feedback_protocols{pr}.show_as,'HorizontalAlignment','left','Tag','Edit name text');
            end
            okay_button = uicontrol('Parent',csp_figure,'Style','pushbutton','Position', [csp_figure.Position(3)*0.75,csp_figure.Position(4)*0.05, csp_figure.Position(3)*0.09,csp_figure.Position(4)*0.12],'String', 'OK','Tag','okay_button','Callback','uiresume','enable','off'); %#ok<NASGU>
            uiwait();
            
            data_sets = {};
            data_names = {};
            
            %%loop through protocols
            for pr = 1:self.next_protocol-1
                %fetch baseline data
                if get(protocol_rb{pr},'Value')
                    idx1 = self.protocol_indices(pr);
                    idx2 = self.protocol_indices(pr+1);
                    baseline_data = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx1:self.derived_signals{1}.collect_buff.fst+idx2-1,:);
                end
                %fetch csp data
                if get(protocol_chb{pr},'Value')
                    idx1 = self.protocol_indices(pr);
                    idx2 = self.protocol_indices(pr+1);
                    %check data length
                    data_sets{end+1} = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx1:self.derived_signals{1}.collect_buff.fst+idx2-1,:);
                    data_names{end+1} = get(protocol_chb{pr},'String');
                    
                end
                
            end
            delete(csp_figure);

            for d = 1:length(data_sets)
                self.CalculateCSP(baseline_data,data_sets{d}, data_names{d});
            end
        end
        function CalculateCSP(self,baseline_data,csp_data,data_name)
            global selected
            
            if ~isfield(self.csp_settings,'iNComp') || isempty(self.csp_settings.iNComp)
                Ncomp = 2;
            else
                Ncomp = self.csp_settings.iNComp;
            end
            
            if ~isfield(self.csp_settings,'dInitBand') || isempty(self.csp_settings.dInitBand)
                init_band = self.init_band;
            else
                init_band = self.csp_settings.dInitBand;
            end
            for ib = 1:4
                band = init_band +ib-1 ;
                [z, p, k] = cheby1(3,1,band/(0.5*self.sampling_frequency),'bandpass');
                [b,a] = zp2tf(z,p,k); 
                filtered_bd = filtfilt(b,a,baseline_data)';
                filtered_csp_d = filtfilt(b,a,csp_data)';
                C10 = filtered_bd * filtered_bd'/fix(size(filtered_bd,2));
                C20 = filtered_csp_d * filtered_csp_d'/fix(size(filtered_csp_d,2));
                
                nchan = size(C10,1);
                
                %%regularize covariances
                Lambda = 0.1;%%%%%%%%%%%%%%%%%%%%%%
                C1 = C10 + Lambda * trace(C10) * eye(nchan) / nchan;
                C2 = C20 + Lambda * trace(C20) * eye(nchan) / nchan;
                %%do generalized eigenvalue decomp
                [V, d] = eig(C1,C2); %#ok<ASGLU>
                iV = inv(V);
                M12{ib} = V(:,[1:Ncomp, end-Ncomp+1:end])';
                G12{ib} = iV([1:Ncomp, end-Ncomp+1:end],:);
                
            end
            %%present the heads
            hh1 = figure;
            StandChannels = load('StandardEEGChannelLocations.mat');
            chan_labels = self.used_ch(:,1)';
            Nbands = length(G12);
            PlotIndex = 1;
            selected = {};
            
            for ib = 1:Nbands
                rearranged_map = rearrange_channels(G12{ib}',chan_labels, StandChannels.channels);
                Nmaps = size(rearranged_map,2);
                for tpm=1:Nmaps
                    sp(PlotIndex) = subplot(Nbands,Nmaps,PlotIndex);
                    
                    %                            topoplot(rearranged_map(:,tpm), StandChannels.channels, 'electrodes', 'labelpoint', 'chaninfo', StandChannels.chaninfo);
                    topoplot(rearranged_map(:,tpm), StandChannels.channels,  'chaninfo', StandChannels.chaninfo);
                    hold on;
                    sibs = get(sp(PlotIndex),'Children');
                    for k = 1:length(sibs)
                        set(sp(PlotIndex).Children(k), 'ButtonDownFcn', @(src,event)toggleplot(src,event));
                    end
                    title(num2str(PlotIndex));
                    %add legend
                    PlotIndex = PlotIndex+1;
                end;
            end
            %setting  figure title
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off'); %#ok<NASGU>
            text(0.5, 1, data_name,'HorizontalAlignment','center','VerticalAlignment', 'top')
            
            okay_btn = uicontrol('Parent',hh1,'Style', 'pushbutton', 'String', 'OK', 'Callback', 'uiresume', 'Position', [230 10 100 35],'Tag','SelectHeadsBtn','enable','off');  %#ok<NASGU>
            uiwait();
            
            MapIndex = cellfun(@str2num,selected);
            BandNumber = fix((MapIndex-1)/(2*Ncomp))+1;
            CompNumber = mod(MapIndex-1,(2*Ncomp))+1;
            close(hh1);
            w_ssd = zeros(length(MapIndex),length(M12{1}));
            for bn = 1:length(MapIndex)
                w_ssd(bn,:) = M12{BandNumber(bn)}(CompNumber(bn),:);
            end
            %%if the band is the same for all inputs, a choice is
            %offered whether to make one combined signal or
            %several; If the bands are different, no such
            %choice is offered
            %determine whether the band is the same
            signals_added = 0;
            if isempty(nonzeros(BandNumber - max(BandNumber))) && length(selected) > 1 %if it is empty, then all the bands are the same
                button = questdlg(['One combined DS or ', num2str(length(MapIndex)), ' different DS?'],'Choose the number of DS', 'One', num2str(length(MapIndex)),'One');
                if strcmp(button, 'One')
                    %make one DS with 2-D sp_filter
                    band = init_band + BandNumber(1)-1;
                    chan_w = cell(length(self.derived_signals{1}.channels),size(w_ssd,1)+1);
                    w_ssd = w_ssd';
                    % 2-d spatial_filter
                    for i=1:length(w_ssd)
                        %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                        chan_w{i,1} = self.derived_signals{1}.channels{i};
                        for ind = 2:size(w_ssd,2)+1
                            chan_w{i,ind} = w_ssd(i,ind-1);
                        end
                        
                    end;
                    clear chs
                    chs(size(chan_w,1)) = struct();
                    for ch = 1:size(chan_w,1)
                        chs(ch).channel_name = chan_w{ch,1};
                        coeff = '';
                        for ind = 2:size(w_ssd,2)+1
                            coeff = [coeff ',' num2str(chan_w{ch,ind})];
                        end
                        
                        chs(ch).coefficient = coeff(2:end);
                    end
                    
                    
                    
                    sn = {''};
                    while isempty(sn{1})
                        sn = inputdlg('Enter derived signal name','Derived signal name',1,{strcat(data_name,'_CombinedDS_',num2str(band(1)),'-',num2str(band(2)))});
                    end
                    
                    %%write spatial_filter to file
                    if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                        mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                    end
                    pth = (strcat(self.path,'\',self.subject_record.subject_name));
                    full_name = [pth '\CSP_' sn{1} '.xml'];
                    
                    if exist(full_name,'file')
                        choice = questdlg('The spatial filter for this person already exists. Rewrite the filter?','Rewrite?','Yes','No, use the old one','No, use the old one');
                        switch choice
                            case 'Yes'
                                chan_structure = struct();
                                chan_structure.channels.channel = chs;
                                x = struct2xml(chan_structure);
                                f = fopen(full_name,'w');
                                fwrite(f,x);
                                fclose(f);
                                spatial_filter = chan_w;
                            case 'No, use the old one'
                                s = xml2struct(full_name);
                                channels_coeff = cell(length(s.channels.channel),2);
                                for i = 1:length(s.channels.channel)
                                    channels_coeff{i,1} = strtrim(s.channels.channel{i}.channel_name.Text);
                                    coeffs = str2num(strtrim(s.channels.channel{i}.coefficient.Text)); %#ok<ST2NM>
                                    for j = 2:length(coeffs)+1
                                        channels_coeff{i,j} = coeffs(j-1);
                                    end
                                end
                                %convert struct to cell!!
                                spatial_filter = channels_coeff;
                        end
                    else
                        spatial_filter = chs;
                        sp_filter_to_write = struct();
                        sp_filter_to_write.channels = chs;
                        x = struct2xml(sp_filter_to_write);
                        f = fopen(full_name,'w');
                        fwrite(f,x);
                        fclose(f);
                    end
                    
                    %%create DS
                    %
                    self.derived_signals{end+1} = CreateNewDS(self,sn{1},chan_w, band,full_name);
                    self.derived_signals{end}.signal_type = 'combined';
                    signals_added = 1;
                    %%set fb
                    
                    %
                    temp_derived_signal = CreateNewDS(self,'Temp',chan_w, band);
                    %
                    temp_derived_signal.signal_type = 'combined';
                    %

                    temp_derived_signal.Apply(baseline_data',1);
                    values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                    %calcuate feedback stats
                    values = abs(values);
                    self.feedback_manager.average(length(self.derived_signals)-1) = mean(values);
                    self.feedback_manager.standard_deviation(length(self.derived_signals)-1) = std(values);
                    
                    
                    self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                    self.fb_manager_set = 1;
                    
                end
            end
            
            %if no signals were added so far
            if ~signals_added
            for ind = 1:length(MapIndex)
                band = init_band + BandNumber(ind)-1;
                for i=1:length(w_ssd)
                    %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                    chan_w{i,1} = self.derived_signals{1}.channels{i};
                    chan_w{i,2} = w_ssd(ind,i);
                end;
                
                %write to file
                try %write
                    %convert cell chan_w to structure
                    clear chs
                    chs(size(chan_w,1)) = struct();
                    for ch = 1:size(chan_w,1)
                        chs(ch).channel_name = chan_w{ch,1};
                        chs(ch).coefficient = chan_w{ch,2};
                    end
                    chan_structure = struct();
                    chan_structure.channels.channel = chs;
                    x = struct2xml(chan_structure);
                    
                    
                    sn = {''};
                    if CompNumber(ind) <= Nmaps/2
                        n = CompNumber(ind);
                    else
                        n = CompNumber(ind)-Nmaps-1;
                    end
                    while isempty(sn{1})
                        sn = inputdlg('Enter derived signal name','Derived signal name',1,{strcat(data_name,'_DS_', num2str(band(1)),'-',num2str(band(2)),'_',num2str(n))});
                    end
                    
                    if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                        mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                    end
                    pth = (strcat(self.path,'\',self.subject_record.subject_name));
                    full_name = [pth '\CSP_' sn{1} '.xml'];
                    if exist(full_name,'file')
                        choice = questdlg('The spatial filter for this person already exists. Rewrite the filter?','Rewrite?','Yes','No, use the old one','No, use the old one');
                        switch choice
                            case 'Yes'
                                f = fopen(full_name,'w');
                                fwrite(f,x);
                                fclose(f);
                                spatial_filter = chan_w;
                            case 'No, use the old one'
                                
                                s = xml2struct(full_name);
                                channels_coeff = cell(length(s.channels.channel),2);
                                for i = 1:length(s.channels.channel)
                                    channels_coeff{i,1} = strtrim(s.channels.channel{i}.channel_name.Text);
                                    channels_coeff{i,2} = str2double(strtrim(s.channels.channel{i}.coefficient.Text));
                                end
                                spatial_filter = channels_coeff;
                        end
                    else
                        spatial_filter = chan_w;
                        f = fopen(full_name,'w');
                        fwrite(f,x);
                        fclose(f);
                    end
                    
                    
                catch
                    'Error while writing to file, function UpdateStatistics' %#ok<NOPRT>
                end
                
                
                try
                    
                    %
                    
                    self.derived_signals{end+1} = CreateNewDS(self,sn{1},spatial_filter, band,full_name);
                    
                    %%temp derived signal to calculate stats
                    temp_derived_signal = CreateNewDS(self,sn{1},spatial_filter, band);
                    
                    %%calculating stats
                    %j = 1;
                    x = zeros(size(baseline_data,2),length(self.derived_signals{1}.spatial_filter));
                    for i = 1:length(self.derived_signals{1}.spatial_filter)
                        if self.derived_signals{1}.spatial_filter(i)
                            x(:,i) = baseline_data(i,:);
                            %j = j+1;
                        else
                            x(:,i) = 0;
                        end
                    end
                    temp_derived_signal.Apply(x',1);
                    values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                    self.feedback_manager.average(length(self.derived_signals)-1) = mean(values);
                    self.feedback_manager.standard_deviation(length(self.derived_signals)-1) = std(values);
                    
                    
                catch
                    'Error while creating a new derived signal, function CalculateCSP' %#ok<NOPRT>
                end
                
            end
            
            
            
            self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
            self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
            self.fb_manager_set = 1;
            end
            
        end
        function CheckIfSelected(self,src,event)  %#ok<INUSD>
            
            okay_button = findobj('Tag','okay_button');
            
            protocol_chbs = findobj('Tag','protocols_chb');
            selected = 0;
            for pr = 1:length(protocol_chbs)
                if get(protocol_chbs(pr),'Value') == 1
                    selected = selected + 1;
                end
            end
            
            if selected < 1
                set(okay_button,'enable','off');
            else
                set(okay_button,'enable','on');
            end
            
        end
        function CSPLearning(self,data_sets,data_names)
            %%data_sets_choices(n,2)
            choices = nchoosek(1:length(data_sets),2);
            
            for ch = 1:size(choices,1)
                if size(choices,1) == 1
                    choice = choices;
                else
                    choice = choices(ch,:);
                end
                
                %%run pairwise csp learning
                first_raw = data_sets{choice(1)};
                second_raw = data_sets{choice(2)};
                first_name = data_names(choice(1));
                second_name = data_names(choice(2));
                try
                    
                    %fetch the data of the choices
                    
                    if ~isempty(self.csp_settings.n_comp)
                        Ncomp = self.csp_settings.n_comp;
                    else
                        Ncomp = 2;
                    end
                catch
                    Ncomp = 2;
                end
                try
                    if ~isempty(self.csp_settings.init_band)
                        init_band = self.csp_settings.init_band;
                    else
                        init_band = self.init_band;
                    end
                catch
                    init_band = self.init_band;
                end
                
                for ib = 1:4
                    band = init_band +ib-1 ;
                    [z, p, k] = cheby1(3,1,band/(0.5*self.sampling_frequency),'bandpass');
                    [b,a] = zp2tf(z,p,k);
                    %                     x = filtfilt(b,a,x_raw)';
                    %                     C10 = x(:,1:fix(end/2))* x(:,1:fix(end/2))'/fix(size(x,2)/2);
                    %                     C20 = x(:,fix(end/2)+1:end)* x(:,fix(end/2)+1:end)'/fix(size(x,2)/2);
                    x1 = filtfilt(b,a, first_raw)';
                    x2 = filtfilt(b,a, second_raw)';
                    C10 = x1 * x1' / size(x1,2);
                    C20 = x2 * x2' / size(x2,2);
                    
                    nchan = size(C10,1);
                    
                    %%regularize covariances
                    Lambda = 0.1;%%%%%%%%%%%%%%%%%%%%%%
                    C1 = C10 + Lambda * trace(C10) * eye(nchan) / nchan;
                    C2 = C20 + Lambda * trace(C20) * eye(nchan) / nchan;
                    %%do generalized eigenvalue decomp
                    [V, d] = eig(C1,C2); %#ok<ASGLU>
                    iV = inv(V);
                    M12{ib} = V(:,[1:Ncomp, end-Ncomp+1:end])';
                    G12{ib} = iV([1:Ncomp, end-Ncomp+1:end],:);
                end
                %%show the heads
                name = [first_name '-' second_name];
                hh1 = figure;
                title(name);
                
                StandChannels = load('StandardEEGChannelLocations.mat');
                chan_labels = self.used_ch(:,1)';
                Nbands = length(G12);
                PlotIndex = 1;
                selected = {};
                for ib = 1:Nbands
                    rearranged_map = rearrange_channels(G12{ib}',chan_labels, StandChannels.channels);
                    Nmaps = size(rearranged_map,2);
                    for tpm=1:Nmaps
                        sp(PlotIndex) = subplot(Nbands,Nmaps,PlotIndex);
                        
                        %                            topoplot(rearranged_map(:,tpm), StandChannels.channels, 'electrodes', 'labelpoint', 'chaninfo', StandChannels.chaninfo);
                        topoplot(rearranged_map(:,tpm), StandChannels.channels,  'chaninfo', StandChannels.chaninfo);
                        hold on;
                        sibs = get(sp(PlotIndex),'Children');
                        for k = 1:length(sibs)
                            set(sp(PlotIndex).Children(k), 'ButtonDownFcn', @(src,event)toggleplot(src,event));
                        end
                        title(num2str(PlotIndex));
                        %add legend
                        PlotIndex = PlotIndex+1;
                    end;
                end
                okay_btn = uicontrol('Parent',hh1,'Style', 'pushbutton', 'String', 'OK', 'Callback', 'uiresume', 'Position', [230 10 100 35],'Tag','SelectHeadsBtn','enable','off');  %#ok<NASGU>
                
                
                uiwait();
                
                MapIndex = cellfun(@str2num,selected);
                BandNumber = fix((MapIndex-1)/(2*Ncomp))+1;
                CompNumber = mod(MapIndex-1,(2*Ncomp))+1;
                close(hh1);
                w_ssd = zeros(length(MapIndex),length(M12{1}));
                for bn = 1:length(MapIndex)
                    w_ssd(bn,:) = M12{BandNumber(bn)}(CompNumber(bn),:);
                end
                
                
                
                
                %%save as DS
            end
            
        end
        function NewDS = CreateNewDS(self,signal_name,spatial_filter, band,filename)
            %%init dummy signal
            dummy_signal = struct();
            if nargin > 1
                dummy_signal.sSignalName = signal_name;
            else
                dummy_signal.sSignalName = '';
            end
            %%if raw signal exists, use its channels to copy to dummy ds
            if ~isempty(self.derived_signals)
                channels = cell(length(self.derived_signals{1}.channels),2);
                for ch = 1:length(channels)
                    channels{ch,1} = self.derived_signals{1}.channels{ch};
                    channels{ch,2} = 1;
                end
                %else use self.channel_labels and the weigth of one
            else
                channels = cell(size(self.channel_labels,2));
                for ch = 1:length(channels)
                    channels{ch,1} = self.channel_labels{ch};
                    channels{ch,2} = 1;
                end
            end
            dummy_signal.filters = cell(0,0);
            dummy_signal.channels = channels;
            
            %%DS to use
            NewDS= DerivedSignal(1,dummy_signal, self.sampling_frequency, self.exp_data_length ,self.channel_labels,self.plot_length);
            if strcmpi(dummy_signal.sSignalName,'raw')
                NewDS.ring_buff = circVBuf(self.plot_size,length(channels),0);
                NewDS.collect_buff = circVBuf(self.exp_data_length,length(channels),0);
            else
                NewDS.ring_buff = circVBuf(self.plot_size,1,0);
                NewDS.collect_buff = circVBuf(self.exp_data_length,1,0);
            end
            if nargin > 4
                NewDS.channels_file = filename;
            end
            if nargin > 3
                NewDS.UpdateTemporalFilter(size(spatial_filter),band);
            end
            if nargin > 2 && ~isempty(self.derived_signals)
                NewDS.UpdateSpatialFilter(spatial_filter,self.derived_signals{1},self.bad_channels);
            elseif nargin > 2
                NewDS.UpdateSpatialFilter(spatial_filter);
            end
            
            
        end
        function NewRP = CreateNewRP(self,protocol_name,actual_size) %#ok<INUSL>
            NewRP = RealtimeProtocol;
            NewRP.protocol_name = protocol_name;
            NewRP.to_update_statistics = true;
            NewRP.actual_protocol_size = actual_size;
            
        end
        function [results,data_names] = Recalculate(self,data_sets,window,av,stddev)
            % >> self = EEGLSL;
            % >> Recalculate(self, 1:10, 10)
            %
            %ans =
            %
            %    5.5000    5.5000    5.5000    5.5000    5.5000    5.5000    5.5000    5.5000    5.5000    5.5000
            %
            % >> Recalculate(self,1:10,2)
            %
            %ans =
            %
            %     1.5000    1.5000    3.5000    3.5000    5.5000    5.5000    7.5000    7.5000    9.5000    9.5000
            %
            % >> sum(Recalculate(self,1:10,10)) == sum(Recalculate(self,1:10,2))
            %
            %ans =
            %
            % 1
            %
            % >> sum(Recalculate(self,1:100,10))
            %
            %ans =
            %
            % 5050
            %
            if nargin < 2 && length(self.derived_signals) > 1
                %open the window
                %select protocol(s)
                f = figure('Tag', 'recalculate_protocol_choice');
                protocols_to_use = uicontrol('Parent', f,'Style', 'text', 'Position',[f.Position(3)*0.12,f.Position(4)*0.85, f.Position(3)*0.2, f.Position(4)*0.1],'Tag','protocols_to_select','String','Data to calculate','FontSize',11,'HorizontalAlignment','center'); %#ok<NASGU>
                max_height = f.Position(4) - f.Position(4)*0.1;
                if self.next_protocol > length(self.feedback_protocols)
                    finish = length(self.feedback_protocols);
                else
                    finish = self.next_protocol-1;
                end
                for pr = 1:finish
                    
                    bgr = 0.94-[0.1 0.1 0.1] * mod(pr-1,2);
                    protocol_chb{pr} = uicontrol('Parent',f,'Style','checkbox','Position',[f.Position(3)*0.1,max_height-f.Position(4)*0.05*pr, f.Position(3)*0.3, f.Position(4)*0.03],'Tag','protocols_chb','BackgroundColor',bgr,'Callback',@self.CheckIfSelected,'String',self.feedback_protocols{pr}.protocol_name);
                    
                    %protocol_count{p} = uicontrol('Parent',csp_figure,'Style','text','Position', [csp_figure.Position(3)*0.03,max_height-csp_figure.Position(4)*0.05*pr, csp_figure.Position(3)*0.05, csp_figure.Position(4)*0.04],'String', num2str(pr),'HorizontalAlignment','left','Tag','Protocol count','BackgroundColor',bgr); %#ok<NASGU>
                    %protocol_name{pr} = uicontrol('Parent',csp_figure,'Style','text','Position', [csp_figure.Position(3)*0.07,max_height-csp_figure.Position(4)*0.05*pr, csp_figure.Position(3)*0.25, csp_figure.Position(4)*0.04],'String', self.feedback_protocols{pr}.protocol_name,'HorizontalAlignment','left','Tag','Protocol name text','BackgroundColor',bgr); %#ok<NASGU>
                    %edit_name{pr} = uicontrol('Parent',csp_figure,'Style','edit','Position', [csp_figure.Position(3)*0.45,max_height-csp_figure.Position(4)*0.09*pr, csp_figure.Position(3)*0.3, csp_figure.Position(4)*0.05],'String', self.feedback_protocols{pr}.protocol_name,'HorizontalAlignment','left','Tag','Edit name text','BackgroundColor',bgr);
                    
                end
                okay_button = uicontrol('Parent',f,'Style','pushbutton','Position', [f.Position(3)*0.75,f.Position(4)*0.05, f.Position(3)*0.09,f.Position(4)*0.12],'String', 'OK','Tag','okay_button','Callback','uiresume','enable','off'); %#ok<NASGU>
                uiwait();
                
                data_sets = {};
                data_names = {};
                window_sizes = {};
                used_pr = 0;
                
                for pr = 1:finish
                    if get(protocol_chb{pr},'Value')
                        used_pr = used_pr + 1;
                        idx1 = self.protocol_indices(pr,2)-self.protocol_indices(self.ssd+1,2);
                        idx2 = self.protocol_indices(pr+1,2)-self.protocol_indices(self.ssd+1,2);
                        %check data length
                        %test whether idx are correct
                        if idx1 > self.derived_signals{2}.collect_buff.lst ||idx2 > self.derived_signals{2}.collect_buff.lst
                            disp(['Wrong data indices, protocol ' num2str(pr) ' ' self.feedback_protocols{pr}.protocol_name])
                        else
                            
                            %grab the data
                            
                            for ds = 2:length(self.derived_signals)
                                data_sets{used_pr, ds-1} = self.derived_signals{ds}.collect_buff.raw(self.derived_signals{ds}.collect_buff.fst+idx1:self.derived_signals{ds}.collect_buff.fst+idx2-1,:);
                            end
                            data_names{end+1} = self.feedback_protocols{pr}.protocol_name;
                            if ~isempty(self.feedback_protocols{pr}.window_size)
                                
                                window_sizes{end+1} = self.feedback_protocols{pr}.window_size;
                            else
                                if isempty(window_sizes)
                                    window_sizes{end+1} = 20;
                                else
                                window_sizes{end+1} = window_sizes{end};
                                end
                            end
                        end
                    end
                end
                close(f);
                
                
                
            elseif nargin < 2 && length(self.derived_signals) < 2
                disp('Nothing to recalculate!')
                return
            end
            
            if iscell(data_sets)
                results = cell(size(data_sets));
                for p = 1:size(data_sets,1)
                    window = window_sizes{p};
                    for ds = 1:size(data_sets,2)
                        res = zeros(length(data_sets{p,ds}),1);
                        av = self.feedback_manager.average(ds);
                        s = self.feedback_manager.standard_deviation(ds);
                        for i = window:window:length(data_sets{p,ds})-window
                            dat = data_sets{p,ds}(i-window+1:i);
                            val = sum(abs(dat))/window;
                            res(i-window+1:i) = (val-av)/s;
                        end
                        results{p,ds} = res;
                    end
                end
                
            elseif isnumeric(data_sets)
                if nargin < 3 || isempty(window)
                    window = self.current_window_size;
                end
                if nargin < 4
                    av = 0;
                end
                if nargin < 5
                    stddev = 1;
                end
                results = zeros(size(data_sets));
                for i = window:window:length(data_sets)
                    dat = data_sets(i-window+1:i);
                    val = sum(abs(dat))/window;
                    results(i-window+1:i) = (val-av)/stddev;
                end
                if nargout > 1
                    
                    data_names = {};
                end
                %takes ds_data, window_size, and recalculates feedback using
                %the same algorithm as in RefreshFb function
            end
        end
        function PlotFB(self)
            [results,protocol_names] = self.Recalculate;
            for pr = 1:size(results,1)
                f_pr(pr) = figure; %#ok<NASGU>
                legend_str = {};
                for ds = 1:size(results,2)
                    plot( 1:length(results{pr,ds}), results{pr,ds});
                    hold on;
                    legend_str{end+1} = self.derived_signals{ds+1}.signal_name;
                end
                title(protocol_names{pr});
                legend(legend_str);
            end
        end
        function TransmitToBluetooth(self)
            try %#ok<TRYNC>
                
                ch1 = self.feedback_manager.feedback_vector(1);
                ch2 = self.feedback_manager.feedback_vector(2);
                
                if ~any([isnan(ch1),isinf(ch1)]) && ~any([isnan(ch2),isinf(ch2)])
                    %-1std -->0
                    % mean --> 32
                    % +7std ~ 255
                    %%%
                    out1 = uint8((ch1+1)*32); %#ok<NASGU>
                    out2 = uint8((ch2+1)*32); %#ok<NASGU>
                    
                    % fwrite(self.bluetooth_connection,out1,out2)
                end
            end
        end
        function DeriveSettingsFromFile(self)
            [self.fnames, self.files_pathname, filterindex] = uigetfile('.bin','Select files to play','MultiSelect','on');
            %                     if any([self.fnames, self.files_pathname, filterindex] )
            %                         %subject_folder = self.streams{1}.source_id;
            %                         [protocols, durations, channels] = GetDataProperties(self.files_pathname,self.fnames);
            %                         self.channel_labels = channels;
            %                         if self.run_protocols
            %                             self.protocol_sequence = protocols;
            %                             for pr = 1:length(protocols)
            %
            %                                 self.feedback_protocols{pr} = RealtimeProtocol;
            %                                 self.feedback_protocols{pr}.protocol_name = protocols{pr};
            %                                 self.feedback_protocols{pr}.protocol_duration = durations(pr);
            %
            %                                 if strfind(protocols{pr},'ssd')
            %                                     self.ssd = 1;
            %                                     self.feedback_protocols{pr}.to_update_statistics = 1;
            %                                     self.feedback_protocols{pr}.band = 1;
            %                                 elseif strfind(protocols{pr},'csp')
            %                                     self.ssd = 1;
            %                                     self.feedback_protocols{pr}.to_update_statistics = 1;
            %                                     self.feedback_protocols{pr}.band = 1;
            %                                 elseif strcmpi(protocols{pr},'baseline')
            %                                     self.feedback_protocols{pr}.to_update_statistics = 1;
            %                                 elseif strfind(lower(protocols{pr}),'feedback')
            %                                     self.feedback_protocols{pr}.fb_type = protocols{pr};
            %                                 end
            %                             end
        end
    end
end






function channels = read_montage_file(fname) %#ok<DEFNU>
montage = xml2struct(fname);
channels = {};
for i = 1:length(montage.neorec.transmission.clogicals.clogical)
    channels{end+1} = montage.neorec.transmission.clogicals.clogical{i}.name.Text;
end
end




function channels = get_channel_labels(input) %#ok<DEFNU> %input = inlet obj
ChS = input.info.desc.child('channels');
ch = ChS.first_child;
channels = {};
try
    
    while ch.PtrHandle
        l = ch.child('label');
        channels{end+1} = l.child_value ;
        ch = ch.next_sibling;
    end
catch
    channels = cell(1,input.channel_count());
    for i = 1:input.channel_count()
        channels{i} = num2str(i);
    end
end
ChS =  input.info.desc.child('channels');
ch = ChS.first_child;
channels = {};
try
    
    % while ch.next_sibling.PtrHandle
    while ch.PtrHandle
        l = ch.child('label');
        channels{end+1} = l.child_value ;
        ch = ch.next_sibling;
    end
catch
    channels = cell(1,input.channel_count());
    for i = 1:input.channel_count()
        channels{i} = num2str(i);
    end
end
end
function channels = read_channel_file(fname)%input = txt file

if nargin<1
    fname = 'channels.txt';
end
t = fileread(fname);
ch = strsplit(t, '\n');
if strcmp(ch(end), '')
    
    ch(end) = [];
end
channel_count = str2double(ch{1});
channels = {};
for c = 2:length(ch)
    channels{end+1} = ch{c};
end
channels = strtrim(channels);
if length(channels) ~= channel_count
    disp('Wrong number of channels')
else
    disp('Channels read successfully')
end

end
%         function Run(self)
%             %draw self.raw_and_ds_figure
%              self.raw_and_ds_figure = figure('Tag','raw_and_ds_figure'); %add Tag
%                 set(self.raw_and_ds_figure,'ResizeFcn',@self.FitFigure);
%
%                 settings_button = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton', ...
%                     'String', 'Settings','Tag','settings_button','Callback',@SetExpSettings);
%                 self.connect_button =  uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton',...
%                     'String', 'Start recording','Tag','connect_button');
%                 self.disconnect_button = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton', ...
%                     'String', 'Disconnect', 'Callback', @self.Disconnect,'Tag','disconnect_button');
%                 self.log_text = uicontrol('Parent', self.raw_and_ds_figure  ,'Style', 'Text','String', {'Log'}, 'Tag','log_text');
%                 self.status_text = uicontrol('Parent', self.raw_and_ds_figure,'Style', 'text', 'String', 'Status: ','HorizontalAlignment','left','Tag','status_text');
%                 self.curr_protocol_text = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'text','String', 'Current protocol: ','Tag','curr_protocol_text');
%                 %self.edit_protocols_button = uicontrol('Parent',self.raw_and_ds_figure,'Style','pushbutton','Callback',@self.EditProtocols,'Tag','edit_protocols_button','String','Edit protocols');
%                 self.raw_subplot = subplot(2,1,1);
%                 self.raw_scale_slider = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'slider', 'String','Raw scale', 'Value', 0,  'Max', 24, 'Min',-24,'SliderStep',[1 1],'Callback',@self.SetYScale,'Tag','raw_slider');
%                 self.raw_line = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '','Tag', 'raw_line');
%                 self.ds_subplot = subplot(2,1,2);
%
%                 self.FitFigure;
%
%
%         end


%         function SetExpSettings(self,obj,event)
%             settings_figure = figure;
%             SetDesignFile
%             SetMontageFile
%             SetSubject
%
%         end
%         function ObtainSettings(self,obj,event)

%         end