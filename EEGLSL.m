classdef EEGLSL < handle
    
    properties
        %%%%%%%%%%%%%%%%%%%%%% plot & user interface related members %%%%%%%%%%
        %figures
        fig_interface %settings
        raw_and_ds_figure %plots
        fig_feedback %fb bar
        %settings
        plot_length %sec
        plot_size %samples
        plot_refresh_rate %sec
        fb_refresh_rate %sec
        %%%plots and axes
        %raw eeg subplot
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
        %derived_signals subplot
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
        %feedback subplot
        fbplot_handle %bar
        feedback_axis_handle
        fb_stub %text
        show_fb %bool
        fb_type %string
        
        %%%uicontrol
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
        %data
        channel_labels
        %%%%%%%%%%%% LSL and data input related objects %%%%%%
        streams
        inlet
        data_receive_rate
        path_text
        settings_file
        settings_file_text
        nd %temporarily collects new data
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
        window
        default_window_size
        buffer_length
        ssd
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
    end
    
    methods
        
        function self = EEGLSL(self) %#ok<INUSD>
            self.plot_length = 4;
            self.sampling_frequency = -1;
            self.channel_count = -1;
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
            
            self.last_to_fb = 0;
            self.to_proceed = [];
            self.samples_acquired = 0;
            self.count = 1;
            %self.montage_fname = 'C:\Users\user1\AppData\Local\MCS\NeoRec\nvx136.nvx136.monopolar-Pz';
            self.montage_fname = 'D:\neurofeedback\settings\nvx136.nvx136.monopolar-Pz.xml';
            self.montage_fname_text = 'nvx136.nvx136.monopolar-Pz';
            self.raw_shift = 1;
            self.ds_shift = 1;
            self.sizes = [0]; %#ok<NBRAK>
            self.window = 0;
            self.fb_type = 'Color';
            self.default_window_size = 0;
            self.buffer_length = 0;
            %self.fb_type= 'Color intensity';
            self.ssd = 0;
            self.raw_yscale_fixed = 0;
            self.fb_manager_set = 0;
            self.protocol_indices = 0;
            self.fb_sigmas = 16;
            self.from_file = 0;
            
            self.settings = struct();
            self.settings.subject = 'Null';
            self.settings.montage_fname = 'D:\neurofeedback\settings\nvx136.nvx136.monopolar-Pz.xml';
            self.settings.settings_file =  'settings\LeftVsRightMu.nss.xml';
            self.bad_channels = {};
            self.raw_data_indices = [];
            
            
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
                
                if self.current_protocol >0 && self.current_protocol <= length(self.feedback_protocols)
                    if self.current_protocol <= length(self.feedback_manager.window_size)
                        n = self.feedback_manager.window_size(self.current_protocol);
                    else
                        n = self.default_window_size;
                    end
                elseif self.default_window_size == 0 %zero protocol at the beginning
                    n = self.feedback_manager.window_size(1);
                else
                    n = self.default_window_size;
                end
                
                for s = 2:length(self.derived_signals)
                    dat = self.derived_signals{s}.ring_buff.raw(self.derived_signals{s}.ring_buff.lst-n+1:self.derived_signals{s}.ring_buff.lst);
                    avg  = self.feedback_manager.average(s-1);
                    sdev = self.feedback_manager.standard_deviation(s-1);
                    val = sum(abs(dat))/n;
                    self.feedback_manager.feedback_vector(s-1)  = (val-avg)/sdev;
                end
                
                if self.recording
                    fb = zeros(n,6);
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
            
            if(self.current_protocol>0 && self.current_protocol <= length(self.feedback_protocols))
                if strcmp(self.feedback_protocols{self.current_protocol}.protocol_name,'SSD')
                    self.inlet.close_stream();
                    stop(self.timer_new_data);
                    stop(self.timer_disp);
                    stop(self.timer_fb);
                    
                    if ~isdir(strcat(self.path,'\',self.subject_record.subject_name))
                        mkdir(strcat(self.path,'\',self.subject_record.subject_name));
                    end
                    pth = (strcat(self.path,'\',self.subject_record.subject_name));
                    
                    
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
                    %while ~peaks_found
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
                    
                    %b_ssd = B(middle_point,:);
                    %a_ssd = A(middle_point,:);
                    for i=1:length(w_ssd)
                        %self.derived_signals{1}: the first DS is ALWAYS RAW signal
                        chan_w{i,1} = self.derived_signals{1}.channels{i};
                        chan_w{i,2} = w_ssd(i);
                    end;
                    
                    
                    %write to file
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
                        dummy_signal = struct();
                        dummy_signal.sSignalName = sn;
                        dummy_signal.channels = self.derived_signals{1}.channels;
                        dummy_signal.filters = cell(0,0);
                        full_name = [pth '\SSD_' sn{1} '.xml'];
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
                        
                        
                        NewDS= DerivedSignal(1,dummy_signal, self.sampling_frequency, self.exp_data_length ,self.channel_labels,self.plot_length);
                        NewDS.signal_name = sn{1};
                        NewDS.ring_buff = circVBuf(self.plot_size,1,0);
                        NewDS.collect_buff = circVBuf(self.exp_data_length,1,0);
                        NewDS.UpdateSpatialFilter(spatial_filter);
                        NewDS.UpdateTemporalFilter(Rng(middle_point,:));
                        NewDS.channels_file = full_name;
                        self.derived_signals{end+1} = NewDS;
                    catch
                        'Error while creating a new derived signal, function UpdateStatistics' %#ok<NOPRT>
                    end
                    % end
                    % end
                    if exist('channel_mapping','var')
                        close(channel_mapping);
                    end
                    
                    %calculate feedback avg, std
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if N > self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                        self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.derived_signals{1}.collect_buff.lst - self.derived_signals{1}.collect_buff.fst + 1;
                    else
                        x_raw = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                    end
                    
                    %filter
                    dummy_signal = struct();
                    dummy_signal.sSignalName = 'Temp';
                    channels = cell(size(self.channel_labels,2));
                    for i = 1:length(channels)
                        channels{i,1} = self.channel_labels{i};
                        channels{i,2} = 1;
                    end
                    dummy_signal.filters = cell(0,0);
                    dummy_signal.channels = channels;
                    
                    
                    
                    temp_derived_signal = DerivedSignal(1,dummy_signal, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                    temp_derived_signal.ring_buff = circVBuf(self.plot_size,1,0);
                    temp_derived_signal.collect_buff = circVBuf(self.exp_data_length,1,0);
                    
                    temp_derived_signal.spatial_filter = w_ssd(self.derived_signals{1}.channels_indices)';
                    temp_derived_signal.UpdateTemporalFilter(Rng(middle_point,:));
                    
                    temp_derived_signal.Apply(x_raw',1);
                    values = temp_derived_signal.collect_buff.raw(temp_derived_signal.collect_buff.fst:temp_derived_signal.collect_buff.lst,:);
                    self.feedback_manager.average(1) = mean(values);
                    self.feedback_manager.standard_deviation(1) = std(values);
                    self.SetRawYTicks;
                    self.SetDSYTicks;
                    self.yscales_fixed = 1;
                    self.raw_yscale_fixed = 1;
                    self.ds_yscale_fixed = 1;
                    self.fb_statistics_set = 1;
                    self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
                    self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 6,0);
                    self.fb_manager_set = 1;
                    
                else
                    N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                    if(N>0)
                        for s = 2:length(self.derived_signals)
                            if self.derived_signals{s}.collect_buff.lst - N+1 < self.derived_signals{s}.collect_buff.fst
                                values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.fst:self.derived_signals{s}.collect_buff.lst,:);
                            else
                                values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.lst - N+1:self.derived_signals{s}.collect_buff.lst,:);
                            end
                            self.feedback_manager.average(s-1) = mean(values);
                            self.feedback_manager.standard_deviation(s-1) = std(values);
                        end
                    end;
                    
                    self.SetRawYTicks;
                    self.SetDSYTicks;
                    self.yscales_fixed = 1;
                    self.raw_yscale_fixed = 1;
                    self.ds_yscale_fixed = 1;
                    self.fb_statistics_set = 1;
                end
            end
            
        end
        function Receive(self,timer_obj, event) %#ok<INUSD>
            
            [sample, timestamp] = self.inlet.pull_chunk(); %#ok<ASGLU>
            self.nd = [self.nd sample];
            sz = size(self.nd,2);
            if self.current_protocol > 0 && self.current_protocol <= length(self.feedback_protocols)
                if self.current_protocol <= length(self.feedback_manager.window_size)
                    window_size = self.feedback_manager.window_size(self.current_protocol);
                    self.default_window_size = window_size;
                else
                    window_size = self.default_window_size;
                end
            elseif self.default_window_size == 0 %zero protocol at the beginning
                window_size = self.feedback_manager.window_size(1);
            else  %zero protocol after some data was recorded
                window_size = self.default_window_size;
            end
            if (sz >= window_size)
                for ds = 1:length(self.derived_signals)
                    self.derived_signals{ds}.Apply(self.nd(:,1:window_size),self.recording);
                end
                self.nd =self.nd(:,window_size+1:end);
                try
                    self.UpdateFeedbackSignal;
                catch
                    'Error while updating Feedback' %#ok<NOPRT>
                end
                
                
                self.samples_acquired = self.samples_acquired+window_size;
                if(self.current_protocol>0 && self.current_protocol <= length(self.feedback_protocols))
                    self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.feedback_protocols{self.current_protocol}.actual_protocol_size + window_size;
                    
                    if self.feedback_protocols{self.current_protocol}.actual_protocol_size + window_size >= self.feedback_protocols{self.current_protocol}.protocol_size
                        self.protocol_indices(self.current_protocol+1,:) = self.derived_signals{1}.collect_buff.lst -self.derived_signals{1}.collect_buff.fst +1;
                        try
                            if self.feedback_protocols{self.current_protocol}.to_update_statistics
                                self.Update_Statistics();
                            end
                            temp_log_str = get(self.log_text,'String');
                            temp_log_str{end+1} = self.feedback_protocols{self.current_protocol}.protocol_name;
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
        function Connect(self,predicate, value)
            RecordedStudyToLSL(self.fnames,self.files_pathname,self.looped);
            lsllib = lsl_loadlib();
            %self.channel_labels = [];
            while isempty(self.streams)
                self.streams = lsl_resolve_byprop(lsllib,predicate, value);
            end
            
            self.sampling_frequency = self.streams{1}.nominal_srate();
            if self.max_ch_count <=0
                self.channel_count = self.streams{1}.channel_count();% else default
            else
                self.channel_count = self.max_ch_count;
            end
            self.inlet = lsl_inlet(self.streams{1});
            
            self.plot_size = self.plot_length * self.sampling_frequency;
            
            %set durations and window size based on sampling frequency
            for pr = 1:length(self.feedback_protocols)
                self.feedback_protocols{pr}.protocol_size = self.feedback_protocols{pr}.protocol_duration * self.sampling_frequency;
            end
            if ~self.from_file
                for i = 1:length(self.feedback_protocols)
                    if strcmp(self.feedback_protocols{i}.protocol_name,'SSD')
                        self.ssd = 1;
                    end
                    self.feedback_protocols{i}.protocol_size = self.feedback_protocols{i}.protocol_duration * self.sampling_frequency;
                    try %#ok<TRYNC>
                        self.feedback_protocols{i}.window_size = self.feedback_protocols{i}.window_duration * self.sampling_frequency/1000;
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
                dummy_signal = struct();
                dummy_signal.sSignalName = 'Raw';
                channels = cell(size(self.channel_labels,2));
                for ch = 1:length(channels)
                    channels{ch,1} = self.channel_labels{ch};
                    channels{ch,2} = 1;
                end
                dummy_signal.channels = channels;
                dummy_signal.filters = [];
                self.derived_signals{1} = DerivedSignal(1,dummy_signal, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                self.derived_signals{1}.UpdateSpatialFilter(ones(length(self.channel_labels)));
            else
                if self.ssd
                    self.derived_signals = cell(1,1);
                else
                    self.derived_signals = cell(1,length(self.signals));
                end
                for i = 1: length(self.derived_signals)
                    self.derived_signals{i} = DerivedSignal(1,self.signals{i}, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.plot_length);
                    self.derived_signals{i}.UpdateSpatialFilter(self.signals{i}.channels);
                end
            end
            for ds = 1:length(self.derived_signals)
                if strcmpi(self.derived_signals{ds}.signal_name, 'raw')
                    raw = self.derived_signals{ds};
                    self.used_ch = raw.channels;
                end
            end
            %self.RunInterface;
            if self.from_file
                self.StartRecording();
            end
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
            self.settings_file_text =uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.settings_file_text,'Position', [120 90 200 35],'HorizontalAlignment','left');
            settings_file_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select exp.design', 'Callback', @self.SetDesignFile, 'Position', [20 100 100 35]); %#ok<NASGU>
            set_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Run the experiment', 'Position', [100 20 200 40],'Callback','uiresume'); %#ok<NASGU>
            montage_file_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select exp. montage', 'Callback', @self.SetMontageFile, 'Position', [20 60 100 35]); %#ok<NASGU>
            self.montage_fname_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.montage_fname_text,'Position', [120 60 200 35],'HorizontalAlignment','left');
            show_feedback = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Show feedback to subject','Position', [20 290 135 20],'HorizontalAlignment','left'); %#ok<NASGU>
            show_fb_check = uicontrol('Parent', self.fig_interface, 'Style', 'checkbox' ,'Position', [160 295 20 20],'HorizontalAlignment','left','Value',1);
            % fb_type_string = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Feedback type','Position', [200 290 100 20],'HorizontalAlignment','left');
            %fb_type_menu = uicontrol('Parent', self.fig_interface, 'Style', 'popupmenu', 'String', {'Bar','Color intensity'},'Position', [310 295 100 20],'HorizontalAlignment','left','Value',1);
            self.subjects_dropmenu = uicontrol('Parent', self.fig_interface,'Style','popupmenu','Position',[170 320 100 20],'String',subjects,'Callback',@self.SetSubject);
            sn_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Choose/Enter subject name', 'Position',[20 315 140 20],'HorizontalAlignment','left'); %#ok<NASGU>
            subj_folder_button = uicontrol('Parent', self.fig_interface,'Style','pushbutton','Position',[285 320 150 20],'String','Or select subject folder','Callback',@self.SetSubjectFolder); %#ok<NASGU>
            from_file_chb =  uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 273 20 20]);
            from_file_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[285 270 50 20],'String','From file'); %#ok<NASGU>
            loop_replay_chb = uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 253 20 20]);
            loop_replay_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[220 250 120 20],'String','Loop the recording?'); %#ok<NASGU>
            use_protocols_chb = uicontrol('Parent', self.fig_interface,'Style','checkbox','Position',[340 233 20 20],'Value',1);
            use_protocols_text = uicontrol('Parent', self.fig_interface,'Style','text','Position',[200 230 140 20],'String','Use protocols from files?'); %#ok<NASGU>
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
                    self.fig_interface.Visible = 'off';
                    self.show_fb = show_fb_check.Value;
                    self.from_file = from_file_chb.Value;
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
                    %self.channel_labels = get_channel_labels(self.inlet);
                    [self.fnames, self.files_pathname, filterindex] = uigetfile('.bin','Select files to play','MultiSelect','on');
                    %subject_folder = self.streams{1}.source_id;
                    [protocols, durations, channels] = GetDataProperties(self.files_pathname,self.fnames);
                    self.channel_labels = channels;
                    if self.run_protocols
                        self.protocol_sequence = protocols;
                        for pr = 1:length(protocols)
                            
                            self.feedback_protocols{pr} = RealtimeProtocol;
                            self.feedback_protocols{pr}.protocol_name = protocols{pr};
                            self.feedback_protocols{pr}.protocol_duration = durations(pr);
                            
                            if strcmpi(protocols{pr},'ssd')
                                self.ssd = 1;
                                self.feedback_protocols{pr}.to_update_statistics = 1;
                                self.feedback_protocols{pr}.band = 1;
                            elseif strcmpi(protocols{pr},'baseline')
                                self.feedback_protocols{pr}.to_update_statistics = 1;
                            elseif strfind(lower(protocols{pr}),'feedback')
                                self.feedback_protocols{pr}.fb_type = protocols{pr};
                            end
                        end
                    end
                    
                    
                else
                    self.channel_labels = read_montage_file(self.montage_fname);
                    nfs = NeurofeedbackSession;
                    nfs.LoadFromFile(self.settings_file);
                    self.protocol_types = nfs.protocol_types;
                    self.feedback_protocols = nfs.feedback_protocols;
                    self.signals = nfs.derived_signals;
                    self.protocol_sequence = nfs.protocol_sequence;
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
                %self.edit_protocols_button = uicontrol('Parent',self.raw_and_ds_figure,'Style','pushbutton','Position',[70 400 80 15],'Callback',@self.EditProtocols,'Tag','edit_protocols_button','String','Edit protocols');
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
                        self.r_ytick_labels{end+1} = self.used_ch{i,1};
                    end
                    self.r_ytick_labels{end+1} = ' ';
                    for i = 1:length(self.used_ch)
                        set(self.raw_plot(i),'DisplayName', self.used_ch{i,1});
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
                        set(self.curr_protocol_text, 'String', {strcat('Current protocol: ',self.feedback_protocols{self.current_protocol}.protocol_name),...
                            strcat('Samples acquired', num2str(self.feedback_protocols{self.current_protocol}.actual_protocol_size),'/', num2str(self.feedback_protocols{self.current_protocol}.protocol_size)),...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            });
                    else
                        self.curr_protocol_text.String = {strcat('Current protocol: ',self.feedback_protocols{self.current_protocol}.protocol_name),...
                            strcat('Samples acquired', num2str(self.feedback_protocols{self.current_protocol}.actual_protocol_size),'/', num2str(self.feedback_protocols{self.current_protocol}.protocol_size)),...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            };
                    end
                else %zero protocol
                    if verLessThan('matlab','8.4.0')
                        set(self.curr_protocol_text, 'String', {'Current protocol: idle, ',...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate), ' s')
                            });
                    else
                        self.curr_protocol_text.String = {'Current protocol: idle, ',...
                            strcat(' avg ', num2str(self.feedback_manager.average(self.signal_to_feedback-1))),...
                            strcat(' std ',num2str(self.feedback_manager.standard_deviation(self.signal_to_feedback-1))),...
                            strcat('feedback vector', num2str(self.feedback_manager.feedback_vector(self.signal_to_feedback-1))),...
                            strcat('Receiving samples every ', num2str(self.data_receive_rate), ' s'),...
                            strcat('Updating plots every ', num2str(self.plot_refresh_rate)), ' s'};
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
                if self.current_protocol
                    try
                        self.fb_type = self.feedback_protocols{self.current_protocol}.fb_type;
                    catch
                        'Error while getting fb_type, function RefreshFB' %#ok<NOPRT>
                    end
                end
                try
                    if (self.current_protocol> 0 && self.current_protocol<=length(self.feedback_protocols))
                        set(self.fb_stub,'String',self.feedback_protocols{self.current_protocol}.string_to_show);
                        if isempty(get(self.fb_stub,'String')) && self.show_fb %feedback
                            set(self.fb_stub, 'Visible', 'off'); %string
                            if strfind(lower(self.fb_type),'bar')
                                set(self.fbplot_handle,'Visible','on'); %feedback if bar
                                self.y_limit = [-1 7];
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
                            set(self.fbplot_handle.Visible,'off');
                            
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
                            self.y_limit = [-1 7];
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
                    %self.current_protocol
                    'Error while setting feedback' %#ok<NOPRT>
                    1/16+1/16*self.feedback_manager.feedback_vector(self.signal_to_feedback-1) %#ok<NOPRT>
                end
            elseif self.current_protocol && strcmp(self.feedback_protocols{self.current_protocol}.protocol_name,'SSD')
                set(self.fb_stub,'String',self.feedback_protocols{self.current_protocol}.string_to_show);
                set(self.fb_stub,'Visible','on');
                set(self.feedback_axis_handle,'Visible','off');
                set(self.fbplot_handle,'Visible','off');
            else%zero protocol before baseline recorded
                set(self.feedback_axis_handle,'Visible','off');
                set(self.fbplot_handle,'Visible','off');
                set(self.fb_stub,'Visible','off'); %string
            end
        end
        function StartRecording(self,obj,event) %#ok<INUSD>
            if self.from_file && ~self.run_protocols
            else
                self.current_protocol = self.next_protocol;
                self.next_protocol = self.next_protocol + 1;
            end
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
                temp_log_text = get(self.log_text,'String');
                temp_log_text{end+1} = 'Finished';
                set(self.log_text,'String',temp_log_text);
                set(self.disconnect_button,'String','Disconnect and write');
                set(self.connect_button, 'String', 'Recording finished');
                set(self.connect_button,'Callback','');
            end
            if  ~self.finished
                if self.feedback_protocols{self.current_protocol}.actual_protocol_size*1.1 < self.feedback_protocols{self.current_protocol}.protocol_size
                    self.feedback_protocols{self.current_protocol}.actual_protocol_size = 0;
                    self.next_protocol = self.current_protocol;
                else
                    self.next_protocol = self.current_protocol + 1;
                end
                self.current_protocol = 0;
                set(self.connect_button, 'String', 'Start recording');
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
            self.subject_record.time_stop = datestr(now,13);
            if self.finished
                self.WriteToFile;
            end
        end
        function PlotErrorBar(self)
            if ~strcmpi(self.feedback_protocols{1}.protocol_name,'SSD')
                pr = length(self.feedback_protocols);
                pr_shift = 0;
            else
                pr = length(self.feedback_protocols)-1;
                pr_shift = 1;
            end
            averages = zeros(pr,1);
            deviations =  zeros(pr,1);
            names = cell(pr,1);
            for i = 1+pr_shift:length(self.feedback_protocols)
                idx11 = self.protocol_indices(i,1);
                idx12= self.protocol_indices(i+1,1);
                idx21 = self.protocol_indices(i,2);
                idx22 = self.protocol_indices(i+1,2);
                if ~strcmpi(self.feedback_protocols{i}.protocol_name,'SSD')
                    if idx22-idx21 ~= idx12-idx11
                        warning('Something went wrong... Function PlotErrorBar')
                    end
                    names{i-pr_shift} = self.feedback_protocols{i}.protocol_name;
                    averages(i-pr_shift) = mean(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst+idx21:self.feedback_manager.feedback_records.fst+idx22-1,2));
                    deviations(i-pr_shift) = std(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst+idx21:self.feedback_manager.feedback_records.fst+idx22-1,2));
                end
            end
            f = figure;
            e = errorbar(averages, deviations);
            
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
                idx11 = self.protocol_indices(i,1);
                idx12= self.protocol_indices(i+1,1);
                idx21 = self.protocol_indices(i,2);
                idx22 = self.protocol_indices(i+1,2);
                
                if strcmpi(self.feedback_protocols{i}.protocol_name,'SSD')
                    %correct second indices
                    start = self.protocol_indices(2,2);
                    for j = 2:length(self.feedback_protocols)+1
                        
                        self.protocol_indices(j,2) = self.protocol_indices(j,2) - start;
                    end
                    
                    raw_data_matrix = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx11:self.derived_signals{1}.collect_buff.fst+idx12-1,:);
                    whole_data = raw_data_matrix;
                    inf_file = fopen('ssd_exp_info.hdr','w');
                    fprintf(inf_file,string); %basically writes channels names only
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
                
                filename = [num2str(i) ' ' self.feedback_protocols{i}.protocol_name ' ' num2str(self.feedback_protocols{i}.actual_protocol_size/self.sampling_frequency) '.bin'];
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
                    set(self.status_text,'String',strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.protocol_name, ': ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration)));
                end
            else
                
                if self.from_file && isempty(self.feedback_protocols)
                    self.status_text.String = 'Playing from file';
                elseif self.current_protocol == 0 || self.current_protocol > length(self.feedback_protocols)
                    self.status_text.String = 'Status: receiving';
                else
                    self.status_text.String = strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.protocol_name, ': ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration));
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
            %  epb = findobj('Tag','edit_protocols_button');
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
            set(cpt, 'Position',[0 0.125*fp(4), 0.12*fp(3), 0.32*fp(4)]);
            set(dm,'Position', [0.45*fp(3), 0.015*fp(4),0.12*fp(3),0.04*fp(4)]);
            set(rl,'Position', [0.8 * fp(3), 0.62 *fp(4), 0.05*fp(3), 0.02*fp(4)]);
            set(dsl,'Position', [0.8 * fp(3), 0.15 *fp(4), 0.05*fp(3), 0.02*fp(4)]);
            %  set(epb,'Position', [0.13*fp(3), 0.95*fp(4), 0.12*fp(3), 0.04*fp(4)]);
            set(sb,'Position', [0.1*fp(3), 0.94*fp(4),0.1*fp(3), 0.05*fp(4)]);
            set(ok,'Position', [0.65*fp(3), 0.94*fp(4), 0.1*fp(3),0.05*fp(4)]);
            set(fin,'Position', [0.75*fp(3), 0.94*fp(4), 0.1*fp(3),0.05*fp(4)]);
            set(bcht, 'Position',[fp(3)*0.2,fp(4)*0.75,fp(3)*0.05,fp(4)*0.20]);
            set(cht,'Position', [fp(3)*0.35,fp(4)*0.92,fp(3)*0.3,fp(4)*0.07]);
            set(sbch,'Position', [fp(3)*0.14,fp(4)*0.945,fp(3)*0.2,fp(4)*0.05]);
            self.SetRawYTicks;
            self.SetDSYTicks;
        end
        function SetRawYTicks(self)
            try %#ok<TRYNC>
                r_sp = get(self.raw_subplot);
                r_yticks = [r_sp.YLim(1):(r_sp.YLim(2)-r_sp.YLim(1))/length(self.raw_data_indices):r_sp.YLim(2)]; %#ok<NBRAK>
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
            [~, b, ~ ] = fileparts(uigetdir);
            if b
                self.subject_record.subject_name = b;
                set(self.subjects_dropmenu,'String',[{b} get(self.subjects_dropmenu,'String')']);
                set(self.subjects_dropmenu,'Value',1);
            end
            
        end
        function EditProtocols(self,obj,event) %#ok<INUSD>
            if min(length(findobj('Tag', 'EditProtocolFigure'))) %if it already exists, bring it to front
                uistack(findobj('Tag', 'EditProtocolFigure'));
            else
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
                delete_protocol_dropmenu = uicontrol('Parent',protocol_figure,'Style','popupmenu','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.61,protocol_figure.Position(3)*0.2 , protocol_figure.Position(4)*0.05],'String', del_pr_dpmenu,'Tag','Delete protocol dropmenu'); %#ok<NASGU>
                delete_protocol_pushbutton = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.62,protocol_figure.Position(4)*0.53, protocol_figure.Position(3)*0.2, protocol_figure.Position(4)*0.06],'String', 'Delete','Tag','delete_protocol_button','Callback',@self.DeleteProtocol); %#ok<NASGU>
                
                okay_button = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.75,delta_y/2, protocol_figure.Position(3)*0.09,protocol_figure.Position(4)*0.12],'String', 'Apply','Tag','okay_button','Callback',@self.ChangeProtocols); %#ok<NASGU>
                cancel_button = uicontrol('Parent',protocol_figure,'Style','pushbutton','Position', [protocol_figure.Position(3)*0.85,delta_y/2, protocol_figure.Position(3)*0.09, protocol_figure.Position(4)*0.12],'String', 'Cancel','Tag','cancel_button','Callback',@self.DoNothing); %#ok<NASGU>
                
            end
        end
        function AddProtocol(self,obj,event) %#ok<INUSD>
            
            add_obj = findobj('Tag','Add protocol dropmenu');
            insert_obj = findobj('Tag','Insert protocol dropmenu');
            protocol_to_add = add_obj.String(add_obj.Value);
            
            %add the protocol
            for p = 1:length(self.protocol_types)
                if strcmp(protocol_to_add, self.protocol_types{p}.sProtocolName)
                    protocol = GetProtocolFields(self.protocol_types{p});
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
            
            protocols_names = [protocols_names(1:insert_obj.Value) protocol.protocol_name protocols_names(insert_obj.Value+1:end)];
            protocols_durations = [protocols_durations(1:insert_obj.Value) num2str(protocol.protocol_duration) protocols_durations(insert_obj.Value+1:end)];
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
            
            protocols_names = [protocols_names(1:delete_obj.Value-1) protocols_names(delete_obj.Value+1:end)];
            protocols_durations = [protocols_durations(1:delete_obj.Value) protocols_durations(delete_obj.Value+1:end)];
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
            %mind the numbers!
            for p = 1:length(protocols_names)
                bgr = 0.94-[0.1 0.1 0.1] * mod(p-1,2);
                protocol_count{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.01,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.05, protocol_figure.Position(4)*0.04],'String', num2str(p),'HorizontalAlignment','left','Tag','Protocol count','BackgroundColor',bgr); %#ok<NASGU>
                protocol_name{p} = uicontrol('Parent',protocol_figure,'Style','text','Position', [protocol_figure.Position(3)*0.04,max_height-protocol_figure.Position(4)*0.05*p, protocol_figure.Position(3)*0.25, protocol_figure.Position(4)*0.04],'String', protocols_names{p},'HorizontalAlignment','left','Tag','Protocol name text','BackgroundColor',bgr); %#ok<NASGU>
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
            insert_protocol_dropmenu.String = [{'0'} new_protocols_names];
            delete_protocol_dropmenu.String = new_protocols_names;
            
            
            
        end
        function ChangeProtocols(self,obj,event) %#ok<INUSD>
            %get the protocols
            self.feedback_protocols(self.next_protocol:end) = [];
            protocols_names_obj = findobj('Tag', 'Protocol name text');
            protocols_durations_obj = findobj('Tag', 'Protocol duration text');
            for j = length(protocols_names_obj):-1:self.next_protocol
                for i = 1:length(self.protocol_types)
                    if strcmp(protocols_names_obj(j).String,self.protocol_types{i}.sProtocolName)
                        rtp = RealtimeProtocol;
                        rtp.protocol_name = self.protocol_types{i}.sProtocolName;
                        rtp.to_update_statistics = self.protocol_types{i}.bUpdateStatistics;
                        rtp.protocol_duration = str2double(protocols_durations_obj(j).String);  %%%%%duration is taken from the figure
                        rtp.protocol_size = rtp.protocol_duration * self.sampling_frequency; %%% size is calculated accordingly
                        rtp.stop_after = self.protocol_types{i}.bStopAfter;
                        rtp.string_to_show = self.protocol_types{i}.cString;
                        try %#ok<TRYNC>
                            
                            rtp.filter_filename = self.protocol_types{i}.sFilterFilename;
                            
                        end
                        try %#ok<TRYNC>
                            rtp.band = self.protocol_types{i}.dBand;
                        end
                        try %#ok<TRYNC>
                            rtp.fb_type = self.protocol_types{i}.sFb_type;
                        end
                        try %#ok<TRYNC>
                            rtp.window_size = self.protocol_types{i}.nMSecondsPerWindow;
                        end
                        
                        self.feedback_protocols{end+1} = rtp;
                        
                    end
                end
            end
            
            %get their durations
            %set new protocol duration and size (s and samples acc.)
            %             for p = self.next_protocol:length(self.feedback_protocols)
            %                 self.feedback_protocols{p}.protocol_duration = str2double(self.protocol_duration_text{p}.String);
            %                 self.feedback_protocols{p}.protocol_size = self.feedback_protocols{p}.protocol_duration * self.sampling_frequency;
            %             end
            
            %check if we reserved enough space
            data_length = 0;
            for p = 1:length(self.feedback_protocols)
                data_length = data_length + self.feedback_protocols{p}.protocol_size;
            end
            if data_length > self.exp_data_length
                self.exp_data_length = fix(data_length*1.1);
                for ds = 1:length(self.derived_signals)
                    new_circbuff = circVBuf(self.exp_data_length, size(self.derived_signals{ds}.collect_buff.raw,2),0);
                    new_circbuff.append(self.derived_signals{ds}.collect_buff.raw(self.derived_signals{ds}.collect_buff.fst:self.derived_signals{ds}.collect_buff.lst,:));
                    self.derived_signals{ds}.collect_buff = new_circbuff;
                end
                new_fb_records = circVBuf(self.exp_data_length,6,0);
                new_fb_records.append(self.feedback_manager.feedback_records.raw(self.feedback_manager.feedback_records.fst:self.feedback_manager.feedback_records.lst,:));
                self.feedback_manager.feedback_records = new_fb_records;
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
                
                vSignals.DerivedSignal(ds).sSignalName = self.derived_signals{ds}.signal_name;
                
                try %#ok<TRYNC>
                    vSignals.DerivedSignal(ds).fBandpassLowHz = self.derived_signals{ds}.temporal_filter{1}.range(1);
                    vSignals.DerivedSignal(ds).fBandpassHighHz = self.derived_signals{ds}.temporal_filter{1}.range(2);
                end
                
                try %#ok<TRYNC>
                    vSignals.DerivedSignal(ds).sType = self.derived_signals{ds}.type;
                end
                %prepare spatial filter matrix
                spatial_filter_matrix_struct = struct();
                for ch = 1:length(self.derived_signals{ds}.channels)
                    spatial_filter_matrix_struct.channels.(self.derived_signals{ds}.channels{ch,1}) = num2str(self.derived_signals{ds}.channels{ch,2} );
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
            a = struct2xml(self.exp_design);
            
            design_file = fopen(strcat(self.path,'\',self.subject_record.subject_name,'\',curr_date,'\',self.subject_record.time_start,'\Exp_design.xml'),'w');
            fwrite(design_file,a);
            fclose(design_file);
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
                    if finished
                        break; %#ok<UNRCH>
                    end
%                     %%%set it gray
%                     if ~isempty(getCursorInfo(dcm_obj))
%                         cursor_info = getCursorInfo(dcm_obj);
%                         
%                         for ch = 1:length(self.derived_signals{1}.channels)
%                             if strcmp(self.derived_signals{1}.channels{ch}, cursor_info.Target.DisplayName)
%                                 set(r_sp.Children(length(r_sp.Children)-ch+1),'Color',[0.75 0.75 0.75]);
%                             end
%                         end
%                     end
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
            %update derived signals
            
            for b_ch = 1:length(self.bad_channels)
                for child = 1:length(r_sp.Children)
                    if strcmp(self.bad_channels{b_ch},r_sp.Children(length(r_sp.Children)-child+1).DisplayName)
                        self.raw_data_indices = self.raw_data_indices(self.raw_data_indices ~= child);
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
            self.SetRawYTicks();

        end
        %         function SetExpSettings(self,obj,event)
        %             settings_figure = figure;
        %             SetDesignFile
        %             SetMontageFile
        %             SetSubject
        %
        %         end
        %         function ObtainSettings(self,obj,event)
        %         end
    end
end


function channels = read_montage_file(fname)
montage = xml2struct(fname);
channels = {};
for i = 1:length(montage.neorec.transmission.clogicals.clogical)
    channels{end+1} = montage.neorec.transmission.clogicals.clogical{i}.name.Text;
end
end
function rtp = GetProtocolFields(protocol_type)
rtp = RealtimeProtocol;
rtp.protocol_name = protocol_type.sProtocolName;
rtp.to_update_statistics = protocol_type.bUpdateStatistics;
rtp.protocol_duration = protocol_type.fDuration;
rtp.stop_after = protocol_type.bStopAfter;
rtp.string_to_show = protocol_type.cString;
try %#ok<TRYNC>
    rtp.filter_filename = protocol_type.sFilterFilename;
end
try %#ok<TRYNC>
    rtp.band = protocol_type.dBand;
end
try %#ok<TRYNC>
    rtp.fb_type = protocol_type.sFb_type;
end
try %#ok<TRYNC>
    rtp.window_size = protocol_type.nMSecondsPerWindow;
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
% function channels = read_channel_file()%input = txt file
%
% fname = 'mitsar_channels.txt';
% t = fileread('mitsar_channels.txt');
% channels = {};
% str = '';
% k = 1;
% j = 1;
% while true
%     if k>length(t)
%         break
%     end
%     if strcmp(t(k), ' ')
%         channels{end+1} = {t(j:k-1)};
%         j = k+1;
%     end
%     k = k + 1;
% end



%                 if self.feedback_manager.standard_deviation
%                     self.y_limit = [self.feedback_manager.average/self.feedback_manager.standard_deviation*0.9,self.feedback_manager.average/self.feedback_manager.standard_deviation*1.1];
%                 else
%                     self.y_limit = [self.feedback_manager.average*0.1,self.feedback_manager.average*10];
%                 end
%self.y_limit = [0 ,self.feedback_manager.average/self.feedback_manager.standard_deviation*1.1];

% axes(self.fb_plot_handle);
%                self.y_limit = [-4 4];
%                ylim(self.feedback_axis_handle,self.y_limit);
%                xlim(self.feedback_axis_handle, [1 3]);
%protocol.feedback_avg = self.feedback_manager.average;
%protocol.feedback_std = self.feedback_manager.standard_deviation;
%self.feedback_manager.samples_acquired = 0;
%             self.ds_plot_min= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 50 100 25],'Tag', 'ds_plot_min');
%             self.ds_plot_max= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 250 100 25],'Tag', 'ds_plot_max');
%             self.ds_plot_shift= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 150 100 25],'Tag', 'ds_plot_shift');
%             self.raw_plot_min= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 350 100 25],'Tag', 'raw_plot_min');
%             self.raw_plot_max= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 550 100 25],'Tag', 'raw_plot_max');
%             self.raw_plot_shift= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [600 450 100 25],'Tag', 'raw_plot_shift');
%             %self.raw_fit_plot =uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [500 280 130 20], ...
%'String', 'Fit raw plot', 'Callback', @self.FitPlot,'Tag','raw_fit_plot');
%self.ds_fit_plot = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [420 50 130 20], ...
%'String', 'Fit DS plot', 'Callback', @self.FitPlot,'Tag','ds_fit_plot');
%f = fopen('cm_136_5.bin','w'), fwrite(f,size(self.composite_montage),'int'); fwrite(f,self.composite_montage,'float'); fclose(f);
%f = fopen('cm_136_5.bin','r'), sz = fread(f,2,'int'); A = fread(f,sz,'float'); fclose(f);


%         function FitPlot(self,obj,event)
%             if strcmp(get(obj, 'Tag'), 'raw_fit_plot')
%                 [signals_averages, signals_stds] = self.CalculateStats('Raw');
%                 m = max(signals_averages);
%                 s = abs(max(signals_stds));
%                 self.raw_shift = abs(m + 6*s);
%                 set(self.raw_subplot,'YLim',[(m - 3*s) (m - 3*s +self.raw_shift*(length(self.used_ch)+1))]);
%                 self.raw_yscale_fixed = 1;
%             elseif strcmp(get(obj,'Tag'),'ds_fit_plot')
%                 [signals_averages, signals_stds] = self.CalculateStats('DS');
%                 m = max(signals_averages);
%                 s = abs(max(signals_stds));
%                 self.ds_shift = m + 6*s;
%                 set(self.ds_subplot,'YLim',[(m - 3*s) (m - 3*s+self.ds_shift*length(self.derived_signals))]);
%                 self.ds_yscale_fixed = 1;
%             end
%         end
%         function [signals_averages, signals_stds] = CalculateStats(self,signals)
%             signals_averages = [];
%             signals_stds = [];
%             if strcmp(signals,'Raw')
%                 signals_averages = zeros(length(self.used_ch),1);
%                 signals_stds = zeros(length(self.used_ch),1);
%                 for signal = 1:length(self.used_ch)
%                     signals_averages(signal) = mean(self.derived_signals{1}.ring_buff.raw(self.derived_signals{1}.ring_buff.fst:self.derived_signals{1}.ring_buff.lst,signal));
%                     signals_stds(signal) = std(self.derived_signals{1}.ring_buff.raw(self.derived_signals{1}.ring_buff.fst:self.derived_signals{1}.ring_buff.lst,signal));
%                 end
%             elseif strcmp(signals,'DS')
%                 signals_averages = zeros(length(self.derived_signals)-1,1);
%                 signals_stds = zeros(length(self.derived_signals)-1,1);
%                 for signal = 2:length(self.derived_signals)
%                     signals_averages(signal-1)= mean(self.derived_signals{signal}.ring_buff.raw(self.derived_signals{signal}.ring_buff.fst:self.derived_signals{signal}.ring_buff.lst));
%                     signals_stds(signal-1) = std(self.derived_signals{signal}.ring_buff.raw(self.derived_signals{signal}.ring_buff.fst:self.derived_signals{signal}.ring_buff.lst));
%                 end
%             end
%
%         end
