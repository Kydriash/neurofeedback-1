classdef EEGLSL < handle
    
    properties
        %%%%%%%%%%%%%%%%%%%%%% plot & user interface related members %%%%%%%%%%
        %figures
        fig_interface
        raw_and_ds_figure
        fig_feedback
        % settings
        plot_length %sec
        plot_size %samples
        plot_refresh_rate
        %plots and axes
        fbplot_handle
        raw_subplot
        raw_shift
        ds_subplot
        feedback_axis_handle
        raw_plot
        ds_plot
        ds_ytick_labels
        r_ytick_labels
        ds_y_text
        r_y_text
        raw_legend
        ds_legend
        fb_stub
        raw_ydata_scale
        ds_ydata_scale
        raw_line
        ds_line

        %data
        ds_shift
        channel_labels
        
        
        
        %uicontrol
        
        text_buff_size
        buff_size
        show_button
        
        connect_button
        disconnect_button
        curr_protocol_text
        sn_to_fb_dropmenu
        sn_to_fb_string
        y_limit
        ds_inc_button
        ds_decr_button
        raw_inc_button
        raw_decr_button

        status_text
        status_field
        log_text
        
        add_notes_window
        add_notes_field
        write_notes
        raw_scale_slider
        ds_scale_slider
        ds_scaling
        raw_scaling
        %%%%%%%%%%%% LSL and data input related objects %%%%%%
        streams
        inlet
        chunksize
        chunkcount
        samples_collected
        data_receive_rate
        path_text
        settings_file
        settings_file_text
        nd
        %%%%%%%%%%%%%%%%%%%%% Data info %%%%%%%%%%%%%%%%%%%%%%
        sampling_frequency
        channel_count
        
        max_ch_count
        %%%%%%%%%%%%%%%%%%%%% Timers and callbacks %%%%%%%%%%%
        timer_new_data
        timer_disp
        timer_new_data_function
        timer_disp_function
        %%%%%%%%%%%%%%%%%%%% Status members
        connected
        recording
        finished
        ylabels_fixed
        yscales_fixed
        fb_statistics_set
        %%%%%% Signal processing and feedback related members%
        composite_montage
        signals %parameters
        derived_signals
        feedback_protocols
        feedback_manager
        
        current_protocol
        next_protocol
        protocol_sequence
        %fbwindow_size
        samples_acquired
        signal_to_feedback
        used_ch
        %%%%%% subject data management related
        subject_record
        path
        exp_data_length
        %other
        tstop
    end
    
    methods
        
        function self = EEGLSL(self)
            self.plot_length = 4;
            self.sampling_frequency = -1;
            self.channel_count = -1;
            self.streams = {};
            self.plot_refresh_rate = 0.075;
            self.data_receive_rate = 0.02;
            self.timer_new_data_function = @self.Receive;
            self.timer_new_data = timer('Name','receive_data','TimerFcn', self.timer_new_data_function,'ExecutionMode','fixedRate',...
                'Period',self.data_receive_rate);
            self.timer_disp_function = @self.PlotEEGData;
            self.timer_disp = timer('Name','plot_data','TimerFcn', self.timer_disp_function,'ExecutionMode','fixedRate',...
                'Period', self.plot_refresh_rate);
            self.max_ch_count = -1; % -1 to get all the channels
            self.connected = false;
            self.fig_feedback = figure('Visible', 'off');
            self.channel_labels = {};
            self.feedback_manager = FeedbackManager;
            self.current_protocol = 0;
            self.feedback_protocols = {};
            self.exp_data_length = 0;
            self.samples_acquired = 0;
            self.y_limit = [0 4];
            self.subject_record = SubjectRecord;
            self.path = 'D:\anna\exps';
            self.signal_to_feedback = 2;
            self.composite_montage = [];
            self.settings_file_text = 'LeftVsRightMu.nss.xml';
            self.settings_file =  'settings\LeftVsRightMu.nss.xml';
            self.recording = 0;
            self.next_protocol = 1;
            self.finished = 0;
            self.ylabels_fixed = 0;
            self.yscales_fixed = 0;
            self.fb_statistics_set = 0;
            self.raw_ydata_scale = 1;
            self.ds_ydata_scale = 1;
            self.nd = [];
            
            
        end
        
        function UpdateFeedbackSignal(self)
            for s = 2:length(self.derived_signals)
                n = self.feedback_manager.window_length;
                dat = self.derived_signals{s}.ring_buff.raw(self.derived_signals{s}.ring_buff.lst-n+1:...
                    self.derived_signals{s}.ring_buff.lst,:);% NaN O_o
                
                avg  = self.feedback_manager.average(s-1);
                sdev = self.feedback_manager.standard_deviation(s-1);
                val = sum(abs(dat))/n;
                self.feedback_manager.feedback_vector(s-1)  = (val-avg)/sdev;
                
            end
        end
        
        function Update_Statistics(self)
            if(self.current_protocol>0)
                N = self.feedback_protocols{self.current_protocol}.actual_protocol_size;
                if(N>0)
                    ds_minimum = Inf;
                    ds_maximum = -Inf;
                    for s = 2:length(self.derived_signals)
                        if self.derived_signals{s}.collect_buff.lst - N+1 < self.derived_signals{s}.collect_buff.fst
                            values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.fst:self.derived_signals{s}.collect_buff.lst,:);
                        else
                            values = self.derived_signals{s}.collect_buff.raw(self.derived_signals{s}.collect_buff.lst - N+1:self.derived_signals{s}.collect_buff.lst,:);
                        end
                        self.feedback_manager.average(s-1) = mean(values);
                        self.feedback_manager.standard_deviation(s-1) = std(values);
                        if ds_minimum > min(values)
                            ds_minimum = min(values);
                        end
                        if ds_maximum < max(values)
                            ds_maximum = max(values);
                        end
                        
                    end
                    
                    
                    self.fb_statistics_set = 1;
                    
                    ds_d = ds_maximum - ds_minimum;
                    set(self.ds_subplot,'YLim',[ds_minimum ds_d*length(self.derived_signals)*2]);
                    self.ds_shift = ds_d*2;
                    raw_minimum = Inf;
                    raw_maximum = -Inf;
                    
                    for ch = 1:length(self.used_ch)
                        if self.derived_signals{1}.collect_buff.lst - N+1 < self.derived_signals{1}.collect_buff.fst
                            v = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst:self.derived_signals{1}.collect_buff.lst,:);
                            
                        else
                            v = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.lst - N+1:self.derived_signals{1}.collect_buff.lst,:);
                        end
                        values = v(:,ch);
                        if raw_minimum > min(values)
                           raw_minimum = min(values);
                        end
                        if raw_maximum < max(values)
                            raw_maximum = max(values);
                        end
                    end
                    raw_d = raw_maximum - raw_minimum;
                     set(self.raw_subplot,'YLim',[raw_minimum raw_minimum+raw_d*length(self.used_ch)*2]);
                    self.raw_shift = raw_d*2;
                     self.SetRawYTicks;
                     self.SetDSYTicks;
                    self.yscales_fixed = 1;
                end;
            end;
        end
        
        function Receive(self,timer_obj, event)
            
            [sample, ~] = self.inlet.pull_chunk();
            self.nd = [self.nd sample];
            if(size(self.nd,2)>5)
                
                %removing DC component
                %sample = sample - repmat(sum(sample,1)/length(self.channel_count),size(sample,1),1);
                
                % apply average reference montage
                self.nd = self.composite_montage*self.nd;
                
                for ds = 1:length(self.derived_signals)
                    self.derived_signals{ds}.Apply(self.nd,self.recording);
                    
                    
                end;
                self.UpdateFeedbackSignal;
                fb = zeros(size(self.nd,2),4);
                fb(:,1) = self.signal_to_feedback;
                fb(:,2) = self.feedback_manager.feedback_vector(self.signal_to_feedback-1);
                fb(:,3) = self.feedback_manager.average(self.signal_to_feedback-1);
                fb(:,4) = self.feedback_manager.standard_deviation(self.signal_to_feedback-1);
                self.feedback_manager.feedback_records.append(fb);
                self.nd =[];
            else
                
            end;
            self.samples_acquired = self.samples_acquired+size(sample,2);
            
            if(self.current_protocol>0 && self.current_protocol <= length(self.feedback_protocols))
                self.feedback_protocols{self.current_protocol}.actual_protocol_size = self.feedback_protocols{self.current_protocol}.actual_protocol_size + size(sample,2);
                try
                    if self.feedback_protocols{self.current_protocol}.actual_protocol_size >= self.feedback_protocols{self.current_protocol}.protocol_size
                        temp_log_str = get(self.log_text,'String');
                        temp_log_str{end+1} = self.feedback_protocols{self.current_protocol}.protocol_name;
                        set(self.log_text,'String', temp_log_str);
                        if self.feedback_protocols{self.current_protocol}.to_update_statistics
                            self.Update_Statistics();
                        end
                        if self.feedback_protocols{self.current_protocol}.stop_after
                            set(self.connect_button, 'String', 'Start recording');
                            set(self.connect_button, 'Callback',@self.StartRecording);%%%%%
                            self.StopRecording();
                            
                        else
                            self.current_protocol = self.next_protocol;
                            self.next_protocol = self.next_protocol + 1;
                            if self.current_protocol > length(self.feedback_protocols)
                                self.StopRecording();
                            end
                            
                            
                        end
                        
                    end
                catch
                end
                
            end
        end
        
        
        function Connect(self,predicate, value)
            lsllib = lsl_loadlib();
            while isempty(self.streams)
                self.streams = lsl_resolve_byprop(lsllib,predicate, value);
            end
            self.sampling_frequency = self.streams{1}.nominal_srate();
            %set protocol size
            if self.max_ch_count <=0
                self.channel_count = self.streams{1}.channel_count();% else default
            else
                self.channel_count = self.max_ch_count;
            end
            self.inlet = lsl_inlet(self.streams{1});
            %self.channel_labels = get_channel_labels(self.inlet);
            self.channel_labels = read_channel_file();
            self.plot_size = self.plot_length * self.sampling_frequency;
            self.RunInterface;
            
            self.InitTimer();
            
        end
        
        function InitTimer(self)
            
            if strcmp(self.timer_new_data.Running,'off')
                start(self.timer_new_data);
            end
            if strcmp(self.timer_disp.Running,'off')
                start(self.timer_disp);
            end
            tic
            set(self.connect_button, 'String', 'Start recording');
            set(self.connect_button, 'Callback',@self.StartRecording);
            
        end
        
        function RunInterface(self)
            self.fig_interface = figure;
            prr_text = uicontrol('Parent',self.fig_interface, 'Style', 'text', 'String', 'Plot refresh rate, s', 'Position',[20 250 120 30],'HorizontalAlignment','left'); %#ok<NASGU>
            prr = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', num2str(self.plot_refresh_rate), 'Position', [125 250 50 30]);
            drr_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Data receive rate, s', 'Position', [20 210 100 30],'HorizontalAlignment','left'); %#ok<NASGU>
            drr = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', num2str(self.data_receive_rate), 'Position', [125 210 50 30]);
            sn_text = uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', 'Subject name', 'Position',[20 170 100 20],'HorizontalAlignment','left'); %#ok<NASGU>
            sn = uicontrol('Parent', self.fig_interface, 'Style', 'edit', 'String', self.subject_record.subject_name, 'Position', [125 170 100 20]);
            self.path_text =uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.path,'Position', [120 125 200 35],'HorizontalAlignment','left');
            path_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select path', 'Callback', @self.SetWorkpath, 'Position', [20 135 100 35]); %#ok<NASGU>
            self.settings_file_text =uicontrol('Parent', self.fig_interface, 'Style', 'text', 'String', self.settings_file_text,'Position', [120 90 200 35],'HorizontalAlignment','left');
            settings_file_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Select exp.design', 'Callback', @self.SetDesignFile, 'Position', [20 100 100 35]); %#ok<NASGU>
            set_button = uicontrol('Parent',self.fig_interface,'Style', 'pushbutton', 'String', 'Run the experiment', 'Position', [100 20 200 40],'Callback','uiresume'); %#ok<NASGU>
            uiwait();
            if verLessThan('matlab','8.4.0')
                self.plot_refresh_rate = str2num(get(prr,'String'));
                self.data_receive_rate = str2num(get(drr,'String'));
                self.subject_record.subject_name = get(sn,'String');
                set(self.fig_interface,'Visible', 'off');
            else
                self.plot_refresh_rate = str2num(prr.String);
                self.data_receive_rate = str2num(drr.String);
                self.subject_record.subject_name = sn.String;
                self.fig_interface.Visible = 'off';
            end
            
            
            nfs = NeurofeedbackSession;
            nfs.LoadFromFile(self.settings_file);
            self.feedback_protocols = nfs.feedback_protocols;
            self.signals = nfs.derived_signals;
            self.protocol_sequence = nfs.protocol_sequence;
            
            for i = 1:length(self.feedback_protocols)
                self.feedback_protocols{i}.protocol_size = self.feedback_protocols{i}.protocol_duration * self.sampling_frequency;
            end
            
            for j = 1:length(self.feedback_protocols)
                self.exp_data_length = self.exp_data_length + self.feedback_protocols{j}.protocol_size;
            end
            self.exp_data_length = fix(self.exp_data_length * 1.1); %just in case
            self.derived_signals = cell(1,length(self.signals));
            for i = 1: length(self.signals)
                self.derived_signals{i} = DerivedSignal(1,self.signals{i}, self.sampling_frequency,self.exp_data_length,self.channel_labels,self.channel_count,self.plot_length);
            end
            self.feedback_manager.window_length = round(self.feedback_manager.window_size*self.sampling_frequency / 1000);
            self.feedback_manager.average = zeros(1,length(self.derived_signals)-1);
            self.feedback_manager.standard_deviation = ones(1,length(self.derived_signals)-1);
            self.feedback_manager.feedback_vector = zeros(1,length(self.derived_signals)-1);
            self.feedback_manager.feedback_records = circVBuf(self.exp_data_length, 4,0);
            self.sn_to_fb_string = '';
            for i = 2:length(self.derived_signals)
                if i == 2
                    self.sn_to_fb_string = self.derived_signals{i}.signal_name;
                else
                    self.sn_to_fb_string = strcat(self.sn_to_fb_string,'|',self.derived_signals{i}.signal_name);
                end
            end
            figure(self.fig_feedback);
            set(self.fig_feedback, 'OuterPosition', [-1280 0 1280 1024]);
            self.feedback_axis_handle = axes;
            self.fbplot_handle = bar(self.feedback_axis_handle,[0 1 0],'FaceColor',[1 1 1]);
            
            self.fb_stub = uicontrol('Parent', self.fig_feedback, 'String', 'Baseline acquisition', 'Style', 'text', 'ForegroundColor',[0 1 0],'Position', [200 500 900 250], 'FontSize', 75, 'BackgroundColor',[1 1 1], 'FontName', 'Courier New', 'Visible', 'off' );
            % assume standardized FB signal
            self.y_limit = [-4 4];
            ylim(self.feedback_axis_handle,self.y_limit);
            xlim(self.feedback_axis_handle, [1 3]);
            %set average reference (composite montage)
            for ds = 1:length(self.derived_signals)
                if strcmpi(self.derived_signals{ds}.signal_name, 'raw')
                    raw = self.derived_signals{ds};
                    self.used_ch = raw.channels;
                end
            end
            all_ch = self.channel_labels;
            self.composite_montage = zeros(length(all_ch),length(all_ch));
            used_ch_labels = self.used_ch(:,1);
            for i = 1:length(used_ch_labels)
                for j = 1:length(all_ch)
                    if strcmp(self.used_ch{i,1},all_ch{j})
                        self.composite_montage(:,j) = 1/length(self.used_ch);
                        self.composite_montage(j,j) = 1-1/length(self.used_ch);
                    end
                end
            end
            %connect
            self.connected = 1;
            self.raw_and_ds_figure = figure;
            set(self.raw_and_ds_figure,'ResizeFcn',@self.FitFigure);
            self.connect_button =  uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [10 10 150 20], ...
                'String', 'Start recording','Tag','connect_button');
            self.disconnect_button = uicontrol('Parent',self.raw_and_ds_figure,'style','pushbutton','Position', [420 10 130 20], ...
                'String', 'Disconnect', 'Callback', @self.Disconnect,'Tag','disconnect_button');
            self.sn_to_fb_dropmenu = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'popupmenu', 'String', self.sn_to_fb_string, 'Position',[300 10 100 20], 'Callback', @self.SelectSignalToFeedback,'Tag','sn_to_fb_dropmenu');
            self.log_text = uicontrol('Parent', self.raw_and_ds_figure  ,'Style', 'Text','String', {'Log'}, 'Position', [0 300 50 100],'Tag','log_text');
            self.status_text = uicontrol('Parent', self.raw_and_ds_figure,'Style', 'text', 'String', 'Status: ', 'Position', [0 210 200 20],'HorizontalAlignment','left','Tag','status_text');
            self.curr_protocol_text = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'text','String', 'Current protocol: ', 'Position', [0 40  190 100],'Tag','curr_protocol_text');
              self.raw_subplot = subplot(2,1,1);
            self.ds_subplot = subplot(2,1,2);
            self.raw_scale_slider = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'slider', 'String','Raw scale', 'Value', 0, 'Position', [520 300 10 100], 'Max', 16, 'Min',-16,'SliderStep',[1 1],'Callback',@self.SetYScale,'Tag','raw_slider');
             
             self.ds_scale_slider= uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'slider', 'String','DS scale', 'Value', 0, 'Position', [520 100 10 100], 'Max', 16, 'Min',-16,'SliderStep',[1 1],'Callback',@self.SetYScale,'Tag','ds_slider');

            ds_temp = zeros(length(self.derived_signals),fix(self.plot_size));
            r_temp = zeros(self.channel_count,fix(self.plot_size));
            self.raw_plot = plot(r_temp', 'Parent', self.raw_subplot);
            self.ds_plot = plot(ds_temp', 'Parent', self.ds_subplot);
            self.raw_line = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [480 320 100 25],'Tag', 'raw_line');
            self.ds_line = uicontrol('Parent', self.raw_and_ds_figure, 'Style', 'Text','String', '', 'Position', [480 120 100 25],'Tag','ds_line');


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
            [fname fpath fspec] = uigetfile('*.*');
            if ~isempty(nonzeros([fname fpath fspec]))
                self.settings_file = strcat(fpath,fname);
                if verLessThan('matlab','8.4.0')
                    set(self.settings_file_text, 'String',self.settings_file);
                else
                    self.settings_file_text.String = self.settings_file;
                end
            end
        end
        function PlotEEGData(self,timer_obj, event)
            r_sp = get(self.raw_subplot);
            ds_sp = get(self.ds_subplot);
            self.raw_shift = (r_sp.YLim(2)-r_sp.YLim(1))/(length(self.used_ch)+1);
            self.ds_shift = (ds_sp.YLim(2)-ds_sp.YLim(1))/(length(self.derived_signals));
            if ~self.ylabels_fixed
                self.ds_ytick_labels = {' '};
                self.r_ytick_labels = {' '};
                for i = 2:length(self.derived_signals)
                    self.ds_ytick_labels{end+1} = self.derived_signals{i}.signal_name;
                end
                self.ds_ytick_labels{end+1} = ' ';
                for i = 1:length(self.used_ch)
                    self.r_ytick_labels{end+1} = self.used_ch{i,1};
                end
                self.r_ytick_labels{end+1} = ' ';
                for i = 1:length(self.derived_signals)
                    set(self.ds_plot(i),'DisplayName', self.derived_signals{i}.signal_name);
                end
                for i = 1:length(self.used_ch)
                    set(self.raw_plot(i),'DisplayName', self.used_ch{i,1});
                end
                self.ylabels_fixed = 1;
            end
            if self.connected
                first_to_show = self.derived_signals{self.signal_to_feedback}.ring_buff.lst-self.plot_size;
                last_to_show = self.derived_signals{self.signal_to_feedback}.ring_buff.lst;
                if last_to_show > first_to_show
                      
                    if ~self.yscales_fixed
                        ds_min = Inf;
                        ds_max = 0;
                        r_min = Inf;
                        r_max = 0;
                        raw_data = self.derived_signals{1}.ring_buff.raw(self.derived_signals{1}.ring_buff.lst-self.plot_size+1:self.derived_signals{1}.ring_buff.lst,:);
                        raw_data = raw_data';
                         for i = 1:size(raw_data,1)

                    
                            if min(raw_data(i:i,:)) < r_min
                                r_min = min(raw_data(i:i,:));
                            end
                            if max(raw_data(i:i,:)) > r_min
                                r_max = max(raw_data(i:i,:));
                            end
                            r_d = r_max - r_min;
                            set(self.raw_plot(i),'YData', raw_data(i:i,:)*self.raw_ydata_scale+self.raw_shift*i);
                        end
                        for i = 2:length(self.derived_signals)
                            pulled = self.derived_signals{i}.ring_buff.raw(first_to_show:last_to_show,:);
                            pulled = pulled';
                            if min(pulled) < ds_min
                                ds_min = min(pulled);
                            end
                            if max(pulled) > ds_max
                                ds_max = max(pulled);
                            end
                            ds_d = ds_max - ds_min;
                            set(self.ds_plot(i-1), 'YData',pulled(1,:)*self.ds_ydata_scale+self.ds_shift*(i-1)); %NaN
                        end
                        if all([isnumeric(ds_d),~isinf(ds_d)])
                            ylim(self.ds_subplot,[ds_d/2 (2*length(self.derived_signals)+1)*ds_d]);
                        elseif all([isnumeric(r_d),~isinf(r_d)])
                            ylim(self.raw_subplot,[r_d/2 2*(length(self.used_ch)+2)*r_d]);
                            
                        end
            self.SetRawYTicks;
            self.SetDSYTicks;
                    else
                        raw_data = self.derived_signals{1}.ring_buff.raw(self.derived_signals{1}.ring_buff.lst-self.plot_size+1:self.derived_signals{1}.ring_buff.lst,:);
                        raw_data = raw_data';
                        for i = 1:size(raw_data,1)
                            set(self.raw_plot(i),'YData', raw_data(i:i,:)*self.raw_ydata_scale+self.raw_shift*i);
                        end
                        for i = 2:length(self.derived_signals)
                             pulled = self.derived_signals{i}.ring_buff.raw(first_to_show:last_to_show,:);
                            pulled = pulled';
                           
                            set(self.ds_plot(i-1), 'YData',pulled(1,:)*self.ds_ydata_scale+self.ds_shift*(i-1)); 
                        end
                    end
                    
                    if self.samples_acquired < self.plot_size
                        xlim(self.ds_subplot, [0 self.plot_size]);
                        xlim(self.raw_subplot, [0 self.plot_size]);
                    else
                        xlim(self.ds_subplot, [0 self.plot_size]);
                        xlim(self.raw_subplot, [0 self.plot_size]);
                        set(self.raw_subplot, 'XTick', [0:self.sampling_frequency:self.plot_size]);
                        set(self.ds_subplot, 'XTick', [0:self.sampling_frequency:self.plot_size]);
                        set(self.ds_subplot, 'XTickLabel', [self.samples_acquired - ds_sp.XTick(end):ds_sp.XTick(2):self.samples_acquired]);
                        set(self.raw_subplot, 'XTickLabel', [self.samples_acquired - r_sp.XTick(end):r_sp.XTick(2):self.samples_acquired]);
                    end

                   
                end
                if(self.current_protocol> 0 && self.current_protocol<=length(self.feedback_protocols))
                    try
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
                    catch
                         
                        
                        
                    end
                else
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
                
                self.fig_feedback;
                try
                if (self.current_protocol> 0 && self.current_protocol<=length(self.feedback_protocols))
                    set(self.fb_stub,'String',self.feedback_protocols{self.current_protocol}.string_to_show);
                    
                    if isempty(get(self.fb_stub,'String'))
                        set(self.fb_stub, 'Visible', 'off');
                        set(self.fbplot_handle,'FaceColor',[1 0 0]);
                        set(self.fbplot_handle,'EdgeColor',[0 0 0]);
                    else
                        set(self.fb_stub,'Visible', 'on');
                        set(self.fbplot_handle,'FaceColor',[1 1 1]);
                        set(self.fbplot_handle,'EdgeColor','none');
                    end
                elseif self.fb_statistics_set %zero protocol after baseline recorded
                        set(self.fbplot_handle,'FaceColor',[1 0 0]);
                        set(self.fbplot_handle,'EdgeColor',[0 0 0]);
                        set(self.fb_stub,'Visible','off');
                else %zero protocol before baseline recorded
                    set(self.fbplot_handle,'FaceColor',[1 1 1],'EdgeColor','none');
                    set(self.fb_stub,'Visible','off');
                end
                catch
                    self.current_protocol
                end
                set(self.fbplot_handle,'YData',[0 self.feedback_manager.feedback_vector(self.signal_to_feedback-1) 0]);
                self.SetRecordingStatus;
            end
        end
        
        
        function WriteToFile(self)
            curr_date = datestr(date,29);
            if ~isdir(strcat(self.path,'\',curr_date))
                mkdir(strcat(self.path,'\',curr_date));
            end
            if ~isdir(strcat(self.path,'\',curr_date,'\',self.subject_record.subject_name))
                mkdir(strcat(self.path,'\',curr_date,'\',self.subject_record.subject_name));
            end
            if ~isdir(strcat(self.path,'\',curr_date,'\',self.subject_record.subject_name,'\',self.subject_record.time_start))
                mkdir(strcat(self.path,'\',curr_date,'\',self.subject_record.subject_name,'\',self.subject_record.time_start));
            end
            cd(strcat(self.path,'\',curr_date,'\',self.subject_record.subject_name,'\',self.subject_record.time_start));
            idx = 0;
            for i = 1:length(self.feedback_protocols)
                filename = strcat(num2str(i),self.feedback_protocols{i}.protocol_name,'.bin');
                string = '';
                
                fb_matrix = self.feedback_manager.feedback_records.raw(self.derived_signals{1}.collect_buff.fst+idx+1:self.derived_signals{1}.collect_buff.fst+idx + self.feedback_protocols{i}.actual_protocol_size, :);
                data_matrix = zeros(self.feedback_protocols{i}.actual_protocol_size, length(self.derived_signals)-1);
                raw_data_matrix = self.derived_signals{1}.collect_buff.raw(self.derived_signals{1}.collect_buff.fst+idx+1:self.derived_signals{1}.collect_buff.fst+idx + self.feedback_protocols{i}.actual_protocol_size, :);
                for c = 1:length(self.used_ch)
                    if c == 1
                        string = self.used_ch{c,1};
                    else
                        string = strcat(string, ',', self.used_ch{c,1});
                    end
                end
                
                for j = 2:length(self.derived_signals)
                    string = strcat(string,',',self.derived_signals{j}.signal_name);
                    data_matrix(:,j) = self.derived_signals{j}.collect_buff.raw(self.derived_signals{j}.collect_buff.fst+idx+1:self.derived_signals{j}.collect_buff.fst+idx + self.feedback_protocols{i}.actual_protocol_size, :);
                end
                
                string = strcat(string,',','Feedbacked signal', ',','Fb values',',','Average',',','Stddev');
                whole_data = [raw_data_matrix data_matrix fb_matrix];
                idx = idx+self.feedback_protocols{i}.actual_protocol_size;
                %write data
                f = fopen(filename,'w');
                fwrite(f,size(whole_data),'int');
                fwrite(f,whole_data, 'double');
                fclose(f);
                %write header
                inf_file = fopen('exp_info.hdr','w');
                fprintf(inf_file,string);
                fclose(inf_file);
            end
            self.AddNotes;
        end
        
        function SelectSignalToFeedback(self,obj,event)
            if verLessThan('matlab','8.4.0')
                self.signal_to_feedback = get(obj,'Value')+1;
            else
                self.signal_to_feedback = obj.Value+1;
            end
        end
        
        function Disconnect(self, obj,event)
            
            self.tstop = toc;
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
        function SetRecordingStatus(self)
            if verLessThan('matlab','8.4.0')
                if self.current_protocol == 0 || self.current_protocol > length(self.feedback_protocols)
                    set(self.status_text, 'String','Status: receiving');
                else
                    set(self.status_text,'String',strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.protocol_name, ': ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration)));
                end
            else
                if self.current_protocol == 0
                    self.status_text.String = 'Status: receiving';
                else
                    self.status_text.String = strcat('Status: Recording  ', self.feedback_protocols{self.current_protocol}.protocol_name, ': ',num2str(round(self.feedback_protocols{self.current_protocol}.actual_protocol_size/self.sampling_frequency)), '/',num2str(self.feedback_protocols{self.current_protocol}.protocol_duration));
                end
            end
        end
        function SetYScale(self,obj,event)
            if strcmp(get(obj,'String'),'DS scale')
                self.ds_ydata_scale =  2^get(obj,'Value');                           
                self.SetDSYTicks;
            
                
            elseif strcmp(get(obj,'String'),'Raw scale')
                self.raw_ydata_scale = 2^get(obj,'Value');
                self.SetRawYTicks;
                
            end
        end
        
        function StopRecording(self, obj,event)

            self.recording = 0;
             if self.current_protocol > length(self.feedback_protocols)
                        self.finished = 1;
                        
                        temp_log_text = get(self.log_text,'String');
                        temp_log_text{end+1} = 'Finished';
                        set(self.log_text,'String',temp_log_text);
                        set(self.disconnect_button,'String','Disconnect and write');
                        set(self.connect_button, 'String', 'Recording finished');
                        toc;
                        
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
        function StartRecording(self,obj,event)
            self.InitTimer();
            self.recording = 1;
            self.current_protocol = self.next_protocol;
            self.next_protocol = self.next_protocol + 1;
            set(self.connect_button, 'String', 'Stop recording');
            set(self.connect_button, 'Callback', @self.StopRecording);
        end
        function AddNotes(self)
            self.add_notes_window = figure;
            self.add_notes_field = uicontrol('Parent', self.add_notes_window, 'Style', 'edit', 'Position', [ 10 30 300 200]);
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
        
        function FitFigure(self,obj, event)
            f = gcbo;
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
            self.SetRawYTicks;
            self.SetDSYTicks;



        end
        function SetRawYTicks(self)
            r_sp = get(self.raw_subplot);
                self.raw_shift = (r_sp.YLim(2)-r_sp.YLim(1))/(length(self.used_ch)+1);
                r_yticks = [r_sp.YLim(1):self.raw_shift:r_sp.YLim(2)];
                set(self.raw_subplot, 'YTick', r_yticks);
                set(self.raw_subplot, 'YTickLabel', self.r_ytick_labels);
                set(self.raw_line,'String',num2str((r_sp.YLim(2)-r_sp.YLim(1))/(length(self.used_ch)+1)/self.raw_ydata_scale));
            
        end
        
        
        function SetDSYTicks(self)
            ds_sp = get(self.ds_subplot);      
                self.ds_shift = (ds_sp.YLim(2)-ds_sp.YLim(1))/(length(self.derived_signals));
                ds_yticks = [ds_sp.YLim(1):self.ds_shift:ds_sp.YLim(2)];
                set(self.ds_subplot, 'YTick', ds_yticks);
                set(self.ds_subplot, 'YTickLabel', self.ds_ytick_labels);
                set(self.ds_line,'String',num2str((ds_sp.YLim(2)-ds_sp.YLim(1))/(length(self.derived_signals))/self.ds_ydata_scale));
        end
        
    end
    
end





function channels = get_channel_labels(input) %input = inlet obj
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
function channels = read_channel_file()

fname = 'settings/mitsar_channels.txt';
t = fileread('settings/mitsar_channels.txt');
channels = {};
str = '';
k = 1;
j = 1;
while true
    if k>length(t)
        break
    end
    if strcmp(t(k), ' ')
        channels{end+1} = {t(j:k-1)};
        j = k+1;
    end
    k = k + 1;
end


end

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

%f = fopen('cm_136_5.bin','w'), fwrite(f,size(self.composite_montage),'int'); fwrite(f,self.composite_montage,'float'); fclose(f);
%f = fopen('cm_136_5.bin','r'), sz = fread(f,2,'int'); A = fread(f,sz,'float'); fclose(f);
