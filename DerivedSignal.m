classdef DerivedSignal < handle
    
    
    properties
        
        spatial_filter %array
        temporal_filter %array
        ring_buff % fb and plot buff
        %current_filename % name of data file to collect data and calculate stdev
        window_size
        window_coefficients
        signal_name
        data
        channels
        collect_buff %collects all the data
        channels_indices
        filtered %temp for debugging
        composite_montage
        spf_times_cm
    end
    
    methods
        
        function self = DerivedSignal(self,signal, sampling_frequency, data_length,channels,channel_count,plot_length)
            self.signal_name = signal.sSignalName;
            
            self.channels = signal.channels;
            self.channels_indices = zeros(1,length(signal.channels));
            if strcmpi(self.signal_name, 'raw')         
                self.collect_buff = circVBuf(data_length, length(signal.channels), 0);
                self.ring_buff = circVBuf(fix(plot_length*sampling_frequency*1.1),length(signal.channels), 0); 
% 
                   
            else
                self.collect_buff = circVBuf(data_length*10,1,0);
                self.ring_buff = circVBuf(fix(plot_length*sampling_frequency*1.1),1, 0);

            end
            for i = 1:channel_count
                    for j = 1:length(self.channels)
                        if strcmp(channels{i},self.channels{j,1})
                            try
                            self.channels_indices(j) = i;
                            catch i,j
                            end
                        end
                    end
                end  
            self.temporal_filter = cell(1,length(signal.filters));
            for c = 1:length(self.temporal_filter)
                    self.temporal_filter{c}.order = signal.filters(c).order;
                    self.temporal_filter{c}.range = signal.filters(c).range;
                    self.temporal_filter{c}.mode = signal.filters(c).mode;
                    [z p k] = butter(self.temporal_filter{c}.order,self.temporal_filter{c}.range/(sampling_frequency/2),self.temporal_filter{c}.mode);
                    [self.temporal_filter{c}.B, self.temporal_filter{c}.A] = zp2tf(z,p,k);
                    self.temporal_filter{c}.Zf = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
                    self.temporal_filter{c}.Zi = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
            end
            
            
            
            self.filtered = 0;
            
        end
        
        function Apply(self, newdata,recording)
            %do projection, i.e. apply spatial filter(s)
            if strcmpi(self.signal_name, 'raw')
                self.data = zeros(length(self.channels), size(newdata,2));
                %select only channels we need
                self.spf_times_cm = self.spatial_filter* self.composite_montage;
                 for i = 1:length(self.channels_indices)
                     sz = newdata(self.channels_indices(i):self.channels_indices(i),:);
                    self.data(i:i,:) = sz*self.spf_times_cm(i);
                 end
                 

                self.ring_buff.append(self.data');
                if recording
                    self.collect_buff.append(self.data');
                end

            else

                sz = self.spatial_filter*self.composite_montage * newdata;
                if size(sz,2) > 5
                    %add selection
                    for i = 1:size(sz,1)
                        for f = 1:length(self.temporal_filter)
                            [sz(i,:), self.temporal_filter{f}.Zf(:,i) ] = filter(self.temporal_filter{f}.B,  self.temporal_filter{f}.A, sz(i,:)', self.temporal_filter{f}.Zi(:,i));
                            self.temporal_filter{f}.Zi(:,i) = self.temporal_filter{f}.Zf(:,i);
                        end
                    end
                    
                    
                    try
                        self.ring_buff.append(sz');
                        if recording
                            self.collect_buff.append(sz');
                        end
                    catch
                        1
                    end
                end
            end
            
        end
        
    end
end

