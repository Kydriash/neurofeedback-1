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
    end
    
    methods
        
        function self = DerivedSignal(self,signal, sampling_frequency, data_length,channels,channel_count,plot_length)
            self.signal_name = signal.sSignalName;
            
            self.channels = signal.channels;
            self.channels_indices = zeros(1,length(signal.channels));
            if strcmpi(self.signal_name, 'raw')
                
                self.collect_buff = circVBuf(data_length, length(signal.channels), 0);
                self.ring_buff = circVBuf(fix(plot_length*sampling_frequency*1.1),length(signal.channels), 0); %why 128
                
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
%                 self.temporal_filter = cell(1,2);
%                 for c = 1:length(self.temporal_filter)
%                     self.temporal_filter{c}.order = signal.filters(c).order;
%                     self.temporal_filter{c}.range = signal.filters(c).range;
%                     self.temporal_filter{c}.mode = signal.filters(c).mode;
%                     [z p k] = butter(self.temporal_filter{c}.order,self.temporal_filter{c}.range/(sampling_frequency/2),self.temporal_filter{c}.mode);
%                     [self.temporal_filter{c}.B, self.temporal_filter{c}.A] = zp2tf(z,p,k);
%                     %[self.temporal_filter{c}.B, self.temporal_filter{c}.A ] = butter(self.temporal_filter{c}.order,self.temporal_filter{c}.range/(sampling_frequency/2),self.temporal_filter{c}.mode);
%                     self.temporal_filter{c}.Zf = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
%                     self.temporal_filter{c}.Zi = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
%                 end
                
            else
                self.collect_buff = circVBuf(data_length,1,0);
                self.ring_buff = circVBuf(fix(plot_length*sampling_frequency*1.1),1, 0);
                self.temporal_filter = cell(1,length(signal.filters));
                % the first temporal filter is [45 55], the second - [0.5
                % 100]
                for c = 3:length(self.temporal_filter)
                    self.temporal_filter{c}.order = signal.filters(c).order;
                    self.temporal_filter{c}.range = signal.filters(c).range;
                    self.temporal_filter{c}.mode = signal.filters(c).mode;
                    [z p k] = butter(self.temporal_filter{c}.order,self.temporal_filter{c}.range/(sampling_frequency/2),self.temporal_filter{c}.mode);
                    [self.temporal_filter{c}.B, self.temporal_filter{c}.A] = zp2tf(z,p,k);
                    %[self.temporal_filter{c}.B, self.temporal_filter{c}.A ] = butter(self.temporal_filter{c}.order,self.temporal_filter{c}.range/(sampling_frequency/2));
                    self.temporal_filter{c}.Zf = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
                    self.temporal_filter{c}.Zi = zeros(max(length(self.temporal_filter{c}.A),length(self.temporal_filter{c}.B))-1,1);
                end
                
            end
            
            
            self.spatial_filter = zeros(1,channel_count);
            for i = 1:size(channels,2)
                for j = 1:length(signal.channels)
                    if strcmp(channels{i},signal.channels{j,1})
                        self.spatial_filter(i) = signal.channels{j,2};
                    end
                end
            end
            
            
            
            self.filtered = 0;
            
        end
        
        function Apply(self, newdata,recording)
            %do projection, i.e. apply spatial filter(s)
            if strcmpi(self.signal_name, 'raw')
                self.data = zeros(length(self.channels), size(newdata,2));
                for i = 1:length(self.channels_indices)
                    sz = newdata(self.channels_indices(i):self.channels_indices(i),:);
%                     for f = 1:2 %temporal filters
%                         [ sz self.temporal_filter{f}.Zf] = filter( self.temporal_filter{f}.B,  self.temporal_filter{f}.A, sz', self.temporal_filter{f}.Zi);
%                         self.temporal_filter{f}.Zi = self.temporal_filter{f}.Zf;
%                     end
%                     
                    self.data(i:i,:) = sz;
                    
                end
                self.data = self.composite_montage*self.data;
                self.ring_buff.append(self.data');
                if recording
                    self.collect_buff.append(self.data');
                end
                
                
            else
                
                
                
                %sz  = self.spatial_filter*self.composite_montage*newdata;
                sz = newdata;
                

                if size(sz,2) > 5
                    %add selection
                    
                    for i = 1:size(sz,1)
                        %add selection
                        if self.spatial_filter(i)
                            for f = 3:length(self.temporal_filter)
                                [sz(i,:), self.temporal_filter{f}.Zf ] = filter( self.temporal_filter{f}.B,  self.temporal_filter{f}.A, sz(i,:)', self.temporal_filter{f}.Zi);
                                self.temporal_filter{f}.Zi = self.temporal_filter{f}.Zf;
                            end
                        end
                    end
                    
                    sz = self.spatial_filter*self.composite_montage * sz;
                    
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

