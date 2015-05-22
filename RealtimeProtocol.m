classdef RealtimeProtocol < handle
    properties
        protocol_name
        window_duration %in sec
        window_size % in samples
        to_update_statistics % (boolean) whether or not to update both global and local to protocol values of avg and tsd
        protocol_duration %sec
        protocol_size %samples
        actual_protocol_size
        stop_after
        string_to_show
        fb_type
        filter_filename
        band %range to calculate spatial filter
        n_comp
        
        
    end 
   
    methods
        function self = RealtimeProtocol(self,protocol_type,sampling_frequency)  %#ok<INUSL>
            
            if nargin == 1
                self.protocol_name           = '';
                self.to_update_statistics    = false;
                self.window_duration         = 0;
                self.protocol_duration       = 0;
                self.stop_after = false;
                self.string_to_show = '';
                self.window_size             = 0;
                self.protocol_size           = 0;
                
            elseif nargin >= 2
                self.protocol_name = protocol_type.sProtocolName;
                self.to_update_statistics    = protocol_type.bUpdateStatistics;
                self.protocol_duration =  protocol_type.fDuration;
                self.stop_after = protocol_type.bStopAfter;
                self.string_to_show = protocol_type.cString;
                try %#ok<TRYNC>
                    self.filter_filename = protocol_type.sFilterFilename;
                end
                try %#ok<TRYNC>
                    self.band = protocol_type.dBand;
                end
                try %#ok<TRYNC>
                    self.fb_type = protocol_type.sFb_type;
                end
                try %#ok<TRYNC>
                    self.window_duration = protocol_type.nMSecondsPerWindow;
                end
                try
                    self.n_comp = protocol_type.NComp;
                end
                self.window_size             = 0;
                self.protocol_size           = 0;
            end
            if nargin == 3
                self.protocol_size = self.protocol_duration*sampling_frequency;
                self.window_size = self.window_duration*sampling_frequency/1000;
            end
            
            self.actual_protocol_size = 0;
            
        end
        
        function Recalculate(self,sampling_frequency)
            self.protocol_size = self.protocol_duration*sampling_frequency;
            self.window_size = self.window_duration*sampling_frequency/1000;
        end
        
%         function set_ds_index(self, protocol_sequence)
%             for j = 1:length(protocol_sequence)
%                 %self.protocol_name
%                 if strcmp(protocol_sequence{j}, self.protocol_name)
%                     self.ds_index = [self.ds_index, j];
%                     
%                 end
%             end
%         end
        % function  set_ds_index(self,eeglsl_obj)
            % self.ds_index = [];
            % for i = 1:length(self.ds_names)
                % for j = 1:length(eeglsl_obj.derived_signals)
                    % if(strcmp(eeglsl_obj.derived_signals{j}.name, self.ds_names(i)))
                        % self.ds_index = [self.ds_index, j];
                    % end;
                % end;
            % end;
        % end;
    
	end
	
    
end

