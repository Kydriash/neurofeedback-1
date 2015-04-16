classdef RealtimeProtocol < handle
    properties
        protocol_name
        window_duration %in sec
        window_size % in samples
        to_update_statistics % (boolean) whether or not to update both global and local to protocol values of avg and tsd
        protocol_duration
        protocol_size
        actual_protocol_size
        ds_names
        ds_index
        rbuff
        stop_after
        string_to_show
        fb_type
        filter_filename
        band %range to calculate spatial filter
        
    end 
   
    methods
        function self = RealtimeProtocol(self)
           self.protocol_name           = '';
           self.to_update_statistics    = false;
           self.window_size             = 0;
           self.window_duration         = 0;
           self.protocol_duration       = 0;
           self.protocol_size           = 0;
           self.actual_protocol_size = 0;
           self.ds_names                = {};
           self.ds_index                = [];
		   self.stop_after = false;
        end
        function set_ds_index(self, protocol_sequence)
            for j = 1:length(protocol_sequence)
                %self.protocol_name
                if strcmp(protocol_sequence{j}, self.protocol_name)
                    self.ds_index = [self.ds_index, j];
                    
                end
            end
        end
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

