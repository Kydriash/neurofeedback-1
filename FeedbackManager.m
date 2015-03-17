classdef FeedbackManager < handle
    properties
        standard_deviation
        average
        sum_of_squares
        sum_of_values
        window_length %ms
        window_size %samples
        current_protocol
        %samples_acquired
        eeglsl_obj
        feedback_vector
        neurofeedback_session
        feedback_records
    end
    
    methods
        function self = FeedbackManager(self)
            self.window_size            = 10;
            self.reset_stats();
             self.standard_deviation     = [1];
             self.average                = [0];
%             self.sum_of_squares         = 0;
%             self.sum_of_values          = 0;
            %self.samples_acquired       = 0;
            self.current_protocol       = {};
            self.feedback_vector        = [];
            self.neurofeedback_session  = [];
            
        end
        
        function reset_stats(self)
            self.standard_deviation = 1;%ones(length(self.current_protocol.ds_index),1);
            self.average            = 0;%zeros(length(self.current_protocol.ds_index),1);
            self.sum_of_squares     = 0;%zeros(length(self.current_protocol.ds_index),1);
            self.sum_of_values      = 0;%zeros(length(self.current_protocol.ds_index),1);
            %self.samples_acquired   = 0;
        end;
        
        function set_eeglsl(self, eeglsl_obj_in)
            self.eeglsl_obj = eeglsl_obj_in;
        end
        
%         function set_protocol(self, curr_protocol)
%             self.current_protocol = curr_protocol;
%             
%             self.standard_deviation = ones(length(self.current_protocol.ds_index),1);
%             self.average            = zeros(length(self.current_protocol.ds_index),1);
%             self.sum_of_squares     = zeros(length(self.current_protocol.ds_index),1);
%             self.sum_of_values      = zeros(length(self.current_protocol.ds_index),1);
%             self.samples_acquired   = 0;
%             self.window_size = self.current_protocol.window_size;
%         end;
        
        
        
%         function update_statistics(self)
%             if(self.samples_acquired>0)
%                 for i = self.current_protocol.ds_index
%                     self.average(i)            = self.sum_of_values(i)/self.samples_acquired;
%                     self.standard_deviation(i) = sqrt(self.sum_of_squares(i)/self.samples_acquired - self.average(i)*self.average(i));
%                     self.samples_acquired    = 0;
%                     self.sum_of_values(i)      = 0;
%                     self.sum_of_squares(i)     = 0;
%                 end;
%             end;
%         end
        
        function update_feedback_signal(self)
            n = self.window_size;
            for i = self.current_protocol.ds_index
                % extract the latest chunk
                if(strcmp(self.eeglsl_obj.derived_signals{i}.type,'plain'))
                    dat = self.eeglsl_obj.derived_signals{i}.ring_buff.raw(self.eeglsl_obj.derived_signals{i}.ring_buff.lst-n+1:...
                        self.eeglsl_obj.derived_signals{i}.ring_buff.lst,:);
                    % calculate z-score like statistics
                    avg  = self.average(i);
                    sdev = self.standard_deviation(i);
                    self.feedback_vector(i)  = (sqrt(sum(dat.^2))/n-avg)/sdev;
                    % update accumulators
                    if(~isempty(self.current_protocol))
                        self.sum_of_values(i)   =   self.sum_of_values(i)   +   self.feedback_vector(i);
                        self.sum_of_squares(i)  =   self.sum_of_squares(i)  +   self.feedback_vector(i);
%                         self.samples_acquired   =   self.samples_acquired   +   1;
                        %self.samples_acquired=self.eeglsl_obj.derived_signals{i}.ring_buff.lst;
                    end;
                elseif (strcmp(self.eeglsl_obj.derived_signals{i}.type,'composite'))
                    %not implemented yet
                    switch self.eeglsl_obj.derived_signals{i}.operation
                        case 'difference'
                            
                            %find d_signals operator1 and operator 2
                            for j = 1:length(self.eeglsl_obj.derived_signals)
                                if strcmp(self.eeglsl_obj.derived_signals{i}.sComponent1, self.eeglsl_obj.derived_signals{j}.name)
                                    component1 = self.eeglsl_obj.derived_signals{j};
                                elseif strcmp(self.eeglsl_obj.derived_signals{i}.sComponent2, self.eeglsl_obj.derived_signals{j}.name)
                                    component2 = self.eeglsl_obj.derived_signals{j};
                                end
                            end
                            dat = component1.ring_buff.raw(component1.ring_buff.lst-n+1:...
                                component1.ring_buff.lst,:) - component2.ring_buff.raw(component2.ring_buff.lst-n+1:...
                                component2.ring_buff.lst,:);
                            
                            
                            avg = self.average(i);
                            sdev = self.standard_deviation(i);
                            self.feedback_vector(i) = sqrt(sum(dat.^2)/n-avg)/sdev;
                            %update accumulators
                            if(~isempty(self.current_protocol))
                                self.sum_of_values(i)   =   self.sum_of_values(i)   +   self.feedback_vector(i);
                                self.sum_of_squares(i)  =   self.sum_of_squares(i)  +   self.feedback_vector(i);
                                %self.samples_acquired   =   self.samples_acquired   +   1;
                                %self.samples_acquired=self.eeglsl_obj.derived_signals{i}.ring_buff.lst;
                            end;
                            
                            
                    end
                    
                    
                end;
            end;
        end
    end
    
end