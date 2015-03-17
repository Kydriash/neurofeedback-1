classdef NeurofeedbackSession < handle
    properties
        derived_signals
        feedback_protocols
        protocol_sequence
    end
    
    methods
        function self = NeurofeedbackSession(self)
            self.derived_signals    = [];
            self.feedback_protocols          = {};
            self.protocol_sequence  = {};
        end
        function self = LoadFromFile(self,fname)
            
            nfs = xml2struct(fname);
            ds = nfs.NeurofeedbackSignalSpecs.vSignals.DerivedSignal;
                       
            for i = 1:length(ds)
                fields = fieldnames(ds{i});
                d = ds{i};
                for j = 1: numel(fields)
                     if strcmp(fields{j},'SpatialFilterMatrix')
                         t = xml2struct(d.SpatialFilterMatrix.Text);
                         chs = fieldnames(t.channels);
                         d.channels = cell(length(chs),2);
                         for ch = 1:numel(chs);
                             d.channels(ch,1) = {chs{ch}};
                             d.channels(ch,2) = {str2num(t.channels.(chs{ch}).Text)};
                         end
                     else
                    try
                        if str2num(d.(fields{j}).Text) || str2num(d.(fields{j}).Text) ==0
                            d.(fields{j}) = str2num(d.(fields{j}).Text);
                        end
                    catch
                        d.(fields{j}) = d.(fields{j}).Text;
                    end
                     end
                end
                d.filters(2,1) = struct();
                
                

                d.filters(1).order = 2;
                d.filters(1).mode = 'stop';
                d.filters(1).range = [45 55];
                d.filters(2).order = 2;
                d.filters(2).mode = 'pass';
                d.filters(2).range = [0.5 100];
                try
                    d.filters(3).range = [d.fBandpassLowHz d.fBandpassHighHz];
                    d.filters(3).order = 5;
                    d.filters(3).mode = 'pass';
                catch 
                end
                
                
                self.derived_signals{end+1} = d;
            end
           
            protocols = nfs.NeurofeedbackSignalSpecs.vProtocols.FeedbackProtocol;
            for i = 1:length(protocols)
                fields = fieldnames(protocols{i});
                pr = protocols{i};
                for j = 1: numel(fields)
                    
                    try
                        if str2num(pr.(fields{j}).Text) || str2num(pr.(fields{j}).Text) ==0
                            pr.(fields{j}) = str2num(pr.(fields{j}).Text);
                        end
                        
                    catch  err
                        switch err.identifier
                            case 'MATLAB:nonLogicalConditional'
                                
                                pr.(fields{j}) = pr.(fields{j}).Text;
                            case  'MATLAB:nonExistentField'
                                pr.(fields{j}) = pr.(fields{j});
                        end
                    end
                    
                    
                end
                protocols{i} = pr;
            end
            ps = nfs.NeurofeedbackSignalSpecs.vPSequence.s;
            if length(ps) == 1
                self.protocol_sequence{end+1} = ps.Text;
            else
                for i = 1:length(ps)
                    self.protocol_sequence{end+1} = ps{i}.Text;
                end
            end
            
            for j = 1: length(self.protocol_sequence)
                for i = 1:length(protocols)
                    if strcmp(self.protocol_sequence{j},protocols{i}.sProtocolName)
                        rtp = RealtimeProtocol;
                        rtp.protocol_name = protocols{i}.sProtocolName;
                        rtp.to_update_statistics = protocols{i}.bUpdateStatistics;
                        rtp.protocol_duration = protocols{i}.fDuration;
                        rtp.stop_after = protocols{i}.bStopAfter;
                        rtp.string_to_show = protocols{i}.cString;
                        try
                            rtp.window_length = protocols{i}.nMSecondsPerWindow;
                        end

                        self.feedback_protocols{end+1} = rtp;
                        rtp.set_ds_index(self.protocol_sequence);
                    end
                end
            end
        end
    end
    
    
    
    
end


