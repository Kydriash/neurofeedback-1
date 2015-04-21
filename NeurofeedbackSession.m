classdef NeurofeedbackSession < handle
    properties
        derived_signals
        feedback_protocols
        protocol_sequence
        protocol_types
    end
    
    methods
        function self = NeurofeedbackSession(self) %#ok<INUSD>
            self.derived_signals    = [];
            self.feedback_protocols          = {};
            self.protocol_sequence  = {};
            self.protocol_types = {};
        end
        function self = LoadFromFile(self,fname)
            
            nfs = xml2struct(fname);
            
            %derived signals
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
                             try
                             d.channels(ch,1) = chs(ch);
                             d.channels(ch,2) = {str2num(t.channels.(chs{ch}).Text)};
                             catch
                                 chs{ch} %#ok<NOPRT>
                             end
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
                if ~strcmp(d.sSignalName, 'Raw')
                    d.filters(1,1) = struct();
                    d.filters(1).range = [d.fBandpassLowHz d.fBandpassHighHz];
                    d.filters(1).order = 3;
                    d.filters(1).mode = 'bandpass';
                else
                    d.filters = cell(1,0);
                end
                
                
                self.derived_signals{end+1} = d;
            end
           %protocols
            self.protocol_types = nfs.NeurofeedbackSignalSpecs.vProtocols.FeedbackProtocol;
            for i = 1:length(self.protocol_types)
                fields = fieldnames(self.protocol_types{i});
                pr = self.protocol_types{i};
                for j = 1: numel(fields)
                    
                    try
                        if str2num(pr.(fields{j}).Text) || str2num(pr.(fields{j}).Text) ==0 %#ok<ST2NM>
                            pr.(fields{j}) = str2num(pr.(fields{j}).Text); %#ok<ST2NM>
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
                self.protocol_types{i} = pr;
            end
            %protocol_sequence
            ps = nfs.NeurofeedbackSignalSpecs.vPSequence.s;
            if length(ps) == 1
                self.protocol_sequence{end+1} = ps.Text;
            else
                for i = 1:length(ps)
                    self.protocol_sequence{end+1} = ps{i}.Text;
                end
            end
            
            for j = 1: length(self.protocol_sequence)
                for i = 1:length(self.protocol_types)
                    if strcmp(self.protocol_sequence{j},self.protocol_types{i}.sProtocolName)
                        rtp = RealtimeProtocol;
                        rtp.protocol_name = self.protocol_types{i}.sProtocolName;
                        rtp.to_update_statistics = self.protocol_types{i}.bUpdateStatistics;
                        rtp.protocol_duration = self.protocol_types{i}.fDuration;
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
                        rtp.set_ds_index(self.protocol_sequence);
                    end
                end
            end
        end
    end
    
    
    
    
end


