classdef NeurofeedbackSession < handle
    properties
        derived_signals
        feedback_protocols
        protocol_sequence
        protocol_types
        csp_settings
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
            [folder, fn, ext] = fileparts(fname); %#ok<ASGLU>
            %% derived signals
            ds = nfs.NeurofeedbackSignalSpecs.vSignals.DerivedSignal;
            for i = 1:length(ds)
                if isstruct(ds)
                    fields = fieldnames(ds(i));
                    d = ds(i);
                elseif iscell(ds)
                    fields = fieldnames(ds{i});
                    d = ds{i};
                end
                
                
                for j = 1: numel(fields)
                    if strcmp(fields{j},'SpatialFilterMatrix')
                        if ~isempty(d.SpatialFilterMatrix.Text)
                            [directory, filename, extension] = fileparts(d.SpatialFilterMatrix.Text); %#ok<ASGLU>
                            if isempty(directory)
                                t = xml2struct(strcat(folder,'\',d.SpatialFilterMatrix.Text));
                            else
                                t = xml2struct(d.SpatialFilterMatrix.Text);
                            end
                            
                            chs = fieldnames(t.channels);
                            d.channels = cell(length(chs),2);
                            for ch = 1:numel(chs);
                                try
                                    d.channels(ch,1) = chs(ch);
                                    coeffs = str2num(t.channels.(chs{ch}).Text); %#ok<ST2NM>
                                    for j = 1:length(coeffs)
                                        
                                        d.channels{ch,j+1} = coeffs(j);
                                    end
                                catch
                                    chs{ch} %#ok<NOPRT>
                                end
                            end
                        else
                            d.channels = cell(0,0);
                        end
                    else
                        try
                            if str2num(d.(fields{j}).Text) || str2num(d.(fields{j}).Text) ==0 %#ok<ST2NM>
                                d.(fields{j}) = str2num(d.(fields{j}).Text); %#ok<ST2NM>
                            end
                        catch
                            d.(fields{j}) = d.(fields{j}).Text;
                        end
                    end
                end
                if ~strcmpi(d.sSignalName, 'Raw')
                    d.filters(1,1) = struct();
                    d.filters(1).range = [d.fBandpassLowHz d.fBandpassHighHz];
                    d.filters(1).order = 3;
                    d.filters(1).mode = 'bandpass';
                else
                    d.filters = cell(1,0);
                end
                
                
                self.derived_signals{end+1} = d;
            end
            %% protocols
            self.protocol_types = nfs.NeurofeedbackSignalSpecs.vProtocols.FeedbackProtocol;
            for i = 1:length(self.protocol_types)
                fields = fieldnames(self.protocol_types{i});
                pr = self.protocol_types{i};
                for j = 1: numel(fields)
                    
                    try
                        % upd on 2015-06-02
                        if any(str2num(pr.(fields{j}).Text)) || all(str2num(pr.(fields{j}).Text)) ==0 %#ok<ST2NM>
                            pr.(fields{j}) = str2num(pr.(fields{j}).Text); %#ok<ST2NM>
                        else
                            pr.(fields{j}) = pr.(fields{j}).Text;
                        end
                        
                    catch  err
                        switch err.identifier
                            %                             case 'MATLAB:nonLogicalConditional'
                            %
                            %                                 pr.(fields{j}) = pr.(fields{j}).Text;
                            case  'MATLAB:nonExistentField'
                                pr.(fields{j}) = pr.(fields{j});
                        end
                    end
                    
                    
                end
                self.protocol_types{i} = pr;
            end
            %% protocol_sequence
            %%% upd on 2015-05-13
            show_as = {};
            try
                seq = nfs.NeurofeedbackSignalSpecs.vPSequence.s;
                %ps = {};
                if length(seq) == 1
                    self.protocol_sequence{end+1} = seq.Text;
                else
                    for s = 1:length(seq)
                        self.protocol_sequence{end+1} = seq{s}.Text;
                    end
                end
                
            catch err
                if strcmp(err.identifier, 'MATLAB:nonExistentField')
                    seq = nfs.NeurofeedbackSignalSpecs.vPSequence.loop;
                    %ps = {};
                    for ss = 1:length(seq)
                        for a = 1:str2double(seq{ss}.Attributes.count)
                            for p = 1:length(seq{ss}.s)
                                if length(seq{ss}.s) == 1
                                    self.protocol_sequence{end+1} = seq{ss}.s(p).Text;
                                else
                                    self.protocol_sequence{end+1} = seq{ss}.s{p}.Text;
                                end
                            end
                        end
                    end
                end
            end
            
            
            
            %%%
            
            
            for j = 1: length(self.protocol_sequence)
                for i = 1:length(self.protocol_types)
                    if strcmp(self.protocol_sequence{j},self.protocol_types{i}.sProtocolName)
                        rtp = RealtimeProtocol(1,self.protocol_types{i});
                        self.feedback_protocols{end+1} = rtp;
                        %rtp.set_ds_index(self.protocol_sequence);
                    end
                end
            end
            try
                csp = nfs.NeurofeedbackSignalSpecs.CSP_settings;
                
                fields = fieldnames(csp);
                if length(fields) == 1 && strcmp(fields,'Text')
                    csp = [];
                else
                    
                    for i = 1:numel(fields)
                        try
                            % upd on 2015-06-02
                            if any(str2num(csp.(fields{i}).Text)) || all(str2num(csp.(fields{i}).Text)) ==0 %#ok<ST2NM>
                                csp.(fields{i}) = str2num(csp.(fields{i}).Text); %#ok<ST2NM>
                            else
                                csp.(fields{i}) = csp.(fields{i}).Text;
                            end
                            
                        catch  err
                            switch err.identifier
                                %                             case 'MATLAB:nonLogicalConditional'
                                %
                                %                                 pr.(fields{j}) = pr.(fields{j}).Text;
                                case  'MATLAB:nonExistentField'
                                    csp.(fields{i}) = csp.(fields{i});
                            end
                        end
                    end
                    
                end
            catch err
                switch err.identifier
                    case 'MATLAB:nonExistentField'
                        csp = [];
                end
                
            end
            
            self.csp_settings = csp;
        end
    end
    
    
    
    
end


