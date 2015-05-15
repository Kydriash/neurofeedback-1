
function RecordedStudyToLSL(fnames,pathname,looped)
global pushed;
%read the files

%check if fnames ets exist
    
    [protocols, durations, channels]  = GetDataProperties(pathname,fnames); %#ok<ASGLU>
    
    %create lsl
    sampling_frequency = 500;
    source_id = pathname;
    lsllib = lsl_loadlib();
    
    eeg_info = lsl_streaminfo(lsllib,'File', 'Data',length(channels), sampling_frequency,'cf_float32',source_id);
    chns = eeg_info.desc().append_child('channels');
    
    for label = channels
        ch = chns.append_child('channel');
        ch.append_child_value('label',label{1});
    end
    outlet = lsl_outlet(eeg_info);
    filenames ={};
    if ischar(fnames)
        filenames{end+1} = strcat(pathname,fnames);
    else
    for f = fnames
        filenames{end+1} = strcat(pathname,f{1}); 
    end
    end
    data = [];
    for fn = filenames
        protocol_data = ReadEEGData(fn{1});
        protocol_data = protocol_data(:,1:length(channels))';
        data = [data protocol_data];
    end
    
    
    pushed = 1;
%     while ~outlet.have_consumers()
%         pause(0.01);
%     end
    timer_push_data = timer('Name','push_data','TimerFcn', {@PushDataToLSL,outlet,data,looped},'ExecutionMode','fixedRate','Period',1/sampling_frequency);

    start(timer_push_data);

%     for fn = filenames
%         protocol_data = ReadEEGData(fn{1});
%         protocol_data = protocol_data(:,1:length(channels))';
%         pushed = 1;
%         
%         while outlet.have_consumers() && pushed <= size(protocol_data,2)
%             try
%                 
%                 outlet.push_chunk(protocol_data(:,pushed));
%             catch
%                 pushed, size(protocol_data,2) %#ok<NOPRT>
%             end
%             pushed = pushed + 1;
%             pause(0.0001);
%             if pushed == size(protocol_data,2)
%                 'Protocol finished. Sent %d samples', pushed
%             end
%         end
%         
%         
%     end
end
    



function PushDataToLSL(timer_obj,event,outlet, data,looped) %#ok<INUSL>
global pushed
if looped
    if outlet.have_consumers()
        try %#ok<TRYNC>
        outlet.push_chunk(data(:,mod(pushed,size(data,2))))
%         catch 
%             mod(pushed,size(data,2))
        end
        pushed = pushed + 1;
       
    end
else
    if outlet.have_consumers() && pushed <= size(data,2)
        outlet.push_chunk(data(:,pushed));
        pushed = pushed + 1;
        
    end
end
end


