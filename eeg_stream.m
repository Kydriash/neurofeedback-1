function eeg_stream()
global connected
global pushed
global is_transmitting
global channel_count
[fnames, pathname, filterindex] = uigetfile('.bin','Select files to play','MultiSelect','on');  %#ok<ASGLU>
if ischar(fnames) || iscell(fnames) %fnames ~= 0
    f = figure('Position', [400, 400, 150,80],'Tag','stream control');
    l = uicontrol(f,'String','loop', 'Position', [10 40 100 20],'Style','checkbox','Tag','loop data'); %#ok<NASGU>
    p = uicontrol(f,'String','Start streaming','Position', [10 10 100 20],'Callback',@ct,'Tag','connect button'); %#ok<NASGU>
    
    
    
    %read the files
    %check if fnames ets exist
    [protocols,protocols_show_as, durations, channels]  = GetDataProperties(pathname,fnames); %#ok<ASGLU>
    
    channel_count = length(channels);
    
    %get data and calculate sampling_frequency
    filenames ={};
    if ischar(fnames)
        filenames{end+1} = strcat(pathname,fnames);
    else
        for f = fnames
            filenames{end+1} = strcat(pathname,f{1});
        end
    end
    data_length = 0;
    duration = sum(durations);
    data = [];
    for fn = filenames
        protocol_data = ReadEEGData(fn{1});
        %size(protocol_data)
        if size(protocol_data,1)
            temp_protocol_data = protocol_data(:,1:length(channels))';
            data = [data temp_protocol_data];
            
        end
        data_length = data_length + size(protocol_data,1);
    end
    
    %create lsl
    sampling_frequency = round(data_length/duration);
    
    source_id = pathname;
    lsllib = lsl_loadlib();
    eeg_info = lsl_streaminfo(lsllib,'NVX136_Data', 'Data',length(channels), sampling_frequency,'cf_float32',source_id);
    chns = eeg_info.desc().append_child('channels');
    
    for label = channels
        ch = chns.append_child('channel');
        ch.append_child_value('label',label{1});
    end
    outlet = lsl_outlet(eeg_info);
    
    pushed = 1;
    timer_push_data = timer('Name','push_data','TimerFcn', {@PushDataToLSL,outlet,data},'ExecutionMode','fixedRate','Period',1/sampling_frequency);
    is_transmitting = 0;
    connected = 0;
    start(timer_push_data);
end
end

function ct(obj,event) %#ok<INUSD>
global connected
global looped
b = findobj('Tag','connect button');
l = findobj('Tag','loop data');
if ~connected
    b.String = 'Stop streaming';
    connected = 1;
    looped = l.Value;
else
    connected = 0;
    b.String = 'Start streaming';
    looped = l.Value;
    
end
end

function PushDataToLSL(timer_obj,event,outlet, data) %#ok<INUSL>
global pushed
global looped
global connected
global channel_count
b = findobj('Tag','connect button');
%eeg_figure = findobj('Tag','raw_and_ds_figure');
if ~connected
    
    outlet.push_chunk(zeros(channel_count,1));
    if pushed > 1 && ~outlet.have_consumers()
        disp('The eeg window has been closed, closing transmission')
        delete(outlet);
        stop(timer_obj);
        f = findobj('Tag','stream control');
        delete(f);
    end
    
else
    
    if looped
        
        
        if outlet.have_consumers()
            outlet.push_chunk(data(:,mod(pushed,size(data,2))))
            pushed = pushed + 1;
        elseif pushed > 1 && ~outlet.have_consumers()
            disp('The eeg window has been closed, closing transmission')
            delete(outlet);
            stop(timer_obj);
            f = findobj('Tag','stream control');
            delete(f);
        elseif connected &&  ~outlet.have_consumers()
            disp('Waiting for consumers');
            connected = 0;
            b.String = 'Start streaming';
            
        end
        
    else
        if outlet.have_consumers() && pushed <= size(data,2) %&& ~isempty(eeg_figure)
            outlet.push_chunk(data(:,pushed));
            pushed = pushed + 1;
        elseif pushed > size(data,2)
            disp('Streaming finished')
            delete(outlet);
            stop(timer_obj);
            f = findobj('Tag','stream control');
            delete(f);
        elseif pushed > 1 && ~outlet.have_consumers()
            disp('The eeg window has been closed, closing transmission')
            delete(outlet);
            stop(timer_obj);
            f = findobj('Tag','stream control');
            delete(f);
        elseif connected && ~outlet.have_consumers()
            %delete(outlet);
            %stop(timer_obj);
            %f = findobj('Tag','stream control');
            %delete(f);
            disp('Waiting for consumers');
            connected = 0;
            b.String = 'Start streaming';
            
        end
    end
end
end




