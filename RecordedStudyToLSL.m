
function RecordedStudyToLSL(fnames,pathname,looped,s_frequency)
global pushed
global is_transmitting
global connected

%read the files

%check if fnames ets exist

[protocols,protocols_show_as, durations, channels]  = GetDataProperties(pathname,fnames); %#ok<ASGLU>


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
if nargin > 3 && s_frequency > 0
    sampling_frequency = s_frequency;
else
    sampling_frequency = round(data_length/duration);
    %sampling_frequency = 500;
end
%sampling_frequency = 500;
source_id = pathname;
lsllib = lsl_loadlib();

eeg_info = lsl_streaminfo(lsllib,'File', 'Data',length(channels), sampling_frequency,'cf_float32',source_id);
chns = eeg_info.desc().append_child('channels');

for label = channels
    ch = chns.append_child('channel');
    ch.append_child_value('label',label{1});
end
outlet = lsl_outlet(eeg_info);




pushed = 1;
%     while ~outlet.have_consumers()
%         pause(0.01);
%     end
timer_push_data = timer('Name','push_data','TimerFcn', {@PushDataToLSL,outlet,data,looped},'ExecutionMode','fixedRate','Period',1/sampling_frequency);
is_transmitting = 0;
connected = 0;
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
%global is_transmitting
global connected
if looped
    if outlet.have_consumers()
        %is_transmitting = 1;
        connected = 1;
        try %#ok<TRYNC>
            outlet.push_chunk(data(:,mod(pushed,size(data,2))))
            %         catch
            %             mod(pushed,size(data,2))
        end
        pushed = pushed + 1;
        %     elseif is_transmitting && ~(outlet.have_consumers())
        %         outlet.wait_for_consumers();
        %         is_transmitting = 0;
    else
        eeg_figure = findobj('Tag','raw_and_ds_figure');
        if connected &&  isempty(eeg_figure)
            delete(outlet);
            stop(timer_obj);
        end
    end
    
else
    if outlet.have_consumers() && pushed <= size(data,2)
        connected = 1;
        % is_transmitting = 1;
        outlet.push_chunk(data(:,pushed));
        pushed = pushed + 1;
        %     elseif is_transmitting && ~(outlet.have_consumers())
        %         outlet.wait_for_consumers();
        %         is_transmitting = 0;
    else
        eeg_figure = findobj('Tag','raw_and_ds_figure');
        if connected && isempty(eeg_figure)
            delete(outlet);
            stop(timer_obj);
        end
    end
    
    
end


end





