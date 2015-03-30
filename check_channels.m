function  check_channels()
lsllib = lsl_loadlib();
streams = [];
channel_count = 0;


while isempty(streams)
    streams = lsl_resolve_byprop(lsllib,'type', 'Data');
end
channel_count = streams{1}.channel_count();
'trying get the inlet'
inlet = lsl_inlet(streams{1});
try
    channel_labels = get_channel_labels(inlet,channel_count);
end
if channel_count ~= length(channel_labels)
    channel_count
    length(channel_labels)
else
    %channel_labels
end
%'finished'
end


function channels = get_channel_labels(input,channel_count) %input = inlet obj
ChS =  input.info.desc.child('channels');
ch = ChS.first_child;
channels = {};
try
    
    % while ch.next_sibling.PtrHandle
    
    'trying retrieve the channels'
    for i = 1:channel_count
        ch.PtrHandle

%     while ch.PtrHandle

        
        if ch.child('label').child_value
            channels{end+1} = ch.child('label').child_value ;
        else
            channels{end+1} = 'NA'; %if ch.PtrHandle == 0
        end
        ch = ch.next_sibling;
    end
catch
    channels = cell(1,channel_count());
    for i = 1:channel_count()
        channels{i} = num2str(i);
    end
end
end