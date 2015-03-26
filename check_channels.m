function  check_channels()
lsllib = lsl_loadlib();
streams = [];
channel_count = 0;


while isempty(streams)
    streams = lsl_resolve_byprop(lsllib,'type', 'Data');
end
channel_count = streams{1}.channel_count();
inlet = lsl_inlet(streams{1});
try
    channel_labels = get_channel_labels(inlet);
end
if channel_count ~= length(channel_labels)
    channel_count
    length(channel_labels)
end
%'finished'
end


function channels = get_channel_labels(input) %input = inlet obj
ChS =  input.info.desc.child('channels');
ch = ChS.first_child;
channels = {};
try
    
    % while ch.next_sibling.PtrHandle
    while ch.PtrHandle
        l = ch.child('label');
        channels{end+1} = l.child_value ;
        ch = ch.next_sibling;
    end
catch
    channels = cell(1,input.channel_count());
    for i = 1:input.channel_count()
        channels{i} = num2str(i);
    end
end
end