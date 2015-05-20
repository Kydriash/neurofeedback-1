function data2 = rearrange_channels(data, chan_labels, chan_locs)
%% matrix data(channels, t),cell chan_labels(channels),struct chan_locs(channels) -->
%% matrix data2(channels2, t)
%% matches the rows in the data to chan_locs labels  using chan_labels 
%% obtained from the stream
%% if the channel was not presented in the stream, it is represented with zeros
data2 = zeros(length(chan_locs), size(data,2));
for j = 1:length(chan_locs)
    index = find(strcmpi(chan_labels,chan_locs(j).labels));
    if ~isempty(index)
        data2(j,:) = data(index,:);
    end
    
end
end


