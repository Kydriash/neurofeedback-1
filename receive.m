function receive()
lsllib = lsl_loadlib();
streams = [];
while isempty(streams)
                streams = lsl_resolve_byprop(lsllib,'type', 'Data');
end
self.inlet = lsl_inlet(streams{1});
   
end