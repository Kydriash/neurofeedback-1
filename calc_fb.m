function res = calc_fb(data, window,av,s)
%%%% vect, int --> vect
%%%% calculates the envelope curve of data
for i = window:window:size(data,1)
     dat = data(i-window+1:i);
val = sum(abs(dat))/window;
res(i-window+1:i) = (val-av)/s;
 end
end