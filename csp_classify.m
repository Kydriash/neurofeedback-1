function class = csp_classify(data, w1,w2)
%%%% returns 1 if the data belongs to the 1st class and 2 if to the 2nd 
DLP = log(w1'*data*data'*w1) - log(w2'*data*data'*w2);  %#ok<MHERM>
DLP %#ok<NOPRT>
class = (DLP < 0)+1;
end