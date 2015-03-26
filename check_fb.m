mulr = ReadEEGData('D:\neurofeedback\results\2015-03-26\Null\12-58-16\2Feedback.bin');
mul = mulr(:,6); % mu from the left side
fb = mulr(:,10); %feedback
fb5 = zeros(5005,1);
fb10 = zeros(5005,1);
fb20 = zeros(5005,1);
av = mulr(1,11); %average
s = mulr(1,12); %std
 for i = 5:5:5000
val = sum(abs(mul(i-4:i)))/5;
fb5(i-4:i) = (val-av)/s;
 end
%  for i = 5:1:5000
% val = sum(abs(mul(i-4:i)))/5;
% fb5(i-4:i) = (val-av)/s;
% end
%  for i = 10:10:5000
% val = sum(abs(mul(i-9:i)))/10;
% fb10(i-9:i) = (val-av)/s;
% end
% for i = 20:20:5000
% val = sum(abs(mul(i-19:i)))/20;
% fb20(i-19:i) = (val-av)/s;
% end
figure;
plot(fb);
hold on;
plot(fb5,'r-');
% hold on;
% plot(fb10, 'k-');
% hold on;
% plot(fb20, 'g-');
XLim([0 5000]);
[R, P] = corrcoef(fb,fb5)
% [R, P] = corrcoef(fb,fb10)
% [R, P] = corrcoef(fb,fb20)
