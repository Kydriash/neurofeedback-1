mulr = ReadEEGData('D:\neurofeedback\results\2015-03-26\Null\15-29-40\2Feedback.bin');
mul = mulr(:,6); % mu from the left side
fb = mulr(:,10); %feedback
fb5 = zeros(size(mulr,1),1);
fb10 = zeros(size(mulr,1),1);
fb20 = zeros(size(mulr,1),1);
av = mulr(1,11); %average
s = mulr(1,12); %std

window = 5;
step = 5;

 for i = window:step:size(mulr,1)

val = sum(abs(mul(i-window+1:i)))/window;
fb5(i-step+1:i) = (val-av)/s;
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
XLim([window size(mulr,1)-step]);
[R, P] = corrcoef(fb,fb5)
% [R, P] = corrcoef(fb,fb10)
% [R, P] = corrcoef(fb,fb20)
