mulr = ReadEEGData('D:\neurofeedback\results\2015-03-27\Null\13-49-44\2Feedback.bin');
mul = mulr(:,6); % mu from the left side
fb = mulr(:,10); %feedback
fb5 = zeros(size(mulr,1),1);
fb10 = zeros(size(mulr,1),1);
fb20 = zeros(size(mulr,1),1);
av = mulr(50,11); %average
s = mulr(50,12); %std

window = 15;
step = 5;
 for i = window:step:size(mulr,1)-step
val = sum(abs(mul(i-window+1:i)))/window;
fb5(i-step+1:i) = (val-av)/s;
 end
%  step = mulr(window,13);
%  i = window;
%  while i < size(mulr,1)
%      i = i + 5;
%      if step ~= fix(mulr(i,13))
%           val = sum(abs(mul(i-window+1:i)))/window;
%          fb5(i-window+1:i) = (val-av)/s;
%          step = fix(mulr(i,13));
%      else
%           fb5(i-window+1:i) = fb5(i-2*window+1);
%      end
%  end
%  for i = window+5:5:size(mulr,1)-5
%      
%     if step ~= fix(mulr(i+5,13))
%         val = sum(abs(mul(i-window+1:i)))/window;
%         fb5(i-5+1:i) = (val-av)/s;
%         step = fix(mulr(i+5,13));
%     else
%         
%         fb5(i-5+1:i) = fb5(i-10+1);
%         
%     end
%      
% %     step = fix(mulr(i-step+1,13));
% %     i = i + step;
% %     [step i]
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
grid on;
% hold on;
% plot(fb10, 'k-');
% hold on;
% plot(fb20, 'g-');
XLim([window size(mulr,1)-step]);
[R, P] = corrcoef(fb,fb5)
% [R, P] = corrcoef(fb,fb10)
% [R, P] = corrcoef(fb,fb20)
