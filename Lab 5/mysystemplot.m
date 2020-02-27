function mysystemplot(t,y)
% Mp
[M,MpIndex] = max(y);
yss = y(end);
Mp = M-yss;
if Mp<0
    Mp = 0;
else
    Mp = 100*Mp/yss;
end
% tr
t_10index = find(  y > .1*yss, 1, 'first');
t_90index  = find(  y > .9*yss, 1, 'first');
tr = t(t_90index)-t(t_10index);
%ts
% ts is the time it takes for the response settle between 95% and 105% of
% the steady-state value.  One way to find ts is to use a while loop,
% initialize a counter (x) to the end of the response array, and move
% forwards through the array until the response is no longer within the
% 95-105% bounds.
x = length(y); %initialize x to the end of the array
while y(x) > 0.95*yss && y(x) < 1.05*yss  %PLACE YOUR CONDITIONS HERE
    x = x-1;
end
ts = t(x)-t(1);
tssIndex = x;
% Plots
figure
plot(t,y)
title('Step Response for KcG(s)Gc(s) for ','FontSize',14,'FontWeight','bold')
%plot the response and the bounds for 10%, 90%, 95% and 105%
ss = 1:1:size(t); % the final value
ss(:) = yss;
per105=1.05*ss;
per95=.95*ss;
per10=.10*ss;
per90=.90*ss;
plot(t,y,'-',t,per10,':r',t,per90,':r',t,per105,'-g',t,per95,'-g',t,ss,'--')
% document Mp
if(Mp > 0)
    text(t(MpIndex),y(MpIndex),'\leftarrow M_p',...
        'HorizontalAlignment','left')
    line([t(MpIndex);t(MpIndex)],[0,y(MpIndex)],...
        'Color','k','LineWidth',0.5,'LineStyle',':')
end
%document tr
text(t(t_10index),y(t_10index),'\leftarrow 10%',...
    'HorizontalAlignment','left')
line([t(t_10index);t(t_10index)],[0,y(t_10index)],...
    'Color','k','LineWidth',0.5,'LineStyle',':')

text(t(t_90index),y(t_90index),'\leftarrow 90%',...
    'HorizontalAlignment','left')
line([t(t_90index);t(t_90index)],[0,y(t_90index)],...
    'Color','k','LineWidth',0.5,'LineStyle',':')
% YOU DOCUMENT tss IN THE SAME WAY AS tr AND Mp
% document tss
text(t(tssIndex),yss,'\leftarrow tss',...
    'HorizontalAlignment','left')
line([t(tssIndex);t(tssIndex)],[0,yss],...
    'Color','k','LineWidth',0.5,'LineStyle',':')
%legend
legend(['Mp = ',num2str(Mp), '%'],...
       '10% (rise time)',...
       ['90% (rise time) = ',num2str(tr), 's'],...
       '105% ',...
       '95%',...
       ['steady-state_{100%} = ',num2str(yss)],...
       ['ts = ',num2str(ts)],...
       'Location','Best')
% Label the axes:
title('Step Response','FontSize',12,'FontWeight','bold')
xlabel('t, time (sec)','FontSize',12,'FontWeight','bold')  
ylabel('y(t)','FontSize',12,'FontWeight','bold')
end