function [yProj, tProj] = projection(t,y)

% To plot the projection, input must be in form of tSample and ySample
%%% tSample = [t1 t2 t3 t4 t5]';
%%% ySample = [y1 y2 y3 y4 y5]';

y = smooth(y);
t = t-2;
trimIdx = length(t);
[ySample,tSample] = findpeaks(y,t); % in case data is not smooth


% use (t1, y1) through (t5, y5) to estimate K, wn and zeta
% the time domain solution of unit response of second order TF
% y(t)=K*(1-1/sqrt(1-zeta^2)*exp(-zeta*wn*t)*sin(wn*sqrt(1-zeta^2)*t+acos(zeta)))
% --- let's say it is eqn.(+)

% denote damped natural frequence as wd = wn*sqrt(1-zeta^2)
% period is the time difference between adjacent peaks, T = 2*pi/wd
Tavg = mean(diff(tSample));
wd = 2*pi/Tavg; % = wn*sqrt(1-zeta^2), so we know wn times some zeta equals a number

% use data points (t1,y1), (t2,y2), (t3,y3)
% denote M = 1/sqrt(1-zeta^2)*exp(-zeta*wn*t1)*sin(wd*t1+acos(zeta))
%        N = exp(-zeta*wn*T)
% by eqn.(+) above, we have
% y1 - y2 = K*M*N - K*M = K*M*(N-1)         --- eqn.(2)
% similarly,
% y1 - y3 = K*M*N^2 - K*M = K*M*(N-1)*(N+1) --- eqn.(3)
% divide eqn.(3) by eqn.(2), we have
% N + 1 = (y1-y3)/(y1-y2)
% so N = (y1-y3)/(y1-y2) - 1                --- eqn.(N cal.)
% note that N by this method, is independent of M, we can repeat this
% using data points (t2,y2), (t3,y3), (t4,y4)
%                   (t3,y3), (t4,y4), (t5,y5) ...
%                   (t(n),y(n)), (t(n+1),y(n+1)), (t(n+2),y(n+2))

N1 = (ySample(1)-ySample(3))/(ySample(1)-ySample(2)) - 1;
N2 = (ySample(2)-ySample(4))/(ySample(2)-ySample(3)) - 1;
N3 = (ySample(3)-ySample(5))/(ySample(3)-ySample(4)) - 1;
Navg = mean([N1 N2 N3]);
% according to the notation of N
% now we have
% zeta*wn = -log(N)/T
% also by notation of wd
% sqrt(1-zeta^2)*wn = 2*pi/T

% solve for wn and zeta by using N and T
% we need to call numerical solver fsolve for a system of nonlinear
% equations
% note: you need Optimisation Toolbox in order to use fsolve
x0 = [1.5; 3];  % init guess for zeta and wn
% call solver; check funcWnZetaVal, it should be almost 0
% parmSolve2ndOrderTF is defined below
if (strcmp(version,'7.14.0.739 (R2012a)'))
    options = optimset('Display','iter'); % for the older matlab
else
    options = optimoptions('fsolve','Display','iter'); % display iteration output
end
[x,funcWnZetaVal] = fsolve(@(x)parmSolve2ndOrderTF(x,Navg,Tavg),x0,options);


zetaEst = x(1); % compare our estimation with the acutal zeta = 0.05
wnEst = x(2); % compare with the actual wn = 10

% how about K, the steady state value - the most important parameter we
% care about? K = y1/(1-M) etc
K1 = ySample(1)/(1-MCal(zetaEst,wnEst,tSample(1))); % Mcal is defined below
K2 = ySample(2)/(1-MCal(zetaEst,wnEst,tSample(2)));
K3 = ySample(3)/(1-MCal(zetaEst,wnEst,tSample(3)));
KEst = mean([K1 K2 K3]); % compare with the actual K = 3.5

% estimation success!

%%%%%%%%%    complete the partial plot with our projection    %%%%%%%%%
% create the projected TF for from t = 4 sec onwards
sysProj = tf(KEst*wnEst^2, [1 2*zetaEst*wnEst wnEst^2]);
% projected response data
[yProj, tProj] = step(sysProj, 0:0.002:4); % UNIT response data
tProj = tProj;
% plot both collected data (blue) and projected data (red)
figure
plot(t(1:trimIdx),y(1:trimIdx),'b','linewidth',1.5)
hold on
scatter(tSample,ySample,'r')
hold on
plot(tProj,yProj,'g')
xlabel('time [s]')
ylabel('step response')
legend('collected data','Peaks','projected data','location','best')
text(9,2,'is this figure quite similar to figure 1?')

end % projection