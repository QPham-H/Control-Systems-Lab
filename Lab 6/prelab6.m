%% Prelab6 part b

clc
clear
close all

K = 18;
Kc = 1;
Kpot = 5/(pi);
Kamp = 2.4;
tm = 0.2;

KT = K*Kc*Kpot*Kamp;
Den = conv([1,0],[tm,1]);
H1 = tf([0,0,KT],Den);

Kc = 8.9;
KT = K*Kc*Kpot*Kamp;
Den = conv([1,0],[tm,1]);
H2 = tf([0,0,KT],Den);

Kc = 2;
KT = K*Kc*Kpot*Kamp;
Den = conv([1,0],[tm,1]);
H3 = tf([0,0,KT],Den);

subplot(131)
margin(H1)
legend('K = 1',...
       'Location','Best')
subplot(132)
margin(H2)
legend('K = 8.9',...
       'Location','Best')
subplot(133)
margin(H3)
legend('K = 2',...
       'Location','Best')

%% Prelab6 part c

clc
clear
close all

K = 10;
NUM = conv(10,[1,10]);
H = tf(NUM,[1,100]);

margin(H)
%% Prelab6 part d
clc
clear
close all

% Given
K = 18;
Kpot = 5/(pi);
Kamp = 2.4;
tm = 0.2;
%Gc Design
Wn = 55.4;
Wn2 = Wn*Wn;
Zeta = 0.5169;
H = tf(Wn2,[1,2*Zeta*Wn,0]);
[GM,PM,Wgm,Wpm] = margin(H);
PM = deg2rad(PM);
safety = deg2rad(10);
PMmotor = deg2rad(10.9);
alpha = (1-sin(PM+safety-PMmotor))/(1+sin(PM+safety-PMmotor)); % add 10 for safety
T = 1/(Wpm*sqrt(alpha));
% Realization of G
Kc = 2;
KT = Kc*K*Kpot*Kamp;
Den = conv([1,0],[tm,1]);
G = tf(KT,Den);
% Realization of Gc
Gc = tf([T,1],[T*alpha,1]);
% Realization of Close Loop System
HCL = tf(Gc*G/(1+Gc*G));
mysystemplot(HCL);



%% Prelab6 part e
clear all
close all
clc
% G
Kc = 2;
KT = Kc*Kpot*Kamp;
Den = conv([1,0],[tm,1]);
G = tf(KT,Den);
% Gc
num = conv(0.4,[1/2.3,1]);
Gc = tf(num,[1/92,1]);
figure(1)
margin(Gc)
% Realization of Close Loop System
HCL = tf(Gc*G/(1+Gc*G));
[y,t] = step(HCL);
mysystemplot(HCL)






