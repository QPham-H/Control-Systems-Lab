%% Prelab 5 part b

clear
close all
clc

Ka = 2.4;
Kp = 1.6;
Kt = 0.03;
% Gain Set 1
P1 = 0.15;
P2 = 0;
a2= 0.2;
a1 = 1+180*Ka*Kt*P2;
a0 = 180*Kp*Ka*P1;
Den = conv(1,[a2,a1,a0]);
SYS = tf(a0,Den);
% Gain Set 2
P1 = 0.25;
P2 = 0.35;
a2= 0.2;
a1 = 1+180*Ka*Kt*P2;
a0 = 180*Kp*Ka*P1;
Den = conv(1,[a2,a1,a0]);
SYS2 = tf(a0,Den);
% Gain Set 3
P1 = 0.1;
P2 = 0.5;
a2= 0.2;
a1 = 1+180*Ka*Kt*P2;
a0 = 180*Kp*Ka*P1;
Den = conv(1,[a2,a1,a0]);
SYS3 = tf(a0,Den);

subplot(131)
rlocus(SYS)
title('Root Locus for P1 = 0.15, P2 = 0','FontSize',12,'FontWeight','bold')
xlabel('Real (\sigma)','FontSize',12,'FontWeight','bold')
ylabel('Imag (\omega)','FontSize',12,'FontWeight','bold')
subplot(132)
rlocus(SYS2)
title('Root Locus for P1 = 0.25, P2 = 0.35','FontSize',12,'FontWeight','bold')
xlabel('Real (\sigma)','FontSize',12,'FontWeight','bold')
ylabel('Imag (\omega)','FontSize',12,'FontWeight','bold')
subplot(133)
rlocus(SYS3)
title('Root Locus for P1 = 0.1, P2 = 0.5','FontSize',12,'FontWeight','bold')
xlabel('Real (\sigma)','FontSize',12,'FontWeight','bold')
ylabel('Imag (\omega)','FontSize',12,'FontWeight','bold')


%% Prelab 5 part d

clear
close all
clc

Ka = 2.4;
Kp = 1.6;
Kt = 0.03;
% Gain Set 1
P1 = 0.8;
P2 = 0.89;
a2= 0.2;
a1 = 1+180*Ka*Kt*P2;
a0 = 180*Kp*Ka*P1;
Den = conv(1,[a2,a1,a0]);
SYS = tf(a0,Den);

rlocus(SYS)
title('Root Locus for P1 = 0.8, P2 = 0.89','FontSize',12,'FontWeight','bold')
xlabel('Real (\sigma)','FontSize',12,'FontWeight','bold')
ylabel('Imag (\omega)','FontSize',12,'FontWeight','bold')


