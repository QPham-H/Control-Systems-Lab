% motor par
tau = 0.2;
Kpot = 5/pi;
Kamp = 2.5;
K = 18;
Kmot = Kamp*Kpot*K;
% controler par
Kc = .4;
tz = 1/2.3;
tp = 1/92;
Z = 1/tz;
P = 1/tp;
Ka = Kc*tz/tp;


sim('DcMotorLeadController')
V50sat = V;
t50sat = tout;
%%
clc
close
clear
load('simulationsPart1_Sat10.mat')
load('simulationsPart1_SatInf.mat')

figure
plot(t05,V05,t05sat,V05sat)
title(' V05 and V05sat')
xlabel('time (sec)')
ylabel('V_{out}(t)')
legend('V05',...
       'V05sat',...
       'Location','Best')
figure
plot(t5,V5,t5sat,V5sat)
title(' V5 and V5sat')
xlabel('time (sec)')
ylabel('V_{out}(t)')
legend('V5',...
       'V5sat',...
       'Location','Best')
figure
plot(t50,V50,t50sat,V50sat)
title(' V50and V50sat')
xlabel('time (sec)')
ylabel('V_{out}(t)')
legend('V50',...
       'V50sat',...
       'Location','Best')
%% Part II

% controler par
Kc = .4;
tz = 1/2.3;
tp = 1/92;
Z = 1/tz;
P = 1/tp;
Ka = Kc*tz/tp;
set_param('rtwin486', 'SimulationCommand', 'start')

%% Part IV

% controler par
Kc = 2;
tz = 1/13;
tp = 1/120;
Z = 1/tz;
P = 1/tp;
Ka = Kc*tz/tp;
set_param('rtwin486PartIV', 'SimulationCommand', 'start')






