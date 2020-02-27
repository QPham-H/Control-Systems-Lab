clear all
close all
clc


Ktach = 0.03;
Kamp = 2.4;
Kpot = 1.6;

%% Part 1
K1 = 1.5;
K2 = 0;
set_param('rtwin486', 'SimulationCommand', 'start')
%sim('rtwin486');
V01 = Scope(:,3);
Vr1 = Scope(:,2);
plot(Scope(:,1),V01)
hold on
plot(Scope(:,1),Vr1)
%% Part 2
K1 = 2.5;
K2 = 3.5;
set_param('rtwin486', 'SimulationCommand', 'start')
%sim('rtwin486');
V02 = Scope(:,3);
Vr2 = Scope(:,2);
plot(Scope(:,1),V02)
hold on
plot(Scope(:,1),Vr2)
%% Part 3
K1 = 1;
K2 = 5;
set_param('rtwin486', 'SimulationCommand', 'start')
%sim('rtwin486');
V03 = Scope(:,3);
Vr3 = Scope(:,2);
plot(Scope(:,1),V03)
hold on
plot(Scope(:,1),Vr3)
%% Part 4
K1 = 6.5;
K2 = 6;
set_param('rtwin486', 'SimulationCommand', 'start')
%sim('rtwin486');
V04 = Scope(:,3);
Vr4 = Scope(:,2);
plot(Scope(:,1),V04)
hold on
plot(Scope(:,1),Vr4)



%% Part III Part 1
K1 = 6.5;
K2 = 6;
%set_param('lab5_friccomp', 'SimulationCommand', 'start')
%F01 = Scope(:,3);
%Fr1 = Scope(:,2);
figure
plot(Scope(:,1),F01)
hold on
plot(Scope(:,1),Fr1)
%% Part III Part 2
K1 = 2.5;
K2 = 3.5;
%set_param('lab5_friccomp', 'SimulationCommand', 'start')
%F02 = Scope(:,3);
%Fr2 = Scope(:,2);
figure
plot(Scope(:,1),F02)
hold on
plot(Scope(:,1),Fr2)
%% Part III Part 3
K1 = 1;
K2 = 5;
%set_param('lab5_friccomp', 'SimulationCommand', 'start')
%F03 = Scope(:,3);
%Fr3 = Scope(:,2);
figure
plot(Scope(:,1),F03)
hold on
plot(Scope(:,1),Fr3)
%% Part III Part 4
K1 = 6.5;
K2 = 6;
%set_param('lab5_friccomp', 'SimulationCommand', 'start')
%F04 = Scope(:,3);
%Fr4 = Scope(:,2);
figure
plot(Scope(:,1),F04)
hold on
plot(Scope(:,1),Fr4)
%% ParI Analog Computer
% obtaining Mp, tr, and ts for the different sets of ps.

clear all
close all
clc

run('/Users/hernandezmarvinalexander/Desktop/Fall2017/ece486/Labs/lab5/lab5Data.m')

myt1 = y1(:,1);
myy1 = y1(:,2);
myy1 = myy1 - myy1(1);
mysystemplot(myt1,myy1)

myt2 = y2(:,1);
myy2 = y2(:,2);
myy2 = myy2 - myy2(1);
mysystemplot(myt2,myy2)

myt3 = y3(:,1);
myy3 = y3(:,2);
myy3 = myy3 - myy3(1);
mysystemplot(myt3,myy3)

myt4 = y4(:,1);
myy4 = y4(:,2);
myy4 = myy4 - myy4(1);
mysystemplot(myt4,myy4)

%% Part II WinCon
% obtaining Mp, tr, and ts for the different sets of ps.

clear all
close all
clc

load('StepResponseFile.mat')

myt1 = Scope(1000:1500,1);
myy1 = V01(1000:1500);
myy1 = myy1-myy1(1);
myr1 = Vr1(1000:1500);
myr1 = myr1 - myr1(1);
mysystemplot(myt1,myy1)

myt2 = Scope(1000:1500,1);
myy2 = V02(1000:1500);
myy2 = myy2-myy2(1);
myr2 = Vr2(1000:1500);
myr2 = myr2 - myr2(1);
mysystemplot(myt2,myy2)

myt3 = Scope(1000:1500,1);
myy3 = V03(1000:1500);
myy3 = myy3-myy3(1);
myr3 = Vr3(1000:1500);
myr3 = myr3 - myr3(1);
mysystemplot(myt3,myy3)

myt4 = Scope(1000:1500,1);
myy4 = V04(1000:1500);
myy4 = myy4-myy4(1);
myr4 = Vr4(1000:1500);
myr4 = myr4 - myr4(1);
mysystemplot(myt4,myy4)

%% Part IV WinCon With Friction
% obtaining Mp, tr, and ts for the different sets of ps.

clear all
close all
clc

load('fric.mat')

myt1 = t(1000:1500);
myy1 = F01(1000:1500);
myy1 = myy1-myy1(1);
myr1 = Fr1(1000:1500);
myr1 = myr1 - myr1(1);
[y1,t1] = projection(myt1,myy1);
mysystemplot(t1,y1)

myt2 = t(1000:1500);
myy2 = F02(1000:1500);
myy2 = myy2-myy2(1);
myr2 = Fr2(1000:1500);
myr2 = myr2 - myr2(1);
[y2,t2] = projection(myt2,myy2);
mysystemplot(t2,y2)

myt3 = t(1000:1500);
myy3 = F03(1000:1500);
myy3 = myy3-myy3(1);
myr3 = Fr3(1000:1500);
myr3 = myr3 - myr3(1);
[y3,t3] = projection(myt3,myy3);
mysystemplot(t3,y3)

myt4 = t(1000:1500);
myy4 = F04(1000:1500);
myy4 = myy4-myy4(1);
myr4 = Fr4(1000:1500);
myr4 = myr4 - myr4(1);
[y4,t4] = projection(myt4,myy4);
mysystemplot(t4,y4)







