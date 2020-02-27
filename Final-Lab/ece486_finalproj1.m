close all
clc

tr = 0.2;
br = 198;
a = 2.2/tr;
K = a/br;

% Ucl1
REF = 100;
set_param('chapter2_P2b', 'SimulationCommand', 'start')
pause(15)
Ucl1n = Ucl(end);


%%
wn = -500:100:-100;
wp = 100:100:500;

figure(1)
plot(wn,[Ucl5n,Ucl4n,Ucl3n,Ucl2n,Ucl1n]);
figure(2)
plot(wp,[Ucl1,Ucl2,Ucl3,Ucl4,Ucl5]);

%%
close all
clc

tr = 0.2;
br = 198;
a = 2.2/tr;
K = a/br;

REF = -200;
set_param('chapter2_P2b', 'SimulationCommand', 'start')

%% Chapter 3: System Identifivacion

clear
clc
close

tr = 0.2;
br = 198;
a = 2.2/tr;
K = a/br;

mp = 0.2164; % mass of pendulum and motor housing/stator, in kg.
mr = 0.0850; % mass of the rotor/wheel, in kg.
m = mp + mr; % combined mass of rotor and pendulum.
Jp = 2.233*10^-4; % moment of inertia of the pendulum about its center of mass, in kg*m^2.
Jr = 2.495*10^-5; % moment of inertia of rotor about its center of mass, in kg*m^2.
lp = 0.1173; % distance from pivot to the center of the pendulum, in m.
lr = 0.1270; % distance of pivot to the center of mass of rotor, in m.
J = Jp + mp*lp*lp + mr*lr*lr; % total inertia when considering a non rotating rotor attached to the end of the pendulum as an object. You treat rotor plus pendulum as an overall object.
l = (mp*lp + mr*lr)/m; % distance from center of mass of pendulum and rotor, in m.
k =     0.0556; % torque constante of the rotor;
g = 9.8; % Earth's gravity, in m/s^2.
Wnp = sqrt( m*g*l/J );  % frequency of small oscillations of the system around the hanging position.
Wnp_prime = sqrt( m*g*l/(J+Jr) ); % frequency of oscillations of the system including Jr around the hanging position.
Wmp_prime_meas = 2*pi*23/16.442;
%br = w_r_dot/(5+.13*w_r+0.47);
%set_param('chapter2_P2', 'SimulationCommand', 'start')
%%
clc
[bpMax, brMax] = DetermineTorqConstants2(theta_p(:,1),theta_p(:,2),theta_r(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chapter 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
bp = .95;
br = 173;
A = [0 1 0 0; a 0 0 0; 0 0 0 1; 0 0 0 0];
B = [ 0; -bp; 0; br];
% Controlable?
R = ctrb(A,B);
Rank = rank(R);
%% Pole Placement + Control Two-State
clc
a = 8.8^2;
bp = .95;
A2 = [0 1; a 0 ];
B2 = [ 0; -bp];
zeta = 0.5/sqrt(2);
Wn = 1.85*Wnp;
sys2 = [1 2*zeta*Wn Wn*Wn];
pls2 = roots(sys2);
Kontrol2 = place(A2,B2,pls2);
k12 = Kontrol2(1);
k22 = Kontrol2(2);
%% Pole Placement + Control Three-State
clc
bp = 1.089;
br = 180.2;
a = 8.8^2;
A = [0 1 0 0; ...
      a 0 0 0; ...
      0 0 0 1; ...
      0 0 0 0];
B = [ 0; -bp; 0; br];
rank(ctrb(A,B));
zeta = 0.5*(1/sqrt(2));
Wn = 1.2*Wnp;
sys = [1 2*zeta*Wn Wn*Wn];
pls = roots(sys);
poles_cl = [pls(1) pls(2) 0 real(pls(1))];
Kontrol3 = place(A,B,poles_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chapter 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
bp = 1.089;
br = 180.2;
a = 8.8^2;
A = [0 1 0 0; a 0 0 0; 0 0 0 1; 0 0 0 0];
B = [ 0; -bp; 0; br];
C = [ 1 0 0 0; 0 0 1 0];
% Controler K
zeta = 0.5*(1/sqrt(2));
Wn = 1.2*Wnp;
sys = [1 2*zeta*Wn Wn*Wn];
pls = roots(sys);
K = place(A,B,[pls(1) pls(2) 0 real(pls(1)) ]); 
% 5a Observer l
L = place(A',C',[ 5*real(pls(1)) 10*real(pls(2)) -14 12*real(pls(1))]); % Very far to the left for fast convergance
L = L';
eig(A-L*C);
pls*20;
% 5.2 
L_hat = [L(1,1) 0; L(2,1) 0; 0 L(3,2); 0 L(4,2)];
A12 = [0 1; a 0];
B1 = [0; -bp];
C12 = [1 0];
L12 = place(A12',C12',[50*real(pls(1)) 100*real(pls(2))]);
A34 = [ 0 1; 0 0];
B2 = [0;br];
C34 = C12;
L34 = place(A34',C34',[-14*10 120*real(pls(1))]);
L_hat_2 = [L12 0 0; 0 0 L34];
L_hat_2 = L_hat_2';
L_hat = L_hat_2;
L = L_hat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chapter 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up Control
clc
bp = 1.089;
br = 180.2;
a = 8.8^2;
A = [0 1 0 0; ...
      a 0 0 0; ...
      0 0 0 1; ...
      0 0 0 0];
B = [ 0; -bp; 0; br];
zeta = 0.8*(1/sqrt(2));
Wn = 1.3*Wnp;
sys = [1 2*zeta*Wn Wn*Wn];
pls = roots(sys);
poles_cl = [pls(1) pls(2) 0 real(pls(1))];
Kup= place(A,B,poles_cl);
% Down control
A_down = A;
A_down(2,1) = -a;
Kdown = place(A_down,B,poles_cl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chapter 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Max_potent = mp*g*lp+mr*g*lr;
%Cur_pot = -mp*g*lp*cos(THETA_P)- mr*g*lr*cos(THETA_P);
%Cur_kin = 0.5*(Jp+mp*lp^2+mr*lr^2)*VELOCITY_PEN^2+0.5*Jr*VELOCITY_ROT^2;

%Tot_ME = Cur_pot + Cur_Kin;











