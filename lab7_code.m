%% Report item 1
clc, close all;
clear all;

N = 512;
del_t = 0.01;
n = (-(N-1):N-1);
t = n*del_t;
psi_1_0 = 2^(1/2)*mother_wavelet(2*t - 0);
psi_0_1 = 2^(0/2)*mother_wavelet(t-1);
psi_1_1 = 2^(1/2)*mother_wavelet(2*t - 1);
psi_2_1 = 2^(2/2)*mother_wavelet(4*t - 1);

figure(1);
subplot(411);
plot(t,psi_1_0);
title('Psi 1,0 in time');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(412);
plot(t,psi_0_1);
title('Psi 0,1 in time');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(413);
plot(t,psi_1_1);
title('Psi 1,1 in time');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(414);
plot(t,psi_2_1);
title('Psi 2,1 in time');
xlabel('Time (sec)');
ylabel('Amplitude');

PSI_1_0 = fftshift(fft(psi_1_0)); 
Mag  = abs(PSI_1_0);
Phase = angle(PSI_1_0);
w = fftshift((0:(2*N-2))*2*pi/(2*N-1));
w(1:N-1) = w(1:N-1) - 2*pi;

figure(2);
subplot(421);
plot(w,Mag);
title('PSI 1,0 in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');
%Plot the phase
subplot(422);
plot(w,Phase);
title('PSI 1,0 Phase');
xlabel('Freq (rad/sec)');
ylabel('Phase')


PSI_0_1 = fftshift(fft(psi_0_1)); 
Mag  = abs(PSI_0_1);
Phase = angle(PSI_0_1);

subplot(423)
plot(w,Mag)
title('PSI 0,1 in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(424)
plot(w,Phase)
title('PSI 0,1 Phase');
xlabel('Freq (rad/sec)');
ylabel('Phase')

PSI_1_1 = fftshift(fft(psi_1_1));
Mag  = abs(PSI_1_1);
Phase = angle(PSI_1_1);

subplot(425)
plot(w,Mag)
title('PSI 1,1 in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(426)
plot(w,Phase)
title('PSI 1,1 Phase');
xlabel('Freq (rad/sec)');
ylabel('Phase')

PSI_2_1 = fftshift(fft(psi_2_1));  %Take DFT and shift over range of (0,2*pi)
Mag  = abs(PSI_2_1);
Phase = angle(PSI_2_1);

subplot(427)
plot(w,Mag)
title('PSI 2,1 in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(428)
plot(w,Phase)
title('PSI 0,1 Phase');
xlabel('Freq (rad/sec)');
ylabel('Phase')

%% Report Item 2 (Types of)
[phi,psi,t] = wavefun('coif1',5);

figure(3);
subplot(321);
plot(t,phi,'linewidth',2);
title('coif1 Father Wavelet');
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);

subplot(322);
plot(t,psi,'linewidth',2);
title('coif1 Mother Wavelet');
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

N = 512;
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

PHI = fftshift(fft(phi,N));
Mag_c = abs(PHI);

PSI = fftshift(fft(psi,N));
Mag_cf = abs(PSI);

[phi,psi,t] = wavefun('db1',5);

subplot(323);
plot(t,phi,'linewidth',2);
title('db1 Father Wavelet');
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);

subplot(324);
plot(t,psi,'linewidth',2);
title('db1 Mother Wavelet');
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

PHI = fftshift(fft(phi,N));
Mag_d = abs(PHI);

PSI = fftshift(fft(psi,N));
Mag_df = abs(PSI);

[phi,psi,t] = wavefun('sym4',5);
subplot(325);
plot(t,phi,'linewidth',2);
title('sym4 Father Wavelet');
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);

subplot(326);
title('sym4 Mother Wavelet');
plot(t,psi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

PHI = fftshift(fft(phi,N));
Mag_s = abs(PHI);

PSI = fftshift(fft(psi,N));
Mag_sf = abs(PSI);

figure(4);
subplot(321);
plot(w,Mag_cf);
title('coif1 Father Wavelet in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(322);
plot(w,Mag_c);
title('coif1 Mother Wavelet in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(323);
plot(w,Mag_df);
title('db1 Father Wavelet in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(324);
plot(w,Mag_d);
title('db1 Mother Wavelet in Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(325);
plot(w,Mag_sf);
title('sym4 in Father Wavelet Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');

subplot(326);
plot(w,Mag_s);
title('sym4 in Mother Wavelet Freq Domain');
xlabel('Freq (rad/sec)');
ylabel('Magnitude');
%% Report item 3a
clear all;
load('signal.mat')
[cA_s,cD_s] = dwt(x, 'sym4');
[cA_c,cD_c] = dwt(x, 'coif1');

figure(5);
subplot(221);
stem(cA_s);
title('sym4 Approximation Coefficients');
xlabel('Coefficients');
ylabel('Value');

subplot(222);
stem(cD_s);
title('sym4 Detail Coefficients');
xlabel('Coefficients');
ylabel('Value');

subplot(223);
stem(cA_c);
title('coif1 Approximation Coefficients');
xlabel('Coefficients');
ylabel('Value');

subplot(224);
stem(cD_c);
title('coif1 Detail Coefficients');
xlabel('Coefficients');
ylabel('Value');

%% Report Item 3b
clc;
x_s = idwt(cA_s, cD_s, 'sym4');
x_c = idwt(cA_c, cD_c, 'coif1');

cA_s_avg = mean(abs(cA_s));
filter = (abs(cA_s) > cA_s_avg);
cA_s_thld = filter .* cA_s;

cD_s_avg = mean(abs(cD_s));
filter = (abs(cD_s) > cD_s_avg);
cD_s_thld = filter .* cD_s;

cA_c_avg = mean(abs(cA_c));
filter = (abs(cA_c) > cA_c_avg);
cA_c_thld = filter .* cA_c;

cD_c_avg = mean(abs(cD_c));
filter = (abs(cD_c) > cD_c_avg);
cD_c_thld = filter .* cD_c;

x_s_thld = idwt(cA_s_thld, cD_s_thld, 'sym4');
x_c_thld = idwt(cA_c_thld, cD_c_thld, 'coif1');

figure(6);
subplot(311);
stem(x);
title('Original sample x');
xlabel('Samples');
ylabel('Value');

subplot(312);
stem(x_s);
title('Inverse DWT of x from sym4');
xlabel('Samples');
ylabel('Value');

subplot(313);
stem(x_s_thld);
title('Inverse DWT of x from sym4 w/ hard threshold');
xlabel('Samples');
ylabel('Value');

figure(7);
subplot(311);
stem(x);
title('Original sample x');
xlabel('Samples');
ylabel('Value');

subplot(312);
stem(x_c);
title('Inverse DWT of x from coif1');
xlabel('Samples');
ylabel('Value');

subplot(313);
stem(x_c_thld);
title('Inverse DWT of x from coif1 w/ hard threshold');
xlabel('Samples');
ylabel('Value');

%% Report Item 4 (De-noise)
thresholds = find_threshold();