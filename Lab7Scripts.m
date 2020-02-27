%% Mother Wavelet
clc;

size = 1001;
time = linspace(-5,5,size);
wavelet_1 = 2^(1/2)*Wavelet_time(2*time, size);
wavelet_2 = Wavelet_time(time-1,size);
wavelet_3 = 2^(1/2)*Wavelet_time(2*time - 1, size);
wavelet_4 = 2^(2/2)*Wavelet_time(4*time-1,size);

figure
subplot(221);
plot(time,wavelet_1);
title('Wavelet_1');
xlabel('Time (seconds)');
ylabel('Amplitude');
subplot(222);
plot(time,wavelet_2);
title('Wavelet_2');
xlabel('Time (seconds)');
ylabel('Amplitude');
subplot(223)
plot(time,wavelet_3);
title('Wavelet_2');
xlabel('Time (seconds)');
ylabel('Amplitude');
subplot(224)
plot(time,wavelet_4);
title('Wavelet_4');
xlabel('Time (seconds)');
ylabel('Amplitude');

Xd = fftshift(fft(wavelet_1));  %Take DFT and shift over range of (0,2*pi)
Mag  = abs(Xd);
Phase = angle(Xd);
N = length(wavelet_1);
%Shift frequency to the right
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

figure;

%Plot the magnitude 
subplot(421)
plot(w,Mag)
title('Wavelet_1 Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%Plot the phase
subplot(422)
plot(w,Phase)
title('Phase of ψ 1,0');
xlabel('Frequency(Radians/sec)');
ylabel('Phase in Hz')

Xd = fftshift(fft(wavelet_2));  %Take DFT and shift over range of (0,2*pi)
Mag  = abs(Xd);
Phase = angle(Xd);
N = length(wavelet_2);
%Shift frequency to the right
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

%Plot the magnitude 
subplot(423)
plot(w,Mag)
title('Wavelet_2 Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%Plot the phase
subplot(424)
plot(w,Phase)
title('Phase of ψ 0,1 ');
xlabel('Frequency(Radians/sec)');
ylabel('Phase in Hz')

Xd = fftshift(fft(wavelet_3));  %Take DFT and shift over range of (0,2*pi)
Mag  = abs(Xd);
Phase = angle(Xd);
N = length(wavelet_3);
%Shift frequency to the right
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

%Plot the magnitude 
subplot(425)
plot(w,Mag)
title('Wavelet_3 Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%Plot the phase
subplot(426)
plot(w,Phase)
title('Phase of ψ 1,1 ');
xlabel('Frequency(Radians/sec)');
ylabel('Phase in Hz)')


Xd = fftshift(fft(wavelet_4));  %Take DFT and shift over range of (0,2*pi)
Mag  = abs(Xd);
Phase = angle(Xd);
N = length(wavelet_4);
%Shift frequency to the right
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

%Plot the magnitude 
subplot(427)
plot(w,Mag)
title('Wavelet_4 Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%Plot the phase
subplot(428)
plot(w,Phase)
title('Phase of ψ 2,1 ');
xlabel('Frequency(Radians/sec)');
ylabel('Phase in Hz')

%% Types of Wavelets
clc
N = 512;    %Zero pad
%Shift frequency to the right
w = fftshift((0:N - 1)/N*2*pi);
w(1:N/2) = w(1:N/2) - 2*pi;

[phi,psi,t] = wavefun('coif1',5);
figure;
subplot(211);
plot(t,phi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);
subplot(212)
plot(t,psi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

DFT_phi = fftshift(fft(phi,N));
Magnitude = abs(DFT_phi);

[phi,psi,t] = wavefun('db1',5);
figure;
subplot(211);
plot(t,phi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);
subplot(212);
plot(t,psi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

DFT_phi = fftshift(fft(phi,N));
Magnitude_2 = abs(DFT_phi);

figure;
subplot(311);
plot(w,Magnitude);
subplot(312);
plot(w,Magnitude_2);

[phi,psi,t] = wavefun('sym4',5);
figure;
subplot(211);
plot(t,phi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\phi(t)','fontsize',16);
set(gca,'fontsize',14);
subplot(212);
plot(t,psi,'linewidth',2);
xlabel('t','fontsize',16);
ylabel('\psi(t)','fontsize',16);
set(gca,'fontsize',14);

DFT_phi = fftshift(fft(phi,N));
Magnitude_3 = abs(DFT_phi);

figure;
subplot(311);
plot(w,Magnitude);
subplot(312);
plot(w,Magnitude_2);
subplot(313);
plot(w,Magnitude_3);
%% Wavelet Decomposition
clc;
load('signal.mat')
[A,D] = dwt(x, 'sym4');
[Lo_D,Hi_D] = wfilters('bior3.5', 'd');
[A,D] = dwt(x,Lo_D,Hi_D);

subplot(221);
stem(A);
title('Approximation Coefficients (sym4)');
subplot(222);
stem(D);
title('Detail Coefficients (sym4)');

[A,D] = dwt(x, 'coif1');
[Lo_D,Hi_D] = wfilters('bior3.5', 'd');
[A,D] = dwt(x,Lo_D,Hi_D);

subplot(223);
stem(A);
title('Approximation Coefficients (coif1)');
subplot(224);
stem(D);
title('Detail Coefficients (coif1)');
















































