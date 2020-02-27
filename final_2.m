% Report Item 2 : Chipotle
clc, clear all, close all;

load samplerate.mat;
X = fftshift(fft(x));
N = length(x);

w1 = fftshift((0:N-1)/N*2*pi);
w1(1:N/2) = w1(1:N/2)-2*pi;
fs = 40;
W1 = fs*w1/(2*pi);

figure(1);
subplot(211);
plot(W1,abs(X));
xlabel('Frequency','fontsize',12); 
title('Chipotle Magnitude Spectrum','fontsize',12);
subplot(212);
stem(x);
xlabel('Time','fontsize',12); 
title('Time Domain','fontsize',12);

%%% Upsample
y = upsample(x,3); % Upsample by 3
Y = fftshift(fft(y));
N2 = length(Y);

w2 = fftshift((0:N2-1)/N2*2*pi);
w2(1:N2/2) = w2(1:N2/2)-2*pi;
W2 = fs*w2/(2*pi);

figure(2);
subplot(211);
plot(W2,abs(Y));
xlabel('Frequency','fontsize',12); 
title('Upsampled Magnitude Spectrum','fontsize',12);
subplot(212);
stem(y);
xlabel('Time','fontsize',12); 
title('Time Domain','fontsize',12);


%%% Ideal Low Pass Filter
for i=1:N2
    if abs(W2(i))>4
        Y(i)=0;
    end
end

y = ifft(ifftshift(Y));

figure(3);
subplot(211);
plot(W2,abs(Y));
xlabel('Frequency','fontsize',12); 
title('LPF-ed Magnitude Spectrum','fontsize',12);
subplot(212);
stem(y);
xlabel('Time','fontsize',12); 
title('Time Domain','fontsize',12);

%%% Downsample
z = downsample(y,2);
Z = fftshift(fft(z));
N3 = length(z);

w3 = fftshift((0:N3-1)/N3*2*pi);
w3(1:N3/2) = w3(1:N3/2)-2*pi;
W3 = fs*w3/(2*pi);

figure(4);
subplot(211);
plot(W3,abs(Z));
xlabel('Frequency','fontsize',12); 
title('Downsampled Magnitude Spectrum','fontsize',12);
subplot(212);
stem(z);
xlabel('Time','fontsize',12); 
title('Time Domain','fontsize',12);