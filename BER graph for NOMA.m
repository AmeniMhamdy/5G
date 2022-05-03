clc; clear variables; close all;
N = 10^5; %Number of monte carlosimulations
SNR = 0:40; %SNR range in dB
snr = db2pow(SNR);%SNR range in linear scale(db to power)% 10*log10(snr)=SNR
%Power asignation
a1 = 0.75; a2 = 0.18; a3 = 0.05; a4 = 0.02;
%random bits
x1 = randi([0 1],1,N);
x2 = randi([0 1],1,N);
x3 = randi([0 1],1,N);
x4 = randi([0 1],1,N);
%MOD BPSK
xmod1 = 2*x1-1;
xmod2 = 2*x2-1;
xmod3 = 2*x3-1;
xmod4 = 2*x4-1;

%fading rayleigh
h1 = raylrnd(1,1,N);
h2 = raylrnd(1,1,N);
h3 = raylrnd(1,1,N);
h4 = raylrnd(1,1,N);

x11 = h1.*(sqrt(a1)*xmod1) + h1.*(sqrt(a2)*xmod2) + h1.*(sqrt(a3)*xmod3) + h1.*(sqrt(a4)*xmod4); %sum of the two signals
x21 = h2.*(sqrt(a1)*xmod1) + h2.*(sqrt(a2)*xmod2) + h2.*(sqrt(a3)*xmod3) + h2.*(sqrt(a4)*xmod4);
x31 = h3.*(sqrt(a1)*xmod1) + h3.*(sqrt(a2)*xmod2) + h3.*(sqrt(a3)*xmod3) + h3.*(sqrt(a4)*xmod4);
x41 = h4.*(sqrt(a1)*xmod1) + h4.*(sqrt(a2)*xmod2) + h4.*(sqrt(a3)*xmod3) + h4.*(sqrt(a4)*xmod4);

for u = 1:length(snr)
y1 = awgn(x11,SNR(u),'measured');%Received signal at user 1 corrupted by AWGN
y2 = awgn(x21,SNR(u),'measured');%Received signal at user 2 corrupted by AWGN
y3 = awgn(x31,SNR(u),'measured');%Received signal at user 2 corrupted by AWGN
y4 = awgn(x41,SNR(u),'measured');%Received signal at user 2 corrupted by AWGN

%AT USER 1
%Direct decoding of x from y1
x1_hat = ones(1,N); %Just a buffer
x1_hat(y1 < 0) = 0; %Final bits for user 1

%AT USER 2
%Direct decoding of x from y2(dd2)
x21_est = ones(1,N); %just a buffer
x21_est(y2 < 0) = -1; %Estimate user 1's signal first
%Subtract remodulated x11_est component from y2
rem2 = y2 - sqrt(a1)*x21_est;%y2 without signal x1

%Decode x2 from rem without signal x1
x2_hat = zeros(1,N);
x2_hat(rem2>0) = 1; %Final bits for user 2

%AT USER 3
%Direct decoding of x from y2(dd3) First of all, estimate the previous
%to subtract from y3
x31_est = ones(1,N); %just a buffer
x31_est(y3 < 0) = -1; %Estimate user 1's signal first

rem3 = y3 - sqrt(a1)*x31_est; %Extract useless signal 1 from y3
x32_est = ones(1,N);
x32_est(rem3<0) = -1;
%Subtract remodulated x32_est component from y3
rem3=rem3 - sqrt(a2)*x32_est;
%Decode x3 from rem
x3_hat = zeros(1,N);
x3_hat(rem3>0) = 1; %Final bits for user 3

%AT USER 4
%Direct decoding of x from y2(dd3) First of all, estimate the previous
%to subtract from y4
x41_est = ones(1,N); %just a buffer
x41_est(y4 < 0) = -1; %Estimate user 1's signal first

rem4 = y4 - sqrt(a1)*x41_est; %Extract useless signal 1 from y4
x42_est = ones(1,N);
x42_est(rem4<0) = -1;
%Subtract remodulated x42_est component from y4
rem4=rem4 - sqrt(a2)*x42_est;

x43_est = ones(1,N);
x43_est(rem4<0) = -1;
rem4 = rem4 - sqrt(a3)*x43_est;

x4_hat = zeros(1,N);
x4_hat(rem4>0) = 1;

%Estimate BER
ber1(u) = biterr(x1,x1_hat)/N;
ber2(u) = biterr(x2,x2_hat)/N;
ber3(u) = biterr(x3,x3_hat)/N;
ber4(u) = biterr(x4,x4_hat)/N;
end

%plot BER curves
semilogy(SNR, ber1,'b-o' ,'linewidth', 1.5); hold on;%Semilog plot
semilogy(SNR, ber2,'r-o', 'linewidth', 1.5); grid on;
semilogy(SNR, ber3, 'y-o','linewidth', 1.5);
semilogy(SNR, ber4,'g-o', 'linewidth', 1.5);
legend('User 1 \alpha_1 = 0.75','User 2 \alpha_2 = 0.18','User 3 \alpha_2 = 0.05','User 4 \alpha_2 = 0.02','Location','southwest');
xlabel('SNR (dB)');
ylabel('BER');
title('BER graph for NOMA in AWGN channel');