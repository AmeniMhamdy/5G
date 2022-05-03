clc; clear variables; close all;

N = 10^5;

d1=1; d2=1;

SNR = 20:2:40;
snr = db2pow(SNR);
eta = 4;    %Path loss exponent

%Rayleigh fading coefficients for both users
h1 = (sqrt(d1^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h2 = (sqrt(d2^-eta))*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

%Generate bits
x1 = randi([0 1],N,1);
x2 = randi([0 1],N,1);

%Do BPSK modulation
x1_bpsk = 2*x1 - 1;
x2_bpsk = 2*x2 - 1;

%Power allocation coefficients
a1 = 0.75; a2 = 0.25;

%Superposition coding
x_bpsk = sqrt(a1).*x1_bpsk + sqrt(a2).*x2_bpsk;

ber1_noma = zeros(1,length(SNR));
ber2_noma = zeros(1,length(SNR));
ber1_oma = zeros(1,length(SNR));
ber2_oma = zeros(1,length(SNR));


for u = 1:length(SNR)

   %For NOMA 
   y1_noma = awgn(h1.*x_bpsk,SNR(u),'measured');
   y2_noma = awgn(h2.*x_bpsk,SNR(u),'measured');
   
   eq1_noma = y1_noma./h1;
   eq2_noma = y2_noma./h2;
   
   %At user 1
   %Direct decoding of x1
   dec1_noma = zeros(N,1);
   dec1_noma(eq1_noma>0) = 1;
   
   %At user 2
   %Direct decoding of x1
   dec12_noma = zeros(N,1);
   dec12_noma(eq2_noma>0)=1;
   dec12_remod_noma = 2*dec12_noma - 1;
   rem_noma = eq2_noma - sqrt(a1).*dec12_remod_noma;
   dec2_noma = zeros(N,1);
   dec2_noma(rem_noma>0)=1;
   
   %BER calculation
   ber1_noma(u) = biterr(dec1_noma,x1)/N;
   ber2_noma(u) = biterr(dec2_noma,x2)/N;
   
   %For OMA 
   y1_oma = awgn(h1.*x1_bpsk,SNR(u)','measured');
   y2_oma = awgn(h2.*x2_bpsk,SNR(u)','measured');
   
   eq1_oma = y1_oma./h1;
   eq2_oma = y2_oma./h2;
   
   dec1_oma = zeros(N,1);
   dec1_oma(eq1_oma>0)=1;
   dec2_oma = zeros(N,1);
   dec2_oma(eq2_oma>0)=1;
   
   %BER calculation
   ber1_oma(u) = biterr(dec1_oma,x1)/N;
   ber2_oma(u) = biterr(dec2_oma,x2)/N;   
end

semilogy(SNR,ber1_noma,'r-o','linewidth',2); hold on; grid on;
semilogy(SNR,ber2_noma,'b-o','linewidth',2);
semilogy(SNR,ber1_oma,'y-o','linewidth',2);
semilogy(SNR,ber2_oma,'g-o','linewidth',2);
legend('NOMA - user 1','NOMA - user 2','OMA - user 1','OMA - user 2');
xlabel('SNR (dB)');
ylabel('BER');





