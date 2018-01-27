clc;
clear all;
close all;
size=64;
fs=2000;
n=1:size;
f1=100;
f2=150;
signala=sin(2*pi*f1*n/fs);
signalb=sin(2*pi*f2*n/fs);
s=(0.5*signala)+(signalb);
s=s';
mean1=0;
% sd= 1.9021;  
sd=sqrt(10^((10*log(var(s))+7)/10));
w=sd*((rand(1,size))')+mean1;
d= s+w;
input=[0 ; d(1:length(d)-1)];
 shiftedsignal=[0  ; s(1:length(s)-1)];
 shiftednoise= [0  ; w(1:length(w)-1)];
R=toeplitz(autom(shiftedsignal+shiftednoise));
 pd=autom(d)';
aWeiner=R\pd;
y=aWeiner.*input;
mu=[0.001 0.01 0.1]; 
MSE=zeros(3,20000);
for nu=1:3
aTemp = aWeiner;
for j=1:20000
        y1=aTemp.*input;
        e = d -y1 ;
        a(:,j) = aTemp + (mu(nu) *( e .* input));
        aTemp = a(:,j);
        mse=e.^2;
        MSE(nu,j)=mean(mse); % end of part 1
end
end
% aFinal(nu) = aTemp;
% yopt=aFinal(nu).*input;

plot(1:20000, MSE(1,:),'r',1:20000, MSE(2,:),'b', 1:20000, MSE(3,:), 'g');
title('MSE');
hleg=legend('mu = 0.001', 'mu = 0.01', 'mu = 0.1'); 
% impulse=[1,zeros(1,size-1)]; 
% b=1;
% h=filter(aFinal,b,impulse);
% 
% dt=(fft(input));
% subplot(4,2,6);
% plot(abs(dt));
% ylabel('dt(k)--------->');
% xlabel('k--------->');
% title('FFT of sequence(Magnitude Plot)');
% subplot(4,2,7);
% plot(angle(dt));
% ylabel('dt(k)--------->');
% xlabel('k--------->');
% title('FFT of sequence(phase Plot)');
% % end of part 3
% SignaldB1 = 10 * log10(var(y));
% NoisedB1 = 10 * log10(var(y-input));
% SNR1 = SignaldB1 - NoisedB1;
% snr1 = var(y)/var(y-input);
% SignaldB2 = 10 * log10(var(yopt));
% NoisedB2 = 10 * log10(var(yopt-input));
% SNR2 = SignaldB2 - NoisedB2; % end of part 4
% snr2 = var(yopt)/var(yopt-input);
% %%% coefficient plot %%%
% figure
% plot(1:10000, a(2,:),'r',1:10000, a(17,:),'b' ,1:10000, a(14,:),'m', 1:10000, a(8,:),'g')