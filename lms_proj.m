clc;
clear all;
close all;
size=64;
fs=2000;
n=1:size;
f1=100;
f2=150;
MSE=zeros(1,10000);
signala=sin(2*pi*f1*n/fs);
signalb=sin(2*pi*f2*n/fs);
s=(0.5*signala)+(signalb);
s=s';
mean1=0;
% sd= 1.9021;  
sd=sqrt(10^((10*log(var(s))+7)/10));
w=sd*((rand(1,size))')+mean1;
d= s+w;
subplot(3,2,1);
plot(n,d);
title('d');
input=[0 ; d(1:length(d)-1)];
subplot(3,2,2);
plot(n,input);   
title('input');
% Fc1=50;
% Fc2=190;
% a=fir1(size-1, [Fc1 Fc2]/(fs), 'bandpass');
% a=a';
 shiftedsignal=[0  ; s(1:length(s)-1)];
 shiftednoise= [0  ; w(1:length(w)-1)];
%  crshiftedsignal= autocorr(shiftedsignal);
%  crshiftednoise= autocorr(shiftednoise);
% crshiftedsignal= autom(shiftedsignal)';
% crshiftednoise= autom(shiftednoise)';
% 
%  R=toeplitz(crshiftedsignal)+toeplitz(crshiftednoise);
R=toeplitz(autom(shiftedsignal+shiftednoise));
%  pd=autocorr(d);
 pd=autom(d)';
aWeiner=R\pd;
y=aWeiner.*input;
subplot(3,2,3);
plot(n,y);
title('FIR output');
mu=0.0005; 
aTemp = aWeiner;
for j=1:10000
        y1=aTemp.*input;
        e = d -y1 ;
        a(:,j) = aTemp + (mu *( e .* input));
        aTemp = a(:,j);
        mse=e.^2;
        MSE(j)=mean(mse); % end of part 1
end
aFinal = aTemp;
yopt=aFinal.*input;
subplot(3,2,4);
plot(1:10000, MSE);
title('MSE');
impulse=[1,zeros(1,size-1)]; 
b=1;
%[h,t]=impz(aFinal,b);
h=filter(aFinal,b,impulse);
subplot(3,2,5);
plot(h);
title('Impulse response');
% dt=zeros(1,size);
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
% end of part 3
SignaldB1 = 10 * log10(var(y));
NoisedB1 = 10 * log10(var(y-input));
SNR1 = SignaldB1 - NoisedB1;
% snr1 = var(y)/var(y-input);
SignaldB2 = 10 * log10(var(yopt));
NoisedB2 = 10 * log10(var(yopt-input));
SNR2 = SignaldB2 - NoisedB2; % end of part 4
% snr2 = var(yopt)/var(yopt-input);
%%% coefficient plot %%%
figure
plot(1:10000, a(2,:),'r',1:10000, a(17,:),'b' ,1:10000, a(14,:),'m', 1:10000, a(8,:),'g')
hleg=legend('coeff = 2 ','coeff = 17','coeff = 14', 'coeff = 8');