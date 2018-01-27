clc;
clear all;
close all;
size=64;
fs=2000;
n=1:3*size;
f1=100;
f2=150;
MSE=zeros(1,30000);
signala=sin(2*pi*f1*n/fs);
signalb=sin(2*pi*f2*n/fs);
s=(0.5*signala)+(signalb);
s=s';
mean1=0;
% sd= 1.9021;    
sd=sqrt(10^((10*log(var(s))+7)/10));
w=sd*((rand(1,3*size))')+mean1;
d= s+w;
subplot(2,2,1);
plot(n,d);
title('d');
input=[0 ; d(1:length(d)-1)];
subplot(2,2,2);
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

shifted_signal=zeros(size,3);
shifted_signal(:,1)=shiftedsignal(1:size,1);
shifted_signal(:,2)=shiftedsignal(size+1:2*size,1);
shifted_signal(:,3)=shiftedsignal(2*size+1:end,1);
shifted_noise=zeros(size,3);
shifted_noise(:,1)=shiftednoise(1:size,1);
shifted_noise(:,2)=shiftednoise(size+1:2*size,1);
shifted_noise(:,3)=shiftednoise(2*size+1:end,1);
input_signal=zeros(size,3);
input_signal(:,1)=input(1:size,1);
input_signal(:,2)=input(size+1:2*size,1);
input_signal(:,3)=input(2*size+1:end,1);
d_signal=zeros(size,3);
d_signal(:,1)=d(1:size,1);
d_signal(:,2)=d(size+1:2*size,1);
d_signal(:,3)=d(2*size+1:end,1);
y=zeros(size,3);

for i=1:3
crshiftedsignal= autom(shifted_signal(:,i))';
crshiftednoise= autom(shifted_noise(:,i))';

 R=toeplitz(crshiftedsignal)+toeplitz(crshiftednoise);
%  pd=autocorr(d);
 pd=autom(d_signal(:,i))';
aWeiner=R\pd;
y(:,i)=aWeiner.*(input_signal(:,i));
end
yC=[y(:,1); y(:,2); y(:,3)];
subplot(2,2,3);
plot(n,yC);
title('FIR output');
mu=0.005; 
aTemp = aWeiner;
% for i=1:3
y1m=zeros(size,3);
em=zeros(size,3);

for j=1:30000
        y1m(:,1)=aTemp.*input_signal(:,1);
        y1m(:,2)=aTemp.*input_signal(:,2);
        y1m(:,3)=aTemp.*input_signal(:,3);
%         y1=[y1m(:,1); y1m(:,2); y1m(:,3)];
        em(:,1) = d_signal(:,1)-y1m(:,1) ;
        em(:,2) = d_signal(:,2)-y1m(:,2) ;
        em(:,3) = d_signal(:,3)-y1m(:,3) ;
        e=[em(:,1); em(:,2); em(:,3)];
        a(:,j)= aTemp + (mu *(( em(:,1) .* input_signal(:,1))+( em(:,2) .* input_signal(:,2))+( em(:,3) .* input_signal(:,3)))/size);
        aTemp = a(:,j);
        mse=e.^2;
        MSE(j)=mean(mse); % end of part 1
end
% end
aFinal = aTemp;
yoptm(:,1)=aFinal.*input_signal(:,1);
yoptm(:,2)=aFinal.*input_signal(:,2);
yoptm(:,3)=aFinal.*input_signal(:,3);
yopt=[yoptm(:,1); yoptm(:,2); yoptm(:,3)];
% MSE=[mMSE(:,1); mMSE(:,2); mMSE(:,3)];
subplot(2,2,4);
plot(1:30000, MSE);
title('MSE');
% impulse=[1,zeros(1,size-1)]; 
% b=1;

%[h,t]=impz(aFinal,b);

% h=filter(aFinal,b,impulse);
% subplot(3,2,5);
% plot(h);
% title('Impulse response');

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
SignaldB1 = 10 * log10(var(yC));
NoisedB1 = 10 * log10(var(yC-input));
% SNR1 = SignaldB1 - NoisedB1
% snr1 = var(yC)/var(yC-input)
SignaldB2 = 10 * log10(var(yopt));
NoisedB2 = 10 * log10(var(yopt-input));
SNR2 = SignaldB2 - NoisedB2; % end of part 4
% snr2 = var(yopt)/var(yopt-input)
%%% coefficient plot %%%
figure
plot(1:30000, a(2,:),'r',1:30000, a(17,:),'b' ,1:30000, a(14,:),'m', 1:30000, a(8,:),'g')
hleg=legend('coeff = 2 ','coeff = 17','coeff = 14', 'coeff = 8');