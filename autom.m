function [Rxx]=autom(x) 
% [Rxx]=autom(x) 
% This function Estimates the autocorrelation of the sequence of 
% random variables given in x as: Rxx(1), Rxx(2),…,Rxx(N), where N is 
% Number of samples in x. 
N=length(x); 
Rxx=zeros(1,N); 
x1=(x-mean(x))/std(x);
for m=1: N+1 
for n=1: N-m+1 
Rxx(m)=Rxx(m)+x1(n)*x1(n+m-1); 
end; 
end;
Rxx=Rxx/Rxx(1);
end