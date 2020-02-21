%% Example6.3
clc;clear;close all;
% syms p(t)
s = dsolve('Dp-p^2-4*p+1=0','p(5)=1');
S = matlabFunction(s);
t=0:0.01:5;
figure,plot(t,real(S(t)));

%%
p_tf = 5;
p_t = 0:0.01:p_tf;

num =(sqrt(5)-2)+(sqrt(5)+2)*(3-sqrt(5))/(3+sqrt(5))*exp(2*sqrt(5)*(p_t-p_tf));
den = 1-((3-sqrt(5))/(3+sqrt(5)))*exp(2*sqrt(5)*(p_t-p_tf));
p = num./den;
figure,plot(p_t,p)


