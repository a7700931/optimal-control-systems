% example 4.3
clc;clear;close all;
syms m(t)% b r a q x0 m real
a=2;
b=1;
q=2;
r=0.1;
eqn = diff(m,t) == a*m+m*a'+m*q*m-b*(r\b');
S=dsolve(eqn,m(10)==0);
t=[0 10];
figure,fplot(S,t);
title('M(t) dsolve solution')
ylim([0 1.5])

beta =sqrt(a^2+q*b^2/r);
tf = 10;
t = 0:0.05:tf;
M=zeros(1,numel(t));
for i=1:numel(t)
f1 = exp(-beta*(t(i)-tf));
f2 = exp(beta*(t(i)-tf));
f3 = (a+beta)*f1;
f4 = (a-beta)*f2;
M(1,i)=(b^2/r)*((f1-f2)/(f3-f4));
end
figure,plot(t,M);title('M(t) 4.3.37 solution');
ylim([0 1.5])

%% Problem 4.2
clc;clear;close all;

A = [0 1;-2 -4]; % Ax+Bu
B = [0 0.5]'; % Ax+Bu
C = [1 0;0 1]; % y=Cx

F = [0 0;0 0]; % Performance index
Q = [8 0;0 12];
R = 0.04;

E = B*(R\B');

tf = 10;
T = 0:0.05:tf;
TFlip = fliplr(T);
% Xic¡Ú0, Xic=X(t0), X is Forward Integral. Xic=X(tf), X is Backward Integral, T=TFlip.
Xic = [2 4];
% Mic=0, Mic=M(t0), M is Forward Integral. Mic=M(tf), M is Backward Integral, T=TFlip.
Mic = [0 0;0 0]; 
Vic = Xic;

[Tm,M] = ode45(@(t,m) miDRE(t,m,A,E,Q),T,Mic(:));
picFEP(Tm,M);

[Tv,V] = ode45(@(t,v) vDRE(t,v,A,Tm,M,Q),T,Vic(:)); % Only Case 3 need V
picFEP(Tv,V);

Xic = [7 3];

% X is stiff in final point so use stiff ode solver
[Tx,X] = ode23s(@(t,x) xDRE(t,x,A,B,E,Tm,M,Tv,V),TFlip,Xic(:));
picFEP(Tx,X);

%% Function
function dmdt = miDRE(~,m,A,E,Q)
M = reshape(m,size(A));
Mdot = A*M+M*A'+M*Q*M-E;
dmdt = Mdot(:);
end

function dvdt = vDRE(t,v,A,Tm,M,Q)
M = interp1(Tm,M,t);
Mr = reshape(M,size(A));
dvdt = Mr*Q*v+A*v;
end

function dxdt = xDRE(t,x,A,B,E,Tm,M,Tv,V)
if nargin <8
    Tv = Tm;
    V = zeros(size(Tv,1),size(B,1));
end
M = interp1(Tm,M,t);
Mr = reshape(M,size(A));
V = interp1(Tv,V,t);
dxdt = A*x-E*pinv(Mr)*(x-V'); % Mr has singular matrix, so use pinv
end

function picU(TData,A,B,R,Tm,M,Tx,X)
u = zeros(size(TData));
Mr = reshape(flipud(M)',size(A,1),size(A,1),[]);
for i=1:length(TData)
    u(:,i) = -(R\B')*pinv(Mr(:,:,i))*X(i,:)';
end
figure,plot(TData,u);
title('Optimal Control U(t)'); xlabel('t'); ylabel('U(t)');
legend('U(t)')  %,'color','none'
end

function picFEP(TData,Data)
figure,plot(TData,Data);
if TData(1) == 0
    title([inputname(2),'(t) ', 'Forward Integral'])
else
    title([inputname(2),'(t) ', 'Backward Integral'])
end
end