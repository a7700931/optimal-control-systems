%% example4.1
clc;clear;close all;
A = [0 1;-2 -3]; % Ax+Bu
B = [0 1]'; % Ax+Bu
C = [1 0;0 1];  % y=Cx

F = [2 0;0 0];	% Performance index
Q = [2 0;0 0];  % Performance index
R = 0.004;   % Performance index

Pf = C'*F*C;
E = B*(R\B');
V = C'*Q*C;
W = C'*Q;
tf = 20;
T = 0:0.05:tf;
z = [1 0]'*ones(size(T)); % z need to become z(T)
Gf = C'*F*z(:,end);
Xic = [-0.5 0]'; % X 

[Tp,P] = ode45(@(t,p) pDREv(t,p,A,E,V),fliplr(T),Pf(:));
pict(A,Tp,P);

[Tg,G] = ode45(@(t,g) gDREv(t,g,A,E,Tp,P,W,z,T),fliplr(T),Gf);
pict(Gf,Tg,G);

[Tx,X] = ode45(@(t,x) xDREv(t,x,A,E,Tp,P,Tg,G),T,Xic);
pict(Xic,Tx,X,z);

pictu(A,B,R,P,X,T,G)

% picerror(C,T,X,z)

% ubook = -250*(flipud(P(:,2)).*X(:,1)+flipud(P(:,4)).*X(:,2)-flipud(G(:,2)));
% figure,plot(T,ubook);title('Book method');ylabel('U(t)');xlabel('t');

%% example4.2
clc;clear;close all;
A = [0 1;-2 -3];
B = [0 1]';
C = [1 0;0 1];

F = zeros(2);
Q = [2 0;0 0];
R = 0.04;

Pf = C'*F*C;
E = B*(R\B');
V = C'*Q*C;
W = C'*Q;
tf = 20;
T = 0:0.05:tf;

z = [2 0]'*T;
Gf = C'*F*z(:,end);
Xic = [-1 0]';

[Tp,P] = ode45(@(t,p) pDREv(t,p,A,E,V),fliplr(T),Pf(:));
pict(A,Tp,P);

[Tg,G] = ode45(@(t,g) gDREv(t,g,A,E,Tp,P,W,z,T),fliplr(T),Gf);
pict(Gf,Tg,G);

[Tx,X] = ode45(@(t,x) xDREv(t,x,A,E,Tp,P,Tg,G),T,Xic);
pict(Xic,Tx,X,z);

pictu(A,B,R,P,X,T,G)

picerror(C,T,X,z)

% ubook = -25*(flipud(P(:,2)).*X(:,1)+flipud(P(:,4)).*X(:,2)-flipud(G(:,2)));
% figure,plot(T,ubook);title('Book method');ylabel('U(t)');xlabel('t');

%% example4.5
clc;clear;close all;
A = [0 1;0 0];
B = [0 1]';
Q = [1 0;0 1];
R = 1;
[K,Pa,EV]=lqr(A,B,Q,R);

%% Problem4.1
clc;clear;close all;
A = [0 1;-2 -3];
B = [0 1]';
C = [1 0;0 1];

F = [2 0;0 0];
Q = [2 0;0 0];
R = 0.0004;

Pf = C'*F*C;
E = B*(R\B');
V = C'*Q*C;
W = C'*Q;
tf = 20;
T = 0:0.05:tf;

z = [2 0]'*T;
Gf = C'*F*z(:,end);
Xic = [-1 0]';

[Tp,P] = ode45(@(t,p) pDREv(t,p,A,E,V),fliplr(T),Pf(:));
pict(A,Tp,P);

[Tg,G] = ode45(@(t,g) gDREv(t,g,A,E,Tp,P,W,z,T),fliplr(T),Gf);
pict(Gf,Tg,G);

[Tx,X] = ode45(@(t,x) xDREv(t,x,A,E,Tp,P,Tg,G),T,Xic);
pict(Xic,Tx,X,z);

pictu(A,B,R,P,X,T,G);

% picerror(C,T,X,z)

%% Problem4.4
clc;clear;close all;
A = [0 1;-2 -3];
B = [1 0;0 1];
Q = [8 0;0 8];
R = [1 0;0 2];
tf = 10;
T = 0:0.05:tf;
Xic = [1 2]';
[K,Pa,EV]=lqr(A,B,Q,R);

[Tx,X] = ode45(@(t,x) p44(t,x,A,B,K),T,Xic);
figure,plot(Tx,X)

%% function for LQV
function dpdt = pDREv(~,p,A,E,V)
P = reshape(p,size(A));
Pdot = -P*A-A'*P+P*E*P-V;
dpdt = Pdot(:);
end

function dgdt = gDREv(t,g,A,E,Tp,P,W,z,T)
%%Gdot = -[A(t)-E(t)*P(t)]'*G(t)-W(t)*z(t)
P = interp1(Tp,P,t);
Pr = reshape(P,size(A));
z = interp1(T,z',t);
dgdt = -(A-E*Pr)'*g-W*z';
end

function dxdt = xDREv(t,x,A,E,Tp,P,Tg,G)
%Xdot = [A(t)-E(t)'*P(t)]*X(t)+E(t)*G(t)
P = interp1(Tp,P,t);
Pr = reshape(P,size(A));
G = interp1(Tg,G,t);
dxdt = (A-E'*Pr)*x+E*G';
end

function picerror(C,TData,X,z)
error = zeros(size(z));
for i=1:length(TData)
    error(:,i) = z(:,i)-C*X(i,:)';
end
figure,plot(TData,error);
title('Error'); xlabel('t'); ylabel('E(t)');
legend('E1(t)','E2(t)');
set(findall(0,'Type','legend'),'FontSize',11) % set all legend font size
end

function dxdt = p44(~,x,A,B,K)
dxdt = (A-B*K)*x;
end
