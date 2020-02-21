%% example3.1
clc;clear;close all;

A = [0 1;-2 1]; % Ax+Bu
B = [0 1]'; % Ax+Bu
F = [1 0.5;0.5 2]; % Performance index
Q = [2 3;3 5];
R = 0.25;
Pf = F; % P Final Condition
x0 = [2 -3]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:)); % P is backward integral
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%% example3.2
clc;clear;close all;

A = [0 1;-2 1]; % Ax+Bu
B = [0 1]'; % Ax+Bu
F = [0 0;0 0]; % Performance index
Q = [2 3;3 5];
R = 0.25;
Pf = F; % P Final Condition
x0 = [2 -3]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:)); % P is backward integral
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%%Algebraic Riccati Equation
[~,P,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*P(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% example3.3
clc;clear;close all;
%%Open Loop
A = -3; % Ax+Bu
B = 1; % Ax+Bu
F = 0; % Performance index
Q = 1;
R = 1;
Pf = F; % P Final Condition
x0 = 1;  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X)

pictu(A,B,R,P,X,T)

olu = -(sqrt(10)-3)*exp(-sqrt(10)*T);
figure,plot(T,olu); title('Open Loop U(t)');

%%Close loop
A = -3; % Ax+Bu
B = 1; % Ax+Bu
F = 0; % Performance index
Q = 2;
R = 2;
Pf = F; % P Final Condition
x0 = 1;  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X)

pictu(A,B,R,P,X,T)

%% Problem3.1
clc;clear;close all;
A = 1; % Ax+Bu
B = 1; % Ax+Bu
F = 0; % Performance index
Q = 4;
R = 0.5;
Pf = F; % P Final Condition
x0 = 2;  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%%Algebraic Riccati Equation
[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% Problem3.2
clc;clear;close all;
A = [0 1;-1 0]; % Ax+Bu
B = [0 1]'; % Ax+Bu
F = zeros(2); % Performance index
Q = [4 0;0 0];
R = 0.5;
Pf = F; % P Final Condition
x0 = [0 1]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:)); % P is backward integral
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%%Algebraic Riccati Equation
[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% Problem3.4
clc;clear;close all;
A = [0 1;-2 4]; % Ax+Bu
B = [0 5]'; % Ax+Bu
F = [1 0;0 1]; % Performance index
Q = [5 0;0 2];
R = 4;
Pf = F; % P Final Condition
x0 = [1 2]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%% Problem3.5
clc;clear;close all;
A = [0 1;-2 -3]; % Ax+Bu
B = [0 1]'; % Ax+Bu
F = [0 0;0 0]; % Performance index
Q = [2 0;0 2];
R = 2;
Pf = F; % P Final Condition
x0 = [1 2]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%%Algebraic Riccati Equation
[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% Problem3.6
clc;clear;close all;
A = [0 1;-1 -1]; % Ax+Bu
B = [1 1]'; % Ax+Bu
F = [0 0;0 0]; % Performance index
Q = [4 0;0 8];
R = 1;
Pf = F; % P Final Condition
x0 = [1 2]';  % X Initial Condition
tf = 10; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T,x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

%%Algebraic Riccati Equation
[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% Problem3.8.1
clc;clear;close all;
A = [0 1 0;0 0 1;-5 -7 -10]; % Ax+Bu
B = [0 0 4]'; % Ax+Bu
F = zeros(3); % F size need as same as A ,no matter it's zero or not
Q = [1 0 0;0 1 0;0 0 1];
R = 1;
Pf = F; % P Final Condition
x0 = [1 2 3]';  % X Initial Condition
tf = 20; % Finial time
T = 0:0.05:tf;
[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P) % Picture of P(t)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T, x0);
pict(x0,Tx,X) % Picture of X(t)

pictu(A,B,R,P,X,T)

[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X) % Picture of X(t)
pictu(A,B,R,Pa,X,T)

%% Problem3.8.2
Q = [10 0 0;0 1 0;0 0 1];

[~,Pa,~]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X)
pictu(A,B,R,Pa,X,T)

%% Problem3.8.3
Q = [1 0 0;0 1 0;0 0 1];
R = 10;
[K,Pa,EV]=lqr(A,B,Q,R);
Pa = ones(size(T'))*Pa(:)';
[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,T,Pa),T,x0);
pict(x0,Tx,X)
pictu(A,B,R,Pa,X,T)

%% Problem3.9
clc;clear;close all;
A = [0 1;1 1]; % Ax+Bu
B = [1 1;0 1]; % Ax+Bu
F = zeros(2);% F siz(); e need as same as A ,no matter it's zero or not
Q = [4 0;0 8];
R = [1 0;0 0.5];
Pf = F; % P Final Condition
x0 = [1 2]';  % X Initial Condition
tf = 5; % Finial time
T = 0:0.05:tf;

[Tp,P] = ode45(@(t,p) pDRE(t,p,A,B,Q,R),fliplr(T),Pf(:));
pict(A,Tp,P)

[Tx,X] = ode45(@(t,x) xDRE(t,x,A,B,R,Tp,P),T, x0);
pict(x0,Tx,X)

pictu(A,B,R,P,X,T)

%% function for LQR
function dpdt = pDRE(~,p,A,B,Q,R)
P = reshape(p,size(A));  % Reshape input p into matrix
Pdot = -P*A-A'*P+P*B*(R\B')*P-Q;
dpdt = Pdot(:);  % Reshape output as a column vector
end

function dxdt = xDRE(t,x,A,B,R,Tp,P)
P = interp1(Tp,P,t);
Pr = reshape(P,[size(A,1) size(A,1)]); % Reshape p like A matrix
dxdt = (A-B*(R\B')*Pr)*x;
end
