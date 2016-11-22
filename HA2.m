% =======================
% === Given constants ===
% =======================

% TRANSMISSION LINE 1
r1 = 0.3; 	% 0.3 ohm/km 
L1 = 1e-3;		% 1 mH/km
c1 = 10e-9;		% 10 nF/km
g1 = 0;		% 
l1 = 100;	% 100 km


% TRANSMISSION LINE 2
r2 = 0.15; 	% 0.15 ohm/km 
L2 = 1e-3;		% 1 mH/km
c2 = 15e-9;		% 15 nF/km
g2 = 0;		% 
l2 = 200;	% 200 km



% =========================
% ===== a)	===============
% =========================

% Calculations for the models.
% TRANSMISSION LINE 1
R1 = r1*l1;
L1 = L1*l1;
C1 = c1*l1/2;

% TRANSMISSION LINE 2
R2 = r2*l2;
L2 = L2*l2;
C2 = c2*l2/2;

% 			R_tli		XL_tli
% *--------=======------0000000--------*
% +		|							|  +
% Vi	= XC_tli			XC_tli	=  V0
% -		|							|  -
% *------------------------------------*
% 		i = 1,2

	
% =========================
% ===== b)	===============
% =========================



% =========================
% ===== c)	===============
% =========================
L_load = 1; 	% 1 H
R_load = 100;	% 100 ohm

V = 10000;
f = 50;

v10 = @(t) sqrt(2)*V.*sin(2*pi*f*t);
v20 = @(t) sqrt(2)*V.*cos(2*pi*f*t);

A = [-R1/L1 0 0 -1/L1;
    0 -R2/L2 0 -1/L2;
    0 0 -R_load/L_load 1/L_load;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0];

% Test
%A = [R1/L1 0 0 1/L1;
%    0 R2/L2 0 1/L2;
%    0 0 -R_load/L_load 1/L_load;
%    1/(C1+C2) 1/(C1+C2) 1/(C1+C2) 0];


B = [1/L1 0;
    0 1/L2;
    0 0;
    0 0];

C = [0 0 0 1];

D = [0, 0];

v00 = @(t) [v10(t);
            v20(t)];

% Steady state
xsteady = @(t) A\(-B*v00(t));

v0steady = xsteady(0);
v0steady = v0steady(4);

t = linspace(0,0.1,10000);

% Test
%t = linspace(-0.01,0.1,10000);

figure;
plot(t,xsteady(t))
title('Steady state','Fontsize',15,'Interpreter','Latex')
xlabel('$t$','Fontsize',15,'Interpreter','Latex')
leg = legend('$i_{11}$','$i_{21}$','$i_0$','$v_0$');
set(leg,'Fontsize',15,'Interpreter','Latex')

v0phase = -0.002503;
v0amp = 8696;

vstar = @(t) v0amp.*cos(2*pi*f*(t+v0phase));
hold on
plot(t,vstar(t),'--')

%%

s = tf('s');

G = C*(s*eye(4)-A)^(-1)*B+D;
G = C/(s*eye(4)-A)*B+D;

[z,p,~] = zpkdata(G);

%%
%%%%% HA2 %%%%%
clc

n = 4;

contr = ctrb(A,B);
obs = obsv(A,C);

rankContr = rank(contr);
rankObs = rank(obs);

if rankContr == n
    disp('Controllable!')
else
    disp('Uncontrollable!')
end
if rankObs == n
    disp('Observable!')
else
    disp('Unobservable')
end

%% Task 1 sweep
clc;
Afunc = @(R,L)[-R1/L1 0 0 -1/L1;
    0 -R2/L2 0 -1/L2;
    0 0 -R/L 1/L;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0]; 

Rloadvec = 10:20:110;
Lloadvec = 0.1:0.3:2;

checkContr = zeros(size(Rloadvec,2),size(Lloadvec,2));
checkObs = zeros(size(Rloadvec,2),size(Lloadvec,2));

checkStab = zeros(size(Rloadvec,2),size(Lloadvec,2));
checkDet = zeros(size(Rloadvec,2),size(Lloadvec,2));

for Ri = 1:size(Rloadvec,2)
    for Li = 1:size(Lloadvec,2)
        Atmp = Afunc(Rloadvec(Ri),Lloadvec(Li));
        
        % Check controllability and observability
        contr = ctrb(A,B);
        obs = obsv(A,C);
        
        % If full rank
        checkContr(Ri,Li) = (n == rank(contr));
        checkObs(Ri,Li) = (n == rank(obs));
        
        % Check stabilizable and detactable
        lambda = eig(A);
        lambdaI = [lambda(1) 0 0 0;
                    0 lambda(2) 0 0;
                    0 0 lambda(3) 0;
                    0 0 0 lambda(4)];
        % [A-lambda*I B] has full row rank for all complex lambda with Re(?) > 0
        checkStab(Ri,Li) = (n == rank([A-lambdaI B]));
        % The (A, C) pair is called detectable, if the non-observable part
        % is asymptotically stable
        checkDet(Ri,Li) = (lambda(1) < 0 && lambda(2) < 0 && lambda(3) < 0 && lambda(4) < 0);
    end
end

%% Our values

clc

n = 4;

contr = ctrb(A,B);
obs = obsv(A,C);

rankContr = rank(contr);
rankObs = rank(obs);

% Check stabilizable and detactable
lambda = eig(A);
% [A- lambda*I B] has full row rank for all complex lambda with Re(lambda) > 0
checkStab = (n == rank([A-repmat(lambda',4,1)*eye(n) B]));
% The (A, C) pair is called detectable, if the non-observable part
% is asymptotically stable
checkDet = (lambda(1) < 0 && lambda(2) < 0 && lambda(3) < 0 && lambda(4) < 0);

if rankContr == n
    disp('Controllable!')
else
    disp('Uncontrollable!')
end
if rankObs == n
    disp('Observable!')
else
    disp('Unobservable')
end
if checkStab == 1
    disp('Stabalizable!')
else
    disp('Unstabalizable!')
end
if checkDet == 1
    disp('Detectable!')
else
    disp('Undetectable!')
end

%% d/e)

% Diagnolize A
[Tinv,VD] = eig(A);
Ac = Tinv*VD/Tinv;

Ts = 5e-3;

% Discretize
Ad = expm(Ac*Ts);
Bd = inv(Ac)*(Ad-eye(4))*B;

%% f)

absEig = abs(eig(Ad));

% Test if this yields the same result as c2d
[b1,a1] = ss2tf(Ad,Bd,C,D,1);
[b2,a2] = ss2tf(Ad,Bd,C,D,2);

% Generate Gd
Gdmatlab = c2d(G,Ts,'zoh');

% Find poles and zeros
[zd, pd, ~] = zpkdata(Gdmatlab);

pd1 = abs(pd{1});
pd2 = abs(pd{2});
