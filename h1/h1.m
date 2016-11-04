% =======================
% === Given constants ===
% =======================

% TRANSMISSION LINE 1
r1 = 0.3; 	% 0.3 ohm/km 
l1 = 1;		% 1 mH/km
c1 = 10;		% 10 nF/km
g1 = 0;		% 
l1 = 100;	% 100 km


% TRANSMISSION LINE 2
r2 = 0.15; 	% 0.15 ohm/km 
l2 = 1;		% 1 mH/km
c2 = 15;		% 15 nF/km
g2 = 0;		% 
l2 = 200;	% 200 km



% =========================
% ===== a)	===============
% =========================

% Calculations for the models.
% TRANSMISSION LINE 1
R1 = r1*l1;
L1 = l1*l1;
C1 = c1*l1/2;

% TRANSMISSION LINE 2
R2 = r2*l2;
L2 = l2*l2;
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

v10 = @(t) sqrt(2)*V.*sind(2*pi*f*t);
v20 = @(t) sqrt(2)*V.*cosd(2*pi*f*t);


% x = [I_l1 I_l2 I_0 v_0]
A = [-R1/L1 0 0 -1/L1;
    0 -R2/L2 0 -1/l2;
    0 0 -R_load/L_load 1/L_load;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0];


% u = [v_1 v_2]
B = [1/L1 0;
    0 1/L2;
    0 0;
    0 0];

v00 = @(t) [v10(t);
            v20(t)];

% Steady state
xsteady = @(t) A\(-B*v00(t));

t = 0:1/(2*f):10;

figure;
plot(t,xsteady(t))
title('Steady state','Fontsize',15,'Interpreter','Latex')
xlabel('$t$','Fontsize',15,'Interpreter','Latex')
leg = legend('$i_{11}$','$i_{21}$','$i_0$','$v_0$');
set(leg,'Fontsize',15,'Interpreter','Latex')


% =========================
% ===== d)	===============
% =========================

Atilde = 1000;
ftilde = 5;
u1 = @(t) Atilde*sind(2*pi*ftilde*t);
u2 = @(t) Atilde*cosd(2*pi*ftilde*t);

u = @(t) []

