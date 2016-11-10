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
Lt1 = L1*l1;
C1 = c1*l1/2;

% TRANSMISSION LINE 2
R2 = r2*l2;
Lt2 = L2*l2;
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

V = 10e3; 		% 10 kV
f = 50;			% 50 Hz


T = 0.5;
dt = 1/10000;
t = 0:dt:T;


e10 = sqrt(2)*V.*sin(2*pi*f*t);
e20 = sqrt(2)*V.*cos(2*pi*f*t);


% x = [I_l1 I_l2 I_0 v_0]
A = [-R1/Lt1 0 0 -1/Lt1;
    0 -R2/Lt2 0 -1/Lt2;
    0 0 -R_load/L_load 1/L_load;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0];


% u = [v_1 v_2]
B = [1/Lt1 0;
    0 1/Lt2;
    0 0;
    0 0];

C = [0 0 0 1];

D = [0 0];

u00 = [e10;
       e20];

% Steady state


xsteady= -A\(B*u00);
% x_s = @(t) A^-1*-B*u00(t);

u_star = [0 0 0 1]*xsteady;

% plot(t,u_star)
figure;
plot(t,xsteady)
title('Steady state','Fontsize',15,'Interpreter','Latex')
xlabel('$t$','Fontsize',15,'Interpreter','Latex')
leg = legend('$i_{11}$','$i_{21}$','$i_0$','$v_0$');
set(leg,'Fontsize',15,'Interpreter','Latex')


% =========================
% ===== d)	===============
% =========================

Atilde = 1000;
ftilde = 5;
u1 = Atilde*sin(2*pi*ftilde*t);
u2 = Atilde*cos(2*pi*ftilde*t);