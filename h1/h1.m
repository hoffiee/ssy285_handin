clear all, clc, clf, close all, format compact


% TRANSMISSION LINE 1
r1 = 0.3; 	% 0.3 ohm/km 
l1 = 1e-3;		% 1 mH/km
c1 = 10e-9;		% 10 nF/km
g1 = 0;		% 
dist1 = 100;	% 100 km


% TRANSMISSION LINE 2
r2 = 0.15; 	% 0.15 ohm/km 
l2 = 1e-3;		% 1 mH/km
c2 = 15e-9;		% 15 nF/km
g2 = 0;		% 
dist2 = 200;	% 200 km


disp('===================')
disp('======== a) =======')
disp('===================')
% Calculations for the models.
% TRANSMISSION LINE 1
R1 = r1*dist1;
L1 = l1*dist1;
C1 = c1*dist1/2;


% TRANSMISSION LINE 2
R2 = r2*dist2;
L2 = l2*dist2;
C2 = c2*dist2/2;

disp('Transmission line 1')
disp(['R1: ', num2str(R1), '[ohm], L1: ', num2str(L1*1e3), '[mH], C1: ', num2str(C1*1e6), '[nF]'])
disp(['R1: ', num2str(R2), '[ohm], L1: ', num2str(L2*1e3), '[mH], C1: ', num2str(C2*1e6), '[nF]'])	

% disp('===================')
% disp('======== c) =======')
% disp('===================')
L_load = 1; 	% 1 H
R_load = 100;	% 100 ohm

V = 10000;
f = 50;

T = 0.5;
dt = 1/1000000;
t = 0:dt:T;

u10 = sqrt(2)*V.*sin(2*pi*f*t);
u20 = sqrt(2)*V.*cos(2*pi*f*t);

% x = [I_l1 I_l2 I_0 v_0]
A = [-R1/L1 0 0 -1/L1;
    0 -R2/L2 0 -1/L2;
    0 0 -R_load/L_load 1/L_load;
    1/(C1+C2) 1/(C1+C2) -1/(C1+C2) 0];

% u = [v_1 v_2]
B = [1/L1 0;
    0 1/L2;
    0 0;
    0 0];

C = [0 0 0 1];

D = [0 0];

u00 = [u10; u20];

xsteady = -A\(B*u00);

figure;
plot(t,xsteady)
title('Steady state','Fontsize',15,'Interpreter','Latex')
xlabel('$t$','Fontsize',15,'Interpreter','Latex')
leg = legend('$i_{11}$','$i_{21}$','$i_0$','$v_0$');
set(leg,'Fontsize',15,'Interpreter','Latex')
print('steadystate','-depsc')


% =========================
% ===== d)	===============
% =========================

Atilde = 1000;
ftilde = 5;
% u1 = Atilde*sin(2*pi*ftilde*t);
% u2 = Atilde*cos(2*pi*ftilde*t);

% v_0 = -A^-1*B*u00;

% u = [u1-u10;
% 	u2-u20];

