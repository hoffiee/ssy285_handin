


% =======================
% === Given constants ===
% =======================

% TRANSMISSION LINE 1
r = 0.3; 	% 0.3 ohm/km 
l = 1;		% 1 mH/km
c = 10;		% 10 nF/km
g = 0;		% 
l1 = 100;	% 100 km


% TRANSMISSION LINE 2
r = 0.15; 	% 0.15 ohm/km 
l = 1;		% 1 mH/km
c = 15;		% 15 nF/km
g = 0;		% 
l2 = 200;	% 200 km



% =========================
% ===== a)	===============
% =========================

% Calculations for the models.
% TRANSMISSION LINE 1
R_tl1 = r*l1
XL_tl1 = l*l1
XC_tl1 = c*l1/2

% TRANSMISSION LINE 2
R_tl2 = r*l2
XL_tl2 = l*l2
XC_tl2 = c*l2/2

% 			R_tli		XL_tli
% *--------=======------0000000--------*
% +		|							|  +
% Vi	= XC_tli			XC_tli	=  V0
% -		|							|  -
% *------------------------------------*
% 		i = 1,2

	

