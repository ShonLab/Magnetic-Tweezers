function [f_low,f_high] = f_cutoff_double(F,R)
% See: 10.1016/j.bpj.2015.04.011
% y-axis is perpendicular to the magnetic field in our setup. 

% Model PSD check
% R = 1400; % bead radius
Lo = 5e3*.338; Lp = 45; Ko = 1000; % parameters for 5 kbp dsDNA
fs = 1200; % sampling rate
kB = 1.38e-2; T = 300;
eta = 8.9e-4*1e-6; % conversion from N*m to pN*nm for Pa*s

L = eWLC_inv(F,Lo,Lp,T,Ko,1); Lr = L/R;
C_par = 1./(1 - 9/16*(1+Lr)^(-1) + 1/8*(1+Lr)^(-3) - 45/256*(1+Lr)^(-4) - 1/16*(1+Lr)^(-5));
C_rot = 1 + 5/16*(1+Lr)^(-3);

gamma_y = 6*pi*eta*R*C_par;
gamma_phi = 8*pi*eta*R^3*C_rot;

f_low = F/L/(2*pi)*( (L+R)*R/(2*gamma_phi) + 1/(2*gamma_y) - 1/2*sqrt(( (L+R)*R/gamma_phi + 1/gamma_y )^2 - 4*L*R/(gamma_y*gamma_phi)) );
f_high = F/L/(2*pi)*( (L+R)*R/(2*gamma_phi) + 1/(2*gamma_y) + 1/2*sqrt(( (L+R)*R/gamma_phi + 1/gamma_y )^2 - 4*L*R/(gamma_y*gamma_phi)) );

end