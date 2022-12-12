function z = mFJC(F,Lo,Lp,T,Ko)
% Lp = Kuhn length/2
kB = 1.38e-2;
z = Lo*(coth(2*F*Lp/kB/T)-(kB*T/2/Lp./F)).*(1+F/Ko);
% [~,idx] = closest(F,5); z_out = z-z(idx);
