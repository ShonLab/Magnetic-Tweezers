function z = WLC_inv(F,Lo,Lp,T,corr)
z = zeros(size(F));
% initial guess
l0 = F*Lp/1.38e-2/T;
l0_high = 1-1./sqrt(4*F*Lp/1.38e-2/T);
l0(l0>.5) = l0_high(l0>.5);
if corr == 0
    % without correction
    for j = 1:numel(F)
        z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-l).^2-1/4+l)-F(j),l0(j))*Lo;
    end
elseif corr == 1
    % with correction
    a = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718];
    for j = 1:numel(F)
        z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-l).^2-1/4+l+sum(a.*l.^(2:7)))-F(j),l0(j))*Lo;
    end
elseif corr == 2
    % with correction by Ogden et al.
    for j = 1:numel(F)
       z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-l).^2-1/4+l-3/4*l.^2)-F(j),l0(j))*Lo;
    end
end
% [~,idx] = closest(F,5); z_out = z-z(idx);