function z = eWLC_inv(F,Lo,Lp,T,Ko,corr)
z = zeros(size(F));
l0 = arrayfun(@(j) WLC_inv(F(j),Lo,Lp,T,false),1:numel(F))/Lo;
if corr == 0
    % without correction
    for j = 1:numel(F)    
        z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko))-F(j),l0(j))*Lo;
    end
elseif corr == 1
    % with correction
    a = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718];
    for j = 1:numel(F)
        z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)+sum(a.*l.^(2:7)))-F(j),l0(j))*Lo;
    end
elseif corr == 2
    % with correction by Ogden et al.
    for j = 1:numel(F)    
        z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)-3/4*(l-F(j)/Ko).^2)-F(j),l0(j))*Lo;
    end
end
% [~,idx] = closest(F,5); z_out = z-z(idx);