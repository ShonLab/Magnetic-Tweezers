function [dzdat_corr,offset] = correctFEC(Fdat,dzdat,dxdat,Fcal,Rbead,Roff_cal,ori_cal)

    % correction
    F_ref = .1:.1:max(Fdat)+.1; % define array that covers tweezing F range 
    dx_ref = zeros(size(F_ref)); % array for storing dx positions at each force
    for Fi = 1:numel(F_ref)
        frange = (abs(Fdat-F_ref(Fi)) < .05) & (abs(dxdat) < 1e3); % indices of Fdat around force near F_ref(Fi)
        dx_ref(Fi) = mean(dxdat(frange)); % mean of dx at the times of Fdat near F_ref(Fi)
    end
    dx_ref = smooth(dx_ref,10);
    [~,idx] = closest(F_ref,Fcal); % get index in F_ref closest to Fcal
    dx_ref = dx_ref-dx_ref(idx);
    
    Roff_ref = Roff_cal - ori_cal*dx_ref; % array of dR's 
    zoff_ref = Rbead*(1 - sqrt(1 - (Roff_ref/Rbead).^2));
    
    dzdat_corr = dzdat + interp1(F_ref,zoff_ref,Fdat);
    offset = interp1(F_ref,zoff_ref,Fcal);
end