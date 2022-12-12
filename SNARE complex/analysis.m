%% load data
% set parameters
Rbead = 1400; % MB radius in nm
T = 300; % temperature
pixel_size = 80; % in nm
correction_factor = 0.878; % correction factor for refractive index mismatch (n_water/n_oil = 0.878)
F0 = 0.05759; z0 = 20.94; A1 = 91.75; d1 = 0.9161; A2 =  53.88; d2 = 2.647;

% read data
folders = dir('raw*');
pths = arrayfun(@(i) {[folders(i).name,'/']}, 1:numel(folders)); npth = numel(pths); %pths =raw*
[finfo,fname,nbead,nframe,fps,Roff,ori,f,t,t2,d,R,P,rz,rx,ry,x,y,z,dx,dy,dz] = deal(cell(npth,1));
nfile = zeros(npth,1);
for p = 1:npth
    disp(p);
    finfo{p} = dir([pths{p},'r*.xls']); % data files
    nfile(p) = numel(finfo{p});
    [fname{p},Roff{p},ori{p},f{p},d{p},R{p},P{p},x{p},y{p},z{p},dx{p},dy{p},dz{p}] = deal(cell(nfile(p),1));
    [nbead{p},fps{p},nframe{p}] = deal(zeros(nfile(p),1));
    for n = 1:nfile(p)
        disp([int2str(n/nfile(p)*100),'% of ',pths{p}(1:end-1),'...']);
        fname{p}{n} = finfo{p}(n).name;
        
        % read calibration info (c###.xls)
        calinfo = dlmread([pths{p},'c',fname{p}{n}(2:4),'.fps']);
        fps{p}(n) = calinfo(1); % sampling rate in Hz (1200 here)
        Roff{p}{n} = calinfo(2,:); % tether offset (see: 10.1126/sciadv.aav1697)
        ori{p}{n} = calinfo(3,:); % +1/-1 for left-/right-tilted MB, respectively
        
        % read motor data (s###-###.xls)
        dat = dlmread([pths{p},regexprep(fname{p}{n},'r','s')]);
        t2{p}{n} = dat(:,1); % timestamp for motor data in seconds
        M{p}{n} = dat(:,2); % magnet distance in mm
        F{p}{n} = F0 + A1*exp(-(z0-M{p}{n})/d1) + A2*exp(-(z0-M{p}{n})/d2);
        R{p}{n} = (dat(:,3)-floor(dat(:,3)))*360; % magnet orientation in degrees
        P{p}{n} = dat(:,4); % piezo position in μm
        
        % read bead data (r###-###.xls)
        dat = dlmread([pths{p},fname{p}{n}]);
        nframe{p}(n) = size(dat,1); % # of data points
        f{p}{n} = dat(:,1); % frame number
        t{p}{n} = f{p}{n}/fps{p}(n); % timestamp for bead data in seconds
        dat = dat(:,2:end); % remove the first column
        dat(:,[1:3:end,2:3:end]) = dat(:,[1:3:end,2:3:end]) - repmat(mean(dat(31:60,[1:3:end,2:3:end]),1),[nframe{p}(n),1]); % subtract xy offset
        nbead{p}(n) = size(dat,2)/3-1; % # of MB (one is for RB)
        
        % bead position
        rx{p}{n} = dat(:,1)*pixel_size;
        ry{p}{n} = dat(:,2)*pixel_size;
        rz{p}{n} = dat(:,3);
        x{p}{n} = dat(:,4:3:end)*pixel_size;
        y{p}{n} = dat(:,5:3:end)*pixel_size;
        z{p}{n} = dat(:,6:3:end);
        
        % MB position relative to RB
        dx{p}{n} = (x{p}{n}-repmat(rx{p}{n},[1,nbead{p}(n)]));
        dy{p}{n} = (y{p}{n}-repmat(ry{p}{n},[1,nbead{p}(n)]));
        dz{p}{n} = (z{p}{n}-repmat(rz{p}{n},[1,nbead{p}(n)]))*correction_factor;
      
        % synchronize motor data to bead data
        F{p}{n} = interp1(t2{p}{n},F{p}{n},t{p}{n});
        R{p}{n} = interp1(t2{p}{n},R{p}{n},t{p}{n});
        P{p}{n} = interp1(t2{p}{n},P{p}{n},t{p}{n});       
    end
end
clear dat;
save('analysis');

%% force ramp
% force-extension model
Fdat_ref = (0:.1:20)';
Lo_dsDNA = (1020)*.338; Lp_dsDNA = 37; Ko_dsDNA = 400; % dsDNA handle
zdat_model_dsDNA = eWLC_inv(Fdat_ref,Lo_dsDNA,Lp_dsDNA,T,Ko_dsDNA,1);
Lo_PEG = 570*(1/80); Lp_PEG = 0.47; % bPEG (1 kD)
zdat_model_PEG = WLC_inv(Fdat_ref,Lo_PEG,Lp_PEG,T,1);
dzdat_model_FZ = zdat_model_PEG + zdat_model_dsDNA + 2;
Lp_PP = .6; nLc_PP = .4;
nLc_helix = .15; % helical rise per AA
nAA_frayed_SB = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX = 14; % 14, 21, 28 for +4, +2, 0 layer
nAA_frayed_SB2 = 53; % complete
nAA_frayed_SX2 = 28; % 14, 21, 28 for +4, +2, 0 layer
dzdat_model_LO = dzdat_model_FZ + WLC_inv(Fdat_ref,nLc_PP*23,Lp_PP,T,1); % 12 and 11 from SB and SX
dzdat_model_HZ = dzdat_model_LO - 2 + WLC_inv(Fdat_ref,nLc_PP*(nAA_frayed_SB+nAA_frayed_SX),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB-nAA_frayed_SX)/nLc_helix));
dzdat_model_UZ = dzdat_model_LO - 2 + WLC_inv(Fdat_ref,nLc_PP*(nAA_frayed_SB2+nAA_frayed_SX2),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB2-nAA_frayed_SX2)/nLc_helix));
dzdat_model_UF = dzdat_model_FZ - 2 + WLC_inv(Fdat_ref,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

h = figure(1); clf; h.WindowState = 'maximized';
set(gcf,'defaultaxesfontsize',12);
p = 1;
for n = 1:4
    frange = 1:nframe{p}(n);
    Fdat = F{p}{n}(frange);
    dxdat = dx{p}{n}(frange);
    dzdat = dz{p}{n}(frange);
    dzdat_corr = correctFEC(Fdat,dzdat,dxdat,5,Rbead,Roff{p}{n},ori{p}{n});
    
    F_clamp = 10;
    frange = abs(Fdat(1:end-1)-F_clamp) < .1 & diff(Fdat) > 0;
    offset = dzdat_model_FZ(Fdat_ref==F_clamp);
    dzdat_corr = dzdat_corr - mean(dzdat_corr(frange)) + offset;
    
    idx1 = find(diff(Fdat) > 1e-4,1,'first');
    idx2 = find(Fdat == max(Fdat));
    idx3 = find(diff(Fdat) < -1e-4,1,'last');
    for pid = 1:2
        subplot(1,2,pid);
        frange_sub = idx2:idx3;
        plot(dzdat_corr(frange_sub),Fdat(frange_sub),'b-','linewidth',1); hold all;
        frange_sub = idx1:idx2;
        plot(dzdat_corr(frange_sub),Fdat(frange_sub),'k-','linewidth',1);
    end
end
subplot(121);
xlim([200,400]); ylim([0,20]); grid on;
xlabel('Extension (nm)'); ylabel('Force (pN)');
plot(dzdat_model_FZ,Fdat_ref,'color',[0,.7,1]);
plot(dzdat_model_UF,Fdat_ref,'color',[1,0,0]);

subplot(122);
xlim([340,390]); ylim([12,17]); grid on;
xlabel('Extension (nm)'); ylabel('Force (pN)');
plot(dzdat_model_FZ,Fdat_ref,'color',[0,.7,1]);
plot(dzdat_model_LO,Fdat_ref,'color',[0,.5,0]);
plot(dzdat_model_HZ,Fdat_ref,'color',[.2,.8,.2]);
plot(dzdat_model_UZ,Fdat_ref,'color',[1,.8,0]);
plot(dzdat_model_UF,Fdat_ref,'color',[1,0,0]);

saveas(gcf,'SNARE complex_force ramp.fig');

%% force jump
h = figure(2); clf; h.WindowState = 'maximized';
set(gcf,'defaultaxesfontsize',12);

p = 2; n = 1; F_meas = 14.3;
frange = 1:nframe{p}(n);
tdat = t{p}{n}(frange); tdat = tdat-tdat(1);
Fdat = F{p}{n}(frange);
dzdat = dz{p}{n}(frange);
dzdat = dzdat - mean(dzdat(end-1200:end)) + dzdat_model_FZ(101);
dzdat_filt = medfilt1(dzdat,4); dzdat_filt([1,end]) = nan;

subplot(411);
plot(tdat,Fdat,'linewidth',1);
ylabel('Force (pN)');

subplot(4,1,2:4);
plot(tdat,dzdat_filt,'linewidth',.1); hold all;
idx_F = find(Fdat_ref==F_meas);
frange = find(abs(Fdat-F_meas) < .1);
plot(tdat(frange([1,end])),dzdat_model_FZ(idx_F)*[1,1],'--','color',[0,.7,1]);
plot(tdat(frange([1,end])),dzdat_model_LO(idx_F)*[1,1],'--','color',[0,.5,0]);
plot(tdat(frange([1,end])),dzdat_model_HZ(idx_F)*[1,1],'--','color',[.2,.8,.2]);
plot(tdat(frange([1,end])),dzdat_model_UZ(idx_F)*[1,1],'--','color',[1,.8,0]);
plot(tdat(frange([1,end])),dzdat_model_UF(idx_F)*[1,1],'--','color',[1,0,0]);
xlabel('Time (s)'); ylabel('Extension (nm)');
plot(3.9+[0,0],333+[0,5],'k','linewidth',2); text(3.9-.1,333+8,'refolding','fontsize',12);
plot(19.5+[0,0],382+[0,5],'k','linewidth',2); text(19.5-.1,382+8,'unzipping','fontsize',12);
plot(20.5+[0,.5],363+[0,0],'k','linewidth',2); text(20.5+.6,363,'rezipping','fontsize',12);

saveas(gcf,'SNARE complex_force clamp.fig');
