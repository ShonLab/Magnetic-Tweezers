%% Force calibraion of magnetic tweezers (Shon Lab @ POSTECH Physics, 20221120)
% Email MJS (mjshon@postech.ac.kr) for any inquiries. 
%
% This code performs force calibration on the sample data in the folder "data". 
%
% r008-001.xls: bead-coordinate data measured at 1.2 kHz
% Rows: data over time
% Column 1: frame number
% Columns 2–4: x,y,z-coordinates of reference bead (RB)
% Columns 5–7: x,y,z-coordinates of magnetic bead (MB)
% x,y values are in pixels and z values are in nm. 
%
% s008-001.xls: motor readings recorded at 50-ms intervals
% Rows: data over time
% Column 1: time stamp (in seconds)
% Column 2: % magnet position (in mm; the larger the closer)
% Column 3: % stepper motor reading for magnet rotation (in deg, usually unchanged)
% Column 4: % piezo position (in μm)
%
% c008.fps: measurement information
% Line 1: frame rate (in Hz)
% Line 2: tether offset (in nm)
% Line 3: bead tilting orientation; +1/-1 for left-/right-tilted MB, respectively
% See: https://advances.sciencemag.org/content/5/6/eaav1697 (doi: 10.1126/sciadv.aav1697)
%
% Please read the original article for more information: 
% Park et al., High-speed magnetic tweezers for nano-mechanical measurements on force-sensitive elements, JoVE (2022)

%% load data
% set parameters
Rbead = 1400; % MB radius in nm
T = 300; % temperature
pixel_size = 80; % in nm
correction_factor = 0.878; % correction factor for refractive index mismatch (n_water/n_oil = 0.878)

% read calibration info
calinfo = dlmread('data\c008.fps');
fps = calinfo(1); % frame rate
Roff = calinfo(2,:); % tether offset
ori = calinfo(3,:); % bead tilting orientation

% read motor data
dat = dlmread('data\s008-001.xls');
t2 = dat(:,1); % timestamp for motor data
d = 20 - dat(:,2); % magnet distance
R = (dat(:,3)-floor(dat(:,3)))*360; % magnet orientation
P = dat(:,4); % piezo position

% read bead data
dat = dlmread('data\r008-001.xls');
nframe = size(dat,1); % number of frames
f = dat(:,1); % frame index
t = f/fps; % timestamp for bead data
dat = dat(:,2:end); % remove the first column
nbead = size(dat,2)/3-1; % # of MB (one is for RB)

% bead position
rx = dat(:,1)*pixel_size;
ry = dat(:,2)*pixel_size;
rz = dat(:,3);
x = dat(:,4:3:end)*pixel_size;
y = dat(:,5:3:end)*pixel_size;
z = dat(:,6:3:end);

% MB position relative to RB
dx = (x-repmat(rx,[1,nbead]));
dy = (y-repmat(ry,[1,nbead]));
dz = (z-repmat(rz,[1,nbead]))*correction_factor;

% synchronize motor data to bead data
d = interp1(t2,d,t);
R = interp1(t2,R,t);
P = interp1(t2,P,t);       
clear dat;
save('analysis');

%% identify regions to analyze
dx = dx - mean(dx(1:1200)); % subtract baseline
dy = dy - mean(dy(1:1200)); % subtract baseline
dz = dz - min(dz(:)); % subtract baseline
nframe_1s = 1*fps; % number of frames in 1 s
dd = medfilt1(diff(d), nframe_1s); %  of Δd
S = regionprops(abs(dd) < 1e-5,'PIxelIdxList','Area'); % regions with minimal motor motion
S = S([S.Area] > 10*fps); % regions longer than 10 s only
f_reg = zeros(numel(S),2); % frame indices (first and last) to analyze
d_reg = zeros(numel(S),1); % average d values in the regions
for i = 1:numel(S)
    f_reg(i,:) = S(i).PixelIdxList([nframe_1s+1,end-nframe_1s]); % exclude some frames (1 s) for stabilization
    d_reg(i) = mean(d(f_reg(i,:)));
end
 
h = figure(1); clf; h.WindowState = 'maximized';
set(gcf,'defaultaxesfontsize',12);
for dir = 1:3
    subplot(3,1,dir);
    yyaxis left;
    plot(t,d); hold all;
    for i = 1:numel(S)
        frange = f_reg(i,1):f_reg(i,2);
        plot(t(frange),d(frange),'r-');
    end
    ylabel('Magnet distance (mm)');

    yyaxis right;
    if dir == 1
        plot(t,dx);
        ylim([-1500,1500]);
        ylabel('x_{MB} (nm)');
    elseif dir == 2
        plot(t,dy);
        ylim([-1500,1500]);
        ylabel('y_{MB} (nm)');
    else
        plot(t,dz);
        ylabel('z_{MB} (nm)');
    end
    xlabel('Time (s)');
    title([num2str(numel(S)),' regions to analyze']);
end
saveas(gcf,'overall data.fig');

%% calculate PSD and fit models (doi: 10.1016/j.bpj.2015.04.011)
h = figure(2); clf; h.WindowState = 'maximized';
set(gcf,'defaultaxesfontsize',12);

PSDfit_coupled = fittype('PSD_coupled_y(F,f)','independent','f','coefficient','F');
nreg = size(f_reg,1);
F_PSD = zeros(nreg,1);
for i = 1:nreg
    frange = f_reg(i,1):f_reg(i,2);
    if mod(numel(frange),2) == 1
        frange = frange(1:end-1);
    end
    N = numel(frange);
    N_sub = floor(N/5);
    if mod(N_sub,2) == 1
        N_sub = N_sub-1;
    end
    fdat = (0:fps/N_sub:fps/2)';
    PSD = zeros(N_sub/2+1,1);
    for j = 1:5
        frange_sub = frange((j-1)*N_sub + (1:N_sub));
        dy_sub = dy(frange_sub);
        dft = fft(dy_sub);
        dft = dft(1:N_sub/2+1);
        PSD_tmp = (1/(fps*N_sub)) * abs(dft).^2;
        PSD_tmp(2:end-1) = 2*PSD_tmp(2:end-1);
        PSD = PSD + PSD_tmp;
    end
    PSD = PSD/5;

    % fitting
    sel = fdat > 1; % ignore noisy, low-frequency data under 1 Hz
    fitres = fit(fdat(sel),PSD(sel),PSDfit_coupled,'start',10,'upper',100,'lower',0);
    [f_low,f_high] = f_cutoff_double(fitres.F,Rbead);
    flim = fdat(sel); flim = flim([1,end]);

    % plot results        
    subplot(5,7,i);
    loglog(fdat,PSD,'linewidth',.1); hold all;
    loglog(fdat,smooth(PSD,7),'linewidth',.1);
    loglog(fdat,fitres(fdat),'y-');
    xlim([1,1e3]); ylim([1e-3,1e6]);
    xline(flim,'k-','linewidth',.5);
    xline(f_low,'b-','linewidth',1);
    xline(f_high,'r-','linewidth',1);
    text(800,1e5,['F = ',num2str(fitres.F,'%.1f'),' pN'],'hori','right');
    F_PSD(i) = fitres.F;
    if i == 1
        xlabel('Frequency (Hz)'); ylabel('PSD_y (nm^2 Hz^{−1})');
    end
end
saveas(gcf,'PSD.fig');

%% obtain calibration curve
% double-exponential fitting
exp2fit = fittype('F0 + A1*exp(-d/d1) + A2*exp(-d/d2)','independent','d','coefficient',{'F0','d1','d2','A1','A2'});

h = figure(3); clf; h.WindowState = 'maximized';
pid_list = [1,2,5,6];
fitres = fit(d_reg,F_PSD,exp2fit,'start',[0,1,5,20,20])
plot(d_reg,F_PSD,'o','markersize',10); hold all;
plot(fitres,'-');
xlim([0,20]); ylim([0,30]);
xlabel('Magnet distance (mm)'); ylabel('F from PSD (pN)');
saveas(gcf,'calibration result.fig');
