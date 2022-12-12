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
[finfo,fname,nbead,nframe,fps,Roff,ori,f,t,t2,d,M,F,R,P,rz,rx,ry,x,y,z,dx,dy,dz] = deal(cell(npth,1));
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
Fdat_model = (0:.1:20)';
Lo_dsDNA = 1e3*.338; Lp_dsDNA = 37; Ko_dsDNA = 400; % dsDNA handle
zdat_model_dsDNA = eWLC_inv(Fdat_model,Lo_dsDNA,Lp_dsDNA,T,Ko_dsDNA,1);
Lo_ssDNA = 1.69*(10+2)*.338; Lp_ssDNA = 0.75; Ko_ssDNA = 800; % ssDNA linker
zdat_model_ssDNA = mFJC(Fdat_model,Lo_ssDNA,Lp_ssDNA,T,Ko_ssDNA);
Lo_PEG = 570*(1/80); Lp_PEG = 0.47; % bPEG (1 kD)
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
Lo_hp_closed = 2; zdat_model_hp_closed = Lo_hp_closed;
Lo_hp_open = 1.69*(2*8+6)*.338; % open hairpin, 34-nt ssDNA
zdat_model_hp_open = mFJC(Fdat_model,Lo_hp_open,Lp_ssDNA,T,Ko_ssDNA);
dzdat_model_closed = zdat_model_PEG + zdat_model_ssDNA + zdat_model_hp_closed + zdat_model_dsDNA;
dzdat_model_open = zdat_model_PEG + zdat_model_ssDNA + zdat_model_hp_open + zdat_model_dsDNA;

% force ramp
p = 1; n = 1; b = 1; frange = 401:8200;
tdat = t{p}{n}(frange,b); tdat = tdat-tdat(1);
Fdat = F{p}{n}(frange,b); 
dxdat = dx{p}{n}(frange,b);
dzdat = dz{p}{n}(frange,b);
dzdat_corr = correctFEC(Fdat,dzdat,dxdat,5,Rbead,Roff{p}{n},ori{p}{n});

F_clamp = 8;
frange = abs(Fdat-F_clamp) < .1;
offset = dzdat_model_open(Fdat_model==F_clamp);
dzdat_corr = dzdat_corr - mean(dzdat_corr(frange)) + offset;

h = figure(1); clf; h.WindowState = 'maximized';
subplot(131);
yyaxis left; set(gca,'ycolor','k');
plot(tdat,dzdat_corr,'k','linewidth',.1);
ylim([0,400]);
ylabel('Extension (nm)');

yyaxis right; set(gca,'ycolor','b');
plot(tdat,Fdat,'b'); grid on;
xlabel('Time (s)'); ylabel('Force (pN)');

subplot(132);
plot(dzdat_corr,Fdat,'k','linewidth',.1); grid on; hold all;
plot(dzdat_model_closed,Fdat_model,'r');
plot(dzdat_model_open,Fdat_model,'r');
xlim([200,400]);
xlabel('Extension (nm)'); ylabel('Force (pN)');

subplot(133);
plot(dzdat_corr,Fdat,'k','linewidth',.5); grid on; hold all;
plot(dzdat_model_closed,Fdat_model,'r');
plot(dzdat_model_open,Fdat_model,'r');
xlim([310,350]); ylim([4,8]);
xlabel('Extension (nm)'); ylabel('Force (pN)');

saveas(gcf,'DNA hairpin_force clamp.fig');

%% force clamp
h = figure(2); clf; h.WindowState = 'maximized';
set(gcf,'defaultaxesfontsize',10);

dzbin = 310:.75:350;
mu_list = [320,323+(0:1.5:6),333]'+[0,8];
sigma = shiftdim([7.5,5],-1);
prop_list = [.01,.15,.3,.5,.7,.85,.99]';
prop_list = [1-prop_list, prop_list];
[fitres,transmat] = deal(cell(7,1));
p = 2;
for n = 1:2
    frange = 1:nframe{p}(n);
    tdat = t{p}{n}(frange);
    Fdat = F{p}{n}(frange); 
    dzdat = dz{p}{n}(frange);
    F_clamp = 4;
    frange = abs(Fdat-F_clamp) < .1;
    offset = dzdat_model_closed(Fdat_model==F_clamp);
    dzdat = dzdat - mean(dzdat(frange)) + offset;
    if n == 1
        [~,idx] = findpeaks(diff(Fdat),'minpeakheight',.05,'minpeakdistance',500);
    else
        [~,idx] = findpeaks(diff(Fdat),'minpeakheight',.05/12,'minpeakdistance',500*12);
        dzdat_filt = medfilt1(dzdat,4);
    end
    idx = [1; idx; numel(frange)];
    
    % plot time series
    axes('position',[.4*(n-1)+.05,.6,.18,.3]);
    yyaxis right; set(gca,'ycolor','b');
    plot(tdat,Fdat,'b'); set(gca,'ygrid','on');
    ylim([2,10]); xlim(tdat([1,end]));
    set(gca,'ytick',2:10);
    xline(tdat(idx),'k');
    ylabel('Force (pN)');

    yyaxis left; set(gca,'ycolor','k');
    if n == 1
        plot(tdat,dzdat,'k-','linewidth',.1);
    else
        plot(tdat,dzdat_filt,'k-','linewidth',.1);
    end
    xlim([0,72]);
    ylim([310,350]);
    xline(tdat(idx),'k');
    ylabel('Extension (nm)');
    xlabel('Time (s)');
    
    % analyze distribution and kinetics
    for Fi = 1:7
        % plot histogram
        ax2(Fi) = axes('position',[.4*(n-1)+.29,.58+(Fi-1)*.05,.08,.05]);
        if n == 1
            frange = idx(Fi)+50:idx(Fi+1)-50;
            dzdat_tmp = dzdat(frange);
        else
            frange = idx(Fi)+600:idx(Fi+1)-600;
            dzdat_tmp = dzdat_filt(frange);
        end
        tdat_tmp = tdat(frange);       
        histogram(dzdat_tmp,dzbin,'norm','prob','edgecolor','none','facea',1);
        ylim([0,.15]);
        if Fi == 1
            xlabel('Extension (nm)'); ylabel('Frequency');
        else
            set(gca,'xticklabel','');
        end
        set(gca,'yticklabel','');
        
        if n == 2
            % kinetic analysis by HMM
            S = struct('mu',mu_list(Fi,:)','Sigma',sigma,'ComponentProportion',prop_list(Fi,:));
            fitres = fitgmdist(dzdat_tmp,2,'start',S,'Options',statset('MaxIter',1000,'TolFun',1e-4));        
            Q = 2; % number of states
            prior0 = fitres.ComponentProportion; % prob distribution of states
            transmat0 = [.99, .01; .01, .99];
            mu0 = reshape(fitres.mu, [1,Q,1]);
            Sigma0 = reshape(fitres.Sigma, [1,1,Q,1]);
            [LL, prior, transmat{Fi}, mu, Sigma, mixmat] = ...
                mhmm_em(dzdat_tmp', prior0, transmat0, mu0, Sigma0, [], 'max_iter', 100, 'adj_mu', 0, 'adj_Sigma', 1, 'adj_trans', 1);
            B = mixgauss_prob(dzdat_tmp', mu, Sigma, mixmat);
            path = viterbi_path(prior, transmat{Fi}, B);
            path_dz = (double(path)-1)*(diff(mu))+mu(1);
        end

        % plot time series_zoom in
        if Fi >= 2 && Fi <= 6
            ax1(Fi) = axes('position',[.4*(n-1)+.03,.2+(Fi-2)*.06,.35,.055]);
            tmid = tdat(round(mean(idx(Fi+(0:1)))));
            length = 1; % segment length in seconds
            frange_sub = abs(tdat_tmp-tmid) < length/2;
            plot(tdat_tmp(frange_sub),dzdat_tmp(frange_sub),'k','linewidth',1,'clip','off'); hold all;
            if n == 2
                plot(tdat_tmp(frange_sub),path_dz(frange_sub),'r','linewidth',1,'clip','off');
                h = get(gca,'Children'); uistack(h(2),'up',1);
            end
            xlim(tmid+[-.5,.5]*length); ylim([320,340]); axis off;
            if Fi == 2
                plot(tmid-length*.525+[0,0,.1],310+[10,0,0],'k','linewidth',2,'clip','off'); % scale: 100 ms, 10 nm
            end
        end
    end
end

[fitres_multi,transmat_multi] = deal(cell(7,5));
for Fi = 1:7
    S = struct('mu',mu_list(Fi,:)','Sigma',sigma,'ComponentProportion',prop_list(Fi,:));
    for j = 1:5
        frange = idx(Fi)+600 + (j-1)*1.5*1200 + (1:1.5*1200);
        dzdat_tmp = dzdat_filt(frange);
        fitres_multi{Fi,j} = fitgmdist(dzdat_tmp,2,'start',S,'Options',statset('MaxIter',1000,'TolFun',1e-4));
        Q = 2; % number of states
        prior0 = fitres_multi{Fi,j}.ComponentProportion; % prob distribution of states
        transmat0 = transmat{Fi};
        mu0 = reshape(fitres_multi{Fi,j}.mu, [1,Q,1]);
        Sigma0 = reshape(fitres_multi{Fi,j}.Sigma, [1,1,Q,1]);
        [LL, prior, transmat_multi{Fi,j}, mu, Sigma, mixmat] = ...
            mhmm_em(dzdat_tmp', prior0, transmat0, mu0, Sigma0, [], 'max_iter', 100, 'adj_mu', 0, 'adj_Sigma', 1, 'adj_trans', 1);
    end
end

% fit open probability
subplot(366);
set(gca,'position',get(gca,'position')+[.03,.04,.035,-.04]);
Fdat = [4,5:.5:7,8];
mudat = cellfun(@(x) x.mu, fitres_multi, 'unif',0);
mudat = cat(2,mudat{:}); mudat = reshape(mudat',[7,5,2]);
mudat_avg = squeeze(mean(mudat,2));
mudat_err = squeeze(std(mudat,[],2));
dmudat_avg = diff(mudat_avg,1,2);
dmudat_err = sqrt(sum(mudat_err.^2,2));
errorbar(Fdat,dmudat_avg,dmudat_err,'o'); hold on;
xlim([3.5,8.5]); ylim([4.5,9.5]);
xlabel('Force (pN)'); ylabel('Opening distance (nm)');
legend off;
for Fi = 1:7
    axes(ax2(Fi));
    xline(mudat_avg(Fi,:),'r');
    if Fi == 4
        disp(diff(mudat_avg(Fi,:)));
    end
end

subplot(3,6,12);
set(gca,'position',get(gca,'position')+[.03,.11,.035,-.04]);
pdat = cellfun(@(x) x.ComponentProportion, fitres_multi, 'unif',0);
pdat = cat(1,pdat{:}); pdat = reshape(pdat(:,2),[7,5]);
pdat_avg = mean(pdat,2);
pdat_err = std(pdat,[],2)/sqrt(5);
errorbar(Fdat,pdat_avg,pdat_err,'o'); hold on;
fit_Boltzmann = fittype('1./(1 + exp((Fmid - F)*dx/4.14))', 'indep','F','coeff',{'Fmid','dx'});
fitres_Boltzmann = fit(Fdat',pdat_avg,fit_Boltzmann,'start',[6,.1]);
Fdat_cal = Fdat(1)-.5:.05:Fdat(end)+.5;
plot(Fdat_cal,fitres_Boltzmann(Fdat_cal));
xline(fitres_Boltzmann.Fmid,'k');
text(fitres_Boltzmann.Fmid+.1,.1,['F1/2 = ',num2str(fitres_Boltzmann.Fmid,'%.1f'),' pN']);
xlim([3.5,8.5]); ylim([-.1,1.1]);
xlabel('Force (pN)'); ylabel('Open probability');
legend off;

% fit kinetics
subplot(3,6,18);
set(gca,'position',get(gca,'position')+[.03,.11,.035,0]);
% transmat_all = cat(3,transmat{:});
% k_open = squeeze(transmat_all(1,2,:))*1200;
% k_close = squeeze(transmat_all(2,1,:))*1200;
% semilogy(Fdat,k_open,'o'); hold on;
% semilogy(Fdat,k_close,'o');
kdat = cat(3,transmat_multi{:})*1200;
k_open_all = reshape(kdat(1,2,:),[7,5]);
k_close_all = reshape(kdat(2,1,:),[7,5]);
k_open_avg = mean(k_open_all,2);
k_open_err = std(k_open_all,[],2)/sqrt(5);
k_close_avg = mean(k_close_all,2);
k_close_err = std(k_close_all,[],2)/sqrt(5);
h(1) = errorbar(Fdat,k_open_avg,k_open_err,'o'); hold on;
h(2) = errorbar(Fdat,k_close_avg,k_close_err,'o');
set(gca,'yscale','log');
fit_Bell = fittype('k0*exp(F*dx/4.14)','indep','F','coeff',{'k0','dx'});
sel = 2:4;
fitres_open = fit(Fdat(sel)',k_open_avg(sel),fit_Bell,'start',[1,3])
sel = 3:7;
fitres_close = fit(-Fdat(sel)',k_close_avg(sel),fit_Bell,'start',[1e4,3])
xlim([3.5,8.5]); ylim([5,1000]);
plot(Fdat_cal,fitres_open(Fdat_cal),'color',get(h(1),'color'));
plot(Fdat_cal,fitres_close(-Fdat_cal),'color',get(h(2),'color'));
xlabel('Force (pN)'); ylabel('Transition rate (s−1)');

saveas(gcf,'DNA hairpin_force clamp.fig');
