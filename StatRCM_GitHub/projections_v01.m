%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will use the estimated Stat-RCM and an SSP scenario to project all climate variables to 2100. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2025)
%
% This code can be used, distributed, and changed freely. 
% Please cite Bennedsen, Hillebrand, and Koopman (2025): 
% "A Statistical Reduced Complexity Climate Model for Probabilistic Analyses and Projections", 
% Journal of Climate, Volume 38, Issue 21, pp. 6329-6350.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
addpath(genpath('Data'));
addpath(genpath('Functions'));
addpath(genpath('Files'));
%% Init
str_param     = 'estparams0_sqrtlog_1959_2022';
str_bootstrap = 'boot_stderrs_sqrtlog_1959_2022';

% Input RCP scenario: 
RCP = 2; % Valid: 1, 2, 3, 4, 5.

MC = 1e4;

q = 5/100; % quantile for confidence bands (q/2 and 1-q/2).

start_year = 1959;
end_year   = 2022;

indx_vec =  [1:8,32:35]; % Physical parameters, subject to pertubations.

n = 2100-end_year; % number of years in projection
%% Load RCP Data
disp(' ');
if RCP == 1
    load('SSP119_output_v2');
    disp('Projections with Stat-RCM using the RCP119 scenario...')
elseif RCP==2
    load('SSP126_output_v2');
    disp('Projections with Stat-RCM using the RCP126 scenario...')
elseif RCP==3
    load('SSP245_output_v2');
    disp('Projections with Stat-RCM using the RCP245 scenario...')
elseif RCP==4
    load('SSP370_output_v2');
    disp('Projections with Stat-RCM using the RCP370 scenario...')
elseif RCP==5
    load('SSP434_output_v2');
    disp('Projections with Stat-RCM using the RCP434 scenario...')
elseif RCP==6
    load('SSP585_output_v2');
    disp('Projections with Stat-RCM using the RCP585 scenario...')
else
    error('RCP scenario not implemented...');
end


%% Load Data
dat_cc   = xlsread('StatRCM_data.xlsx',1);
dat_ebm  = xlsread('StatRCM_data.xlsx',2);

%%% Difference wrt first obs in EBM
% 1850 baseline
dat_ebm(:,5) = dat_ebm(:,5) -  mean(dat_ebm(1:50,5)); % TAS 
% 1940 baseline
dat_ebm(:,6) = dat_ebm(:,6) - dat_ebm(dat_ebm(:,1)==1940,6); % OTemp
dat_ebm(:,7) = dat_ebm(:,7) - dat_ebm(dat_ebm(:,1)==1940,7); % OHC

%%% Align with start and end years:
dat_cc(dat_cc(:,1)<start_year,:) = [];
dat_cc(dat_cc(:,1)>end_year,:) = [];
dat_ebm(dat_ebm(:,1)<start_year,:) = [];
dat_ebm(dat_ebm(:,1)>end_year,:) = [];

%%% Assign data values
y_FF = dat_cc(:,2);
y_LUC = dat_cc(:,3);
y_ATM = dat_cc(:,4);
y_OCN = dat_cc(:,5);
y_LND = dat_cc(:,6);

y_FCO2    = dat_ebm(:,2);  
x_FnonCO20 = dat_ebm(:,3);
x_Fnat0    = dat_ebm(:,4);
y_TAS     = dat_ebm(:,5);
y_OcT     = dat_ebm(:,6);
y_OHC     = dat_ebm(:,7);

%%% Other data
conc_1750 =  278; % ppm
conc_1959 = 315.39; %ppm
C0  = conc_1959*2.127;
C00 = conc_1750*2.127;

tmp = y_ATM; tmp(isnan(tmp)) = 0;
y_C = C0 + cumsum(tmp);
y_C(isnan(y_ATM)) = nan;
x_E0 = y_FF+y_LUC;
%% Setup obs.+time vector
y_orig = [y_C,y_OCN,y_LND,y_FCO2,y_TAS,y_OcT,y_OHC];
t0 = dat_cc(:,1);
t  = dat_cc(end,1)+(1:n);

T = n;

%% Make experiment data
x_Fnat = x_Fnat0(end)*ones(n,1); % Keep natural forcings to last in-sample value

% Use CO2 emissions and non-CO2 forcings as output from MAGICC:
x_E       = E_magicc(1+end_year-1750+1:end)+1.1918; %=1.1918 = Diff in E_magicc(2022) and actual FF emissions 2022.
x_FnonCO2 = FnonCO2_magicc(1+end_year-1995+1:end); 

%% Load parameters
load(str_param);

%% Load bootstrap params
load(str_bootstrap);
p_boot_trans = nan(size(p_boot));
for iP = 1:size(p_boot,1)
    p_boot_trans(iP,:) = pInvTrans_sqrtlog(p_boot(iP,:));
end
boot_SIG = cov(p_boot_trans);

boot_SIG2 = boot_SIG(indx_vec,indx_vec);
cholCov2 = chol(boot_SIG2);

numB = size(p_boot,1);
indx = randi(numB, MC, 1);

P00 = blkdiag(zeros(6),0,0*eye(7));
x00 = x_end_smooth';
%% Run MC simulation
y_deterministic = nan(n,7);
y_theta = nan(MC,n,7);
y_theta_eta = nan(MC,n,7);
y_theta_eta_eps = nan(MC,n,7);
for iMC = 1:MC
    %disp(iMC/MC); 
    
    if iMC == 1
        %% Deterministic run
        [T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams0);
        
        H = H*0;
        Q = Q*0;

        x0 = x00 + sqrt(P00)*randn(length(x00),1);

        y = nan(T,7);
        eta = randn(size(Q,2),T);
        eps = randn(7,T);
        for i = 1:T
            if i == 1
                x = T_fct{i}(x0) + R*Q*eta(:,i);
            else
                x = T_fct{i}(x) + R*Q*eta(:,i);
            end

            y(i,:) = Z_fct(x) + H*eps(:,i);
        end

        for i=1:7
            y_deterministic(:,i) = y(:,i)';
        end
    end
    
    




    %% Epistemic uncertainty (uncertainty from parameters)
    estParams00 = estParams0;
    estParams00(indx_vec) = estParams00(indx_vec) + cholCov2'*randn(length(indx_vec),1);
    
    [T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams00);

    H = H*0;
    Q = Q*0;
    
    x0 = x00 + sqrt(P00)*randn(length(x00),1);
    
    y = nan(T,7);
    eta = randn(size(Q,2),T);
    eps = randn(7,T);
    for i = 1:T
        if i == 1
            x = T_fct{i}(x0) + R*Q*eta(:,i);
        else
            x = T_fct{i}(x) + R*Q*eta(:,i);
        end
        
        y(i,:) = Z_fct(x) + H*eps(:,i);
    end
    
    for i=1:7
        y_theta(iMC,:,i) = y(:,i)';
    end


    %% Epistemic uncertainty (uncertainty from parameters) + eta (aleatoric uncertainty from state evolution)
    estParams00 = estParams0;
    estParams00(indx_vec) = estParams00(indx_vec) + cholCov2'*randn(length(indx_vec),1);
    
    [T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams00);

    H = H*0;
    Q(6:end,6:end) = zeros(7); % M-errs
    
    x0 = x00 + sqrt(P00)*randn(length(x00),1);
    
    y = nan(T,7);
    eta = randn(size(Q,2),T);
    eps = randn(7,T);
    for i = 1:T
        if i == 1
            x = T_fct{i}(x0) + R*Q*eta(:,i);
        else
            x = T_fct{i}(x) + R*Q*eta(:,i);
        end
        
        y(i,:) = Z_fct(x) + H*eps(:,i);
    end
    
    for i=1:7
        y_theta_eta(iMC,:,i) = y(:,i)';
    end


    %% Epistemic uncertainty (uncertainty from parameters) + eta (aleatoric uncertainty from state evolution) + eps (aleatoric uncertainty from measurement errors and other transient sources)
    estParams00 = estParams0;
    estParams00(indx_vec) = estParams00(indx_vec) + cholCov2'*randn(length(indx_vec),1);
    
    [T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams00);

    x0 = x00 + sqrt(P00)*randn(length(x00),1);
    
    y = nan(T,7);
    eta = randn(size(Q,2),T);
    eps = randn(7,T);
    for i = 1:T
        if i == 1
            x = T_fct{i}(x0) + R*Q*eta(:,i);
        else
            x = T_fct{i}(x) + R*Q*eta(:,i);
        end
        
        y(i,:) = Z_fct(x) + H*eps(:,i);
    end
    
    for i=1:7
        y_theta_eta_eps(iMC,:,i) = y(:,i)';
    end


end



%% plot
fig1 = figure;
x_total0 = x_E0 + x_FnonCO20 + x_Fnat0;
x_total  = x_E + x_FnonCO2 + x_Fnat;
subplot(4,2,1);
plot(t0,x_E0,'k-','LineWidth',2), hold on
plot(t0,x_FnonCO20,'b-','LineWidth',2), hold on
plot(t0,x_Fnat0,'g-','LineWidth',2), hold on

plot(t0,x_E0,'k-','LineWidth',2), hold on
plot(t0,x_FnonCO20,'b-','LineWidth',2), hold on
plot(t0,x_Fnat0,'g-','LineWidth',2), hold on

plot(t,x_E,'k-.','LineWidth',2), hold on
plot(t,x_FnonCO2,'b-.','LineWidth',2), hold on
plot(t,x_Fnat,'g-.','LineWidth',2), hold on


xlim([start_year,2100])

ylim([-12,12])
lgd = legend('CO2 emissions (GtC/yr)','Non-CO2 forcings (W/m$^2$)','Natural forcings (W/m$^2$)','Interpreter','latex','Location','SouthEast');
lgd.FontSize = 8;
legend('boxoff');
title('Input: Emissions and exogeneous forcings','FontSize',10);
set(gca,'fontsize',8)

for i=1:7
    subplot(4,2,1+i);
    plot(t0,y_orig(:,i),'r-','LineWidth',2), hold on
     
    %% Deterministic run
    plot(t,y_deterministic(:,i),'k-','LineWidth',1), hold on
    
    %% Uncertainty: theta+eta+eps
    yy = [quantile(y_theta_eta_eps(:,:,i),q/2),fliplr(quantile(y_theta_eta_eps(:,:,i),1-q/2))];

    xx = [t,fliplr(t)];
    fill(xx,yy,[0.3010, 0.7450, 0.9330]), hold on

    %% Uncertainty: theta+eta
    yy = [quantile(y_theta_eta(:,:,i),q/2),fliplr(quantile(y_theta_eta(:,:,i),1-q/2))];

    xx = [t,fliplr(t)];
    fill(xx,yy,[0, 0.4470, 0.7410]), hold on
    
    %% Uncertainty: theta
    yy = [quantile(y_theta(:,:,i),q/2),fliplr(quantile(y_theta(:,:,i),1-q/2))];

    xx = [t,fliplr(t)];
    fill(xx,yy,[0.3010, 0.7450, 0.9330]/2), hold on
    
    %% Deterministic run
    plot(t,y_deterministic(:,i),'k-','LineWidth',1), hold on
    
    if i == 1
        title('Output: Atmospheric concentrations (C)','FontSize',10);
        ylabel('GtC','FontSize',10,'Interpreter','latex');
    elseif i == 2
        title('Output: Oceanic carbon sink (OCN)','FontSize',10);
        ylabel('GtC/yr','FontSize',10,'Interpreter','latex');
    elseif i == 3
        title('Output: Terrestrial carbon sink (LND)','FontSize',10);
        ylabel('GtC/yr','FontSize',10,'Interpreter','latex');
    elseif i == 4
        title('Output: Forcings from CO2','FontSize',10);
        ylabel('W/m$^2$','FontSize',10,'Interpreter','latex');
    elseif i == 5
        title('Output: Surface temperature','FontSize',10);
        ylabel('K','FontSize',10,'Interpreter','latex');
     elseif i == 6
        title('Output: Ocean temperature','FontSize',10);
        ylabel('K','FontSize',10,'Interpreter','latex');
    elseif i == 7
        title('Output: Ocean heat content','FontSize',10);
        ylabel('W/m$^2$','FontSize',10,'Interpreter','latex');
    end
    set(gca,'fontsize',8)
    xlim([start_year,2100])

end



%% plot only TAS - individual runs
thres = 1.5; % Temperature threshold
numMCplot = 50; 

fig3 = figure;
i = 5; % TAS index

subplot(1,3,1);
plot(t0,y_orig(:,i),'r-','LineWidth',1.5), hold on

%%% Uncertainty: theta
plot(t,y_theta(1:numMCplot,:,i),'LineWidth',1), hold on

plot([t0;t'],thres*ones(length([t0;t']),1),'k-.','LineWidth',1.5), hold on
title('Surface temperature (uncertainty: $\theta$)','FontSize',10,'Interpreter','latex');
ylabel('K','FontSize',10,'Interpreter','latex');
set(gca,'fontsize',8)
xlim([start_year,2100])

subplot(1,3,2);
plot(t0,y_orig(:,i),'r-','LineWidth',1.5), hold on

%%% Uncertainty: theta+eta
plot(t,y_theta_eta(1:numMCplot,:,i),'LineWidth',1), hold on

plot([t0;t'],thres*ones(length([t0;t']),1),'k-.','LineWidth',1.5), hold on
title('Surface temperature (uncertainty: $\theta, \eta$)','FontSize',10,'Interpreter','latex');
ylabel('K','FontSize',10,'Interpreter','latex');

set(gca,'fontsize',8)
xlim([start_year,2100])

subplot(1,3,3);
plot(t0,y_orig(:,i),'r-','LineWidth',1.5), hold on

%%% Uncertainty: theta+eta+eps
plot(t,y_theta_eta_eps(1:numMCplot,:,i),'LineWidth',1), hold on

plot([t0;t'],thres*ones(length([t0;t']),1),'k-.','LineWidth',1.5), hold on
title('Surface temperature (uncertainty: $\theta, \eta, \epsilon$)','FontSize',10,'Interpreter','latex');
ylabel('K','FontSize',10,'Interpreter','latex');

set(gca,'fontsize',8)
xlim([start_year,2100])

%% display
disp(['Probability of crossing of ',num2str(thres),' degrees:  theta   theta+eta   theta+eta+eps'])
disp(['                                        ', num2str([mean(sum(y_theta(:,:,5)>thres,2)>0), mean(sum(y_theta_eta(:,:,5)>thres,2)>0), mean(sum(y_theta_eta_eps(:,:,5)>thres,2)>0)])]);

