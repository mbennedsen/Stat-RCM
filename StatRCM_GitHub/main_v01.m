%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will estimate and report main findings from Stat-RCM applied to the historical data period. 
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
addpath(genpath('Files'));
addpath(genpath('Functions'));
addpath(genpath('Data'));
%% Init 
rng_seed = 666; % Random seed.

% Number of random intializations used for estimation:
numIte = 0;  % If =0, then use pre-run parameters (from paper). In the paper, numIte=100 was used.

% Number of bootstrap replications to get std errs:
B = 0;       % Only relevant if numIte>0.

% Start and end year of data:
start_year = 1959; % Vary these for sub-sample analyses.
end_year   = 2022;

% Options for optimizer:
Opts = optimset('Display','off','TolFun',1e-12,'MaxFunEvals',1e6,'MaxIter',1e6);

%% Load Data
dat_cc   = xlsread('StatRCM_data.xlsx',1);
dat_ebm  = xlsread('StatRCM_data.xlsx',2);

% Normalize to 1850 baseline
dat_ebm(:,5) = dat_ebm(:,5) -  mean(dat_ebm(1:50,5));%dat_ebm(1,5); % TAS 
% Normalize to 1940 baseline
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
x_FnonCO2 = dat_ebm(:,3);
x_Fnat    = dat_ebm(:,4);
y_TAS     = dat_ebm(:,5);
y_OcT     = dat_ebm(:,6);
y_OHC     = dat_ebm(:,7);

%%% Other data
conc_1750 =  278;   %ppm
conc_1959 = 315.39; %ppm
C0  = conc_1959*2.127;
C00 = conc_1750*2.127;

y_C = C0 + cumsum(y_ATM); % Concentrations
x_E = y_FF+y_LUC; % Total anthro CO2 emissions

%% Setup obs.+time vector
y = [y_C,y_OCN,y_LND,y_FCO2,y_TAS,y_OcT,y_OHC];
x = [x_E,x_FnonCO2,x_Fnat];
t = dat_cc(:,1);


%% Retrieve output
output = estimate_sqrtlog(t,y,x,start_year,end_year,numIte,B,Opts,rng_seed);

%% Retrieve output
ML_est = output.ML_est;
ML_std_boot = output.ML_std_boot;
ML_tstat_boot = output.ML_tstat_boot;

ML_est_sig = output.ML_est_sig;
ML_std_sig_boot = output.ML_std_sig_boot;
ML_tstat_sig_boot = output.ML_tstat_sig_boot;

ML_est_ARMA = output.ML_est_ARMA;
ML_std_ARMA_boot = output.ML_std_ARMA_boot;
ML_tstat_ARMA_boot = output.ML_tstat_ARMA_boot;

x_smooth = output.x_smooth;

Cm = output.Cm;
Cd = output.Cd;
mud = output.mud;
muT = output.muT;
muO = output.muO;

h = output.h;
h_std = output.h_std;

ECS = output.ECS;
std_ECS = output.std_ECS;
ECS_boot = output.ECS_boot;

loglik = output.loglik;
BIC = output.BIC;

GoF = output.GoF;
std_prediction_residuals = output.std_prediction_residuals;

%% Print parameters to screen
disp(' ');
disp(['log-likelihood = ',num2str(loglik),', BIC = ',num2str(BIC)]);
disp(' ');
disp( 'Parameters estimated by ML :    b1        b2      delta1    delta2       f1       f2       f3       gamma     lambda      Cm        Cd       effi      mu_T      mu_d      mu_O      hh');
disp(['Estimates:                 :   ',num2str(ML_est',3)]);
disp(['(Std. Errs, boot)          :   ',num2str(ML_std_boot',3)]);
disp(['t-stats (boot)             :   ',num2str(ML_tstat_boot',3)]);

disp(' ');
disp( 'Variance params estimated by ML:  s_eta_OCN    s_eta_LND   s_eta_F   s_eta_m   s_eta_d   s_sig_C   s_sig_OCN   s_sig_LND   s_sig_F   s_sig_m   s_sig_d   s_sig_O   rho');
disp(['Estimates:                     :   ',num2str(ML_est_sig',3)]);
disp(['(Std. Errs, boot)              :   ',num2str(ML_std_sig_boot',3)]);
disp(['t-stats (boot)                 :   ',num2str(ML_tstat_sig_boot',3)]);


disp(' ');
disp( 'ARMA params estimated by ML:    phi1      phi2     phi3      phi4      phi5      phi6      phi7');
disp(['Estimates:                 :   ',num2str(ML_est_ARMA',3)]);
disp(['(Std. Errs, boot)          :   ',num2str(ML_std_ARMA_boot',3)]);
disp(['t-stats (boot)             :   ',num2str(ML_tstat_ARMA_boot',3)]);

%% Print GoF to screen
disp('Standardized prediction residuals (diagnostics):');
disp('     nObs      Mean       std      Skew      Kurt       AR1       JB        DW       LB1       LB5      LB10      ARCH');
disp(round(GoF,2))


disp(' ')
disp(['Estimated ECS: ',num2str(mean(ECS))]);
disp(['5% and 95% quantiles of ECS (using bootstrap sample): ',num2str(quantile(ECS_boot,0.05)),', ',num2str(quantile(ECS_boot,1-0.05)),'.']);    
disp(['90% CI (Gaussian approx): [',num2str(ECS - 1.6449*std_ECS),', ',num2str(ECS + 1.6449*std_ECS),'].'])
disp(' ')
disp(['Estimated height of mixed layer, h = ',num2str(h),' (',num2str(h_std),')'] );

disp(' ')
disp(['Filtered temperature estimate in ',num2str(t(end)),': ',num2str(x_smooth(end,5)+muT)]);
%% plot smoothed states
stVal = 2;
for i = 1:7
    fig3 = figure(3);
    subplot(4,2,i)
    plot(t,y(:,i),'k-','LineWidth',2), hold on
    if i == 4
        plot(t(stVal:end),x_smooth(stVal:end,i),'r-.','LineWidth',2), hold on
    elseif i == 5
        plot(t(stVal:end),x_smooth(stVal:end,i)+muT,'r-.','LineWidth',2), hold on  
    elseif i == 6
        plot(t(stVal:end),x_smooth(stVal:end,i-1)*h/2000+ x_smooth(stVal:end,i)*(1-h/2000) +  mud,'r-.','LineWidth',2), hold on 
    elseif i == 7
        plot(t(stVal:end),x_smooth(stVal:end,i-1)*Cd + x_smooth(stVal:end,i-2)*Cm + muO  ,'r-.','LineWidth',2), hold on
    else
        plot(t(stVal:end),x_smooth(stVal:end,i),'r-.','LineWidth',2), hold on
    end
    
    if i == 1
        title('Atmospheric concentrations (C)','FontSize',10);
        ylabel('GtC','FontSize',10,'Interpreter','latex');
        %legend('Data','Smoothed state','95p CI','Location','NorthWest');
    elseif i == 2
        title('Oceanic carbon sink (OCN)','FontSize',10);
        ylabel('GtC/yr','FontSize',10,'Interpreter','latex');
    elseif i == 3
        title('Terrestrial carbon sink (LND)','FontSize',10);
        ylabel('GtC/yr','FontSize',10,'Interpreter','latex');
    elseif i == 4
        title('Forcings from CO2','FontSize',10);
        ylabel('W/m$^2$','FontSize',10,'Interpreter','latex');
    elseif i == 5
        title('Surface temperature','FontSize',10);
        ylabel('K','FontSize',10,'Interpreter','latex');
     elseif i == 6
        title('Ocean temperature','FontSize',10);
        ylabel('K','FontSize',10,'Interpreter','latex');
    elseif i == 7
        title('Ocean heat content','FontSize',10);
        ylabel('W/m$^2$','FontSize',10,'Interpreter','latex');
        set(gca,'fontsize',8)
        xlim([start_year,end_year])
        
        subplot(4,2,8);
        plot(nan,nan,'k-','LineWidth',2), hold on
        plot(nan,nan,'r-.','LineWidth',2), hold on
        plot(nan,nan,'r--','LineWidth',1), hold on
        lgd = legend('Data','Smoothed state','Interpreter','latex','Location','NorthEast');
        lgd.FontSize = 8;
        set(gca, 'XTick', [], 'YTick', [])
        axis off
    end
    set(gca,'fontsize',8)
    xlim([start_year,end_year])


    fig4 = figure(4);
    subplot(4,2,i)
    plot(t(stVal:end),std_prediction_residuals(stVal:end,i),'r-','LineWidth',2), hold on
    plot(t,zeros(length(t),1),'k-','LineWidth',1), hold on
     plot(t(stVal:end),std_prediction_residuals(stVal:end,i),'r-','LineWidth',2), hold on
    ylabel('Residuals','FontSize',10,'Interpreter','latex');
    if i == 1
        title('Atmospheric concentrations (C)','FontSize',10);
    elseif i == 2
        title('Oceanic carbon sink (OCN)','FontSize',10);
    elseif i == 3

        title('Terrestrial carbon sink (LND)','FontSize',10);
    elseif i == 4
        title('Forcings from CO2','FontSize',10);
    elseif i == 5
        title('Surface temperature','FontSize',10);
     elseif i == 6
        title('Ocean temperature','FontSize',10);
    elseif i == 7
        title('Ocean heat content','FontSize',10);
        set(gca,'fontsize',8)
        xlim([start_year,end_year])
        ylim([-3.25,3.25]);
        
        subplot(4,2,8);
        plot(nan,nan,'r-','LineWidth',2), hold on
        legend('Standardized prediction residuals','Location','NorthEast');
        lgd = legend('Standardized prediction residuals','Interpreter','latex','Location','NorthEast');
        lgd.FontSize = 10;
        set(gca, 'XTick', [], 'YTick', [])
        axis off
    end
    set(gca,'fontsize',8)
        xlim([start_year,end_year])
        ylim([-3.25,3.25]);
end

