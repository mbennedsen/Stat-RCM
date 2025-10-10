function output = estimate_sqrtlog(t,y,x,start_year,end_year,numIte,B,Opts,rng_seed)
if nargin < 6
    numIte = 0;
end
if nargin < 7
    B = 1000;
end
if nargin < 8
    Opts = optimset('Display','off','TolFun',1e-12,'MaxFunEvals',1e6,'MaxIter',1e6);
end
if nargin < 9
    rng_seed = 666;
end

%% Load in data
y_C    = y(:,1);
y_OCN  = y(:,2);
y_LND  = y(:,3);
y_FCO2 = y(:,4);
y_TAS  = y(:,5);
y_OcT  = y(:,6);
y_OHC  = y(:,7);

x_E       = x(:,1);
x_FnonCO2 = x(:,2);
x_Fnat    = x(:,3);

if end_year>2019
    y_FCO2(end-(end_year-2019-1):end) = y_FCO2(end-(end_year-2019)); % For use in OLS regression below
end
%%% Other data
conc_1750 =  278;   %ppm
conc_1959 = 315.39; %ppm
C00 = conc_1750*2.127;

%% Estimate or retrieve parameters
if numIte == 0 
    if exist(['Files/estimated_params_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'])>0
        disp(' ')
        disp(' --- Loading pre-saved parameter estimates of Stat-RCM --- ')
        load(['Files/estimated_params_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat']);
    else
        error('No parameters pre-estimated for this setting. Please run with numIte>0.')
    end
else
    disp(' ')
    disp([' --- Estimating parameters of Stat-RCM (best fit from ',num2str(numIte),' initializations) --- '])
    rng(rng_seed);

    save_loglik = nan(numIte,1);
    save_params = nan(numIte,35);
    save_StParams = nan(numIte,35);
    for iP = 1:numIte
         disp(['Estimating parameters, run ',num2str(iP), ' of ',num2str(numIte),'.']);
    
        %% Setup param0: Randomized
        params0     = [abs(randn(31,1));sqrt(.1);sqrt(.1);3;-3];
        
        %% Setup param0: Initial estimates via OLS
        gammO = 0.1;
        gammL = 0.1;

        XX = [y_C(1:end-1).*exp(-gammO*y_TAS(1:end-1)) - C00];
        yy = y_OCN(2:end);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_O_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        params0(1) = sqrt(b_tmp(1));
          
        
        XX = [y_C(1:end-1).*exp(-gammL*y_TAS(1:end-1))-C00];
        yy = y_LND(2:end);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_L_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        params0(2) = sqrt(b_tmp(1));
        
        
        XX = [log(y_C(1:end-1))-log(C00),sqrt(y_C(1:end-1))-sqrt(C00)];%[log(y_C/C00),sqrt(y_C-C00)];
        yy = y_FCO2(2:end);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_F_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        muF_tmp = 0;
        params0(3) = sqrt(abs(b_tmp(1)));
        params0(4) = sqrt(abs(b_tmp(2)));
        
        XX = [ones(length(t),1),y_OcT];
        yy = y_OHC;
        b_tmp = (XX'*XX)\XX'*yy;
        eta_D_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t) - length(b_tmp));
        Cd_tmp = b_tmp(2);
        
        XX = [ones(length(t)-1,1),y_TAS(1:end-1)-y_OcT(1:end-1)];
        yy = diff(y_OcT);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_OcT_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        gam_tmp = b_tmp(2)*Cd_tmp;
        a2_tmp = b_tmp(1);
        
        XX = [ones(length(t)-1,1),y_TAS(1:end-1),y_OcT(1:end-1),y_FCO2(1:end-1)+x_Fnat(1:end-1) + x_FnonCO2(1:end-1)];
        yy = y_TAS(2:end);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_T_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        Cm_tmp = 1/b_tmp(4);
        lam_tmp = Cm_tmp*(1-b_tmp(2)) - gam_tmp;
        a1_tmp = b_tmp(1);
         
        
        muT_tmp = Cm_tmp/lam_tmp*(a1_tmp + muF_tmp/Cm_tmp) + Cd_tmp/lam_tmp*a2_tmp;
        mud_tmp = Cm_tmp/lam_tmp*(a1_tmp + muF_tmp/Cm_tmp) + (gam_tmp+lam_tmp)/lam_tmp/gam_tmp*Cd_tmp*a2_tmp;
        
        cc = 4;
        params0(cc+1) = sqrt(gam_tmp);
        params0(cc+2) = sqrt(lam_tmp);
        params0(cc+3) = sqrt(Cm_tmp);
        params0(cc+4) = sqrt(Cd_tmp);
        
        cc = cc+4;
        params0(cc+1:cc+5) = sqrt(sqrt([eta_O_tmp;eta_L_tmp;eta_F_tmp;eta_T_tmp;eta_OcT_tmp]));
        cc = cc+5;
        
        cc = cc+7; % m-errors.
        
        params0(cc+1) = mud_tmp;
        params0(cc+2) = muT_tmp;

        %% Estimate parameters with EKF:
        P00 = blkdiag(1e5*eye(6),0,1e2*eye(7));
        x00 = [672;zeros(5,1);1;zeros(7,1)];
        obj_fct = @(params)( -1*LogLik_EKF_sqrtlog(y',x_E,x_FnonCO2,x_Fnat,params,x00,P00) );
        
        [est_params,fval] = fminunc(obj_fct,params0, Opts);
        
        fval_diff = 1e4; 
        while fval_diff>1e-3
            fval2 = fval;
            [est_params,fval] = fminunc(obj_fct,est_params, Opts);
            fval_diff = abs(fval-fval2);
        end
        [est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
        [est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
        [est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
        [est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
        [est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
        
        loglik = -1*obj_fct(est_params);


        save_StParams(iP,:) = params0;
        save_params(iP,:) = est_params;
        save_loglik(iP) = loglik;
    end
    if exist(['Files/estimated_params_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'])==0
        save(['Files/estimated_params_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat']);
    end
 
end

%% Get best params
[val,indx] = max(save_loglik);
params0 = save_params(indx,:)';

%% Estimate parameters with EKF once more
P00 = blkdiag(1e5*eye(6),0,1e2*eye(7));
x00 = [672;zeros(5,1);1;zeros(7,1)];
obj_fct = @(params)( -1*LogLik_EKF_sqrtlog(y',x_E,x_FnonCO2,x_Fnat,params,x00,P00) );

[est_params,fval] = fminunc(obj_fct,params0, Opts);

fval_diff = 1e4; 
while fval_diff>1e-3
    fval2 = fval;
    [est_params,fval] = fminunc(obj_fct,est_params, Opts);
    fval_diff = abs(fval-fval2);
end
[est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
[est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
[est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
[est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);
[est_params,fval,exitflag,output,grad,hessian] = fminunc(obj_fct,est_params, Opts);

loglik = -1*obj_fct(est_params);
AIC = 2*length(est_params) - 2*loglik;
AICc = AIC + (2*length(est_params)^2 + 2*length(est_params))/(length(t)-length(est_params)-1);
BIC = length(est_params)*log(length(t)) - 2*loglik;

hessian = -1*hessian;

estParams0 = est_params;

%% Load bootstrap std errs
if numIte == 0 % Load pre-saved bootstrap standard errors
    if exist(['Files/boot_stderrs_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'])>0
        disp(' ')
        disp(' --- Loading pre-saved bootstrap sample --- ')
        load(['Files/boot_stderrs_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat']);
    else
        error('No bootstrap standard errors pre-estimated for this setting. Please run with numIte>0 and B>0.')
    end
else % Run bootstrap to get std errs.
    disp(' ')
    disp(' --- Running bootstrap procedure (time consuming; consider parallelizing) --- ')
    [p_boot,ECS_boot] = bootstrap_stderrs_sqrtlog(t,B,x,estParams0,numIte,x00,P00,Opts);

    if exist(['Files/boot_stderrs_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'])==0
        B_boot = B;
        numIte_boot = numIte;
        save(['Files/boot_stderrs_sqrtlog_',num2str(start_year),'_',num2str(end_year)],'p_boot','ECS_boot','B_boot','numIte_boot');
    end
end
boot_stderrs = std(p_boot)';


%% Re-parametrize estimated parameters
b1 = (estParams0(1))^2;
b1_std_boot = boot_stderrs(1);

b2 = (estParams0(2))^2;
b2_std_boot = boot_stderrs(2);

f1 = (estParams0(3))^2;
f1_std_boot = boot_stderrs(3);


f2 = 0;
f2_std_boot = nan;

f3 = (estParams0(4))^2;
f3_std_boot = boot_stderrs(4);

cc = 4;
gamma = (estParams0(cc+1))^2;
gamma_std_boot = boot_stderrs(cc+1);


lambda = (estParams0(cc+2))^2;
lambda_std_boot = boot_stderrs(cc+2);


Cm = (estParams0(cc+3))^2;
Cm_std_boot = boot_stderrs(cc+3);


Cd = (estParams0(cc+4))^2;
Cd_std_boot = boot_stderrs(cc+4);


cc=cc+4;
s_eta = (estParams0(cc+1:cc+5)).^2;
s_eta_std_boot = boot_stderrs(cc+1:cc+5);


cc = cc+5;
s_sig = (estParams0(cc+1:cc+7)).^2;
s_sig_std_boot = boot_stderrs(cc+1:cc+7);


cc = cc+7;
mud = estParams0(cc+1);
mud_std_boot = boot_stderrs(cc+1);


muT = estParams0(cc+2);
muT_std_boot = boot_stderrs(cc+2);

muO = estParams0(cc+3);
muO_std_boot = boot_stderrs(cc+3);


cc = cc+3;
sigmoid = @(x)( exp(x)./(1+exp(x)) );

rho = -1+2*sigmoid(estParams0(cc+1:cc+7));
rho_std_boot = boot_stderrs(cc+1:cc+7);

cc = cc+7;

rho_d = -1+2*sigmoid(estParams0(cc+1:cc+1));
rho_d_std_boot = boot_stderrs(cc+1);


cc = cc+1;
gammO = (estParams0(cc+1))^2;
gammO_std_boot = boot_stderrs(cc+1);

gammL = (estParams0(cc+2))^2;
gammL_std_boot = boot_stderrs(cc+2);

ww = sigmoid(estParams0(cc+3));
ww_std_boot = boot_stderrs(cc+3);


efficacy = 1+sigmoid(estParams0(cc+4));
efficacy_std_boot = boot_stderrs(cc+4);

hh = 2000*(1-ww);
hh_std_boot = 2000*ww_std_boot;

ML_est = [b1;b2;gammO;gammL;f1;f2;f3;gamma;lambda;Cm;Cd;efficacy;muT;mud;muO;hh];
ML_std_boot = [b1_std_boot;b2_std_boot;gammO_std_boot;gammL_std_boot;f1_std_boot;f2_std_boot;f3_std_boot;gamma_std_boot;lambda_std_boot;Cm_std_boot;Cd_std_boot;efficacy_std_boot;muT_std_boot;mud_std_boot;muO_std_boot;hh_std_boot];%;s_E_std];
ML_tstat_boot = ML_est./ML_std_boot;

ML_est_ARMA = rho;
ML_std_ARMA_boot = rho_std_boot;
ML_tstat_ARMA_boot = ML_est_ARMA./ML_std_ARMA_boot;

ML_est_sig = [s_eta;s_sig;rho_d];
ML_std_sig_boot = [s_eta_std_boot;s_sig_std_boot;rho_d_std_boot];
ML_tstat_sig_boot = ML_est_sig./ML_std_sig_boot;


ECS = 3.93/lambda; % 3.93 is from AR6
std_ECS = std(ECS_boot);
    

%% Ocean layer height
h = 2000*(1-ww);
h_std = 2000*ww_std_boot;

%% Goodness-of-fit
[T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams0,x00,P00);
[x_pred, x_filter, P, x_smooth] = EKF_Model1_v01(y',T_fct,Tp_fct,Q*Q',Z_fct,Zp_fct,H*H',R,x00,P00);
%[~, ~, ~, ~,~,V_smooth] = EKF_Model1_wV_v01(y',T_fct,Tp_fct,Q*Q',Z_fct,Zp_fct,H*H',R,x00,P00);
stVal = 2;

y_pred     = Z_fct(x_pred);
residuals_filter = y_pred'-y;

Z_tmp = Zp_fct(x_pred(:,1));
std_prediction_residuals = nan(size(residuals_filter));
for ii = 1:size(y,1)
    cov_pred = Z_tmp*P(:,:,ii)*Z_tmp' + H*H';
    std_y_pred = sqrt(diag(cov_pred));
    std_prediction_residuals(ii,:) = residuals_filter(ii,:)./std_y_pred';
end

if exist(['Files/estparams0_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'])==0
    x_end_smooth = x_smooth(end,:);
    save(['Files/estparams0_sqrtlog_',num2str(start_year),'_',num2str(end_year),'.mat'],'estParams0','x_end_smooth');
end

%%% Run Goodness-of-fit
GoF = nan(size(y,2),12);
for iY = 1:size(y,2)
    e_hat = std_prediction_residuals(stVal:end,iY);
    if iY == 4
        e_hat(end-2:end) = []; % Forcings NaN in last 3 entries
    end
    
    v_t = e_hat;
    m    = nan(4,1);
    m(1) = mean(v_t);
    for i = 2:4
        m(i) = mean( (v_t-m(1)).^i );
    end
    
    S = m(3)/sqrt(m(2)^3);
    K = m(4)/m(2)^2;
    NN = length(e_hat)*(S^2/6 + (K-3)^2/24);
    
    DW = sum( diff(v_t).^2 )/sum(v_t.^2);
    
    c_ar1 = mvregress(e_hat(1:end-1),e_hat(2:end));
    
    %[h1,pValue1,LB,cValue1] = lbqtest(e_hat);
    [h1,pValue1,LB1,cValue1] = lbqtest(e_hat,'Lags',1);
    [h2,pValue2,LB5,cValue2] = lbqtest(e_hat,'Lags',5);
    [h3,pValue3,LB10,cValue3] = lbqtest(e_hat,'Lags',10);
    
    [h4,pValue4,ARCH,cValue4] = archtest(e_hat);

    GoF(iY,:) = [length(e_hat),mean(e_hat),std(e_hat),skewness(e_hat),kurtosis(e_hat),c_ar1,NN,DW,LB1,LB5,LB10,ARCH];
    
end


%% Save relevant output to struct
output.ML_est = ML_est;
output.ML_std_boot = ML_std_boot;
output.ML_tstat_boot = ML_tstat_boot;

output.ML_est_sig = ML_est_sig;
output.ML_std_sig_boot = ML_std_sig_boot;
output.ML_tstat_sig_boot = ML_tstat_sig_boot;

output.ML_est_ARMA = ML_est_ARMA;
output.ML_std_ARMA_boot = ML_std_ARMA_boot;
output.ML_tstat_ARMA_boot = ML_tstat_ARMA_boot;

output.x_smooth = x_smooth';

output.Cm = Cm;
output.Cd = Cd;
output.mud = mud;
output.muT = muT;
output.muO = muO;

output.h = h;
output.h_std = h_std;

output.ECS = ECS;
output.std_ECS = std_ECS;
output.ECS_boot = ECS_boot;

output.loglik = loglik;
output.BIC = BIC;

output.GoF = GoF;

output.std_prediction_residuals = std_prediction_residuals;