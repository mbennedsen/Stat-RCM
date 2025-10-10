function [p_boot,ECS] = bootstrap_stderrs_sqrtlog(t,B,x,estParams0,numIte,x00,P00,Opts)

conc_1750 =  278; % ppm
conc_1959 = 315.39; %ppm
C00 = conc_1750*2.127;
gammO = 0.1;
gammL = 0.1;

x_E       = x(:,1);
x_FnonCO2 = x(:,2);
x_Fnat    = x(:,3);

T = length(x_FnonCO2);
%% Set-up transition matrices
[T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,estParams0,x00,P00);

P000 = blkdiag(1e2*eye(6),0,1e2*eye(7));

sigmoid = @(x)( exp(x)./(1+exp(x)) );
%% Run bootstrap simulation
p_boot = nan(B,length(estParams0));
ECS = nan(B,1);
for iB = 1:B
    disp(['Bootstrap run ',num2str(iB),' out of ',num2str(B),'.']);
    
    warning('off');
    
    x0 = x00 + 0*sqrt(P000)*randn(length(x00),1);
    
    y = nan(T,size(H,2));
    eta = randn(size(Q,2),T);
    eps = randn(size(H,2),T);
    for i = 1:T
        if i == 1
            x = T_fct{i}(x0) + R*Q*eta(:,i);
        else
            x = T_fct{i}(x) + R*Q*eta(:,i);
        end
        
        y(i,:) = Z_fct(x) + H*eps(:,i);
    end


    save_loglik = nan(numIte,1);
    save_params = nan(numIte,35);
    save_StParams = nan(numIte,35);
    for iP = 1:numIte

        %% Setup param0: Randomized
        params0     = [abs(randn(31,1));sqrt(.1);sqrt(.1);3;-3];
        
        %%% Setup param0: Initial OLS estimates
        XX = y(1:end-1,1).*exp(-gammO*y(1:end-1,5))-C00;
        yy = y(2:end,2);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_O_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        params0(1) = sqrt(b_tmp(1));
    
        XX = y(1:end-1,1).*exp(-gammL*y(1:end-1,5))-C00;
        yy = y(2:end,3);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_L_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        params0(2) = sqrt(b_tmp(1));
    
        XX = [log(y(1:end-1,1))-log(C00),sqrt(y(1:end-1,1))-log(C00)];
        yy = y(2:end,4);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_F_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        muF_tmp = 0;
        params0(3) = sqrt(abs(b_tmp(1)));
        params0(4) = sqrt(abs(b_tmp(2)));
    
        XX = [ones(length(t),1),y(:,6)];
        yy = y(:,7);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_D_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t) - length(b_tmp));
        Cd_tmp = b_tmp(2);
      
        XX = [ones(length(t)-1,1),y(1:end-1,5)-y(1:end-1,6)];
        yy = diff(y(:,6));
        b_tmp = (XX'*XX)\XX'*yy;
        eta_OcT_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        gam_tmp = b_tmp(2)*Cd_tmp;
        a2_tmp = b_tmp(1);
    
    
        XX = [ones(length(t)-1,1),y(1:end-1,5),y(1:end-1,6),y(1:end-1,4)+x_Fnat(1:end-1) + x_FnonCO2(1:end-1)];
        yy = y(2:end,5);
        b_tmp = (XX'*XX)\XX'*yy;
        eta_T_tmp = sum( (yy - XX*b_tmp).^2 )/(length(t)-1 - length(b_tmp));
        Cm_tmp = 1/b_tmp(4);
        lam_tmp = Cm_tmp*(1-b_tmp(2)) - gam_tmp;
        a1_tmp = b_tmp(1);
    
        muT_tmp = Cm_tmp/lam_tmp*(a1_tmp + muF_tmp/Cm_tmp) + Cd_tmp/lam_tmp*a2_tmp;
        mud_tmp = Cm_tmp/lam_tmp*(a1_tmp + muF_tmp/Cm_tmp) + (gam_tmp+lam_tmp)/lam_tmp/gam_tmp*Cd_tmp*a2_tmp;
    
        cc = 4;
        params0(cc+1) = sqrt(abs(gam_tmp));
        params0(cc+2) = sqrt(abs(lam_tmp));
        params0(cc+3) = sqrt(abs(Cm_tmp));
        params0(cc+4) = sqrt(abs(Cd_tmp));
    
        cc = cc+4;
        params0(cc+1:cc+5) = sqrt(sqrt([eta_O_tmp;eta_L_tmp;eta_F_tmp;eta_T_tmp;eta_OcT_tmp]));
        cc = cc+5;
    
        cc = cc+7; % m-errors.
        
        params0(cc+1) = mud_tmp;
        params0(cc+2) = muT_tmp;
    
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

    %% Get best params
    [val,indx] = max(save_loglik);
    params0 = save_params(indx,:)';

    %% Estimate
    obj_fct = @(params)( -1*LogLik_EKF_Model15_sqrtlog_effi_sq_wARall_v45(y',x_E,x_FnonCO2,x_Fnat,params,x00,P00) );
    
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


    p_adjust = [est_params(1:20).^2;
                est_params(21:23);
                -1+2*sigmoid(est_params(24:31));
                est_params(32:33).^2;
                sigmoid(est_params(34));
                est_params(35).^2];

    p_boot(iB,:) = p_adjust;

    lambda = est_params(6)^2;
    ECS(iB) = 3.93/lambda;
end
