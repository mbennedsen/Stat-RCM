function LogLik = EKF_Model1_getLogLik_v01(y,T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R,x0,P0)
% Function to compute Log-lik in EKF model.

%%
if nargin < 8
    R = eye( size(Q,1) );
end
if nargin < 9
    x0 = zeros(size(Q,1),1);
end
if nargin < 10
    P0 = 1e7*eye(size(Q));
end

N = length(y);

x_pred = nan(size(R,1),N);

LogLik = 0;
for iN = 1:N
    if iN == 1
        x_pred(:,iN) = T_fct{iN}(x0); 
        P_pred = P0;  
    else       
        x_pred(:,iN) = T_fct{iN}(a_tt);
        
        Tp = Tp_fct{iN}(a_tt);
        P_pred = Tp*P_tt*Tp' + R*Q*R';
    end
    
    Zp = Zp_fct(x_pred(:,iN)); 
    
    if isnan( sum(y(:,iN)) )   
        W = eye(length(y(:,iN)));
        W(isnan(y(:,iN)),:) = [];
        
        Zpstar = W*Zp;
        Hstar  = W*H*W';
        
        v_t = y(~isnan(y(:,iN)),iN) - W*Z_fct(x_pred(:,iN));
        F = Zpstar*P_pred*Zpstar' + Hstar;

        a_tt = x_pred(:,iN) + P_pred*Zpstar'*(F\v_t);
        P_tt = P_pred - P_pred*Zpstar'*(F\(Zpstar*P_pred));
    else
        v_t = y(:,iN) - Z_fct(x_pred(:,iN));
        F = Zp*P_pred*Zp' + H;

        a_tt = x_pred(:,iN) + P_pred*Zp'*(F\v_t);
        P_tt = P_pred - P_pred*Zp'*(F\(Zp*P_pred));
        
    end
    LogLik = LogLik - 0.5*log(det(F)) - 0.5*v_t'*(F\v_t);

end

