function [x_pred, x_filter, P, x_smooth, FF] = EKF_Model1_v01(y,T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R,x0,P0)
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
if nargout > 1
    x_filter = nan(size(R,1),N);
end
if nargout > 2
    P = nan(size(R,1),size(R,1),N);
end

if nargout > 4
    FF = nan(size(y,1),size(y,1),N);
end

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

        % if iN == 1
        %     W
        %     Zpstar
        %     v_t
        %     P_pred
        % end
    else
        v_t = y(:,iN) - Z_fct(x_pred(:,iN));
        F = Zp*P_pred*Zp' + H;

        a_tt = x_pred(:,iN) + P_pred*Zp'*(F\v_t);
        P_tt = P_pred - P_pred*Zp'*(F\(Zp*P_pred));
        
    end
    LogLik = LogLik - 0.5*log(det(F)) - 0.5*v_t'*(F\v_t);
    
    if nargout > 1
        x_filter(:,iN) = a_tt;
    end
    if nargout > 2
        P(:,:,iN) = P_pred;
    end
    if nargout > 4
        FF(:,:,iN) = F;
    end
end



%% Smoothing using DK p 91 (only smooth state vector)

if nargout > 3
    x_smooth = nan(size(R,1),N);
    r = zeros(size(R,1),1);
    for iN = N:-1:1
        
        Zp = Zp_fct(x_pred(:,iN));      
        y_pred = Z_fct(x_pred(:,iN));
        T = T_fct{iN}(x_pred(:,iN));
        Tp = Tp_fct{iN}(x_filter(:,iN));
        
        if isnan( sum(y(:,iN)) )   
            W = eye(length(y(:,iN)));
            W(isnan(y(:,iN)),:) = [];

            Zpstar = W*Zp;
            Hstar  = W*H*W';  
            
            

            %K = T*P(:,:,iN)*Z'*inv(F);
            F = Zpstar*P(:,:,iN)*Zpstar' + Hstar;
            
            %size(Zpstar)
            %size(F)
            %size(T)
            
            K = (F\(Zpstar*P(:,:,iN)'*Tp'))';
            %K = (F\(Zpstar*P(:,:,iN)'*T))';
            v_t = y(~isnan(y(:,iN)),iN) - W*y_pred;
            L = Tp-K*Zpstar;

            r   = Zpstar'*(F\v_t) + L'*r; %r_{t-1}
            x_smooth(:,iN) = x_pred(:,iN) + P(:,:,iN)*r;
        else
            %K = T*P(:,:,iN)*Z'*inv(F);
            F = Zp*P(:,:,iN)*Zp' + H;
            K = (F\(Zp*P(:,:,iN)'*Tp'))';
            v_t = y(:,iN) - y_pred;
            L = Tp-K*Zp;
           
            r   = Zp'*(F\v_t) + L'*r; %r_{t-1}
            x_smooth(:,iN) = x_pred(:,iN) + P(:,:,iN)*r;
            

        end    
    end
    %x_smooth = x_smooth';
end

% x_pred = x_pred';
% if nargout > 1
%     x_filter = x_filter';
% end
