function [T_fct,Tp_fct,Q,Z_fct,Zp_fct,H,R,x0,P0] = getMat_EKF_sqrtlog(x_E,x_FnonCO2,x_Fnat,params,x0,P0)
% Model 14: OHC = Hd*Td + Hm*Tm


conc_1750 =  278; % ppm
%conc_1850 =  284.7; % ppm
conc_1959 = 315.39; %ppm
C0  = conc_1959*2.127;
C00 = conc_1750*2.127;
%% Load in parameters
a1 = 0;%-exp(params(1))*log(C00);
a2 = 0;%-exp(params(2))*log(C00);
b1 = (params(1))^2;
b2 = (params(2))^2;

f1 = (params(3))^2;%(params(3))^2;%/log(C00);
f2 = 0;%(params(4))^2;%/sqrt(C00);
f3 = (params(4))^2;%/sqrt(C00);

cc = 4;
gamma = (params(cc+1))^2;
lambda = (params(cc+2))^2;
Cm = (params(cc+3))^2;
Cd = (params(cc+4))^2;

cc = cc+4;
s_eta1 = (params(cc+1))^2;
s_eta2 = (params(cc+2))^2;
s_eta3 = (params(cc+3))^2;
s_eta4 = (params(cc+4))^2;
s_eta5 = (params(cc+5))^2;

cc = cc+5;
s_sig1 = (params(cc+1))^2;
s_sig2 = (params(cc+2))^2;
s_sig3 = (params(cc+3))^2;
s_sig4 = (params(cc+4))^2;
s_sig5 = params(cc+5)^2;%(params(cc+5))^2;
s_sig6 = params(cc+6)^2;%(params(cc+6))^2;
s_sig7 = params(cc+7)^2;%(params(cc+7))^2;

cc = cc+7;
mud  = params(cc+1);
muF  = 0;%params(cc+2);
muT  = params(cc+2);
muO  = params(cc+3);

cc = cc+3;

sigmoid = @(x)( 1./(1+exp(-x)) );
rho = -1+2*sigmoid(params(cc+1:cc+7));

RHO = diag(rho);

cc = cc+7;
rho_d = -1+2*sigmoid(params(cc+1)); % OHC corr with Ocean Temp

cc = cc+1;
g1 = params(cc+1)^2;%0.0423;
g2 = params(cc+2)^2;

ww = sigmoid(params(cc+3));

efficacy = 1+sigmoid(params(cc+4));
%% Construct matrices and functions
Z = [[eye(5),zeros(5,2)];
     zeros(1,4),1-ww,ww,0;
     zeros(1,4),Cm,Cd,muO];
Z(4,end) = muF;
Z(5,end) = muT;
Z(6,end) = mud;
%Z(7,end) = muO;

Z = [Z,eye(7)];

% Z = [Z,...
%     [1,zeros(1,5);
%      0,1,zeros(1,4);
%      0,0,1,zeros(1,3);
%      0,0,0,1,zeros(1,2);
%      zeros(1,6);
%      0,0,0,0,1,0;
%      zeros(1,5),1]];
for iT = 1:length(x_E)
    T = [1      ,0,0,0,0,0,x_E(iT)-a1-a2+(b1+b2)*C00 ;
             0      ,0,0,0,0,0,a1 - b1*C00 ;
             0      ,0,0,0,0,0,a2 - b2*C00 ;
             0      ,0,0,0,0,0,- f1*log(C00 + f2*C00^2) - f3*sqrt(C00);%f1     ,0,0,0,0,0,-C00*f1;
             0      ,0,0,1/Cm,   1-(lambda+efficacy*gamma)/Cm,  efficacy*gamma/Cm,   (x_FnonCO2(iT)+x_Fnat(iT))/Cm ;
             0      ,0,0,0,gamma/Cd, 1-gamma/Cd,0;
             zeros(1,6),1];
         
         T = blkdiag(T,RHO);
         
    T_fct{iT} = @(x)( [-b1*x(1)*exp(-g1*x(5))-b2*x(1)*exp(-g2*x(5)) + T(1,:)*x;...
                       b1*x(1)*exp(-g1*x(5)) + T(2,:)*x;...
                       b2*x(1)*exp(-g2*x(5)) + T(3,:)*x;...
                       T(4,:)*x + f1*log(x(1) + f2*x(1)^2) + f3*sqrt(x(1));...
                       T(5:end,:)*x] );
           
    Tp_fct{iT} = @(x)( [T(1,1) - b1*exp(-g1*x(5))-b2*exp(-g2*x(5)),T(1,2:4),g1*b1*x(1)*exp(-g1*x(5))+g2*b2*x(1)*exp(-g2*x(5)),T(1,6:end);...
                       b1*exp(-g1*x(5)),T(2,2:4),-g1*b1*x(1)*exp(-g1*x(5)),T(2,6:end);...
                       b2*exp(-g2*x(5)),T(3,2:4),-g2*b2*x(1)*exp(-g2*x(5)),T(3,6:end);...
                       f1/x(1)*((1+2*f2*x(1))/(1+f2*x(1))) + 0.5*f3/sqrt(x(1)),T(4,2:end);...
                       T(5:end,:)] );
         
         
   % A = [A;{A_tmp}];
end

H = zeros(7,7);%diag([s_sig1;s_sig2;s_sig3;s_sig4;s_sig5;s_sig6;s_sig7]);
%H(5,5) = s_sig5;

R = [-1,-1,0,0,0;
     eye(5);
     zeros(1,5)];
%R_tmp = blkdiag([1;1],[1;1],[1;1],[1;1],[1;1],[1;1],[1;1]); 
R = blkdiag(R,eye(7));

Qd = [s_sig6,0;
     s_sig7*rho_d, s_sig7*sqrt(1-rho_d^2)];
Q = blkdiag(diag([s_eta1;s_eta2;s_eta3;s_eta4;s_eta5;s_sig1;s_sig2;s_sig3;s_sig4;s_sig5]),Qd);

Z_fct  = @(x)( Z*x );
Zp_fct = @(x)( Z );
%R = eye( size(T,1) );

%%
if nargin < 4
    x0 = zeros(size(Q,1),1);
end
if nargin < 5
    P0 = 1e7*eye(size(Q));
end
