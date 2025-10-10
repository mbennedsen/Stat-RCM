function p_out = pInvTrans_sqrtlog(params)
% Inverse transformation of parameter vector

p_out = nan(length(params),1);
%% Load in parameters
p_out(1) = sqrt(params(1));
p_out(2) = sqrt(params(2));

p_out(3) = sqrt(params(3));
f2 = 0;
p_out(4) = sqrt(params(4)); 

cc = 4;
p_out(cc+1) = sqrt(params(cc+1));
p_out(cc+2) = sqrt(params(cc+2));
p_out(cc+3) = sqrt(params(cc+3));
p_out(cc+4) = sqrt(params(cc+4));

cc = cc+4;
p_out(cc+1) = sqrt(params(cc+1));
p_out(cc+2) = sqrt(params(cc+2));
p_out(cc+3) = sqrt(params(cc+3));
p_out(cc+4) = sqrt(params(cc+4));
p_out(cc+5) = sqrt(params(cc+5));

cc = cc+5;
p_out(cc+1) = sqrt(params(cc+1));
p_out(cc+2) = sqrt(params(cc+2));
p_out(cc+3) = sqrt(params(cc+3));
p_out(cc+4) = sqrt(params(cc+4));
p_out(cc+5) = sqrt(params(cc+5)); 
p_out(cc+6) = sqrt(params(cc+6)); 
p_out(cc+7) = sqrt(params(cc+7)); 

cc = cc+7;
p_out(cc+1)  = params(cc+1);
p_out(cc+2)  = params(cc+2);
p_out(cc+3)  = params(cc+3);

cc = cc+3;

sigmoid = @(x)( 1./(1+exp(-x)) );
invsigmoid = @(y)( log(y) - log(1-y) );
p_out(cc+1:cc+7) = invsigmoid( (params(cc+1:cc+7) +1)/2);

cc = cc+7;
p_out(cc+1) = invsigmoid( (params(cc+1) +1)/2); 

cc = cc+1;
p_out(cc+1) = sqrt(params(cc+1)); 
p_out(cc+2) = sqrt(params(cc+2)^2);

p_out(cc+3) = invsigmoid(params(cc+3));

p_out(cc+4) = invsigmoid(params(cc+4)-1); 
