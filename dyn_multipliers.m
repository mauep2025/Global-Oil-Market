function PHI = dyn_multipliers(K,maxlag,F,N)

%**************************************************************************
% PURPOSE: Computing dynamic multipliers' matrices of the VAR
%--------------------------------------------------------------------------
% INPUT:
% - K: total number of endogenous variables in the VAR
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
% - F: three-dim matrix (K x K x maxlag) containing estimated 
% coefficients of the VAR (i.e. of reduced form representation)
% - N: forecast horizon
%--------------------------------------------------------------------------
% OUTPUT:
% - PHI: three-dim matrix (K x K x (N+1)) containing dynamic multipliers of
% the VAR
%--------------------------------------------------------------------------

%**************************************************************************

for i=1:maxlag
PHIx(:,:,i) = zeros(K); %#ok
end
PHIx(:,:,maxlag+1) = eye(K);

for t=maxlag+2:maxlag+N+1
    acc = 0;
    for j=1:maxlag
        acc = acc + F(:,:,j)*PHIx(:,:,t-j);
    end
    PHIx(:,:,t) = acc;
end

% reindicize
PHI = PHIx(:,:,maxlag+1:end);

