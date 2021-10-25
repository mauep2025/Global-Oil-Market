function  SIRFres = Sirf(K,N,PHI,Impact,eslct)

%**************************************************************************
% PURPOSE: Computing STRUCTURAL impulse response functions 
%--------------------------------------------------------------------------
% INPUT:
% - K: total number of endogenous variables in the VAR
% - N: forecast horizon
% - PHI: three-dim matrix (K x K x (N+1)) containing dynamic multipliers of
% the VAR
% - Impact: K x K matrix containing impact coefficients of the VAR
% (contemporaneous IRF to all variable/shocks)
% - eslct: K x 1 array, selection vector for the shock to simulate
%--------------------------------------------------------------------------
% OUTPUT:
% - SIRFres: K x (N+1) matrix containing the STRUCTURAL impulse response
% functions
%--------------------------------------------------------------------------
%**************************************************************************

% preallocate matrix containing GIRFs
SIRFres = zeros(K,N+1);


for i=1:N+1
    SIRFres(:,i) = (PHI(:,:,i)*Impact*eslct);
end

