function [AR_3d,Chol_Var,ee] = VAR_OLS(series,p,const,ex_var)
%
%  VAR_OLS_EST estimates a Vector Autoregression of order p [VAR(p)]
%  
%   [b,Sig] = VAR_OLS_EST(series,p,const,ex_var)
%   b: VAR coefficient vector [A1,...,Ap]' 
%   consider the VAR(p) without any constant/trend: 
%       Y(t)=A1*Y(t-1)+...+Ap*Y(t-p)+e(t)
%   	Y=X*b+e 
%   Sig: Variance-Covariance Matrix of Reduced form error
%   theta: Cholesky decomposition of Cov Matrix,
%          i.e. Sigma=theta*theta';
%   
%   p: VAR order
%   const: 0 - NO const/ 1 - constant/ 2 - linear trend
%          
%   ex_var: include possibly exogenous variables
%           [] no exogenous variables
  



[Tnobs,nvar] = size(series); 

    CC=[];
    if const==1
        CC=[CC ones(Tnobs,1)];
    elseif const==2 
        trend=[1:1:Tnobs]';
        CC=[CC ones(Tnobs,1) trend];
    end
    CC=[CC ex_var];
    
[Xlag] = mlag(series,p);
YY = series(p+1:end,:);
XX = [CC(p+1:end,:) Xlag(p+1:end,:)];

% b = XX\YY;
b = inv(XX'*XX)*(XX'*YY);
e = YY - XX*b;
% Sig = (e'*e)/(Tnobs-p);
Sig = cov(e);

AR = b(size(CC,2)+1:end,:)';            % AR coefficient [A1,A2,...,Ap]

AR_3d=NaN*zeros(nvar,nvar,p);
for ii=1:1:p
    AR_3d(:,:,ii)=AR(:,(ii-1)*nvar+1:ii*nvar);
end
Chol_Var = chol(Sig)';

ee = e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xlag] = mlag(X,p)
%MLAG Summary of this function goes here
%   Detailed explanation goes here
[Traw,N]=size(X);
Xlag=zeros(Traw,N*p);
for ii=1:p
    Xlag(p+1:Traw,(N*(ii-1)+1):N*ii)=X(p+1-ii:Traw-ii,1:N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






