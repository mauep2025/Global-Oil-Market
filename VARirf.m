%estimate VAR impulse response function
%takes estimtes of slope coefficients and structural innovation variance matrix
%h is horizon
function [IRM, K]=VARirf(BETAnc,SIGMA,h)
  
  [K, n]=size(BETAnc);
  p=n/K;  %determine number of lags used in original estimation
      
  A=  [[BETAnc; eye(K*(p-1),K*(p-1)), zeros(K*(p-1),K)]];
  J=[eye(K,K) zeros(K,K*(p-1))];
  
  Theta=J*A^(0)*J'*chol(SIGMA)';
  IRM =reshape(Theta,K^2,1);
  for i=1:h;
      Theta=J*A^(i)*J'*chol(SIGMA)';
      eval(['IRM' num2str(i) '=reshape(Theta,K^2,1);']);
      eval(['IRM=[IRM, IRM' num2str(i) '];' ]); 
  end;

