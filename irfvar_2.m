% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

function [IRF]=irfvar(A,B0inv,p,K,H)

J=[eye(K,K) zeros(K,K*(p-1))];
IRF=reshape(J*A^0*J'*B0inv,K^2,1);
for i=1:H
	IRF=([IRF reshape(J*A^i*J'*B0inv,K^2,1)]);
end;


