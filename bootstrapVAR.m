function yb = bootstrapVAR(y,p,k)

% y : series
% p : order of lags
% k : # of initial trimming observation


[nobs,nvar] = size(y);

YY = y(p+1:nobs,:);
XX = ones(nobs-p,1+nvar*p);

for i = 1:p
   XX(:,2+nvar*(i-1):1+nvar*i) = y(p+1-i:nobs-i,:);
end

b = inv(XX'*XX)*XX'*YY;
e = YY - XX*b;

% generating Pseudo-Sample

% Pseudo-Disturbance
T = size(YY,1);
segment = (1:T)/T;
eb = zeros(T+k+p,nvar);

for i=1:T+k+p
 u=rand(1,1);
 eb(i,:) = e(min(find(segment>=u)),:);
end

% Pseudo-Sample
yb = zeros(T+k+p,nvar);
r = XX(1,:);
for i=1:T+k+p
  yb(i,:) = r*b + eb(i,:);
  r = [1,yb(i,:),r(2:end-nvar)];
end

% trim dataset
yb = yb(k+1:end,:);


