%////////////////////////////////////////////////////////////////%
%//////- Structural-VAR (Assessment Idetificatiojn) -     ///////%
%////////////////////////////////////////////////////////////////%

clc;
clear;

%--- 1.- Shock identification with oil prices, oil production and an index of global economic activity
z1=xlsread('Data_oil_1.xlsx');
 time=(1973+1/12:1/12:2019)';  
RAC   = [z1(:,1)]; % US Refiner Acquisition Cost of Crude Oil
WTI   = [z1(:,2)]; % West Texas Intermediate Oil Price
Oil_p = [z1(:,3)]; % Global Oil production 
Oil_i = [z1(:,4)]; % OECD oil invetories
Kil_i = [z1(:,5)]; % A proxy of Global economic activity 
Ham_i = [z1(:,6)]; % 2 years growth in global industry production 

%Montly Percentage change in global crude oil production
  [T,~]=size(Oil_p);
    for i =1: size(Oil_p,2)
        for ii=2 : size(Oil_p,1)
        goil_p(ii-1,i)=((Oil_p(ii,i)-Oil_p(ii-1,i)))*100;
        end
    end
% Anual Percentage change in global crude oil production
   gaoil_p=100*((Oil_p(13:T,1)./Oil_p(1:T-12,1))-1);  
  
% Information assemble. The sample goes from 1973m2 to  2006m12
   z2= [goil_p Kil_i(2:T,:) RAC(2:T,:)];
     
%z3= table(goil_p Kil_i[2:T,:] WTI[2:T,:]); % Information assemble

% VAR EStimation 
   numseries=3;
   seriesnam={'Global oil production', 'Global Economic Activity', 'Oil price' }; 
   VAR1f=varm(numseries,24); % the number is the lag 
   VAR1f.SeriesNames=seriesnam;
   [EstMdl1,EstSE1,logL1,E1]= estimate(VAR1f, z2);
   EstMdl1.Description
   results1=summarize(EstMdl1);
   np1=results1.NumEstimatedParameters;
   summarize(EstMdl1);
   Sig=cov(E1)

%--------Fiting the VAR model (Var Estimation using a toolbox (reduced-Form
%model))

  %  a) We estimate the reduce-form VAR parameters and compute the residual
  %  variance-covariance matrix (Sigma)
  %  b) Then we estimate the structural impact multiplier matrix inv(B) based
  %  on the lower-triangular Cholesky decomposition of the residual
  %  variance-covaraince matrix (sigma) 
pp = 24;
hh = 18;
 
[T, N] = size(z2);
Horiz=24;
[AR_3d,Chol_Var] = VAR_OLS(z2,pp,1,[]); 
mm=size(z2,2);
Ai_mat = dyn_multipliers(mm,pp,AR_3d,Horiz);
%Generating the companion matrix and the residuals
%[t,K]=size(z2);
[A_2,SIGMA_2,Uhat_1,V_2,X_2]=olsvarc(z2,pp);
 BETAnc_2=A_2(1:N,:)
%[BETAnc,B1,XZ, SIGMA, U, V]=lsvarc(z2,pp);
global h q 
B0inv_1=chol(Sig(1:q,1:q))';


%------- 1.i.- VAR Identification (Structural Impulse Response)
%nvar=size(Ai_mat,1);
%Q=eye(nvar);
Shock = zeros(size(z2,2),1); Shock(1,1) = 1; %Shock for oil supply
Shock_1 = zeros(size(z2,2),1); Shock_1(2,1) = 1; %Shock for global demand
Shock_2 = zeros(size(z2,2),1); Shock_2(3,1) = 1; %Shock for oil price

%Computing the structural Impulse Response 
SIRF = Sirf(N,hh,Ai_mat,B0inv_1,Shock)'; % Oil supply
SIRF_1 = Sirf(N,hh,Ai_mat,B0inv_1,Shock_1)'; % Global demand
SIRF_2 = Sirf(N,hh,Ai_mat,B0inv_1,Shock_2)'; %Oil price 
SIRFaer=[SIRF(1:18,:) SIRF_1(1:18,:) SIRF_2(1:18,:)];

no_boot=1000;  
%- Error bands corresponding ro the Oil supply shocks
zirf_b=NaN*zeros(size(SIRF,1),size(SIRF,2),no_boot);
for boot=1:1:no_boot
    yb = bootstrapVAR(z2,pp,floor(.25*size(z2,1)));
    [AR_3d_b,Chol_Var_b] = VAR_OLS(yb,pp,1,[]); 
    Ai_mat_b = dyn_multipliers(N,pp,AR_3d_b,hh);
    zirf_b(:,:,boot) = Sirf(N,hh,Ai_mat_b,Chol_Var_b,Shock)';    
end

low_perc=.05; 
up_perc=.95;

for ii=1:size(zirf_b,1)
    for jj=1:size(zirf_b,2)        
        zirf_low_oils(ii,jj)=quantile(zirf_b(ii,jj,:),low_perc);
        zirf_up_oils(ii,jj)=quantile(zirf_b(ii,jj,:),up_perc);
    end
end
%--
%- Error bands corresponding ro the Demand shocks
zirf_b_d=NaN*zeros(size(SIRF_1,1),size(SIRF_1,2),no_boot);
for boot=1:1:no_boot
    yb = bootstrapVAR(z2,pp,floor(.25*size(z2,1)));
    [AR_3d_b,Chol_Var_b] = VAR_OLS(yb,pp,1,[]); 
    Ai_mat_b = dyn_multipliers(N,pp,AR_3d_b,hh);
    zirf_b_d(:,:,boot) = Sirf(N,hh,Ai_mat_b,Chol_Var_b,Shock_1)';    
end

low_perc=.05; 
up_perc=.95;

for ii=1:size(zirf_b_d,1)
    for jj=1:size(zirf_b_d,2)        
        zirf_low_demand(ii,jj)=quantile(zirf_b_d(ii,jj,:),low_perc);
        zirf_up_demand(ii,jj)=quantile(zirf_b_d(ii,jj,:),up_perc);
    end
end
%--
zirf_b_p=NaN*zeros(size(SIRF_2,1),size(SIRF_2,2),no_boot);
for boot=1:1:no_boot
    yb = bootstrapVAR(z2,pp,floor(.25*size(z2,1)));
    [AR_3d_b,Chol_Var_b] = VAR_OLS(yb,pp,1,[]); 
    Ai_mat_b = dyn_multipliers(N,pp,AR_3d_b,hh);
    zirf_b_p(:,:,boot) = Sirf(N,hh,Ai_mat_b,Chol_Var_b,Shock_2)';    
end

low_perc=.05; 
up_perc=.95;

for ii=1:size(zirf_b_p,1)
    for jj=1:size(zirf_b_p,2)        
        zirf_low_p(ii,jj)=quantile(zirf_b_p(ii,jj,:),low_perc);
        zirf_up_p(ii,jj)=quantile(zirf_b_p(ii,jj,:),up_perc);
    end
end
%--

figure(1) 
subplot(3,3,1); 
plot([0:hh]',[SIRF(:,1)],'linewidth',2)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_oils(:,1) zirf_up_oils(:,1)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-1  1.5]);
ylabel('Oil production','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,4); 
plot([0:hh]',[SIRF(:,2)],'linewidth',2)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_oils(:,2) zirf_up_oils(:,2)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-5  4]);
ylabel('Oil production','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,7); 
plot([0:hh]',[SIRF(:,3)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_oils(:,3) zirf_up_oils(:,3)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-.05  .05]);
ylabel('Oil production','fontsize',11)
title('Oil-specific demand shock','fontsize',11)
subplot(3,3,2); 
plot([0:hh]',[SIRF_1(:,1)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_demand(:,1) zirf_up_demand(:,1)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-.4  .4]);
ylabel('Real activity','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,5); 
plot([0:hh]',[SIRF_1(:,2)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_demand(:,2) zirf_up_demand(:,2)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-5  17]);
ylabel('Real activity','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,8); 
plot([0:hh]',[SIRF_1(:,3)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_demand(:,3) zirf_up_demand(:,3)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-.05  .07]);
ylabel('Real activity','fontsize',11)
title('Oil-specific demand shock','fontsize',11)
subplot(3,3,3); 
plot([0:hh]',[SIRF_2(:,1)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_p(:,1) zirf_up_p(:,1)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-.5  .5]);
ylabel('Real oil price','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,6); 
plot([0:hh]',[SIRF_2(:,2)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_p(:,2) zirf_up_p(:,2)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-7  7]);
ylabel('Real oil price','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,9); 
plot([0:hh]',[SIRF_2(:,3)],'linewidth',3)
hold; plot([0:hh]',[zeros(hh+1,1) zirf_low_p(:,3) zirf_up_p(:,3)],'linewidth',1,'LineStyle','--','color','k'); 
xlim([0 hh])
ylim([-.1  .2]);
title('Oil-specific demand shock','fontsize',11)
ylabel('Real oil price','fontsize',11)
sgtitle('Graph: Structural Impulse Responses' )


%--- Impulse Responses (acording with the paper of Lutz kilian (2009)
xmax=18; 
HO=(0:1:xmax-1)'; 
[IRFaer, K1]=VARirf(BETAnc_2,Sig,xmax-1);
% Saving the structural Impulse responses
ir_aer=reshape(IRFaer,3,3,18);

%Computing the graphs 
 %this are the results that the paper presents
zirf_1_1=-cumsum(squeeze(ir_aer(1,1,:)));
zirf_2_1=cumsum(squeeze(ir_aer(1,2,:)));
zirf_3_1=cumsum(squeeze(ir_aer(1,3,:)));
zirf_4_1=-squeeze(ir_aer(2,1,:));
zirf_5_1=squeeze(ir_aer(2,2,:));
zirf_6_1=squeeze(ir_aer(2,3,:));
zirf_7_1=-squeeze(ir_aer(3,1,:));
zirf_8_1=squeeze(ir_aer(3,2,:));
zirf_9_1=squeeze(ir_aer(3,3,:));

zirf_1=[zirf_1_1 zirf_2_1 zirf_3_1];
zirf_2=[zirf_4_1 zirf_5_1 zirf_6_1];
zirf_3=[zirf_7_1 zirf_8_1 zirf_9_1];
 

figure(2)
subplot(3,3,1) 
box on; 
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,-cumsum(squeeze(ir_aer(1,1,:))),'r:','linewidth',2)
axis([0 xmax-1 -2 1.5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:1:1.5)
ylabel('Oil production','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,4)
box on; 
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,cumsum(squeeze(ir_aer(1,2,:))),'r:','linewidth',2)
axis([0 xmax-1 -.5 .5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:1:1.5)
ylabel('Oil production','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,7)
box on;
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,cumsum(squeeze(ir_aer(1,3,:))),'r:','linewidth',2)
axis([0 xmax-1 -.3 .3])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:1:1.5)
ylabel('Oil production','fontsize',11)
title('Oil-specific demand shock','fontsize',11)
xlabel('Months','fontsize',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,2)
box on;
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,-squeeze(ir_aer(2,1,:)),'r:','linewidth',2)
axis([0 xmax-1 -2 2])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:10)
ylabel('Real activity','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,5) 
box on;
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,squeeze(ir_aer(2,2,:)),'r:','linewidth',2) 
axis([0 xmax-1 -5 17])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:10)
ylabel('Real activity','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,8) 
box on; 
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,squeeze(ir_aer(2,3,:)),'r:','linewidth',2) 
axis([0 xmax-1 -5 5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:35)
ylabel('Real activity','fontsize',11)
title('Oil-specific demand shock','fontsize',11)
xlabel('Months','fontsize',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,3)
box on; 
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,-squeeze(ir_aer(3,1,:)),'r:','linewidth',2)
axis([0 xmax-1 -.05 .05])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:10)
ylabel('Real oil price','fontsize',11)
title('Oil supply shock','fontsize',11)
subplot(3,3,6) 
box on;
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,squeeze(ir_aer(3,2,:)),'r:','linewidth',2)
axis([0 xmax-1 -.05 .05])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:10)
ylabel('Real oil price','fontsize',11)
title('Aggregate demand shock','fontsize',11)
subplot(3,3,9) 
box on; 
plot(HO,zeros(xmax,1),'k:','linewidth',1), hold on, plot(HO,squeeze(ir_aer(3,3,:)),'r:','linewidth',2)
axis([0 xmax-1 -.1 .1])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-5:5:10)
title('Oil-specific demand shock','fontsize',11)
ylabel('Real oil price','fontsize',11)
xlabel('Months','fontsize',10)
sgtitle('Graph: Structural Impulse Responses' )

%--- 1.ii.- VAR Identification (Structural Forecast Error Decomposition)
% How much of the forecast error variance or prediction mean squared error
% (MSPE) of Y_t+h horizon h=0,1....,H is accounted for by each structural
% shock w_kt, where k=1,..,K 

%We need to compute the structural response matrices
%The contribution of the shock j to the MSPE of y_kt


[t,K]=size(z2);
SIGMA_1=SIGMA_2(1:K,1:K); %Reshaping the SIGMA from the impulse response function

sfed=[1 2 3 10 15 20 600]' %This are the number of forecast horizon for the Error Decomposition
[r_1,r]=size(sfed);
sfed_suply=zeros(r_1,K) %This will store the Error decomposition of supply shocks
sfed_demand=zeros(r_1,K) %This will store the Error decomposition of demand
sfed_price=zeros(r_1,K) %This will store the Error decomposition of prices

J=[eye(K,K) zeros(K,K*(pp-1))];
TH1=J*A_2^0*J'; 
TH=TH1*chol(SIGMA_1)'; 
TH=TH'; TH2=(TH.*TH);
TH3=TH2;
%for r=1=size(sfed,1);
VC=zeros(K,K);

for i=2:sfed(1,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; 
    TH=TH'; 
    TH2=(TH.*TH);
    TH3=TH3+TH2;
end;
TH4=sum(TH3);

for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;

% VDC in percentage terms at horizon h
% Columns refer to shocks j=1,...,q that explain any given variable
% Rows refer to variables whose variation is to be explained
disp(VC'*100)
% Forecast error variance decomposition at horizon h for real dividen
% growth
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)
sfed_suply(1,:)=VC(:,1)*100
sfed_demand(1,:)=VC(:,2)*100
sfed_price(1,:)=VC(:,end)*100

%-----------------------------
for i=2:sfed(3,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;

disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)

sfed_suply(3,:)=VC(:,1)*100
sfed_demand(3,:)=VC(:,2)*100
sfed_price(3,:)=VC(:,end)*100
%---------------------------------
for i=2:sfed(2,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;

disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)
sfed_suply(2,:)=VC(:,1)*100
sfed_demand(2,:)=VC(:,2)*100
sfed_price(2,:)=VC(:,end)*100
%-------------------------------------------
for i=2:sfed(4,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;
disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)
sfed_suply(4,:)=VC(:,1)*100
sfed_demand(4,:)=VC(:,2)*100
sfed_price(4,:)=VC(:,end)*100
%---------------------------------------

for i=2:sfed(5,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;
disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)
sfed_suply(5,:)=VC(:,1)*100
sfed_demand(5,:)=VC(:,2)*100
sfed_price(5,:)=VC(:,end)*100
%---------------------------------------------------------

for i=2:sfed(6,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;
disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)
sfed_suply(6,:)=VC(:,1)*100
sfed_demand(6,:)=VC(:,2)*100
sfed_price(6,:)=VC(:,end)*100
%----------------------------------------------

for i=2:sfed(7,1)
    TH=J*A_2^(i-1)*J'*chol(SIGMA_1)'; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
end;
TH4=sum(TH3);
VC=zeros(K,K);
for j=1:K
    VC(j,:)=TH3(j,:)./TH4;
end;
disp(VC'*100)
disp('Row h of Table 4.1')
disp(VC(:,end)'*100)

%Storing the results from the Forecast error variance decomposition 
sfed_suply(7,:)=VC(:,1)*100
sfed_demand(7,:)=VC(:,2)*100
sfed_price(7,:)=VC(:,end)*100

disp(' Forecast Error Variance Decomposition for Oil Production');
disp('Percent of h-step Ahead Forecast Error Variance Explained by');
disp('Oil Supply Shock   Aggregate Demand Shock   Oil-specific Demand Shock');
disp (sfed_suply);

disp(' Forecast Error Variance Decomposition for Global Economic Activity');
disp('Percent of h-step Ahead Forecast Error Variance Explained by');
disp('Oil Supply Shock   Aggregate Demand Shock   Oil-specific Demand Shock');
disp (sfed_demand);

disp(' Forecast Error Variance Decomposition for Real Price of Oil ');
disp('Percent of h-step Ahead Forecast Error Variance Explained by');
disp('Oil Supply Shock   Aggregate Demand Shock   Oil-specific Demand Shock')
disp (sfed_price);


%--- 1.iii.- VAR Identification (Historical decomposition pf the oil price for the period 2000-2019)
% we are interested in quantifying how much a given structural shock
% explains of the historically observed fluctations in VAr variables. In
% other words, we would like to know the cumulative effect of a given
% structural shock on each variable at every given point in time. 

%The historical decomposition: 1) Compute the structural MA coefficientes
%matrices. 2) Compute the structural shocks. 3) Match up each structural
%shock, say shock j, with the appropiate impulse response weight, as
%required by the structural moving average representation, to form
%T*1 vectors of fitted values for variables k. 

%1)Computing the structural MA coeficients matrices. Also Cpunting the
%companion matrix A and invB(structural impact multiplier). 

% Where A denotes the matrix of slope parametrers obtained by expresing the
% VAR(24) (in its VAR(1) companion format 

horizon=0:h;
[IRF_1]=irfvar_2(A_2,B0inv_1,Horiz,N,T-Horiz-1);
% Compute structural shocks What from reduced form shocks Uhat
What=inv(B0inv_1)*E1';

% Historical evolution of the structural shocks(Folowing the Paper of
% Kilian(2009))
time=1973+2/12+Horiz/12:1/12:2019+6/12; 
time_1=1974+2/12+Horiz/12:1/12:2019+6/12;
figure(3)
subplot(3,1,1)
plot(time,What(1,:),'k-','linewidth',3);
title('Oil Supply shock','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -5 5])
grid on

subplot(3,1,2)
plot(time,What(2,:),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -5 5])
grid on

subplot(3,1,3)
plot(time,What(3,:),'k-','linewidth',3);
title('Oil-Specific demand shock','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -5 5])
grid on
sgtitle('Graph: Historical Evolution of the Structural Shocks (2000-2019)' )
% Cross-multiply the weights for the effect of a given shock on the real
% oil price (given by the relevant row of IRF) with the structural shock
% in question

yhat1=zeros(T-Horiz,1);
yhat2=zeros(T-Horiz,1);
yhat3=zeros(T-Horiz,1);
yhat4=zeros(T-Horiz,1);
yhat5=zeros(T-Horiz,1);
yhat6=zeros(T-Horiz,1);
yhat7=zeros(T-Horiz,1);
yhat8=zeros(T-Horiz,1);
yhat9=zeros(T-Horiz,1);

for i=1:T-Horiz
    yhat1(i,:)=dot(IRF_1(1,1:i),What(1,i:-1:1));
    yhat2(i,:)=dot(IRF_1(2,1:i),What(1,i:-1:1));
    yhat3(i,:)=dot(IRF_1(3,1:i),What(1,i:-1:1));  
    yhat4(i,:)=dot(IRF_1(4,1:i),What(2,i:-1:1)); 
    yhat5(i,:)=dot(IRF_1(5,1:i),What(2,i:-1:1)); 
    yhat6(i,:)=dot(IRF_1(6,1:i),What(2,i:-1:1)); 
    yhat7(i,:)=dot(IRF_1(7,1:i),What(3,i:-1:1)); 
    yhat8(i,:)=dot(IRF_1(8,1:i),What(3,i:-1:1));
    yhat9(i,:)=dot(IRF_1(9,1:i),What(3,i:-1:1));
    
end;


zhd=-cumsum(squeeze(yhat3))
zhd=squeeze((yhat3))

% Ploting the historical decomposition corresponding to the oil supply Shock 
%here the sample will goes from 1975m1 to 2006m12

figure(6)
subplot(3,1,1)
plot(time,yhat1,'k-','linewidth',3);
title('Cumulative effect of Oil Supply Shock on Oil Production','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -10 10])
grid on

subplot(3,1,2)
plot(time,yhat2,'k-','linewidth',3);
title('Cumulative effect of  Oil Supply Shock on Real Economic Activity','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -10 10])
grid on

subplot(3,1,3)
plot(time,yhat3,'k-','linewidth',3);
title('Cumulative effect of  Oil Supply Shock on Real price of Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -.2 .2])
grid on
sgtitle('Graph: Historical Decomposition (2000-2019)' )

% Ploting the historical decomposition corresponding to the Aggregate Demand Shock 
%here the sample will goes from 1975m1 to 2006m12

figure(7)
subplot(3,1,1)
plot(time,yhat4,'k-','linewidth',3);
title('Cumulative effect of Aggregate Demand Shock on Oil Production','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -1.5 1.5])
grid on

subplot(3,1,2)
plot(time,yhat5,'k-','linewidth',3);
title('Cumulative effect of Aggregate Demand Shock on Real Econommic Activity','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -200 230])
grid on

subplot(3,1,3)
plot(time,yhat6,'k-','linewidth',3);
title('Cumulative effect of Aggregate Demand Shock on Real price of Crude Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -1 1])
grid on
sgtitle('Graph: Historical Decomposition (2000-2019)' )

% Ploting the historical decomposition corresponding to Oil supply Shock 
%here the sample will goes from 1975m1 to 2006m12

figure(8)
subplot(3,1,1)
plot(time,yhat7,'k-','linewidth',3);
title('Cumulative effect of Oil-market specific demand Shock on Oil Production','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -1.5 2])
grid on

subplot(3,1,2)
plot(time,yhat8,'k-','linewidth',3);
title('Cumulative effect of Oil-market Specific Demand Shock on Real Economic Acivity','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -70 80])
grid on

subplot(3,1,3)
plot(time,yhat9,'k-','linewidth',3);
title('Cumulative effect of Oil-market specific demand Shock on Real price of Crude Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -1 1])
grid on
sgtitle('Graph: Historical Decomposition (2000-2019)' )



%Acording with the paper of Kilian (2009)
figure(9)
subplot(3,1,1)
plot(time,yhat3,'k-','linewidth',3);
title('Cumulative effect of  Oil Supply Shock on Real price of Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -.3 .3])
grid on

subplot(3,1,2)
plot(time,yhat6,'k-','linewidth',3);
title('Cumulative effect of Aggregate Demand Shock on Real price of Crude Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+2/12 2019+6/12 -1 1])
grid on

subplot(3,1,3)
plot(time,yhat9,'k-','linewidth',3);
title('Cumulative effect of Oil-market specific demand Shock on Real price of Crude Oil','fontsize',18)
%ylabel('Percent','fontsize',18)
axis([2000+1/12 2019+6/12 -1 1])
grid on
sgtitle('Graph: Historical Decomposition (2000-2019)' )
    


