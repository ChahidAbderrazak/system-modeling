% %% ################################################################
% %%                      PDE model based on Covid-19 
% %% ################################################################
% close all; clear all; clc;
% addpath function
% global Xx
% %%
% [time, table_data]=retreive_data_china();
% % position
% % R, D, I
R=table_data.Recovered;
D=table_data.Deaths;
I=table_data.Confirmed;
position=table_data.r;

%%
Nt=30e6;
S=Nt-(R+I+D);
% use small portion of point 
% 
% S=S(1:4,:);
% I=I(1:4,:);
% R=R(1:4,:);
% D=D(1:4,:);
%% 

S0=S(:,1);
I0=I(:,1);
R0=R(:,1);
D0=D(:,1);
X0=[S0 I0 R0 D0];

%% buid the mesh 
Nx=size(position,1);
Lx=1.0;
dx=Lx/(Nx-1);
Xx=0:dx:Lx;

N=size(time,2);
Tf=360.0;
dt=Tf/(N-1);
t = [0:N-1].*dt;

%%  fitting of the simulated data   
DiffCoef0=[0, 0, 0, 0];
guess = [DiffCoef0, 0.5, 0.5, 0.5, 0.5]; % initial guess

[DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit] = fit_PDE(I,R,D,S,t,guess);

fit_param=[DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit];
%% 
[S_fit,I_fit,R_fit,D_fit, N_fit] = SEIQRDP_pde(DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit,X0,t);


data=[fit_param];


%% save summaries
T = array2table(data,'VariableNames',{'DiffCoefS','DiffCoefI','DiffCoefR','DiffCoefD','beta','gamma','delta','kappa'});%,'ReMSE','RMSE','computationTime'});
 
T.param=["fitted"]'

%% Plot the results 

Tt=0:dt:Tf;

 [X,Y]=meshgrid(Tt,Xx);
 
 
 
    figure;
    subplot(211)
    mesh(X,Y,S) 
    colormap(jet)
    legend('S(t) Susceptible individials')

    subplot(212)
    mesh(X,Y,S_fit) 
    colormap(jet)
    legend('S(t) Susceptible individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

  
             
             
              
    figure;
    subplot(211)
    mesh(X,Y,I) 
    colormap(jet)
    legend('I(t) infected individials')
    
    subplot(212)
    mesh(X,Y,I_fit) 
    colormap(jet)
    legend('I(t) infected individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

             
              
    figure;
    subplot(211)
    mesh(X,Y,R) 
    colormap(jet)
    legend('R(t) Recovered individials')
    
    subplot(212)
    mesh(X,Y,R_fit) 
    colormap(jet)
    legend('R(t) Recovered individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))


             
    figure;
    subplot(211)
    mesh(X,Y,D) 
    colormap(jet)
    legend('D(t) Death individials')
    
    subplot(212)
    mesh(X,Y,D_fit) 
    colormap(jet)
    legend('D(t) Death individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

             
             
%      figure;
%     subplot(211)
%     mesh(X,Y,Nt) 
%     colormap(jet)
%     legend('N(t) total number of population')
%     
%     subplot(212)
%     mesh(X,Y,N_fit) 
%     colormap(jet)
%     legend('N(t) total number of population [fitted]')
%     title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
%                  ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))
% 
%   