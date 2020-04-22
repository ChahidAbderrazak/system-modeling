%% ################################################################
%%                                Model based on Covid-19 
%% ################################################################
close all; clear all; clc;
addpath function
global Xx


beta=0.25;
gamma=0.03;
delta=0.00175;
kappa=0.95;          %% Lockdown parameter
DiffCoefS=0.0005;
DiffCoefI=0.0005;
DiffCoefR=0.0005;
DiffCoefD=0.0005;
DiffCoef=[DiffCoefS, DiffCoefI, DiffCoefR, DiffCoefD];

%% ################################################################
%%                 Space Discretization
%% ################################################################
Nx=4;
Lx=1.0;
dx=Lx/(Nx-1);
Xx=0:dx:Lx;
%% ################################################################
%%                 Time Discretization
%% ################################################################
N=10000;
Tf=360.0;
dt=Tf/(N-1);
t = [0:N-1].*dt;
%% ################################################################
%%                 Initializations 
%% ################################################################
S0=105529845*abs(cos(pi*Xx))';%ones(Nx,1); 
I0=13.6275*abs(cos(pi*Xx))';%ones(Nx,1);
R0=0.0*ones(Nx,1);
D0=0.0*ones(Nx,1);
X0=[S0 I0 R0 D0];

%% Generate the data from the model 
param=[DiffCoef,beta,gamma,delta,kappa];
[S,I,R,D, N] = SEIQRDP_pde(DiffCoef,beta,gamma,delta,kappa,X0,t);

%%  fitting of the simulated data   
DiffCoef0=[0, 0, 0, 0];
guess = [DiffCoef0, 0.5, 0.5, 0.5, 0.5]; % initial guess

[DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit] = fit_PDE(I,R,D,S,t,guess);

fit_param=[DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit];
%% 
[S_fit,I_fit,R_fit,D_fit, N_fit] = SEIQRDP_pde(DiffCoef_fit,beta_fit,gamma_fit,delta_fit,kappa_fit,X0,t);
data=[param;fit_param];


%% save summaries
T = array2table(data,'VariableNames',{'DiffCoefS','DiffCoefI','DiffCoefR','DiffCoefD','beta','gamma','delta','kappa'});%,'ReMSE','RMSE','computationTime'});
 
T.param=["real","fitted"]'

%% Plot the results 

Tt=0:dt:Tf;

 [X,Y]=meshgrid(Tt,Xx);
 
 
 
    figure;
    subplot(211)
    mesh(X,Y,S) 
    colormap(jet)
    legend('S(t) Susceptible individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
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
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
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
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
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
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,D_fit) 
    colormap(jet)
    legend('D(t) Death individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

             
             
     figure;
    subplot(211)
    mesh(X,Y,N) 
    colormap(jet)
    legend('N(t) total number of population')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,N_fit) 
    colormap(jet)
    legend('N(t) total number of population [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

  