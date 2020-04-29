%% ################################################################
%%                                Model based on Covid-19 
%% ################################################################
close all; clear all; clc;
addpath function
global Xx
alpha = 0.08; %  protection rate
beta = 0.9; %  infection rate
gamma= 1/2; % inverse of average latent time
delta= 1/8; % rate at which infectious people enter in quarantine
lambda = 0.03;%[0.03 0.05]; % cure rate (time dependant)
kappa =  0.03;%[0.03 0.05]; % mortality rate (time dependant)

% DiffCoefS=0.0005;
% DiffCoefI=0.0005;
% DiffCoefR=0.0005;
% DiffCoefD=0.0005;
% DiffCoef=[DiffCoefS, DiffCoefI, DiffCoefR, DiffCoefD];
DiffCoef=0%.0005
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
N=100;
Tf=360.0;
dt=Tf/(N-1);
t = [0:N-1].*dt;
%% ################################################################
%%                 Initializations 
%% ################################################################
S0=105529845*ones(Nx,1); %abs(cos(pi*Xx))';%ones(Nx,1); %
E0=13.6275*ones(Nx,1); %abs(cos(pi*Xx))';
I0=10.6275*ones(Nx,1); %abs(cos(pi*Xx))';
Q0=0.0*ones(Nx,1);
R0=0.0*ones(Nx,1);
D0=0.0*ones(Nx,1);
P0=0.0*ones(Nx,1);
X0=[S0 E0 I0 Q0 R0 D0 P0];
Npop=303e6;
%% Generate the data from the model 
param=[DiffCoef,alpha,beta,gamma,delta,lambda, kappa];
[S,E,I,Q,R,D,P] = SEIQRDP_pde(DiffCoef, alpha,beta,gamma,delta,lambda,kappa,Npop,X0,t)

%%  fitting of the simulated data   
DiffCoef0=[0];
guess = [DiffCoef0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; % initial guess

[DiffCoef_fit, alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit] = fit_PDE(I,Q,R,D,Npop,E,t,guess);

fit_param=[DiffCoef_fit,alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit];


[S_fit,E_fit,I_fit,Q_fit,R_fit,D_fit,P_fit] = SEIQRDP_pde(DiffCoef_fit, alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit,Npop,X0,t);

data=[param;fit_param];


% save summaries
T = array2table(data,'VariableNames',{'DiffCoef','alpha', 'beta','gamma','delta','lambda', 'kappa'});%,'ReMSE','RMSE','computationTime'});
 
T.param=["real","fitted"]'

% Plot the results 

Tt=0:dt:Tf;

 [X,Y]=meshgrid(Tt,Xx);
 

    figure;
    subplot(211)
    mesh(X,Y,S) 
    colormap(jet)
    legend('S(t) Susceptible individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),...'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,S_fit) 
    colormap(jet)
    legend('S(t) Susceptible individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),...'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

  
    figure;
    subplot(211)
    mesh(X,Y,E) 
    colormap(jet)
    legend('E(t) Susceptible individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),...'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,E_fit) 
    colormap(jet)
    legend('E(t) Susceptible individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),...'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

          
             
              
    figure;
    subplot(211)
    mesh(X,Y,I) 
    colormap(jet)
    legend('I(t) infected individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),...'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,I_fit) 
    colormap(jet)
    legend('I(t) infected individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),...'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

        
             
               
             
              
    figure;
    subplot(211)
    mesh(X,Y,Q) 
    colormap(jet)
    legend('Q(t) Recovered individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),...'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,Q_fit) 
    colormap(jet)
    legend('Q(t) Recovered individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),...'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))


             
             
             
    figure;
    subplot(211)
    mesh(X,Y,D) 
    colormap(jet)
    legend('D(t) Death individials')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef(1)),...'-',num2str(DiffCoef(2)),'-',num2str(DiffCoef(3)),'-',num2str(DiffCoef(4)),']',...
                 ', beta=',num2str(beta),', gamma=',num2str(gamma),', delta=',num2str(delta),', kappa=',num2str(kappa)))
    subplot(212)
    mesh(X,Y,D_fit) 
    colormap(jet)
    legend('D(t) Death individials [fitted]')
    title(strcat('parameters:DiffCoef=[',num2str(DiffCoef_fit(1)),...'-',num2str(DiffCoef_fit(2)),'-',num2str(DiffCoef_fit(3)),'-',num2str(DiffCoef_fit(4)),']',...
                 ', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

 