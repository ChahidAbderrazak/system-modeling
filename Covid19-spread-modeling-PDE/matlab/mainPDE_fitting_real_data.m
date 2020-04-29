% %################################################################
% %                      PDE model based on Covid-19 
% % ################################################################
close all; clear all; clc;
addpath function
global Xx
%%
Location='China'
[time, table_data]=retreive_data(Location);

% Initial conditions
E = table_data.Confirmed; % Initial number of exposed cases. Unknown but unlikely to be zero.
I = table_data.Confirmed; % Initial number of infectious cases. Unknown but unlikely to be zero.
Q = table_data.Confirmed-table_data.Recovered-table_data.Deaths;
R = table_data.Recovered;
D = table_data.Deaths;

        
Npop= 30e6; % population (It affects the values of the parameters)
S= Npop-Q-R-D-E-I;
      
%%
% use small portion of point 
% % 
% nn=10
% position=position(1:nn,:);
% S=S(1:nn,:);
% I=I(1:nn,:);
% R=R(1:nn,:);
% D=D(1:nn,:);
%% 



E0=E(:,1);
I0=I(:,1);
Q0=Q(:,1);
R0=R(:,1);
D0=D(:,1);
P0=0*D0;
X0=[E0 I0 Q0 R0 D0];

%% buid the mesh 
Nx=size(X0,1);
Lx=1;
dx=Lx/(Nx-1);
Xx=0:dx:Lx;

N=size(time,2);
Tf=size(time,2);%360.0;
dt=Tf/(N-1);
t = [0:N-1].*dt;
Tt=0:dt:Tf;

%%  fitting of the real data   
DiffCoef_guess=[0];%,0,0,0,0,0,0];
% 
alpha_guess = 0.1; % protection rate
beta_guess = 1; % Infection rate
gamma_guess = 1; % 1/latent time in days
delta_guess = 0.1; % rate at which infectious people enter in quarantine
lambda_guess = 0.1; % recovery rate
kappa_guess = 0.1; % death rate

guess = [DiffCoef_guess, alpha_guess,beta_guess,gamma_guess, delta_guess,lambda_guess,kappa_guess];
    
[DiffCoef_fit, alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit] = fit_PDE(I,Q,R,D,Npop,E,t,guess);

fit_param=[DiffCoef_fit,alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit];


N=1*size(time,2);
Tf=size(time,2);%360.0;
dt=Tf/(N-1);
Tt2 = [0:N-1].*dt;

[S_fit,E_fit,I_fit,Q_fit,R_fit,D_fit,P_fit] = SEIQRDP_pde(DiffCoef_fit, alpha_fit,beta_fit,gamma_fit,delta_fit,lambda_fit,kappa_fit,Npop,X0,Tt2);

data=[fit_param];
II_fit=cumsum(I_fit,2);
EE_fit=cumsum(E_fit,2);
SS_fit=Npop-Q_fit-R_fit-D_fit-EE_fit-II_fit;


% save summaries
% T = array2table(data,'VariableNames',{'DiffCoef','alpha', 'beta','gamma','delta','lambda', 'kappa'});%,'ReMSE','RMSE','computationTime'});
% T.param=["fitted"]'

% Plot the results 



 [X,Y]=meshgrid(Tt,Xx);
 [X2,Y2]=meshgrid(Tt2,Xx);
 
 
    figure;
    subplot(211)
    mesh(X,Y,S) 
    colormap(jet)
    legend('S(t) Susceptible individials')
    
             
    subplot(212)
    mesh(X2,Y2,SS_fit) 
    colormap(jet)
    legend('S(t) Susceptible individials [fitted]')
    title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

  
    figure;
    subplot(211)
    mesh(X,Y,E) 
    colormap(jet)
    legend('E(t) Susceptible individials')

        subplot(212)
    mesh(X2,Y2,EE_fit) 
    colormap(jet)
    legend('E(t) Susceptible individials [fitted]')
    title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))
      
%              
%               
%     figure;
%     subplot(211)
%     mesh(X,Y,I) 
%     colormap(jet)
%     legend('I(t) infected individials')
%     
%     
%     subplot(212)
%     mesh(X2,Y2,I_fit) 
%     colormap(jet)
%     legend('I(t) infected individials [fitted]')
%     title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))
%             
%    


    figure;
    subplot(211)
    mesh(X,Y,I) 
    colormap(jet)
    legend('I(t) confirmed individials')

    
    subplot(212)
    mesh(X2,Y2,II_fit) 
    colormap(jet)
    legend('I(t) confirmed individials [fitted]')
    title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

    


              
    figure;
    subplot(211)
    mesh(X,Y,Q) 
    colormap(jet)
    legend('Q(t) Infected individials')

    
    subplot(212)
    mesh(X2,Y2,Q_fit) 
    colormap(jet)
    legend('Q(t) Infected individials [fitted]')
    title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))


    figure;
    subplot(211)
  
    mesh(X,Y,R) 
    colormap(jet)
    legend('R(t) Recovered individials')

    
    subplot(212)
    mesh(X2,Y2,R_fit) 
    colormap(jet)
    legend('R(t) Recovered individials [fitted]')
    title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))


                        
             
             
    figure;
    subplot(211)
    mesh(X,Y,D) 
    colormap(jet)
    legend('D(t) Death individials')

    
    
    subplot(212)
    mesh(X2,Y2,D_fit) 
    colormap(jet)
    legend('D(t) Death individials [fitted]')
       title(strcat('parameters:DiffCoef_fit=[',num2str(DiffCoef_fit),', beta=',num2str(beta_fit),', gamma=',num2str(gamma_fit),', delta=',num2str(delta_fit),', kappa=',num2str(kappa_fit)))

 