addpath function
% clear all;close all;clc;

% Time definition
dt = 0.1; % time step
time1 = datetime(2010,01,01,0,0,0):dt:datetime(2010,09,01,0,0,0);
N = numel(time1);
t = [0:N-1].*dt;

%%
Npop= 60e6; % population (60 millions)
Q0 = 200; % Initial number of infectious that have bee quanrantined
I0 = Q0; % Initial number of infectious cases non-quarantined
E0 = 0; % Initial number of exposed
R0 = 0; % Initial number of recovereds
D0 = 1; % Initial number of deads
alpha = 0.08; %  protection rate
beta = 0.9; %  infection rate
gamma= 1/2; % inverse of average latent time
delta= 1/8; % rate at which infectious people enter in quarantine
Lambda01 = [0.03 0.05]; % cure rate (time dependant)
Kappa01 =  [0.03 0.05]; % mortality rate (time dependant)


%% Generate the data from the model 
[S,E,I,Q,R,D,P] = Generate_SEIQRDP(alpha,beta,gamma,delta,Lambda01,Kappa01,Npop,E0,I0,Q0,R0,D0,t);

%% fitting of the simulated data 
guess = [0.05,0.9,1/4,1/10,0.03,0.03,0.02,0.06]; % initial guess

[alpha_fit, beta_fit, gamma_fit, delta_fit, Lambda_fit, Kappa_fit] = fit_ODE(Q,R,D,Npop,E0,I0,time1,guess);

fit_param=[alpha_fit, beta_fit, gamma_fit, delta_fit, Lambda_fit, Kappa_fit];

[S_fit,E_fit,I_fit,Q_fit,R_fit,D_fit,P_fit] = SEIQRDP_ode(alpha_fit, beta_fit, gamma_fit, delta_fit, Lambda_fit, Kappa_fit, Npop,E0,I0,Q0,R0,D0,t);



%% Compare all states

%%
Tt=time;
figure(2)
subplot(2,2,1)
plot(Tt,S,'-g',Tt,S_fit,'-.k', 'linewidth',3)
legend('S(t) Susceptible individials','S(t) fitted' )
title(Location)

subplot(2,2,2)
plot(Tt,I,'-g',Tt,I_fit,'-.k', 'linewidth',3)
legend('I(t) infected individials','I(t) fitted')
title(Location)

subplot(2,2,3)
plot(Tt,R,'-r',Tt,R_fit,'-.k', 'linewidth',3)
legend('R(t) Recovered individials','R(t) fitted')
title(Location)






 %% Comapre the obtained results
% figure
% clf;close all;
% plot(time1,Q+R,'r',time1,D,'g','linewidth',2);
% hold on
% plot(time1,Q_fit+R_fit,'r-.',time1,D_fit,'g--','linewidth',2);
% 
% % ylim([0,1.1*Npop])
% ylabel('Number of cases')
% xlabel('Time (days)')
% leg = {'Tested positive minus the deceased cases','Deceased cases',...
%        'ODE Fitted Tested positive minus the deceased cases','Fitted deceased'};
%        
% legend(leg{:},'location','southoutside')
% set(gcf,'color','w')
% axis tight


