% %% ################################################################
% %%                      PDE model based on Covid-19 
% %% ################################################################
% close all; clear all; clc;
% addpath function
% global Xx
% %%
% Location = 'China';%'Netherlands';%'France';%
% Npop= 300e6; % population (It affects the values of the parameters)
% 
% [time, table_data]=retreive_data(Location);
% % position
% % R, D, I

% %%Initial conditions Sum
% E = sum(table_data.Confirmed); % Initial number of exposed cases. Unknown but unlikely to be zero.
% I = sum(table_data.Confirmed); % Initial number of infectious cases. Unknown but unlikely to be zero.
% Q = sum(table_data.Confirmed-table_data.Recovered-table_data.Deaths);
% R = sum(table_data.Recovered);
% D = sum(table_data.Deaths);       
% S= Npop-Q-R-D-E-I;
% %       
%% Initial One city 
cityn=9;
E = table_data.Confirmed(cityn,:); % Initial number of exposed cases. Unknown but unlikely to be zero.
I = table_data.Confirmed(cityn,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q = table_data.Confirmed(cityn,:)-table_data.Recovered(cityn,:)-table_data.Deaths(cityn,:);
R = table_data.Recovered(cityn,:);
D = table_data.Deaths(cityn,:);      
S= Npop-Q-R-D-E-I;
      
% use small portion of point 
% plot
figure(1)
subplot(2,2,1)
plot(S,'b', 'linewidth',3)
legend('S(t) Susceptible individials','S(t) fitted' )

subplot(2,2,2)
plot(I,'-g', 'linewidth',3)
legend('I(t) infected individials','I(t) fitted')

subplot(2,2,3)
plot(R,'-r', 'linewidth',3)
legend('R(t) Recovered individials','R(t) fitted')

%% 


E0=E(1);
I0=I(1);
Q0=Q(1);
R0=R(1);
D0=D(1);
P0=0*D0;
X0=[E0 I0 Q0 R0 D0];

%% buid the mesh 
N=size(time,2);
Tf=size(time,2);%360.0;
dt=Tf/(N-1);
t = [0:N-1].*dt;

%%  fitting of the simulated data   
guess = 0*[ 1, 1, 1, 1, 1 ,1]; % initial guess

[alpha_fit,beta_fit,gamma_fit, delta_fit, lambda_fit, kappa_fit] = fit_ODE(Q,R,D,Npop,E,I,t,guess);

fit_param=[alpha_fit,beta_fit,gamma_fit, delta_fit, lambda_fit, kappa_fit];
% 
[S_fit,E_fit,I_fit,Q_fit,R_fit,D_fit,P_fit] = SEIQRDP_ode(alpha_fit,beta_fit,gamma_fit, delta_fit, lambda_fit, kappa_fit,Npop,X0,t);

data=[fit_param];

% save summaries
% T = array2table(data,'VariableNames',{'DiffCoefS','DiffCoefI','DiffCoefR','DiffCoefD','beta','gamma','delta','kappa'});%,'ReMSE','RMSE','computationTime'});
%  
% T.param=["fitted"]'

% Plot the results 

time1=0:dt:Tf;
close all; figure;
plot(time1,Q)%,'r',time1,R,'c',time1,D,'g','linewidth',2);
hold on
plot(time1,Q_fit)%,'k-.',time1,R_fit,'k:',time1,D_fit,'k--','linewidth',2);
% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('Time (days)')
leg = {'Quarantined','Recovered','Dead',...
       'ODE fitted Quarantined','ODE fitted Recovered','ODE fitted Dead'}

   
legend(leg{:},'location','eastoutside')

% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('Time (days)')
leg = {'Tested positive minus the deceased cases','Deceased cases',...
       'ODE Fitted Tested positive minus the deceased cases','Fitted deceased'};
       
legend(leg{:},'location','southoutside')
set(gcf,'color','w')
axis tight