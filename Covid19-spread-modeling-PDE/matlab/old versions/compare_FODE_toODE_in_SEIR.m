addpath ../function
clear all;close all;clc;
global Rf Ef If Qf Df Pf Sf qs qe qi qq qr qd qp

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

[S,E,I,Q,R,D,P] = FSEIQRDP(alpha,beta,gamma,delta,Lambda01,Kappa01,Npop,E0,I0,Q0,R0,D0,t);

%% Compare fractionl
plot(time1,Q,'r',time1,R,'c',time1,D,'g','linewidth',2);
hold on
plot(time1,Qf,'k-.',time1,Rf,'k:',time1,Df,'k--','linewidth',2);
% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('Time (days)')
leg = {'Quarantined','Recovered','Dead',strcat(' fractional quarantined \alpha=' , num2str(qq))',...
    strcat(' fractional recovered \alpha=' , num2str(qr))',...
    strcat(' fractional Dead \alpha=' , num2str(qd))'};
legend(leg{:},'location','eastoutside')
set(gcf,'color','w')
axis tight



