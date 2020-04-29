function [S_,E_,I_,Q_,R_,D_,P_] = SEIQRDP_pde(DiffCoef, alpha,beta,gamma,delta,lambda,kappa,Npop,X0,t)
global Xx


%% Initial conditions
N = numel(t);
dt = median(diff(t));
Nx = numel(Xx);
dx = median(diff(Xx));

%%
Lambda=1/(dx*dx);
Ix=-2*Lambda*ones(Nx,1);
Ix1=Lambda*ones(Nx-1,1);
Ix2=Lambda*ones(Nx-1,1);
Matr=diag(Ix)+diag(Ix1,1)+diag(Ix2,-1);
%% Boundary Term
Matr(1,1)=Matr(1,1)+Lambda;
Matr(Nx,Nx)=Matr(Nx,Nx)+Lambda;


% lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
% kappa = kappa0(1)*exp(-kappa0(2).*t); % I use these functions for illustrative purpose only

E_=X0(:,1);
I_=X0(:,2);
Q_=X0(:,3);
R_=X0(:,4);
D_=X0(:,5);
P_=0*D_;
S_= Npop-Q_-R_-D_-I_-E_;

% DiffCoefS=DiffCoef(1);
% DiffCoefE=DiffCoef(2);
% DiffCoefI=DiffCoef(3);
% DiffCoefQ=DiffCoef(4);
% DiffCoefR=DiffCoef(5);
% DiffCoefD=DiffCoef(6);
% DiffCoefP=DiffCoef(7);

%%  Simulate the model 
for i=2:N
     S_(:,i) = S_(:,i-1) + dt*( -alpha*S_(:,i-1) - (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoef*Matr*S_(:,i-1);
     E_(:,i) = E_(:,i-1) + dt*( -gamma*E_(:,i-1) + (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoef*Matr*E_(:,i-1);
     I_(:,i) = I_(:,i-1) + dt*( gamma*E_(:,i-1) - delta.*I_(:,i-1))                     + dt*DiffCoef*Matr*I_(:,i-1);
     Q_(:,i) = Q_(:,i-1) + dt*(  delta.*I_(:,i-1) -lambda*Q_(:,i-1) - kappa.*Q_(:,i-1)) + dt*DiffCoef*Matr*Q_(:,i-1);
     R_(:,i) = R_(:,i-1) + dt*(lambda)*Q_(:,i-1)                                        + dt*DiffCoef*Matr*R_(:,i-1);
     D_(:,i) = D_(:,i-1) + dt*(kappa)*Q_(:,i-1)                                         + dt*DiffCoef*Matr*D_(:,i-1);
     P_(:,i) = P_(:,i-1) + dt*(alpha)*S_(:,i-1)                                         + dt*DiffCoef*Matr*P_(:,i-1);
            
            
                          
%      S_(:,i) = S_(:,i-1) + dt*( -alpha*S_(:,i-1) - (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoefS*Matr*S_(:,i-1);
%      E_(:,i) = E_(:,i-1) + dt*( -gamma*E_(:,i-1) + (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoefE*Matr*E_(:,i-1);
%      I_(:,i) = I_(:,i-1) + dt*( gamma*E_(:,i-1) - delta.*I_(:,i-1))                     + dt*DiffCoefI*Matr*I_(:,i-1);
%      Q_(:,i) = Q_(:,i-1) + dt*(  delta.*I_(:,i-1) -lambda*Q_(:,i-1) - kappa.*Q_(:,i-1)) + dt*DiffCoefQ*Matr*Q_(:,i-1);
%      R_(:,i) = R_(:,i-1) + dt*(lambda)*Q_(:,i-1)                                        + dt*DiffCoefR*Matr*R_(:,i-1);
%      D_(:,i) = D_(:,i-1) + dt*(kappa)*Q_(:,i-1)                                         + dt*DiffCoefD*Matr*D_(:,i-1);
%      P_(:,i) = P_(:,i-1) + dt*(alpha)*S_(:,i-1)                                         + dt*DiffCoefP*Matr*P_(:,i-1);
 end 
      
 
end