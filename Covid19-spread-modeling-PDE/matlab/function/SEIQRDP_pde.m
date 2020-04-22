function [Sn,In,Rn,Dn, Nn] = SEIQRDP_pde(DiffCoef,beta,gamma,delta,kappa,X0,t)
global Xx

DiffCoefS=DiffCoef(1);
DiffCoefI=DiffCoef(2);
DiffCoefR=DiffCoef(3);
DiffCoefD=DiffCoef(4);

% DiffCoefS=0;
% DiffCoefI=DiffCoefS;
% DiffCoefR=DiffCoefS;
% DiffCoefD=DiffCoefS;

% beta=0.5;
% gamma=0.03;
% delta=0.0175;
% kappa=0.4;


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

Sn=X0(:,1);
In=X0(:,2);
Rn=X0(:,3);
Dn=X0(:,4);

Nn=Sn+In+Rn+Dn;

for i=2:N
    %%  ================================================================================== 
    Sn(:,i)=Sn(:,i-1) + dt*DiffCoefS*Matr*Sn(:,i-1)  - dt*(1./Nn(:,i-1))*beta*(1-kappa).*Sn(:,i-1).*In(:,i-1);

    %% ================================================================ 

    In(:,i)=In(:,i-1)  + dt*DiffCoefI*Matr*In(:,i-1) + dt*(1./Nn(:,i-1))*beta*(1-kappa).*Sn(:,i-1).*In(:,i-1) -dt*gamma*(1-delta)*In(:,i-1)-dt*delta*In(:,i-1) ;


    %% ================================================================ 

    Rn(:,i)=Rn(:,i-1)  +dt*DiffCoefR*Matr*Rn(:,i-1) + dt*gamma*(1-delta)*In(:,i-1);

    %% ================================================================ 

    Dn(:,i)=Dn(:,i-1) +dt*DiffCoefD*Matr*Dn(:,i-1) + dt*delta*In(:,i-1);
    
    
    Nn(:,i)=Sn(:,i)+In(:,i)+Rn(:,i)+Dn(:,i);

    
end


end