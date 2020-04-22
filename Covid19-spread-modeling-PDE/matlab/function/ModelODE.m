function [Sn1,In1,Rn1,Dn1]= ModelODE(dt,X0,Nn,DiffCoef1,DiffCoef2,DiffCoef3,DiffCoef4,kappa,beta,gamma,delta)
%% ================================================================ 

Sn=X0(1);
In=X0(2);
Rn=X0(3);
Dn=X0(4);

Sn1= Sn +dt*DiffCoef1*Sn - dt*(1/Nn)*beta*(1-kappa)*Sn*In;

%% ================================================================ 

In1=In +dt*DiffCoef2*In + dt*(1/Nn)*beta*(1-kappa)*Sn*In -dt*gamma*(1-delta)*In-dt*delta*In;


%% ================================================================ 

Rn1=Rn +dt*DiffCoef3*Rn + dt*gamma*(1-delta)*In;

%% ================================================================ 

Dn1=Dn +dt*DiffCoef4*Dn + dt*delta*In;

%Dn1=Dn +dt*DiffCoef4*Dn + dt*(1-gamma+delta)*In;


%Dn1=Dn +dt*DiffCoef4*Dn + dt*(1/Nn)*beta*Sn*In -dt*delta*Rn;



end