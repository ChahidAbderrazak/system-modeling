function [S_,E_,I_,Q_,R_,D_,P_] = SEIQRDP_ode(alpha,beta,gamma,delta,lambda,kappa,Npop,X0,t)

%% Initial conditions
N = numel(t);
dt = median(diff(t));

%lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
%kappa = kappa0(1)*exp(-kappa0(2).*t); % I use these functions for illustrative purpose only


E_=X0(1);
I_=X0(2);
Q_=X0(3);
R_=X0(4);
D_=X0(5);
P_=0*D_;
S_= Npop-Q_-R_-D_-I_-E_;


for i=2:N

            S_(i) = S_(i-1)  + dt*( -alpha*S_(i-1) - (beta./Npop).*I_(i-1).*S_(i-1)) ;
            E_(i) = E_(i-1)  + dt*( -gamma*E_(i-1) + (beta./Npop).*I_(i-1).*S_(i-1));
            I_(i) = I_(i-1)  + dt*( gamma*E_(i-1) - delta.*I_(i-1))  ;
            Q_(i) = Q_(i-1)  + dt*(  delta.*I_(i-1) -lambda*Q_(i-1) - kappa.*Q_(i-1)) ;
            R_(i) = R_(i-1)  + dt*(lambda)*Q_(i-1);
            D_(i) = D_(i-1)  + dt*(kappa)*Q_(i-1);
            P_(i) = P_(i-1)  + dt*(alpha)*S_(i-1);

        

end

end


