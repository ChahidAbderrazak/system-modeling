function [S,E,I,Q,R,D,P] = SEIQRDP_ode(alpha,beta,gamma,delta,lambda0,kappa0,Npop,E0,I0,Q0,R0,D0,t)

%% Initial conditions
N = numel(t);
dt = median(diff(t));

lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
kappa = kappa0(1)*exp(-kappa0(2).*t); % I use these functions for illustrative purpose only

S=Npop-Q0-E0-R0-D0-I0;
E(1)=E0;
I(1)=I0;
Q(1)=Q0;
R(1)=R0;
D(1)=D0;
P(1)=0;
        
for i=2:N

     S(i) = S(i-1)  + dt*( -alpha*S(i-1) - (beta./Npop).*I(i-1)*S(i-1)) ;
     E(i) = E(i-1) + dt*( -gamma*E(i-1) + (beta./Npop).*I(i-1)*S(i-1));
     I(i) = I(i-1)  + dt*( gamma*E(i-1) - delta.*I(i-1))  ;
     Q(i) = Q(i-1) + dt*(  delta.*I(i-1) -lambda(i-1)*Q(i-1) - kappa(i-1).*Q(i-1)) ;

     R(i)= R(i-1) + (dt)*(lambda(i-1))*Q(i-1);
     D(i)= D(i-1) + (dt)*(kappa(i-1))*Q(i-1);
     P(i)= P(i-1) + (dt)*(alpha)*S(i-1);


end

end


