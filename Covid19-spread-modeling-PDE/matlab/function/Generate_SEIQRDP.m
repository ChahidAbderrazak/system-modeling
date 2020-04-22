function [S,E,I,Q,R,D,P] = Generate_SEIQRDP(alpha,beta,gamma,delta,lambda0,kappa0,Npop,E0,I0,Q0,R0,D0,t)


%% Initial conditions
N = numel(t);
Y = zeros(7,N);
Y(1,1) = Npop-Q0-E0-R0-D0-I0;
Y(2,1) = E0;
Y(3,1) = I0;
Y(4,1) = Q0;
Y(5,1) = R0;
Y(6,1) = D0;

if round(sum(Y(:,1))-Npop)~=0
    error('the sum must be zero because the total population (including the deads) is assumed constant');
end
%%
modelFun = @(Y,A,F) A*Y + F;
dt = median(diff(t));
% ODE resolution

lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
kappa = kappa0(1)*exp(-kappa0(2).*t); % I use these functions for illustrative purpose only


for ii=1:N-1
    A = getA(alpha,gamma,delta,lambda(ii),kappa(ii));
    SI = Y(1,ii)*Y(3,ii); %  S*I
    F = zeros(7,1);
    F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
end


S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
Q = Y(4,1:N);
R = Y(5,1:N);
D = Y(6,1:N);
P = Y(7,1:N);



    function [A] = getA(alpha,gamma,delta,lambda,kappa)
        A = zeros(7);
        % S
        A(1,1) = -alpha;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-delta];
        % Q
        A(4,3:4) = [delta,-kappa-lambda];
        % R
        A(5,4) = lambda;
        % D
        A(6,4) = kappa;
        % P
        A(7,1) = alpha;
    end


    function [Y] = RK4(Fun,Y,A,F,h) % fun = A*Y +F
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
%         k_2 = Fun(Y+0.5*h*k_1,A,F);
%         k_3 = Fun(Y+0.5*h*k_2,A,F);
%         k_4 = Fun(Y+k_3*h,A,F);
        % output
        Y = Y + k_1*h;
    end






end


