function [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1, varargout] = fit_ODE(Q,R,D,Npop,E0,I0,time,guess,varargin)

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  %  option for optimset
p.addOptional('tolFun',1e-4);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;
dt  = p.Results.dt ;

%% Options for lsqcurvfit

options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',800,'Display',Display);
%% Fitting the data

% Write the target input into a matrix
Q(Q<0)=0; % negative values are not possible
R(R<0)=0; % negative values are not possible
D(D<0)=0; % negative values are not possible

if isempty(R)
    warning(' No data available for "Recovered" ')
    input = [Q;D];
else
    input = [Q;R;D];
end

if size(time,1)>size(time,2) && size(time,2)==1,    time = time';end
if size(time,1)>1 && size(time,2)>1,  error('Time should be a vector');end

fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal 

t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges



modelFun1 = @SEIQRDP_for_fitting; % transform a nested function into anonymous function

% call Lsqcurvefit
[Coeff,~,residual,~,~,~,jacobian] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess,tTarget(:)',input,zeros(1,numel(guess)),[1 3 1 1 2 3 2 2],options);


if nargout ==7
    varargout{1} = residual;
elseif nargout==8
    varargout{1} = residual;
    varargout{2} = jacobian;
elseif nargout==9
    varargout{1} = residual;
    varargout{2} = jacobian;
    varargout{3} = modelFun1;
elseif nargout>9
    error('Too many output specified')
end




%% Write the fitted coeff in the outputs
alpha1 = abs(Coeff(1));
beta1 = abs(Coeff(2));
gamma1 = abs(Coeff(3));
delta1 = abs(Coeff(4));
Lambda1 = abs(Coeff(5:6));
Kappa1 = abs(Coeff(7:8));


%% nested functions


      
 function [output] = SEIQRDP_for_fitting(para,t0)

        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5:6));    lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
        kappa0 = abs(para(7:8));     kappa = kappa0(1)*exp(-kappa0(2).*t);
        %% Initial conditions
        N = numel(t);
        
        E_(1)=E0;
        I_(1)=I0;
        Q_(1)=Q(1);
        
        if ~isempty(R)
            R_(1)=R(1);
            S_(1)= Npop-Q(1)-R(1)-D(1)-E0-I0;
        else
            S_(1)= Npop-Q(1)-D(1)-E0-I0;
              
        end
        D_(1)=D(1);
        P_(1)=0;
        

        dt = median(diff(t));


        
        for i=2:N
           
             S_(i) = S_(i-1)  + dt*( -alpha*S_(i-1) - (beta./Npop).*I_(i-1)*S_(i-1)) ;
             E_(i) = E_(i-1) + dt*( -gamma*E_(i-1) + (beta./Npop).*I_(i-1)*S_(i-1));
             I_(i) = I_(i-1)  + dt*( gamma*E_(i-1) - delta.*I_(i-1))  ;
             Q_(i) = Q_(i-1) + dt*(  delta.*I_(i-1) -lambda(i-1)*Q_(i-1) - kappa(i-1).*Q_(i-1)) ;

             R_(i)= R_(i-1) + (dt)*(lambda(i-1))*Q_(i-1);
             D_(i)= D_(i-1) + (dt)*(kappa(i-1))*Q_(i-1);
             P_(i)= P_(i-1) + (dt)*(alpha)*S_(i-1);


        end
        
                
                
%         I1 = Y(3,1:N);
        Q1 = Q_;
        R1 = R_;
        D1 = D_;
        
        Q1 = interp1(t,Q1,t0);
        R1 = interp1(t,R1,t0);
        D1 = interp1(t,D1,t0);
        
        
        if ~isempty(R)
            output = [Q1;R1;D1];
        else
            output = [Q1+R1;D1];
        end

 end

end



