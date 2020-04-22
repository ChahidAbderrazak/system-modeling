function [DiffCoef1,beta1,gamma1,delta1,kappa1, varargout] = fit_PDE(I,R,D,S0,time,guess,varargin)

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  %  option for optimset
p.addOptional('tolFun',1e-4);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting
% p.addOptional('dt',0.0360); % time step for the fitting


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
I(I<0)=0; % negative values are not possible
R(R<0)=0; % negative values are not possible
D(D<0)=0; % negative values are not possible

if isempty(R)
    warning(' No data available for "Recovered" ')
    input = [I;D];
else
    input = [I;R;D];
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
DiffCoef1=abs(Coeff(1:4));
beta1 = abs(Coeff(5));
gamma1 = abs(Coeff(6));
delta1 = abs(Coeff(7));
kappa1 = abs(Coeff(8));

%% nested functions


      
 function [output] = SEIQRDP_for_fitting(para,t0)
     
         global Xx
        DiffCoefS=abs(para(1));
        DiffCoefI=abs(para(1));
        DiffCoefR=abs(para(3));
        DiffCoefD=abs(para(4));
        beta  = abs(para(5));
        gamma = abs(para(6));
        delta = abs(para(7));
        kappa = abs(para(8));


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

        In(:,1)=I(:,1);
        Rn(:,1)=R(:,1);
        Dn(:,1)=D(:,1);
        Sn(:,1)=S0(:,1);


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
        
               
        
        [Xq,Yq] = meshgrid(time,Xx);
        [X,Y] = meshgrid(t,Xx);

        I1 = interp2(X,Y,In,Xq,Yq);
        R1 = interp2(X,Y,Rn,Xq,Yq);
        D1 = interp2(X,Y,Dn,Xq,Yq);
        
        if ~isempty(R)
            output = [I1;R1;D1];
        else
            output = [I1+R1;D1];
        end

 end

end

