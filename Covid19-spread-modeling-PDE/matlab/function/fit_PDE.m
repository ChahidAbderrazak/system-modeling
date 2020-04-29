function [DiffCoef1, alpha1,beta1,gamma1,delta1,lambda1,kappa1, varargout] = fit_PDE(I,Q,R,D,Npop,E,time,guess,varargin)

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
% I(I<0)=0; % negative values are not possible
Q(Q<0)=0; % negative values are not possible
R(R<0)=0; % negative values are not possible
D(D<0)=0; % negative values are not possible

if isempty(R)
    warning(' No data available for "Recovered" ')
%     input = [I;Q;D];
    input = [Q;D];
else
%     input = [I;Q;R;D];
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
DiffCoef1=abs(Coeff(1));
alpha1 = abs(Coeff(2));
beta1 = abs(Coeff(3));
gamma1 = abs(Coeff(4));
delta1 = abs(Coeff(5));
lambda1 = abs(Coeff(6));
kappa1 = abs(Coeff(7));

% DiffCoef1=[ abs(Coeff(1)) abs(Coeff(8:13))] ;



%% nested functions


      
 function [output] = SEIQRDP_for_fitting(para,t0)
     
        
        global Xx
        DiffCoef=abs(para(1));
        alpha = abs(para(2));
        beta = abs(para(3));
        gamma = abs(para(4));
        delta = abs(para(5));
        lambda = abs(para(6));
        kappa = abs(para(7));
                
%         DiffCoefS = abs(para(1));
%         DiffCoefE = abs(para(8));
%         DiffCoefI = abs(para(9));
%         DiffCoefQ = abs(para(10));
%         DiffCoefR = abs(para(11));
%         DiffCoefD = abs(para(12));
%         DiffCoefP = abs(para(13));
        
        %% Initial conditions
        N = numel(t);
        dt = median(diff(t));
        Nx = numel(Xx);
        dx = median(diff(Xx));

        
        
        E_(:,1)=E(:,1);
        I_(:,1)=I(:,1);
        Q_(:,1)=Q(:,1);
        
        if ~isempty(R)
            R_(:,1)=R(:,1);
            S_(:,1)= Npop-Q(:,1)-R(:,1)-D(:,1)-E(:,1)-I(:,1);
        else
            S_(:,1)= Npop-Q(:,1)-D(:,1)-E(:,1)-I(:,1);
              
        end
        D_(:,1)=D(:,1);
        P_(:,1)=0*D_(:,1);


        %%
        Lambda=1/(dx*dx);
        Ix=-2*Lambda*ones(Nx,1);
        Ix1=Lambda*ones(Nx-1,1);
        Ix2=Lambda*ones(Nx-1,1);
        Matr=diag(Ix)+diag(Ix1,1)+diag(Ix2,-1);
        %% Boundary Term
        Matr(1,1)=Matr(1,1)+Lambda;
        Matr(Nx,Nx)=Matr(Nx,Nx)+Lambda;

%%  Simulate the model 
for i=2:N
     S_(:,i) = S_(:,i-1) + dt*( -alpha*S_(:,i-1) - (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoef*Matr*S_(:,i-1);
     E_(:,i) = E_(:,i-1) + dt*( -gamma*E_(:,i-1) + (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoef*Matr*E_(:,i-1);
     I_(:,i) = I_(:,i-1) + dt*( gamma*E_(:,i-1) - delta.*I_(:,i-1))                     + dt*DiffCoef*Matr*I_(:,i-1);
     Q_(:,i) = Q_(:,i-1) + dt*(  delta.*I_(:,i-1) -lambda*Q_(:,i-1) - kappa.*Q_(:,i-1)) + dt*DiffCoef*Matr*Q_(:,i-1);
     R_(:,i) = R_(:,i-1) + dt*(lambda)*Q_(:,i-1)                                        + dt*DiffCoef*Matr*R_(:,i-1);
     D_(:,i) = D_(:,i-1) + dt*(kappa)*Q_(:,i-1)                                         + dt*DiffCoef*Matr*D_(:,i-1);
     P_(:,i) = P_(:,i-1) + dt*(alpha)*S_(:,i-1)                                         + dt*DiffCoef*Matr*P_(:,i-1);
            
            
            
%             
%             
%             
%      S_(:,i) = S_(:,i-1) + dt*( -alpha*S_(:,i-1) - (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoefS*Matr*S_(:,i-1);
%      E_(:,i) = E_(:,i-1) + dt*( -gamma*E_(:,i-1) + (beta./Npop).*I_(:,i-1).*S_(:,i-1))  + dt*DiffCoefE*Matr*E_(:,i-1);
%      I_(:,i) = I_(:,i-1) + dt*( gamma*E_(:,i-1) - delta.*I_(:,i-1))                     + dt*DiffCoefI*Matr*I_(:,i-1);
%      Q_(:,i) = Q_(:,i-1) + dt*(  delta.*I_(:,i-1) -lambda*Q_(:,i-1) - kappa.*Q_(:,i-1)) + dt*DiffCoefQ*Matr*Q_(:,i-1);
%      R_(:,i) = R_(:,i-1) + dt*(lambda)*Q_(:,i-1)                                        + dt*DiffCoefR*Matr*R_(:,i-1);
%      D_(:,i) = D_(:,i-1) + dt*(kappa)*Q_(:,i-1)                                         + dt*DiffCoefD*Matr*D_(:,i-1);
%      P_(:,i) = P_(:,i-1) + dt*(alpha)*S_(:,i-1)                                         + dt*DiffCoefP*Matr*P_(:,i-1);
 end 
                                   

        
        [Xq,Yq] = meshgrid(time,Xx);
        [X,Y] = meshgrid(t,Xx);
        
%         I1 = interp2(X,Y,I_,Xq,Yq);
        Q1 = interp2(X,Y,Q_,Xq,Yq);
        R1 = interp2(X,Y,R_,Xq,Yq);
        D1 = interp2(X,Y,D_,Xq,Yq);
        
        if ~isempty(R)
%             output = [I1;Q1;R1;D1];
            output = [Q1;R1;D1];

        else
%             output = [I1;Q1+R1;D1];
            output = [Q1+R1;D1];

        end

 end

end

