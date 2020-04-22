clear all;close all;clc;warning ('off','all')
addpath function2

% Download the data from ref [1] and read them with the function getDataCOVID
[tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID();
% time = time(1:end-1);
fprintf(['Most recent update: ',datestr(time(end)),'\n'])


Location = 'China';


try
    indR = find(contains(tableRecovered.CountryRegion,Location)==1);
    indC = find(contains(tableConfirmed.CountryRegion,Location)==1);
    indD = find(contains(tableDeaths.CountryRegion,Location)==1);
catch exception
    searchLoc = strfind(tableRecovered.CountryRegion,Location);
    indR = find(~cellfun(@isempty,searchLoc)) ; 
    
    searchLoc = strfind(tableConfirmed.CountryRegion,Location);
    indC = find(~cellfun(@isempty,searchLoc)) ; 
    
    searchLoc = strfind(tableDeaths.CountryRegion,Location);
    indD = find(~cellfun(@isempty,searchLoc));    
end

% disp(tableRecovered(indR,1:2))
disp(tableConfirmed(indC,1:2))
% disp(tableDeaths(indD,1:2))

% If the number of confirmed Confirmed cases is small, it is difficult to know whether
% the quarantine has been rigorously applied or not. In addition, this
% suggests that the number of infectious is much larger than the number of
% confirmed cases


cities=tableConfirmed(indC,1:2);

 %% Iterative application fiiting for each province    
        
timeRef = time; % Used in the loop only

disp('Fitting in process. Please wait...!!')        
        
for ii = 1:min([numel(indR),numel(indC),numel(indD)])
    Recovered = table2array(tableRecovered(indR(ii),5:end));
    Deaths = table2array(tableDeaths(indD(ii),5:end));
    Confirmed = table2array(tableConfirmed(indC(ii),5:end));
    minNum= max(100,round(0.1*max(Confirmed)));
    % Warning: a dummy value of Npop is used here. 
     Npop= 30e6; % population (It affects the values of the parameters)
     
     
    % Remove case where only few infectious are recorded (to avoid bad
    % initial conditions)
    Recovered(Confirmed<=minNum)=[];
    Deaths(Confirmed<=minNum)=[];
    time = timeRef; % trick to avoid reloading the variable "time" at each new loop
    time(Confirmed<=minNum)= [];
    Confirmed(Confirmed<=minNum)=[];
    
    % The fitting is only applied if enough data is collected (that is why
    % I use the case of China)
    if numel(Confirmed)>30 % If more than 30 days of data, run the fit

        fprintf('%s...',cities.ProvinceState(ii))
        % Definition of the first estimates for the parameters
        alpha_guess = 0.6; % protection rate
        beta_guess = 0.8; % Infection rate
        LT_guess = 0.5; % latent time in days
        Q_guess = 0.5; % rate at which infectious people enter in quarantine
        lambda_guess = [0.1,0.05]; % recovery rate
        kappa_guess = [0.1,0.05]; % death rate
        
        guess = [alpha_guess,beta_guess,1/LT_guess, Q_guess,lambda_guess,kappa_guess];
        
        % Initial conditions
        E0 = Confirmed(1); % Initial number of exposed cases. Unknown but unlikely to be zero.
        I0 = Confirmed(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
        Q0 = Confirmed(1)-Recovered(1)-Deaths(1);
        R0 = Recovered(1);
        D0 = Deaths(1);
        
        Active = Confirmed-Recovered-Deaths;
        Active(Active<0) = 0; % No negative number possible
        
        tic
        [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] = fit_ODE(Active,Recovered,Deaths,Npop,E0,I0,time,guess,'Display','off');
        time_fit_ode=toc;

        dt = 0.1; % time step
        time1 = datetime(time(1)):dt:datetime(datestr(floor(datenum(now))+datenum(10)));
        N = numel(time1);
        t = [0:N-1].*dt;
        
        
        % Call of the function SEIQRDP.m with the fitted parameters
        [S,E,I,Q,R,D,P] = SEIQRDP_ode(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t);
        

        
        %% plot figure
%         clf;close all;
        figure
        semilogy(time,Active,'ro',time,Recovered,'bo',time,Deaths,'ko');
        hold on
        semilogy(time1,Q,'r--',time1,R,'b--',time1,D,'k--');
        hold off
        % ylim([0,1.1*Npop])
        ylabel('Number of cases')
        xlabel('time (days)')
        leg = {'Confirmed (reported)','Recovered (reported)','Deceased  (reported)',...
                'ODE Confirmed (fitted)','ODE Recovered (fitted)','ODE Deceased (fitted)'}
               
        lh =legend(leg{:});%,'location','southoutside');
                
%         lh.Position(1) = 0.5 - lh.Position(3)/2; 
        lh.Position(2) = 0.5 - lh.Position(4)/2;

        set(gcf,'color','w')

        %%% title %%%
        subLoc = char(table2array(tableRecovered(indR(ii),1)));
        Loc = char(table2array(tableRecovered(indR(ii),2)));
        title(['Location: ',subLoc,' (',Loc,')']);
        %%%%%%%%%%%%%
        
        %%%%%%%%%%%%%

        grid on
        axis tight
        set(gca,'yscale','lin')
        pause(1)

        
        %% save performance as table
        % metrics
        Q_=downsample_data(Q, dt, 11);      rmse_q=error_vectors(Active,Q_);    

    
        
        % store data
        data_ode(:,ii)=[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,rmse_q,  time_fit_ode];
        path_reslt=strcat('result/',cities.CountryRegion(ii),'/');
        mkdir(path_reslt)
        saveas(gcf,strcat(path_reslt,'/ODE-fitting-',cities.ProvinceState(ii),'.fig'))

    end

end

%% save summaries
T1 = array2table(data_ode','VariableNames',{'alpha','beta','gamma','delta','lambda0','lambda1','kappa0','kappa1','ReMSE','RMSE','computationTime'});
T1.city=cities.ProvinceState
writetable(T1, strcat(path_reslt,'/ODE-fitting-',Location,'_summary.csv'))
    
        