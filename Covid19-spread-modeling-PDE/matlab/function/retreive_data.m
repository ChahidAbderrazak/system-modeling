function   [time, Sort_Table]=retreive_data()

warning ('off','all')
disp('Data colelcion online. Please wait...!!')        

% Download the data from ref [1] and read them with the function getDataCOVID
[tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID();
% time = time(1:end-1);
fprintf(['Most recent update: ',datestr(time(end)),'\n'])
Location = 'China';%'Netherlands';%'France';%

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

Recovered_ = table2array(tableRecovered(indR,3:end));
Deaths_ = table2array(tableDeaths(indD,3:end));
Confirmed_ = table2array(tableConfirmed(indC,3:end));
% minNum= max(100,round(0.1*max(Confirmed_)));

% Warning: a dummy value of Npop is used here. 
% Npop= 30e6; % population (It affects the values of the parameters)


%% get position/records

%% save summaries
position_=Recovered_(:,1:2);
Recovered = Recovered_(:,3:end);
Deaths = Deaths_(:,3:end);
Confirmed = Confirmed_(:,3:end);

T = array2table(position_,'VariableNames',{'Lat','Long'})
T.r=T.Lat.^2 + T.Long.^2; 
T.ProvinceState=tableRecovered(indR,1);
T.CountryRegion=tableRecovered(indR,2);

T.Deaths=Deaths;
T.Confirmed=Confirmed;
T.Recovered=Recovered;
T.s=sum(T.Deaths+T.Recovered+T.Confirmed,2)
Sort_Table = sortrows(T,'s');
end
    