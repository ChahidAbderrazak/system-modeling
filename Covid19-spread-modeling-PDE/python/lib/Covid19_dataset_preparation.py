import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta


def load_data(filename):
       
    #-- filename='data/time_series_covid19_confirmed_global_narrow.csv'
    data= pd.read_csv(filename)
    data=data[{'Country/Region','Lat','Long','Date','Value'}]
    data=data.drop([0])
    data.dropna(inplace = True) 
    
    


    return data



def load_dataset(country):
    
    
    file_confirmed='data/time_series_covid19_confirmed_global_narrow.csv'
    file_recovered='data/time_series_covid19_recovered_global_narrow.csv'
    file_deaths='data/time_series_covid19_deaths_global_narrow.csv'
    
    data_confirmed_=load_data(file_confirmed); 
    data_confirmed=Build_data_for_countries(country, data_confirmed_)

    data_recovered_=load_data(file_recovered); 
    data_recovered=Build_data_for_countries(country, data_recovered_)
    
    data_deaths_=load_data(file_deaths);       
    data_deaths=Build_data_for_countries(country, data_deaths_)


    return data_confirmed, data_recovered, data_deaths


def Build_data_for_countries(countries, data_):
    
    #-- country='US'; data_=data
    #-- plt.plot(day,Value, label = country)
    
    Dict_data_= dict()

    
    for country in countries:
        data_n=data_.iloc[ data_.index[ data_['Country/Region'] == country ] ]
        data_n.drop(data_n.tail(1).index,inplace=True)
        country_name=np.unique(data_n['Country/Region'])
        
        if np.max(country_name.shape) > 1:
            print(' more than one contry is retreived ==>', country_name)
            
        else :
        
            country_name=country_name[0]
            
        #% update colectde data for each country
        data_n = data_n.astype({"Lat": float, "Long": float ,"Value": float})
            
        # convert the graduation date column to datetime objects
        # data_n['Date0'] = pd.to_datetime(data_n['Date'])
        Nrecord=data_n.shape[0]
        data_n['Date0'] = data_n['Date']

        data_n['Date'] = [Nrecord - i for i in range(0,Nrecord)]

        
        Country, Long, Lat, Value, day, day0 = get_sample(data_n) 
        
        
        Dict_data_.update({Country+'-p' : [Long, Lat]} ) # update the position 
        Dict_data_.update( {Country+'-t' : day} ) # update the time 
        Dict_data_.update( {Country+'-t0' : day0} ) # update the date 
        Dict_data_.update( {Country+'-v' : Value } ) # update the values 

    
    return Dict_data_




def get_sample(data_m):
    
    country = data_m['Country/Region'].iloc[0]
    Lat = data_m['Lat'].iloc[0]
    Long = data_m['Long'].iloc[0]
    Value = data_m['Value']
    day = data_m['Date']
    day0 = data_m['Date0']


    return country, Long, Lat, Value, day, day0


def Get_data_of_country(country, dict_confirmed, dict_recovered, dict_deaths, cut_time_end, cut_time_start):

    Long, Lat = dict_confirmed[country+'-p']
    t = dict_confirmed[country+'-t'].values
    t0 = dict_confirmed[country+'-t0'].values


    
    C_ = dict_confirmed[country+'-v'].values
    R_ = dict_recovered[country+'-v'].values
    D_ = dict_deaths[country+'-v'].values

    # idx0=np.min(np.where( C_<1 ) )+time_offset
    # cut_time_end, cut_time_start
    
    return Long, Lat, t[cut_time_start:-cut_time_end-1], t0[cut_time_start:-cut_time_end-1],\
        C_[cut_time_start:-cut_time_end-1], R_[cut_time_start:-cut_time_end-1], D_[cut_time_start:-cut_time_end-1]



