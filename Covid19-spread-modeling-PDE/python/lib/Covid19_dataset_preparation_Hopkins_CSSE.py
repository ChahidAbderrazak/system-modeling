
     

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta


def data_location():
    
    import glob
    
    root='data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/'
    place='global' #'US'
    
    # get the path of the differnt data
    file_confirmed=glob.glob(root+'*confirmed*_'+place+'.csv')[0]
    file_recovered=glob.glob(root+'*recovered*_'+place+'.csv')[0]
    file_deaths=glob.glob(root+'*deaths*_'+place+'.csv')[0]

    
    return file_confirmed, file_recovered,file_deaths



def load_data(filename):
       
    #-- filename=file_confirmed
    data= pd.read_csv(filename)
    
    
    # data=data[{'Country/Region','Lat','Long','Date','Value'}]
    # data=data.drop([0])
    # data.dropna(inplace = True) 
    
    # cases ODE
    data=data.drop(columns=['Province/State'])


    # case PDE
    #data.dropna(inplace = True) 

    return data



def load_dataset(country):
    
    file_confirmed, file_recovered,file_deaths=data_location()
    data_confirmed_=load_data(file_confirmed); 
    
    data_confirmed=Build_data_for_countries(country, data_confirmed_)

    data_recovered_=load_data(file_recovered); 
    data_recovered=Build_data_for_countries(country, data_recovered_)
    
    data_deaths_=load_data(file_deaths);       
    data_deaths=Build_data_for_countries(country, data_deaths_)


    return data_confirmed, data_recovered, data_deaths


def Build_data_for_countries(countries, data_):
    
    #-- country='France'; data_=data_confirmed_
    
    
    
    Dict_data_= dict()

    
    for country in countries:
        data_n=data_.iloc[ data_.index[ data_['Country/Region'] == country ] ]
        num_regions= data_n.shape[0]       
        
        # if num_regions>1:
        data_n=data_n.groupby(['Country/Region']).sum()
        
        data_n['Lat']=data_n['Lat'].values/num_regions
        data_n['Long']=data_n['Long'].values/num_regions
        
        
    
            
        #% update colectde data for each country
        # data_n = data_n.astype({"Lat": float, "Long": float ,"Value": float})
            
        # convert the graduation date column to datetime objects
        # data_n['Date0'] = pd.to_datetime(data_n['Date'])

    

        Lat=data_n['Lat'].values[0];      data_n=data_n.drop(columns=['Lat'])
        Long=data_n['Long'].values[0];    data_n=data_n.drop(columns=['Long'])
        
        country_name=country ;            
        # data_n=data_n.drop(columns=['Country/Region'])
        day0=data_n.columns
        Value=data_n.values[0]
        
        day=np.array([n for n in range(0,len(day0))])
        
        # day=np.array([n for n in range(len(day0)-1,-1,-1)])

        #-- plt.plot(day,Value) ; plt.show()
        Dict_data_.update({country+'-p' : [Long, Lat]} ) # update the position 
        Dict_data_.update( {country+'-t' : day} ) # update the time 
        Dict_data_.update( {country+'-t0' : day0} ) # update the date 
        Dict_data_.update( {country+'-v' : Value } ) # update the values 

    
    return Dict_data_



def Get_data_of_country(country, dict_confirmed, dict_recovered, dict_deaths, cut_time_end=0, cut_time_start=0):

    Long, Lat = dict_confirmed[country+'-p']
    t = dict_confirmed[country+'-t']
    t0 = dict_confirmed[country+'-t0']


    
    C_ = dict_confirmed[country+'-v']
    R_ = dict_recovered[country+'-v']
    D_ = dict_deaths[country+'-v']

    cut_time_start=np.min(np.where( C_<1 ) )
    # cut_time_end, cut_time_start
    
    return Long, Lat, t[cut_time_start:-cut_time_end-1], t0[cut_time_start:-cut_time_end-1],\
        C_[cut_time_start:-cut_time_end-1], R_[cut_time_start:-cut_time_end-1], D_[cut_time_start:-cut_time_end-1]




