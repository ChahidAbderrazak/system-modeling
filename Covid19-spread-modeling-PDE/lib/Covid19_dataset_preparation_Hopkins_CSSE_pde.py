
     

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
    #data=data.drop(columns=['Province/State'])

    return data



def load_dataset(countries):
    
    file_confirmed, file_recovered,file_deaths=data_location()
    data_confirmed_=load_data(file_confirmed); 
    data_confirmed=Build_data_for_countries(countries, data_confirmed_)

    data_recovered_=load_data(file_recovered); 
    data_recovered=Build_data_for_countries(countries, data_recovered_)
    
    data_deaths_=load_data(file_deaths);       
    data_deaths=Build_data_for_countries(countries, data_deaths_)

    return data_confirmed, data_recovered, data_deaths


def Build_data_for_countries(countries, data_):
    
    #-- country='France'; data_=data_confirmed_

    Dict_data_= dict()
    cnt=0
    for country in countries:
        cnt=cnt+1
        data_n=data_.iloc[ data_.index[ data_['Country/Region'] == country ] ]
    
        # data_n=data_n.drop(columns=['Country/Region'])
        num_regions= data_n.shape[0] 
        data_n=data_n.groupby(['Country/Region']).sum()
        data_n['Lat']=data_n['Lat'].values/num_regions
        data_n['Long']=data_n['Long'].values/num_regions
        
        
        
        
        
        if cnt ==1:
            # data_pos=data_n[{'Province/State', 'Country/Region', 'Lat', 'Long'}]
            # data_n=data_n.drop(columns={'Province/State', 'Country/Region', 'Lat', 'Long'})
            
            
            data_pos=data_n[{'Lat', 'Long'}]
            data_n=data_n.drop(columns={'Lat', 'Long'})
            
            day0=data_n.columns            
            day=np.array([n for n in range(0,len(day0))])
            
            value=data_n.values
            
        else:
            
            # data_pos=data_n[{'Province/State', 'Country/Region', 'Lat', 'Long'}]
            # data_n=data_n.drop(columns={'Province/State', 'Country/Region', 'Lat', 'Long'})

            data_pos_=data_n[{'Lat', 'Long'}]
            data_n=data_n.drop(columns={'Lat', 'Long'})
            
            day0_=data_n.columns            
            day_=np.array([n for n in range(0,len(day0))])
            value_=data_n.values
                        
            # update the  dictionary 
            data_pos = pd.concat([data_pos, data_pos_ ], ignore_index=True)
            value = np.vstack([value, value_])

    
    data_pos=data_pos.reset_index()
    data_pos=data_pos.drop(columns=['index'])

    
    #-- plt.plot(day,Value) ; plt.show()
    Dict_data_.update({'position' : data_pos} ) # update the position 
    Dict_data_.update( {'t' : day} ) # update the time 
    Dict_data_.update( {'t0' : day0} ) # update the date 
    Dict_data_.update( {'value' : value } ) # update the values #-- Dict_data_[country+'-v'][0]

    return Dict_data_



def Get_data_of_country_pde(dict_confirmed, dict_recovered, dict_deaths):

    position= dict_confirmed['position']

    t = dict_confirmed['t']
    t0 = dict_confirmed['t0']


    
    C_ = dict_confirmed['value']  ; C_=C_.T
    R_ = dict_recovered['value']  ; R_=R_.T
    D_ = dict_deaths['value']     ; D_=D_.T

    # cut_time_start=np.min(np.min(np.where( C_<1 ) ))
    # cut_time_end, cut_time_start
    
    return position, t, t0, C_, R_, D_




