
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta



def Get_data_of_country(country, dict_confirmed, dict_recovered, dict_deaths, cut_time_end=0, cut_time_start=0):

    Long, Lat = dict_confirmed[country+'-p']
    t = dict_confirmed[country+'-t']
    t0 = dict_confirmed[country+'-t0']


    
    C_ = dict_confirmed[country+'-v']
    R_ = dict_recovered[country+'-v']
    D_ = dict_deaths[country+'-v']
    
    patient0=np.where( C_<1 )
    
    if len(patient0[0]) >0: 
        cut_time_start=np.min(patient0)

    # cut_time_start=np.min( )
    # cut_time_end, cut_time_start
    
    return Long, Lat, t[cut_time_start:-cut_time_end-1], t0[cut_time_start:-cut_time_end-1],\
        C_[cut_time_start:-cut_time_end-1], R_[cut_time_start:-cut_time_end-1], D_[cut_time_start:-cut_time_end-1]





def plot_data_of_country(country, dict_confirmed, dict_recovered, dict_deaths, Nt=6e7):
    
    import matplotlib.ticker as ticker
    Long, Lat, t_exp, t_exp0, Value_C, Value_R, Value_D = Get_data_of_country(country, dict_confirmed,dict_recovered, dict_deaths)
    

    my_xticks0=[str(s)  for s in t_exp0.tolist()]
    #
    Value_S0= Value_C +Value_R+Value_D
    Value_N= Nt*np.ones(Value_S0.shape[0])
    Value_S=Nt-Value_S0
    
    plt.subplot(511)
    plt.plot(t_exp, Value_S, label = country)
    ax = plt.gca()
    plt.title('The CVID19 data of '+ country)
    plt.xticks(rotation=45)
    plt.ylabel('S(t)')
    # plt.legend()
    plt.xticks([])
    
    
    plt.subplot(512)
    plt.plot(t_exp, Value_C, label = country)
    ax = plt.gca()
    plt.xticks(rotation=45)
    plt.ylabel('I(t)')
    # plt.legend()
    plt.xticks([])


    plt.subplot(513)
    plt.plot(t_exp,  Value_R, label = country)
    plt.xticks(rotation=45)
    plt.ylabel('R(t)')
    plt.xticks([])
    # plt.legend()
    
    plt.subplot(514)
    plt.plot(t_exp, Value_D,  label = country)
    plt.xticks(rotation=45)
    plt.ylabel('D(t)')
    # plt.legend()
    plt.xticks([])

    
    plt.subplot(515)
    plt.xticks(t_exp, my_xticks0)

    plt.plot(t_exp, Value_N,  label = country)
    plt.xticks(rotation=45)
    plt.ylabel('N(t)')
    # plt.legend()
    # ax = plt.gca()

    plt.locator_params(axis='x', nbins=8)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    
    
        
    plt.show()

    
        # plt.savefig('results/data-country/Covid19-'+country+'.png', format='png', dpi=1000)




def plot_data_of_locations(position, dict_confirmed, dict_recovered, dict_deaths, Nt=6e7):
    
    import matplotlib.ticker as ticker
    Long, Lat, t_exp, t_exp0, Value_C, Value_R, Value_D = Get_data_of_country(country, dict_confirmed,dict_recovered, dict_deaths)
    

    my_xticks0=[str(s)  for s in t_exp0.tolist()]
    #
    Value_S0= Value_C +Value_R+Value_D
    Value_N= Nt*np.ones(Value_S0.shape[0])
    Value_S=Nt-Value_S0
    
    plt.subplot(511)
    plt.plot(t_exp, Value_S, label = country)
    ax = plt.gca()
    plt.title('The CVID19 data of '+ country)
    plt.xticks(rotation=45)
    plt.ylabel('S(t)')
    # plt.legend()
    plt.xticks([])
    
    
    plt.subplot(512)
    plt.plot(t_exp, Value_C, label = country)
    ax = plt.gca()
    plt.xticks(rotation=45)
    plt.ylabel('I(t)')
    # plt.legend()
    plt.xticks([])


    plt.subplot(513)
    plt.plot(t_exp,  Value_R, label = country)
    plt.xticks(rotation=45)
    plt.ylabel('R(t)')
    plt.xticks([])
    # plt.legend()
    
    plt.subplot(514)
    plt.plot(t_exp, Value_D,  label = country)
    plt.xticks(rotation=45)
    plt.ylabel('D(t)')
    # plt.legend()
    plt.xticks([])

    
    plt.subplot(515)
    plt.xticks(t_exp, my_xticks0)

    plt.plot(t_exp, Value_N,  label = country)
    plt.xticks(rotation=45)
    plt.ylabel('N(t)')
    # plt.legend()
    # ax = plt.gca()

    plt.locator_params(axis='x', nbins=8)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    
    
        
    plt.show()
    
    
    

def plot_fitting_results_country_4s(country, t_exp, S_nz, I_nz, R_nz, D_nz, t, S_hat, I_hat, R_hat, D_hat, t_exp0, filename_rslt):
    
    import matplotlib.ticker as ticker
    t=t_exp
    Nticks=4
    stp=int(t_exp0.shape[0]/Nticks)
    
    s0=t_exp0.shape[0]-stp*Nticks-1
    
    my_xticks0=[str(t_exp0[s])  for s in range(s0,t_exp0.shape[0]+1,stp)]
    
    # plot figure
    plt.subplot(221)
    plt.scatter(t_exp, S_nz, c='r', s=8, label = 'S(t) Data')
    plt.plot(t, S_hat, '--g', label = 'S(t) fitted')
    plt.ylabel('S(t)')
    plt.legend()
    plt.xticks([])
    
    
    
    plt.subplot(222)
    plt.scatter(t_exp, I_nz, c='r', s=8, label = 'I(t) Data')
    plt.plot(t, I_hat, '--g', label = 'I(t) fitted')
    plt.ylabel('I(t)')
    plt.legend()
    plt.xticks([])
    
    plt.subplot(223)
    plt.scatter(t_exp, R_nz, c='r', s=8, label = 'R(t) Data')
    plt.xticks(t, my_xticks0)
    plt.plot(t, R_hat, '--g', label = 'R(t) fitted')
    plt.ylabel('R(t)')
    plt.legend()
    plt.locator_params(axis='x', nbins=Nticks)
    plt.xticks(rotation=45)
    
    plt.subplot(224)
    plt.scatter(t_exp, D_nz, c='r', s=8, label = 'D(t) Data')
    plt.xticks(t, my_xticks0)
    plt.plot(t, D_hat, '--g', label = 'D(t) fitted')
    plt.ylabel('D(t)')
    plt.legend()
    plt.locator_params(axis='x', nbins=Nticks)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(filename_rslt +'.png', format='png', dpi=1000)
    plt.show()

 
    
def draw_data_on_map(day, dict_, data_state, countries, tune_area=50, show_nmber=0):

    from mpl_toolkits.basemap import Basemap
    from matplotlib import animation, rc
    from IPython.display import HTML
    
    
    
    
    
    fig = plt.figure(figsize=(10, 10))
    cmap = plt.get_cmap('coolwarm')
    
    map = Basemap(projection='cyl')
    map.drawmapboundary()
    map.fillcontinents(color='lightgray', zorder=1)
    
    START_YEAR = 1950
    LAST_YEAR = 2013
    
    
    idx_p=[s + '-p' for s in countries]
    idx_t=[s + '-t' for s in countries]
    idx_v=[s + '-v' for s in countries]
    
    Long=[]
    Lat=[]
    value=[]
    k=-1
    
    
    # tag_xy=
    # ["{}{:0}".format(b_, a_) for a_, b_ in zip(country, b)]

    for country in countries:
        k=k+1
    
        Long_, Lat_ = dict_[idx_p[k]]
        value_ = dict_[idx_v[k]]
        value_=value_.iloc[day]
            

        Long.append(Long_)
        Lat.append(Lat_)
        value.append(value_)
        
        
        # add anotstion to the cases
        if show_nmber==1:
            
            plt.annotate(country+'-'+str(int(value_)), (Long_, Lat_))
            
        else:
                
            plt.annotate(country, (Long_, Lat_))
            

     
    plt.title('Number of '+ data_state +' '+str(day)+' days ago')
   
    year_text = plt.text(-170, 80, str(day),fontsize=15)
    

    M_value=max(value)

    Number=[s / tune_area for s in value]
    
    degree_colr=Number/M_value
    
    xs, ys = map( Long , Lat)
    
    print( ' the cases are: ', value)
    
    # year_text = plt.text(-170, 80, str(START_YEAR),fontsize=15)
    
    scat = map.scatter(xs, ys, s=Number, c=np.log10(value), cmap=cmap, marker='o', alpha=0.5, zorder=10)
    
    


    # make legend with dummy points
    for a in [1000, 5000,10000]:
        plt.scatter([], [], c='k', alpha=0.5, s=a/tune_area,
                    label=str(a) + ' cases')
    plt.legend(scatterpoints=1, frameon=False,
               labelspacing=1, loc='lower left');
    
