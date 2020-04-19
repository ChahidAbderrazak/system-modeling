#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from lib.Covid19_dataset_preparation_Hopkins_CSSE import *
from lib.Display_and_plotting import *

# import os
# os.chdir('../')
  
def Covid19System_4s(state, t, *args):
    
    DiffCoef1=args[0]
    DiffCoef2=args[1]
    DiffCoef3=args[2]
    DiffCoef4=args[3]
    kappa=args[4]        
    beta=args[5]
    gamma=args[6]
    delta=args[7]
    zeta=args[8]
 
    Sn, In, Rn, Dn = state
    
    
    Nn=Sn + In + Rn + Dn;

    dSn= DiffCoef1*Sn - (1/Nn)*beta*(1-kappa)*Sn*In;
    dIn= DiffCoef2*In + (1/Nn)*beta*(1-kappa)*Sn*In -gamma*(1-delta)*In -delta*In ;
    dRn= DiffCoef3*Rn + gamma*(1-delta)*In;
    dDn= DiffCoef4*Dn + delta*In;


    return [dSn, dIn, dRn, dDn]

# ###############################################################
# Input parameter of the ODE

# ##################################################################################################
# input parameters
countries={'Spain', 'France','Qatar', 'Tunisia','Switzerland', 'Italy',\
           'Germany','United Kingdom','Turkey','US','Saudi Arabia','Greece', 'Korea, South',\
          'China','Iran','Morocco','Algeria','United Arab Emirates','Japan', 'Australia'}
         
population_dict={'US':328.2e6}  # 'Spain':46.94e6, 'France':66.99e6,'US':328.2e6,'Saudi Arabia':33.7e6

cut_time_start=0        # cut the data from the end    
cut_time_end=0          # cut the data from the end
Nup=1                   # data upsampling factor       
Ndown=1                 # data downsampling factor 
NL=20                    # #interation of optimal initial guess loop for solve the fitting 
params_name = ['DiffCoef1', 'DiffCoef2', 'DiffCoef3', 'DiffCoef4', 'kappa', 'beta', 'gamma', 'delta','zeta']
#-- params = (DiffCoef1,DiffCoef2,DiffCoef3,DiffCoef4,kappa,beta,gamma,delta,zeta)


#%% load data of a contry
dict_confirmed, dict_recovered, dict_deaths = load_dataset(countries)


#%% Fitting a  a specific contry 
#-- country='Spain'; Nt=46.94e6

for country, Nt in population_dict.items():
    print('\n\n--> fitting the COVID10 data of ',country, '\n     -> Total population =',Nt,'\n\n Please wait!....... ')
    
    
    #%% load data of a contry
    dict_confirmed, dict_recovered, dict_deaths = load_dataset(countries)
    
    # plot the data of the country
    # plot_data_of_country(country, dict_confirmed, dict_recovered, dict_deaths)

    
    # #% ODE parameters
    # N=85;
    # Tf=85.0;
    # dt=Tf/(N-1);
    # t = np.linspace(0, Tf, N)
    
    # load the data / conpute the relative states S
    Long, Lat, t_exp, t_exp0, I_nz, R_nz, D_nz = Get_data_of_country(country, dict_confirmed,dict_recovered, dict_deaths, cut_time_end, cut_time_start)    
    S_nz0= I_nz +R_nz+D_nz ; 
    S_nz=Nt-S_nz0
    
    t_exp = np.linspace(0 , max(t_exp), I_nz.shape[0])
    
    #-- plt.plot(t_exp,I_nz) ; plt.show()
    
    # build the fiiting dataset
    Noisy_data=np.array([S_nz, I_nz, R_nz, D_nz]).T
    
    #% very important dant  chnge
    #-----------------------------------------!!!
    Noisy_data=np.flipud(Noisy_data)
    t_exp0=np.flipud(t_exp0)
    #-----------------------------------------!!!


    #% add more sample to imporve fitting accuracy
    from scipy import signal
    # Noisy_data = Noisy_data[1500:2000]
    # t_exp = t_exp[1500:2000]
    # t_exp0 = t_exp0[1500:2000]
    

    Nup=Nup*t_exp.shape[0]
    Noisy_data = signal.resample(Noisy_data, Nup)
    t_exp = np.linspace(t_exp[-1],0 , Nup)
        
    if Ndown>1: # Downsampling
        Noisy_data = Noisy_data[::Ndown]
        t_exp = t_exp[::Ndown]
        t_exp0 = t_exp0[::Ndown]
           
        

    # store the states
    S_nz= Noisy_data[:,0]
    I_nz= Noisy_data[:,1]
    R_nz= Noisy_data[:,2]
    D_nz= Noisy_data[:,3]
    NT_nz=S_nz+I_nz+R_nz+D_nz
    
    # initial states
    S0=S_nz[0]
    I0=I_nz[0]
    R0=R_nz[0]
    D0=D_nz[0]
    init_state = [S0, I0, R0, D0]

    t=t_exp
    
    # % Define the fitting models 
    def residuals(p):
        p = tuple(p)
        sim_P = odeint(Covid19System_4s, init_state, t_exp, args = p)
        res = sim_P - Noisy_data
        return res.flatten()
    
    
    #  Fit the genrate data to the ODE
    from scipy.optimize import leastsq
    # initial_guess = [1, 1, 1, 1, 1, 1, 1, 1]
    initial_guess = [0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5]
    init_state = [S0, I0, R0, D0]
    
    
    # LS method fitting
    for cnt in range(0,NL):
        fitted_params = leastsq(residuals, initial_guess)[0]
        initial_guess=list(fitted_params)
    
    
    # simulate the fitted model
    fitted_params = tuple(fitted_params)
    sol_fitted =odeint(Covid19System_4s, init_state, t, args = fitted_params)
    
    # store the states
    S_hat= sol_fitted[:,0]
    I_hat= sol_fitted[:,1]
    R_hat= sol_fitted[:,2]
    D_hat= sol_fitted[:,3]
    
    NT_hat=S_hat+I_hat+R_hat+D_hat
    
    #%%  plot figure  with name
    from datetime import datetime
    import os
    now = datetime.now() # current date and time
    tag_time = now.strftime("%m-%d-%Y")
    my_path='results/ode4s_fitted_hpkins'+tag_time + '/' ; 
    if not os.path.exists(my_path):
        os.mkdir(my_path)
        
    filename_rslt=my_path+'NL-'+str(NL)+'_Ns-stice'+str(cut_time_start)+'_'+str(cut_time_end)+'_of'+country
    
    plot_fitting_results_country_4s(country, t_exp, S_nz, I_nz, R_nz, D_nz, t, S_hat, I_hat, R_hat, D_hat, t_exp0, filename_rslt)
    
    #%
    import numpy as np
    import pandas as pd
    Err=np.asarray(fitted_params)#np.linalg.norm()
    param_array=np.array ([np.asarray(fitted_params),Err])
    output = pd.DataFrame(param_array,index=['Fitted ','Error'],columns=params_name)
    output.head(-1)
    print(output)
    
    output.to_csv(filename_rslt +'.csv')
    


Noisy_data_hpks=Noisy_data
t_exp0_hpks=t_exp0
t_exp_hpks=t_exp
