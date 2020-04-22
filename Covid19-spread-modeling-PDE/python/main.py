#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from lib.Covid19_dataset_preparation_Hopkins_CSSE_pde import *
from lib.Display_and_plotting import *

# import os
# os.chdir('../')
  
def Covid19System_pde(state, t, *args):
    
    DiffCoefS=args[0]
    DiffCoefI=args[1]
    DiffCoefR=args[2]
    DiffCoefD=args[3]
    kappa=args[4]        
    beta=args[5]
    gamma=args[6]
    delta=args[7]
    # PDE defusion laplacian
    Nx=4
    Lx=1
    dx=Lx/(Nx-1)
    Lambda=1/(dx*dx)
    
    
    Ix=[-2*Lambda  for i in range(1,Nx+1)]
    Ix1=Lambda*np.eye(Nx, k=1)
    Ix2=Lambda*np.eye(Nx, k=-1)
    
    Matr=np.diag(np.array([Ix])[0]) + Ix1  + Ix2
    # Boundary Term
    Matr[0,0]=Matr[0,0]+Lambda;
    Matr[Nx-1,Nx-1]=Matr[Nx-1,Nx-1]+Lambda;


    #-- state=init_state
    # t the states
    Sn= state[0:Nx]
    In= state[1*Nx:2*Nx]
    Rn= state[2*Nx:3*Nx]
    Dn= state[3*Nx:4*Nx]

    Nn=Sn + In + Rn + Dn

    dSn= DiffCoefS*Matr.dot(Sn) - (1/Nn)*beta*(1-kappa)*Sn*In
    dIn= DiffCoefI*Matr.dot(In) + (1/Nn)*beta*(1-kappa)*Sn*In -gamma*(1-delta)*In -delta*In ;
    dRn= DiffCoefR*Matr.dot(Rn) + gamma*(1-delta)*In;
    dDn= DiffCoefD*Matr.dot(Dn) + delta*In;

    return np.concatenate((dSn, dIn, dRn, dDn)) 

# ###############################################################
# Input parameter of the ODE
Nx=4
Lx=1
    
# ##################################################################################################
# input parameters
countries={'Tunisia', 'Algeria', 'Morocco','Spain','Portugal', 'France','Italy',\
           'Switzerland','Germany','Austria','Belgium','Netherlands'}

 
countries= {'Spain', 'Austria', 'Germany','Belgium'}
        
population_dict={'France':66.99e6,}  # 'Spain':46.94e6, 'France':66.99e6,'US':328.2e6,'Saudi Arabia':33.7e6
cut_time_start=0        # cut the data from the end    
cut_time_end=0          # cut the data from the end
Nup=1                   # data upsampling factor       
Ndown=1                 # data downsampling factor 
NL=100                    # #interation of optimal initial guess loop for solve the fitting 
params_name = ['DiffCoef1', 'DiffCoef2', 'DiffCoef3', 'DiffCoef4', 'kappa', 'beta', 'gamma', 'delta','zeta']
#-- params = (DiffCoef1,DiffCoef2,DiffCoef3,DiffCoef4,kappa,beta,gamma,delta,zeta)


#%% load data of a contry
dict_confirmed, dict_recovered, dict_deaths = load_dataset(countries)

# plot the data of the country
# cities_location=[1,3,5]
# plot_data_of_locations(cities_location, dict_confirmed, dict_recovered, dict_deaths)

# print('\n\n--> fitting the COVID10 data of ',country, '\n     -> Total population =',Nt,'\n\n Please wait!....... ')


# #% ODE parameters
# N=85;
# Tf=85.0;
# dt=Tf/(N-1);
# t = np.linspace(0, Tf, N)


#%%

# load the data / conpute the relative states S
Nt=66.99e6
position, t_exp, t_exp0, I_nz, R_nz, D_nz = Get_data_of_country_pde(dict_confirmed,dict_recovered, dict_deaths)    

S_nz0= I_nz +R_nz+D_nz ;  S_nz=Nt-S_nz0


position['r']=position['Lat']**2  + position['Long']**2
position['country']=countries
position_ = position.sort_values(by=['r'], ascending=True)

x=np.linspace(0 , Lx, Nx)
position_['x']=x
position_=position_.reset_index()



#-- plt.plot(t_exp,I_nz) ; plt.show()

# build the fiiting dataset
Noisy_data=np.concatenate((S_nz,I_nz, R_nz, D_nz ),axis=1) 
t_exp = np.linspace(0 , max(t_exp), I_nz.shape[0])

#%%   Fitting 
#% very important dant  chnge
#-----------------------------------------!!!
Noisy_data=np.flipud(Noisy_data)
t_exp0=np.flipud(t_exp0)
#-----------------------------------------!!!


# #% add more sample to imporve fitting accuracy
# from scipy import signal
# # Noisy_data = Noisy_data[1500:2000]
# # t_exp = t_exp[1500:2000]
# # t_exp0 = t_exp0[1500:2000]


# Nup=Nup*t_exp.shape[0]
# Noisy_data = signal.resample(Noisy_data, Nup)
# t_exp = np.linspace(t_exp[-1],0 , Nup)
    
# if Ndown>1: # Downsampling
#     Noisy_data = Noisy_data[::Ndown]
#     t_exp = t_exp[::Ndown]
#     t_exp0 = t_exp0[::Ndown]
       
    
nb_loc=position.shape[0]


# store the states
S_nz_= Noisy_data[:,0:nb_loc]
I_nz= Noisy_data[:,1*nb_loc:2*nb_loc]
R_nz= Noisy_data[:,2*nb_loc:3*nb_loc]
D_nz= Noisy_data[:,3*nb_loc:4*nb_loc]


NT_nz=S_nz+I_nz+R_nz+D_nz

# initial states
S0=S_nz[0,:]
I0=I_nz[0,:]
R0=R_nz[0,:]
D0=D_nz[0,:]

init_state=np.concatenate((S0,I0, R0, D0 )) 


t=t_exp

# % Define the fitting models 
def residuals(p):
    p = tuple(p)
    sim_P = odeint(Covid19System_pde, init_state, t_exp, args = p)
    res = sim_P - Noisy_data
    return res.flatten()


#  Fit the genrate data to the ODE
from scipy.optimize import leastsq
# initial_guess = [1, 1, 1, 1, 1, 1, 1, 1]
initial_guess = [0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5]


# LS method fitting
for cnt in range(0,NL):
    fitted_params = leastsq(residuals, initial_guess)[0]
    initial_guess=list(fitted_params)


# simulate the fitted model
fitted_params = tuple(fitted_params)
sol_fitted =odeint(Covid19System_pde, init_state, t, args = fitted_params)

# store the states

S_nz_= Noisy_data[:,0:nb_loc]
I_nz= Noisy_data[:,1*nb_loc:2*nb_loc]
R_nz= Noisy_data[:,2*nb_loc:3*nb_loc]
D_nz= Noisy_data[:,3*nb_loc:4*nb_loc]


S_hat= sol_fitted[:,0:nb_loc]
I_hat= sol_fitted[:,0:nb_loc]
R_hat= sol_fitted[:,0:nb_loc]
D_hat= sol_fitted[:,0:nb_loc]

NT_hat=S_hat+I_hat+R_hat+D_hat

#%%  plot figure  with name
from datetime import datetime
import os
now = datetime.now() # current date and time
tag_time = now.strftime("%m-%d-%Y")
my_path='results/ode4s_fitted_hpkins'+tag_time + '/' ; 
if not os.path.exists(my_path):
    os.mkdir(my_path)
    
filename_rslt=my_path+'NL-'+str(NL)+'_Ns-stice'+str(cut_time_start)+'_'+str(cut_time_end)+'_of'


country='four contries'

# plot_fitting_results_country_4s(country, t_exp, S_nz[:,1], I_nz[:,1], R_nz[:,1], D_nz[:,1], t, S_hat[:,1], I_hat[:,1], R_hat[:,1], D_hat[:,1], t_exp0, filename_rslt)

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





#%%

import plotly.graph_objects as go

fig = go.Figure(go.Surface(
    contours = {
        "x": {"show": True, "start": x[0], "end": x[-1], "size": 0.04, "color":"white"},
        "z": {"show": True, "start": 0.5, "end": 0.8, "size": 0.05}
    },
    
    y = t_exp,
    z = I_hat
    
    ))
fig.update_layout(
        scene = {
            "xaxis": {"nticks": 20},
            "zaxis": {"nticks": 4},
            'camera_eye': {"x": 0, "y": -1, "z": 0.5},
            "aspectratio": {"x": 1, "y": 1, "z": 0.2}
        })
fig.show()

#%%   
fig = go.Figure(go.Surface(
    contours = {
        "x": {"show": True, "start": x[0], "end": x[-1], "size": 0.04, "color":"white"},
        "z": {"show": True, "start": 0.5, "end": 0.8, "size": 0.05}
    },
    
    y = t_exp,
    z = I_nz
    
    ))
fig.update_layout(
        scene = {
            "xaxis": {"nticks": 20},
            "zaxis": {"nticks": 4},
            'camera_eye': {"x": 0, "y": -1, "z": 0.5},
            "aspectratio": {"x": 1, "y": 1, "z": 0.2}
        })
fig.show()


  
#%% 
fig = go.Figure(go.Surface(
    contours = {
        "x": {"show": True, "start": x[0], "end": x[-1], "size": 0.04, "color":"white"},
        "z": {"show": True, "start": 0.5, "end": 0.8, "size": 0.05}
    },
    
    y = t_exp,
    z = np.abs(I_nz-I_hat)
    
    ))
fig.update_layout(
        scene = {
            "xaxis": {"nticks": 20},
            "zaxis": {"nticks": 4},
            'camera_eye': {"x": 0, "y": -1, "z": 0.5},
            "aspectratio": {"x": 1, "y": 1, "z": 0.2}
        })
fig.show()



