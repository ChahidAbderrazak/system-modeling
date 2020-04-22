import numpy as np
import matplotlib.pyplot as plt



def laplacian(Z):
    Ztop = Z[0:-2, 1:-1]
    Zleft = Z[1:-1, 0:-2]
    Zbottom = Z[2:, 1:-1]
    Zright = Z[1:-1, 2:]
    Zcenter = Z[1:-1, 1:-1]
    
    return (Ztop + Zleft + Zbottom + Zright -
            4 * Zcenter) / dx**2


def show_patterns(S, ax=None):
    ax.imshow(S, cmap=plt.cm.copper,
              interpolation='bilinear',
              extent=[-1, 1, -1, 1])
    ax.set_axis_off()


  
def SIRD_pde(init_state, Nx, t, *args):
    

    DiffCoefS=args[0]
    DiffCoefI=args[1]
    DiffCoefR=args[2]
    DiffCoefD=args[3]
    kappa=args[4]        
    beta=args[5]
    gamma=args[6]
    delta=args[7]


    #-- state=init_state
    # t the states
    Nx=int(init_state.shape[0]/4)
    


    S= init_state[0:Nx]
    I= init_state[1*Nx:2*Nx]
    R= init_state[2*Nx:3*Nx]
    D= init_state[3*Nx:4*Nx]
    N=S + I + R + D
    

    
    # We take the Values of S and I inside the grid.
    Sn = S[1:-1, 1:-1]
    In = I[1:-1, 1:-1]
    Rn = R[1:-1, 1:-1]
    Dn = D[1:-1, 1:-1]
    Nn=Sn + In + Rn + Dn
    
    
    
     # We Spdate the Iariables.
    S[1:-1, 1:-1]= Sn + dt * (DiffCoefS * deltaS - (1/Nn)*beta*(1-kappa)*Sn*In)
    I[1:-1, 1:-1]= In + dt * (DiffCoefI * deltaI + (1/Nn)*beta*(1-kappa)*Sn*In -gamma*(1-delta)*In -delta*In)
    R[1:-1, 1:-1]= In + dt * (DiffCoefR * deltaR + gamma*(1-delta)*In) 
    D[1:-1, 1:-1]= In + dt * (DiffCoefD * deltaD + delta*In) 

    # Neumann conditions: derivatives at the edges
    # are null.
    for Z in (S, I, R, D):
        Z[0, :] = Z[1, :]
        Z[-1, :] = Z[-2, :]
        Z[:, 0] = Z[:, 1]
        Z[:, -1] = Z[:, -2]
             
            
    return [S, I, R, D]


# #%% anmation defScion
# print('\n\n --> preparing the animation please wait..:):)')
# # skip_sz=10
# # sol_S=sol_S[::skip_sz,:,:]
# # sol_I=sol_I[::skip_sz,:,:]
# # sol_R=sol_R[::skip_sz,:,:]
# # sol_D=sol_D[::skip_sz,:,:]
       
# # plt.plot(sol_I[:,0,0])
# # plt.show()
    

# xs = np.linspace(0, 1, size)
# ys = np.linspace(0, 1, size)
# X, Y = np.meshgrid(xs, ys)

# Z=sol_I.T

# #
# from mpl_toolkits import mplot3d
# from matplotlib import cm
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# T=len(Z) # frame nSmber of the animation



# N = X.shape[0] # Meshsize
# fps = 10 # frame per sec

# zmin= np.min(Z); #print(zmin)
# zmax= np.max(Z); #print(zmax)

# def update_plot(frame_number, zarray, plot):
#     plot[0].remove()
#     plot[0] = ax.plot_surface(X, Y, Z[frame_number,:,:], cmap="magma")

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# plot = [ax.plot_surface(Y, X, Z[0,:,:], color='0.75', rstride=1, cstride=1)]
# ax.set_zlim(zmin,1.2*zmax)


# ani = animation.FuncAnimation(fig, update_plot, T, fargs=(Z, plot), interval=1000/fps)

# fn = '../animation/PDE-2D'
# ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)
# ani.save(fn+'.gif',writer='imagemagick',fps=fps)


