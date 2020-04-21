import numpy as np
import matplotlib.pyplot as plt

# parameters of the model:
DiffCoefS=0.0005;
DiffCoefI=0.0005;
DiffCoefR=0.0005;
DiffCoefD=0.0005;
kappa=0.95; 
beta=0.25;
gamma=0.03;
delta=0.00175;

# We discretize time and space.
size = 20  # size of the 2D grid
Lx=1.0;
dx=Lx/(size-1);  # space step


n=1000;      # number of iterations
T=360.0;      # total time
dt=T/(n-1)   # time step


# We initialize the Cariables  
S = np.ones([size, size])*105529845
I = np.ones([size, size])*13.6275
R = np.zeros([size, size])
D = np.zeros([size, size])


def laplacian(Z):
    Ztop = Z[0:-2, 1:-1]
    Zleft = Z[1:-1, 0:-2]
    Zbottom = Z[2:, 1:-1]
    Zright = Z[1:-1, 2:]
    Zcenter = Z[1:-1, 1:-1]
    
    # return (Ztop + Zleft + Zbottom + Zright -
    #         4 * Zcenter) / dx**2

    return ( Zleft + Zright -
            2 * Zcenter) / dx**2


def show_patterns(S, ax=None):
    ax.imshow(S, cmap=plt.cm.copper,
              interpolation='bilinear',
              extent=[-1, 1, -1, 1])
    ax.set_axis_off()
    
    
    
    
# we simSlate the system 

fig, axes = plt.subplots(3, 3, figsize=(8, 8))

shot_t=4
step_plot = n // 9

# We simSlate the PDE with the finite difference
# method.

sol_S=np.zeros((n,size,size))
sol_I=np.zeros((n,size,size))
sol_R=np.zeros((n,size,size))
sol_D=np.zeros((n,size,size))




for i in range(n):
    # We compSte the Laplacian of S and I.
    deltaS = laplacian(S)
    deltaI = laplacian(I)
    deltaR = laplacian(R)
    deltaD = laplacian(D)
    
    
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

    sol_S[i,:,:]=S
    sol_I[i,:,:]=I
    sol_R[i,:,:]=R
    sol_D[i,:,:]=D

    # NeSmann conditions: deriIatiIes at the edges
    # are nSll.

    # for Z in (S, I):
    #     Z[0, :] = Z[1, :]
    #     Z[-1, :] = Z[-2, :]
    #     Z[:, 0] = Z[:, 1]
    #     Z[:, -1] = Z[:, -2]
        
        
    # for Z in (R, D):
    #     Z[0, :] = Z[1, :]
    #     Z[-1, :] = Z[-2, :]
    #     Z[:, 0] = Z[:, 1]
    #     Z[:, -1] = Z[:, -2]


    
  
    # We plot the state of the system at
    # 9 different times.
    if i % step_plot == 0 and i < 9 * step_plot:
        ax = axes.flat[i // step_plot]
        show_patterns(S, ax=ax)
        ax.set_title(f'$t={i * dt:.2f}$')

            

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
show_patterns(S, ax=ax)


#%% anmation defScion
print('\n\n --> preparing the animation please wait..:):)')
skip_sz=10
sol_S=sol_S[::skip_sz,:,:]
sol_I=sol_I[::skip_sz,:,:]
sol_R=sol_R[::skip_sz,:,:]
sol_D=sol_D[::skip_sz,:,:]
        
xs = np.linspace(-1, 1, size)
ys = np.linspace(-1, 1, size)
X, Y = np.meshgrid(xs, ys)

Z=sol_I

#%%
from mpl_toolkits import mplot3d
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

T=len(Z) # frame nSmber of the animation



N = X.shape[0] # Meshsize
fps = 10 # frame per sec

zmin= np.min(Z); #print(zmin)
zmax= np.max(Z); #print(zmax)

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, Z[frame_number,:,:], cmap="magma")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plot = [ax.plot_surface(X, Y, Z[0,:,:], color='0.75', rstride=1, cstride=1)]
ax.set_zlim(zmin,1.2*zmax)
ani = animation.FuncAnimation(fig, update_plot, T, fargs=(Z, plot), interval=1000/fps)

fn = 'animation/PDE-2D'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)
ani.save(fn+'.gif',writer='imagemagick',fps=fps)





