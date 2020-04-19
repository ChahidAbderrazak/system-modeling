
plt.plot(t_exp_nrw, Noisy_data_nrw[:,1], label='data0')
plt.plot(t_exp_hpks, Noisy_data_hpks[:,1], label='Hpks')
plt.legend()


Noisy_data_hpks2=np.flipud(Noisy_data_hpks)
