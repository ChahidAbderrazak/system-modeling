function x=downsample_data(x, dt, shitf)
x=x(1:1/dt:end);x=x(1:end-shitf);
end