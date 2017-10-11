import numpy as np

nx = 2
nt = 4

field3d = np.random.rand(nt, nx, nx)
field2d = np.random.rand(nx, nx)

print field3d
print np.nanmean(field3d.flatten())
print field2d

#print np.nanmean(field3d.flatten())

field3d_mask = np.zeros(field3d.shape, dtype=bool)

for t in range(nt):
    field3d_mask[t,:,:] = field2d > 0.3



field3d = np.ma.array(field3d, mask=field3d_mask)

print np.nanmean(field3d.flatten())
print field2d
print field3d