using DynamicalSystems

# Basic tests if modules were loaded properly
ds = Systems.towel()
traj = trajectory(ds, 100)
R = Reconstruction(traj[:, 1], 2, 2)

ls = lyapunovs(ds, 1000)
l = numericallyapunov(R, 1:5)
