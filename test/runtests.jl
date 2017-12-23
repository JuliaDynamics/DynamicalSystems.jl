using DynamicalSystems

# Basic tests if modules were loaded properly
ds = Systems.towel()
traj = trajectory(ds, 1000)
R = Reconstruction(traj[:, 1], 2, 2)
a = lyapunov(ds, 100)
