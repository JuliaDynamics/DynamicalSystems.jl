using DynamicalSystems, Makie

ds = Systems.nosehoover()

chaotic = get_state(ds)
stable = [0., 0.1, 0.5, 0.]

plane = (3, 0.0)
tf = 5000.0

psos = poincaresos(ds, plane, tf; u0 = [0, 0.1, 0.1])

scene = Makie.Scene(resolution = (750, 750))
Makie.scatter!(scene, psos[:, 1], psos[:, 2], markersize = 0.05, color=RGBf0(rand(3)...))

clicks = Node(Point2f0[get_state(ds)[1:2]])
on(scene.events.mousebuttons) do buttons
   if ispressed(scene, Mouse.left)
       pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
       psos = poincaresos(ds, plane, tf; u0 = [pos..., 0.1], warning = false)
       Makie.scatter!(scene, psos[:, 1], psos[:, 2], markersize = 0.05,
       color=RGBf0(rand(3)...))
   end
   display(scene)
   return
end
