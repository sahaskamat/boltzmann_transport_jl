using Symbolics
using LinearAlgebra
using Plots

include("dispersion_2deg.jl")
include("createinitialpoints.jl")
include("createorbits.jl")
include("conductivity.jl")

B=[0,0,1]  #a field of 45T at 45 deg to Z axis


disp = make_2deg_dispersion()
initialpoints = createinitialpoints(disp["c_eff"],disp["E"],50)
extended_intialpoints = extendedzonemultiply(initialpoints,2,disp["c_eff"])
interpolatedcurves = interpolate3d(initialpoints,disp["c_eff"])
intersectionpoints,dkz = makeorbitpoints(extended_intialpoints,interpolatedcurves,B,5,disp["c_eff"])
intersectionpoints_temp = deepcopy(intersectionpoints) #this is because createorbits mutates intersectionpoints
orbits = createorbits!(intersectionpoints_temp,B,disp["gradE"],15,0.5,1000000,1e-11,1e-10,0.1)
new_orbits = orbitCleanUp(orbits,0.2)
sigma,area = createSigma(new_orbits,disp["invtau"],disp["gradE"],B,dkz)
rho = inv(sigma)
println("RHO XX = ",rho[1,1]," uOhm-mm")
println("RHO XY = ",rho[1,2]," uOhm-mm")

#diagnostics, this plots a x-z plane projection of calculated initial points
plot([0],[0],[0])
for points_along_phi in extended_intialpoints
    matrix_of_points = hcat(points_along_phi...)
    plot!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
end
savefig("last_initial_points_list.png")

#diagnostics to make sure interpolated curves are correct
plot([0],[0],[0])
for curve_along_phi in interpolatedcurves
    z_coords = LinRange(-5*pi/(disp["c_eff"]),5*pi/(disp["c_eff"]),110)
    matrix_of_points = hcat(broadcast(curve_along_phi,z_coords)...)
    plot!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
end

#diagnostics to make sure intersectionpoints are correct
for set_of_intersections_with_a_plane in intersectionpoints
    matrix_of_points = hcat(set_of_intersections_with_a_plane...)
    scatter!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
end

#diagnostics to make sure orbits are being made correctly
for orbits_in_plane in orbits
    for curve in orbits_in_plane
        matrix_of_points = hcat(curve...)
        plot!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
    end
end

#diagnostics to make sure equally spaced orbits are being made correctlyfor orbits_in_plane in orbits
for orbits_in_plane in new_orbits
    for curve in orbits_in_plane
        matrix_of_points = hcat(curve...)
        scatter!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
    end
end
savefig("interpolatedpoints.png")

