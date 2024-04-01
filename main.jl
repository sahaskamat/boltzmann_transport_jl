using Symbolics
using LinearAlgebra
using Plots

include("dispersion.jl")
include("createorbits.jl")

disp = make_NdLSCO_dispersion(160f-3,-0.1364,0.0682,0.0651,0.8243,true)
initialpoints = createinitialpoints(disp["c_eff"],disp["E"],90)
extended_intialpoints = extendedzonemultiply(initialpoints,2,disp["c_eff"])
interpolatedcurves = interpolate3d(initialpoints,disp["c_eff"])

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
savefig("interpolatedpoints.png")