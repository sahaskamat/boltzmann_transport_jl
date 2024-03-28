using Symbolics
using LinearAlgebra
using Plots

include("dispersion.jl")
include("createorbits.jl")

disp = make_NdLSCO_dispersion(160f-3,-0.1364,0.0682,0.0651,0.8243)
initialpoints = createinitialpoints(disp["c"],disp["E"],10,true)
dense_initialpoints = interpolate3d(initialpoints,10)
extended_intialpoints = extendedzonemultiply(dense_initialpoints,2,disp["c"],true)

#diagnostics, this plots a x-z plane projection of calculated initial points
for points_along_phi in extended_intialpoints
    matrix_of_points = hcat(points_along_phi...)
    plot!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
end
savefig("last_initial_points_list.png")