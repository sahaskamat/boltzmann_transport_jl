using Symbolics
using LinearAlgebra
using Plots

include("dispersion.jl")
include("createinitialpoints.jl")
include("createorbits.jl")
include("conductivity.jl")

B=[0,0,45]  #a field of 45T at 0 deg to Z axis


disp = make_NdLSCO_dispersion(160f-3,-0.1364,0.0682,0.0651,0.8243,true)
initialpoints = createinitialpoints(disp["c_eff"],disp["E"],200)
extended_intialpoints = extendedzonemultiply(initialpoints,5,disp["c_eff"])
interpolatedcurves = interpolate3d(initialpoints,disp["c_eff"])
intersectionpoints,dkz = makeorbitpoints(extended_intialpoints,interpolatedcurves,B,5,disp["c_eff"])
intersectionpoints_temp = deepcopy(intersectionpoints) #this is because createorbits mutates intersectionpoints
orbits = createorbits!(intersectionpoints_temp,B,disp["gradE"],10)
new_orbits = orbitCleanUp(orbits,0.1)
sigma,area = createSigma(new_orbits,disp["invtau"],disp["gradE"],B,dkz)

print(sigma)

println("rhozz = ",inv(sigma)[3,3])
println("Area = ",area)

open("orbitpoints_julia.dat", "w") do f
    for orbitsinplane in new_orbits
        for orbit in orbitsinplane
            for point in orbit
                println(f,point[1]," ",point[2]," ",point[3])
            end
        end
    end
end