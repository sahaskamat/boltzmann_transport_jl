using Symbolics
using LinearAlgebra
using Plots

include("dispersion.jl")
include("createinitialpoints.jl")
include("createorbits.jl")
include("conductivity.jl")

#functions that DONT depend on B
disp = make_NdLSCO_dispersion(160e-3,-0.1364,0.0682,0.0651,0.8243,true)
initialpoints = createinitialpoints(disp["c_eff"],disp["E"],1000)
extended_intialpoints = extendedzonemultiply(initialpoints,5,disp["c_eff"])
interpolatedcurves = interpolate3d(initialpoints,disp["c_eff"])

sigmalist = []
rholist = []
thetalist = []
arealist = []

#functions that DO, swept over B
for theta in LinRange(0,80,40)
    B = [45*sin(deg2rad(theta)),0,45*cos(deg2rad(theta))]

    intersectionpoints,dkz = makeorbitpoints(extended_intialpoints,interpolatedcurves,B,300,disp["c_eff"])
    intersectionpoints_temp = deepcopy(intersectionpoints) #this is because createorbits mutates intersectionpoints
    orbits = createorbits!(intersectionpoints_temp,B,disp["gradE"],15)
    new_orbits = orbitCleanUp(orbits,0.01)
    sigma,area = createSigma(new_orbits,disp["invtau"],disp["gradE"],B,dkz)

    push!(sigmalist,sigma)
    push!(arealist,area)
    push!(thetalist,theta)
    push!(rholist,inv(sigma))

    println("Calculation completed for theta = ",theta)
end

rhozzlist= [rho[3,3]/10000 for rho in rholist]
plot(thetalist,rhozzlist)
savefig("rhovstheta.png")

plot(thetalist,arealist)
savefig("areavstheta.png")

open("rhoZZvsT_julia.dat", "w") do f
    for (id,rhozz) in enumerate(rhozzlist)
        println(f,thetalist[id]," ",rhozz)
    end
end