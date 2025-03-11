using Symbolics
using LinearAlgebra
using Plots

include("dispersion.jl")
include("createinitialpoints.jl")
include("createorbits.jl")
include("conductivity.jl")


for H in [10,20,50,100]
    sigmalist = []
    rholist = []
    Tlist = []
    arealist = []

    for T in LinRange(0,300,30)
        invtau_aniso = 63.823 #ps-1
        invtau_iso  = 1.2*0.13*T + 8.7 #ps-1, agrees for the 25K value with the correct slope
        #functions that DONT depend on B
        disp = make_NdLSCO_dispersion(160e-3,-0.1364,0.0682,0.0651,0.8243,true,invtau_aniso,invtau_iso)
        initialpoints = createinitialpoints(disp["c_eff"],disp["E"],1000)
        extended_intialpoints = extendedzonemultiply(initialpoints,5,disp["c_eff"])
        interpolatedcurves = interpolate3d(initialpoints,disp["c_eff"])

        B = [0,0,H]

        intersectionpoints,dkz = makeorbitpoints(extended_intialpoints,interpolatedcurves,B,300,disp["c_eff"])
        intersectionpoints_temp = deepcopy(intersectionpoints) #this is because createorbits mutates intersectionpoints
        orbits = createorbits!(intersectionpoints_temp,B,disp["gradE"],15,4,10000,1e-10,1e-9,0.02)
        new_orbits = orbitCleanUp(orbits,0.01)
        sigma,area = createSigma(new_orbits,disp["invtau"],disp["gradE"],B,dkz)

        push!(sigmalist,sigma)
        push!(arealist,area)
        push!(Tlist,T)
        push!(rholist,inv(sigma))

        println("Calculation completed for T = ",T,", H = ",H)
    end


    rhoxylist= [rho[1,2]/10000 for rho in rholist]
    plot(Tlist,rhoxylist)
    savefig("rhoXYvsT"*string(H)*".png")

    plot(Tlist,arealist)
    savefig("areavstheta.png")

    open("rhoXYvsT_julia"*string(H)*".dat", "w") do f
        for (id,rhoxy) in enumerate(rhoxylist)
            println(f,Tlist[id]," ",rhoxy)
        end
    end
end