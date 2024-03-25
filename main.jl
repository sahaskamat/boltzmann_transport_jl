using Symbolics
using LinearAlgebra

include("dispersion.jl")
include("createorbits.jl")

disp = make_NdLSCO_dispersion(160f-3,-0.1364,0.0682,0.0651,0.8243)
createinitialpoints(disp["c"],disp["E"],100,true)

#plot disp["E"] as a diagnostic
#function energyAlongPhi(r0)
#    return disp["E"]([r0,0,0])
#end

#r_range = LinRange(0,5,1000)
#E_range = broadcast(energyAlongPhi,r_range)
#plot(r_range,E_range)
#savefig("myplot.png")