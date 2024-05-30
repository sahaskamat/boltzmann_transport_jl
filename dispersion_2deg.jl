using Symbolics
using LinearAlgebra

function make_2deg_dispersion(doublefermisurface=false)
    #No T parameters to fit in 2deg case, this is just for checking if the code is consistent
    #These are known lattice parameters that will not be fit
    a=1 #angstroms
    b=a #angstroms
    c=2 #angstroms

    #makes c = c/2 if doublefermisurface is true
    if doublefermisurface
        c_eff = c/2
    else
        c_eff = c
    end 

    mu = 7

    @syms kx ky kz #variables in terms of which energies will be defined

    E = -mu + 3.8099820794*(kx^2 + ky^2 + 0.3*cos(kz*(c_eff))) #2d free electron dispersion, E is in eV and k is in (angstrom-1)
    #hbar^2/2m (kx^2 + ky^2 + arbitrary warping) - mu
    #the arbitrary warping is defined such that total fermi surface area remains constant

     
    Dx,Dy,Dz = Differential(kx),Differential(ky),Differential(kz) #define derivative operators to take gradient
    gradE = [expand_derivatives(Dx(E)),expand_derivatives(Dy(E)),expand_derivatives(Dz(E))] #gradient of dispersion

    @syms Bx By Bz
    force = cross(gradE,[Bx,By,Bz]) #Force due to magentic field on wave packet f = de/dk x B
    force_built = build_function(force,[kx,ky,kz,Bx,By,Bz],expression=Val{false})[1] #force_built takes a 6-element kx ky kz Bx By Bz and returns de/dk x B

    E_built = build_function(E,[kx,ky,kz], expression=Val{false}) #E_built is a function that takes a k vector as an input and returns the Energy (a scalar)
    gradE_built = build_function(gradE,[kx,ky,kz],expression=Val{false})[1] #grad_E is a function that takes a k vector as an input and returns dE/dk (a 3-vector)

    """
        invtau(k(3-vector))

    takes a k vector as an input and outputs the scattering rate (inverse scattering time)
    units of tau are ps, invtau are ps-1
    """
    function invtau(k)
        invtau_iso  = (1/100) #ps-1

        return invtau_iso #no angle dependence for 2deg case
    end

    """
        dkperp(B(3-vector),dkz(3-vector),dedk(function))

    the length element lying on the Fermi surface, perpendicular to the k-propogation direction under a magnetic field B

    B is a 3-vector pointing in the direction of the magnetic field
    dkz is a 3-vector that goes from any point on the plane of the current orbit to the plane of the next orbit
    dedk is a 3-vector that points perpendicular to the FS (gradient dE/dk)
    """
    function dkperp(B,dkz,dedk)
        nvec = cross(dedk,cross(dedk,B)) #nvec = dedk x (dedk x B)
        scalar_term = dot(dkz,B)/dot(nvec,B)
        dkperp = scalar_term*nvec

        return dkperp
    end

    return_dict = Dict(
        "a" => a,
        "b" => b,
        "c" => c,
        "c_eff" => c_eff,
        "mu" => mu,
        "E" => E_built,
        "gradE" => gradE_built,
        "force" => force_built,
        "invtau" => invtau,
        "dkperp" => dkperp
    )

    return return_dict
end