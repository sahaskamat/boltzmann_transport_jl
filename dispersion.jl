using Symbolics
using LinearAlgebra

function make_NdLSCO_dispersion(T,T1multvalue,T11multvalue,Tzmultvalue,mumultvalue)
    #These are known lattice parameters that will not be fit
    a=3.75 #angstroms
    b=a #angstroms
    c=2*6.6 #angstroms

    #Unknown tight binding parameters that may have to be fit
    T1 = T1multvalue*T
    T11 = T11multvalue*T
    Tz = Tzmultvalue*T
    mu = -mumultvalue*T

    @syms kx ky kz #variables in terms of which energies will be defined

    #Dispersion restricted by tetragonal symmetry for NdLSCO
    E = -mu - 2*T*(cos(kx*a) + cos(ky*a)) -4*T1*cos(kx*a)*cos(ky*a) - 2*T11*(cos(2*kx*a)+cos(2*ky*a)) - 2*Tz*cos((kx*a)/2)*cos((ky*b)/2)*cos((kz*c)/2)*((cos(kx*a) - cos(ky*b))^2)
     
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
        invtau_iso  = 12.595 #ps-1
        invtau_aniso = 63.823 #ps-1
        nu=12

        angledependence=abs((k[1]^2 - k[2]^2)/(k[1]^2 + k[2]^2))^nu

        return invtau_iso + invtau_aniso*angledependence
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
        "T1" => T1,
        "T11" => T11,
        "Tz" => Tz,
        "mu" => mu,
        "E" => E_built,
        "gradE" => gradE_built,
        "force" => force_built,
        "invtau" => invtau,
        "dkperp" => dkperp
    )

    return return_dict
end

