using LinearAlgebra
using DifferentialEquations

"""
    createorbits(intersectionpoints::vector of vector of 3-vectors,B::3-vector,gradE::function)

    returns: orbits(vector of vector of vector of 3-vectors, with each vector of 3-vectors corresponding to  ONE orbit in a plane)
    i.e. orbits[i][j] is a vector of points lying along one orbit
         orbits[i] is a vector of vectors, with each vector containing curves lying on a single plane
"""
function createorbits(intersectionpoints,B,gradE)
    orbits = []
    B_normalized = B/norm(B)

    for intersectionpoints_in_plane in intersectionpoints
        orbits_in_plane = []

        for point in intersectionpoints_in_plane
            function rhs!(dk,k,p,t)
                dk .= cross(gradE(k),B_normalized) #dk/dt = v x B
            end

            k0 = point
            timespan = (0,8)
            prob = ODEProblem(rhs!,k0,timespan)

            sol = solve(prob)
            push!(orbits_in_plane,sol.u)
        end

        push!(orbits,orbits_in_plane)
    end
    return orbits
end