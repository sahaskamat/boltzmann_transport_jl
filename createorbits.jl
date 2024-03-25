using NonlinearSolve
using Plots

"""
    createinitialpoints(c::Float64,E::function,n::Int64,doublefermisurface::Bool)

returns a vector of equally spaced points on the Fermi surface that can be used as starting points to create orbits
"""
function createinitialpoints(c,E,n,doublefermisurface)

    #makes c = c/2 if doublefermisurface is true
    if doublefermisurface
        c_eff = c/2
    else
        c_eff = c
    end 

    starting_z_coords = LinRange(-pi/c_eff,pi/c_eff,n+1) #create zcoordinates, each defining a plane on which points used for interpolation will be found. Exclude endpoint so that zone can be multiplied easily
    #pop!(starting_z_coords) #this deletes the last element of starting_z_coords, since it is the same as the first upto a reciprocal lattice vector. now starting_z_coords has exactly n elements

    philist = [0,pi] #list of phis along which to create the points. n points are created along each curve defined by phi and lying along the fermi surface

    vector_of_r0s = [] #empty vector, ith element of this vector corresponds to the list of 3-d points lying on philist[i]

    for phi in philist
        r0_list = [] #for each z_coord and phi, one element [x,y,z] will be found

        #a function that returns the energy at (r=r0,phi=phi,z=z0) in cylindrical coordinates
        function energyAlongPhi(r0,z0)
            return E([r0*cos(phi),r0*sin(phi),z0])
        end

        #find an r0 for each z0 by solving energyAlongPhi==0
        for i in eachindex(starting_z_coords)
            z0=starting_z_coords[i]
            u0= 0.5 #starting value for r0
            p=z0 #z0 given as parameter to NewtonRaphson method

            prob = NonlinearProblem(energyAlongPhi,u0,p)
            r0 = solve(prob,NewtonRaphson())
            push!(r0_list,[r0.u*cos(phi),r0.u*sin(phi),z0])
        end
         
        #once r0_list is made, append it to vector_of_r0s
        push!(vector_of_r0s,Array(r0_list))
    end 

    #diagnostics, this plots a x-z plane projection of calculated initial points
    for points_along_phi in vector_of_r0s
        matrix_of_points = hcat(points_along_phi...)
        plot!(matrix_of_points[1,:],matrix_of_points[2,:],matrix_of_points[3,:])
    end
    savefig("last_initial_points_list.png")
    
    #return [v1,v2,v3...] where each v_n is a vector 
    return vector_of_r0s
end

