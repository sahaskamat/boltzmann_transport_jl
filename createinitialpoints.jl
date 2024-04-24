using NonlinearSolve
using DataInterpolations
using LinearAlgebra

"""
    createinitialpoints(c::Float64,E::function,n::Int64,doublefermisurface::Bool)

returns a vector of equally spaced points on the Fermi surface that can be used as starting points to create orbits
"""
function createinitialpoints(c_eff,E,n)

    starting_z_coords = collect(LinRange(-pi/c_eff,pi/c_eff,n)) #create zcoordinates, each defining a plane on which points used for interpolation will be found. Exclude endpoint so that zone can be multiplied easily

    philist = [0,pi] #list of phis along which to create the point s. n points are created along each curve defined by phi and lying along the fermi surface

    vector_of_r0s = [] #empty vector, ith element of this vector corresponds to the list of 3-d points lying on philist[i]

    for phi in philist
        r0_list = Vector{Float64}[] #for each z_coord and phi, one element [x,y,z] will be found

        #a function that returns the energy at (r=r0,phi=phi,z=z0) in cylindrical coordinates
        function energyAlongPhi(r0,z0)
            return E([r0*cos(phi),r0*sin(phi),z0])
        end

        #find an r0 for each z0 by solving energyAlongPhi==0
        for i in eachindex(starting_z_coords)
            z0=starting_z_coords[i]
            u0= 0.5 #starting  value for r0
            p=z0 #z0 given as parameter to NewtonRaphson method

            prob = NonlinearProblem(energyAlongPhi,u0,p)
            r0 = solve(prob,NewtonRaphson())
            push!(r0_list,[r0.u*cos(phi),r0.u*sin(phi),z0])
        end
         
        #once r0_list is made, append it to vector_of_r0s
        push!(vector_of_r0s,Array(r0_list))
    end 
    
    #return [v1,v2,v3...] where each v_n is a vector 
    return vector_of_r0s
end

"""
    interpolate3d(vector_of_r0s::vector of vector of 3-vectors,c_eff::Float64)

returns a vector of vector functions that returns the [x,y,z] point for a z input
"""
function interpolate3d(vector_of_r0s,c_eff)

    vector_of_interpolated_curves=[] #ith element of this vector is an interpolating function giving [x,y,z] along phi_i for an input z

    for points_along_phi in vector_of_r0s
        z_points = hcat(points_along_phi...)[3,:] #vector of z coordinates corresponding to each point in points_along_phi
        itp = CubicSpline(points_along_phi,z_points,extrapolate=true) #interpolate all the points in points_along_phi, each parameterised by its z coordinate. extrapolate has been set to true since function ignores last point in interpolation bounds
        itp_extended(z) = Float64[itp((mod(z+pi/c_eff,2*pi/c_eff) - pi/c_eff))[1],itp((mod(z+pi/c_eff,2*pi/c_eff) - pi/c_eff))[2],z] #take the x and y coordinates from the interpolating function acting on z modulo G
        push!(vector_of_interpolated_curves,itp_extended)
    end

    return vector_of_interpolated_curves
end


"""
    extendedzonemultiply(vector_of_r0s::vector of vector of 3-vectors,n::Int,c::Float64,doublefermisurface::Bool)

returns a vector of vector of 3-vectors, now replicated 2n+1 times along the c axis direction, n below and n above (extended zone scheme)
this should ideally be replaced by an algorithm that respects the periodic nature of the fermi surface better in the future
"""
function extendedzonemultiply(vector_of_r0s,n,c_eff)

    extended_vector_of_r0s=[] #same as vector_of_r0s, but with extended n times in either c axis direction

    for points_along_phi in vector_of_r0s
        extended_points_along_phi=points_along_phi
        pop!(extended_points_along_phi) #this deletes the last element of extended_points_along_phi, since it is the same as the first upto a reciprocal lattice vector
        for i = 1:n
            points_before = points_along_phi .- [[0,0,(2*pi*i)/c_eff]]
            points_after = points_along_phi .+ [[0,0,(2*pi*i)/c_eff]]
            extended_points_along_phi = [points_before;extended_points_along_phi;points_after]
        end

        push!(extended_vector_of_r0s,extended_points_along_phi)
    end 

    return extended_vector_of_r0s
end

"""
    makeorbitpoints(extended_vector_of_r0s::vector of vector of 3-vectors,interpolatedcurves::vector of interpolating functions,B::3-vector,n::Int,c_eff::Float64)

returns vector_of_intersectionpoints,dkz
vector_of_intersectionpoints: a n-length vector. Each element contains the points that lie on the ith orbital plane, which can now be used as an initial point to solve the orbtial differential equation.
dkz: vector connecting a given plane to the plane above it
"""
function makeorbitpoints(extended_vector_of_r0s,interpolatedcurves,B,n,c_eff)
    B_normalized = B/norm(B)

    vector_of_points_on_plane= [[0,0,z] for z in LinRange(-pi/c_eff,pi/c_eff,n+1)]
    pop!(vector_of_points_on_plane)

    planeequation(r,normalvector,pointonplane) = dot(r,normalvector) - dot(pointonplane,normalvector) #planeequation = 0 is the equation of the plane passing through pointonplane and perpendicular to normalvector

    vector_of_intersectionpoints = [] #each element of this vector contains intersectionpoints of the curves in interpolatedcurves and the ith point in vector_of_points_on_plane

    for pointonplane in vector_of_points_on_plane
        intersectionpoints = [] #vector of points that lie on the plane of pointonplane and intersect with interpolatedcurves
        for (id,discrete_curve) in enumerate(extended_vector_of_r0s)
            value_of_planeequation = broadcast(planeequation,discrete_curve,[B_normalized],[pointonplane]) #value of planeequation(x,B,pointonplane) for each x in discrete_curve
            twos_where_intersections = diff(broadcast(sign,value_of_planeequation)) #this looks like [0,0,0,2,0,0,1,1,0,0,0] where the 2 corresponds to the plane passing between two points in discrete_curve and the 1,1 is where the plane exactly hits a point in discrete_curve
            crossing_indices = findall(x -> (abs(x)==2),twos_where_intersections) #this counts all points where the plane passes between points in discrete_curve 
            hitting_indices = findall(x -> (abs(x)==1),twos_where_intersections) #this counts all points where the plane hits points in discrete_curve 

            for index in crossing_indices #use a binary search to find where planeequation goes to zero
                lowerbound = discrete_curve[index]
                upperbound = discrete_curve[index+1]

                intersection_point = binaryfindintersection(lowerbound,upperbound,planeequation,B,pointonplane,interpolatedcurves[id])
                push!(intersectionpoints,intersection_point)
            end

            for id in eachindex(hitting_indices) #this adds every even indexed element in hitting_indices to intersectionpoints
                if isodd(id)
                    index = hitting_indices[id]
                    intersection_point = discrete_curve[index+1]
                    push!(intersectionpoints,intersection_point)
                end
            end

        end
        push!(vector_of_intersectionpoints,intersectionpoints)
    end

    dkz = vector_of_points_on_plane[2] - vector_of_points_on_plane[1]
    return vector_of_intersectionpoints,dkz
end

"""
    binaryfindintersection(lowerbound::3-vector,upperbound::3-vector,planeequation::function,B::3-vector,pointonplane::3-vector,interpolatingcurve::function)

returns the 3-vector x where planeequation(x,B,pointonplane) = 0 on the interpolatingcurve connecting lowerbound and upperbound
"""
function binaryfindintersection(lowerbound,upperbound,planeequation,B,pointonplane,interpolatingcurve)

    z_middle  = (lowerbound[3] + upperbound[3])/2
    middlepoint = interpolatingcurve(z_middle)

    for i in 1:10
        if planeequation(lowerbound,B,pointonplane)*planeequation(middlepoint,B,pointonplane)<0 upperbound=middlepoint
        elseif planeequation(upperbound,B,pointonplane)*planeequation(middlepoint,B,pointonplane)<0 lowerbound=middlepoint
        end

        z_middle  = (lowerbound[3] + upperbound[3])/2
        middlepoint = interpolatingcurve(z_middle)
    end
    
    return middlepoint
end
