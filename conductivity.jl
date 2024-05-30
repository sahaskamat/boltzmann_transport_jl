using LinearAlgebra
using SparseArrays


"""
    createSigma(orbits::vector of vector of vector of 3-vectors,invtau::function,dedk::function,B::3-vector,dkz::3-vector)
    sigma is returned in units of 10^9 S/m (SI units of conductivity), totalarea is returned in units of Angstrom^-2

returns: sigma::(3x3 matrix), totalarea::Float64
"""
function createSigma(orbits,invtau,dedk,B,dkz)
    dedk_list,unitvec_dedk_list = create_dedk_list(orbits,dedk)
    alpha = solveAlphaSparse(orbits,invtau,dedk,dedk_list,B)

    perptermlist = [dkperp(B,dkz,dedk_vec) for dedk_vec in eachrow(dedk_list)] #perptermlist[i] is a vector that lies along the fermi surface pointing from the ith point to the orbit above it

    #create a list of vectors pointing from one point on an orbit to the next
    pointingvectors = Vector{Float64}[]
    for orbitsinplane in orbits
        for orbit in orbitsinplane
            m=length(orbit)
            for id in eachindex(orbit)
                push!(pointingvectors,orbit[mod1(id+1,m)] - orbit[id])
            end
        end
    end

    patcharealist = broadcast(norm,broadcast(cross,perptermlist,pointingvectors)) #patcharealist[i] is the integration patch area corresponding to the ith point

    sigma = zeros(3,3) #conductivity tensor
    for mu in 1:3
        for nu in 1:3
            #create the (mu,nu) element of sigma
            sigma[mu,nu] = (3.699/(4*pi^3))*sum(patcharealist .* unitvec_dedk_list[:,nu] .* alpha[:,mu]) #e^2/4pi^3hbar^2 int(dk3 unitvec_dedk_list^b alpha^a), this is in units of 10^9 S/m (SI units)
        end
    end

    totalarea = sum(patcharealist) #total area of the fermi surface that is integrated over

    return sigma,totalarea
end

function getorbitsize(orbits) #finds the number of 3 vectors (i.e. points) in orbits
    n=0
    for orbitsinplane in orbits
        for orbit in orbitsinplane
            n+=length(orbit)
        end
    end
    return n
end

"""
    createAmatrix(orbits::vector of vector of vector of 3-vectors,invtau::function,dedk::function,B::3-vector)

returns: Adata,Aposition_i,Aposition_j,n
Adata[k] is the kth non zero element of A at position Aposition_i[k],Aposition_j[k] with size nxn
"""
function createAmatrix(orbits,invtau,dedk,B)
    n = getorbitsize(orbits)

    Adata = Float64[] #A[n] element corresponds to non zero matrix element in A
    Aposition_i = Int64[] # Aposition_i[n] element corresponds to ith coordinate of element in A
    Aposition_j = Int64[] #Aposition_j[n] element corresonds to jth coordinate of element in A

    submatrixindex = 1 #denotes the position of the beginning of the current submatrix

    i=1 #iterates over the hilbert space

    for orbitsinplane in orbits
        for orbit in orbitsinplane
            m=length(orbit) #mxm is the size of the current submatrix

            for (id,point) in enumerate(orbit)
                #diagonal term coming from scattering out
                push!(Adata,invtau(point))
                push!(Aposition_i,i)
                push!(Aposition_j,i)

                #positions of the index corresponding to the next and previous states on this orbit in the hilbert space
                i_next = mod((i + 1) - submatrixindex,m) + submatrixindex
                i_prev = mod((i - 1) - submatrixindex,m) + submatrixindex

                graddata = norm(cross(dedk(point),B))/norm(orbit[mod1(id+1,m)]-orbit[mod1(id-1,m)]) #graddata = norm(dedk(orbit[i]) x B)/norm(orbit[i+1] - orbit[i-1]) indices are upto modulo array size
                graddata_with_units = graddata/(6.582119569^2)

                #off diagonal terms that simulate the derivative
                push!(Adata,graddata_with_units)
                push!(Aposition_i,i)
                push!(Aposition_j,i_next)

                push!(Adata,-graddata_with_units)
                push!(Aposition_i,i)
                push!(Aposition_j,i_prev)

                i+=1
            end

            submatrixindex+= m
        end
    end

    return Adata,Aposition_i,Aposition_j,n
end

"""
creates dedk_list, where dedk_list[i] is the ith point in orbits (ordered by iterating over each orbit, then a plane, then planes of orbits i.e. in to out)
"""
function create_dedk_list(orbits,dedk)
    dedk_list = Vector{Float64}[] #dedk_list[i] is the value of de/dk for the ith point in the hilbert space

    for orbits_in_plane in orbits
        for orbit in orbits_in_plane
            for point in orbit
                push!(dedk_list,dedk(point))
            end
        end
    end

    unitvec_dedk_list = dedk_list./[norm(point) for point in dedk_list] #unitvec_dedk_list[i] is the unitvector pointing along dedk_list[i]
    unitvec_dedk_list = transpose(hcat(unitvec_dedk_list...)) #rearranges to a form where the ith column of unitvec_dedk_list is the ith component (in xyz) of unitvec_dedk_list

    dedk_list = transpose(hcat(dedk_list...)) #rearranges dedk_list such that dedk_list[:,i] is now dedk_list[:][i] 

    return dedk_list,unitvec_dedk_list
end

"""
    solveAlphaSparse(orbits::vector of vector of vector of 3-vectors,invtau::function,dedk::function,dedk_list::nx3 matrix,B::3-vector)

solver that returns dedk_list/A using relaxation time approximation and exploiting the resulting sparesity of A
returns: alpha = dedk_list/A
"""
function solveAlphaSparse(orbits,invtau,dedk,dedk_list,B)
    Adata,Aposition_i,Aposition_j,n = createAmatrix(orbits,invtau,dedk,B)
    A_sparse = sparse(Aposition_i,Aposition_j,Adata,n,n,+)
    
    alpha = zero(dedk_list) #dedk_list/A
    ldiv!(alpha,factorize(A_sparse),dedk_list)

    println(size(A_sparse))
    return alpha
end

function dkperp(B,dkz,dedk) #calculates vector perpendicular to a given point to tile the fermi surface
    #dkz is a vector connecting this orbital plane to the next
    nvec = cross(dedk,cross(dedk,B)) #nvec = dedk x (dedk x B)
    scalar_term = dot(dkz,B)/dot(nvec,B)
    return scalar_term*nvec
end