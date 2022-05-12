module LinesCurvesNodes

using LinearAlgebra, Rotations


function transform_vector(L, start_node, Θ)

    vector = [L; 0.0] 

    #rotate
    vector = Angle2d(Θ) * vector

    #shift
    vector = vector + start_node

    return vector

end


function generate_fillet(A, B, C, r, n)

    BC = C-B
    BA = A-B

    #https://www.youtube.com/watch?v=64jtpAV3XYY
    BR = BA + BC/norm(BC) * norm(BA)

    BR_unit = BR/norm(BR)

    #https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
    Θ = -atan(det([BA'
BC']), BA⋅BC)

    #https://www.cpp.edu/~hturner/ce220/circular_curves.pdf

    Δ = π - abs(Θ)

    T = tan(Δ/2) * r

    M = hypot(T, r)

    R = B + BR_unit * M

    a = B + BA/norm(BA) * T  

    Ra = a-R

    dΔ = range(0.0, sign(Θ) * Δ, n+1)

    fillet = [Angle2d(dΔ[i]) * Ra + R for i in eachindex(dΔ)]  #Angle2d is 2D rotation matrix 

    fillet = permutedims(reshape(hcat(fillet...), (length(fillet[1]), length(fillet))))

    return fillet

end


function discretize_vector(A, B, n)

    AB = B - A
    AB_unit = AB / norm(AB)
    spacing = range(0.0, norm(AB), n+1)
    segments = [A .+ AB_unit .* spacing[i] for i in eachindex(spacing)] 

    # segments = permutedims(reshape(hcat(segments...), (length(segments[1]), length(segments))))

    return segments

end


function rotate_nodes(nodes; rotation_axis, rotation_center, θ)

    num_points = size(nodes)[1]

    R = Angle2d(θ)

    rotated_nodes = nodes

    for i = 1:num_points

        x_centered = nodes[i,1] - rotation_center[1]
        y_centered = nodes[i, 2] - rotation_center[2]
        z_centered = nodes[i, 3] - rotation_center[3]

        if rotation_axis == "x"

            rotated_nodes[i,[2,3]] = R * [y_centered; z_centered] .+ [rotation_center[2]; rotation_center[3]]
            
        end

        if rotation_axis == "y"

            rotated_nodes[i,[1,3]] = R * [x_centered; z_centered] .+ [rotation_center[1]; rotation_center[3]]
            
        end

        if rotation_axis == "z"

            rotated_nodes[i,[1,2]] = R * [x_centered; y_centered] .+ [rotation_center[1]; rotation_center[2]]
            
        end

    end

    return rotated_nodes

end

function shift_nodes(nodes; Δx, Δy, Δz)

    nodes[:, 1] = nodes[:, 1] .+ Δx
    nodes[:, 2] = nodes[:, 2] .+ Δy
    nodes[:, 3] = nodes[:, 3] .+ Δz

    return nodes

end

function find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)

    x_index = findall(x->isapprox(x, xloc, atol=atol_x), nodes[:,1])
    y_index = findall(y->isapprox(y, yloc, atol=atol_y), nodes[:,2])
    z_index = findall(z->isapprox(z, zloc, atol=atol_z), nodes[:,3])

    index = intersect(intersect(x_index, y_index), z_index)

    return index

end


#use for subdividing line elements
function subdivide_line_segments(coords, num_sub_segments)

    num_nodes = size(coords)[1]
    num_segments = num_nodes - 1

    sub_segments = Array{Vector{Vector{Float64}}}(undef, 0)

    for i = 1:num_segments

            A = coords[i]
            B = coords[i+1]

            AB = B-A

            AB_unit = AB/norm(AB)

            length_AB = norm(AB)

            range_AB = range(0.0, length_AB, num_sub_segments+1)

            sub_segments = [sub_segments; [A .+ AB_unit .* range_AB[j] for j in eachindex(range_AB)]] 

            sub_segments = [round.(sub_segments[i], sigdigits=8) for i in eachindex(sub_segments)]  #round to help unique function work

    end 

    coords = unique(sub_segments)  
    
    return coords

end

end # module

