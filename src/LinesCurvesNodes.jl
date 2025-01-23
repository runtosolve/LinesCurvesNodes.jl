module LinesCurvesNodes

using LinearAlgebra, Rotations, StaticArrays


mutable struct Line

    point_start::Vector{Float64}
    point_end::Vector{Float64}
    
end

mutable struct Arc

    center::Vector{Float64}
    radius::Float64
    angle_start::Float64
    angle_end::Float64
    
end



function transform_vector(L, start_node, Θ)

    vector = [L; 0.0*unit(L[1])] 

    #rotate
    vector = Angle2d(Θ) * vector

    #shift
    vector = vector + start_node

    #convert mutable form of StaticArray
    vector = MVector(vector)

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

    fillet = convert(Vector{Vector{Any}}, [Angle2d(dΔ[i]) * Ra + R for i in eachindex(dΔ)])  #Angle2d is 2D rotation matrix 

    # fillet = permutedims(reshape(hcat(fillet...), (length(fillet[1]), length(fillet))))

    return fillet

end


function discretize_vector(A, B, n)

    AB = B - A
    AB_unit = AB / norm(AB)
    spacing = range(0.0 * unit(A[1]), norm(AB), n+1)
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


function find_nodes(nodes, xloc, yloc, zloc)

   
    x_index = findall(x-> (x >= xloc[1]) & (x <= xloc[2]), nodes[:,1])
    y_index = findall(y-> (y >= yloc[1]) & (y <= yloc[2]), nodes[:,2])
    z_index = findall(z-> (z >= zloc[1]) & (z <= zloc[2]), nodes[:,3])

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

            range_AB = range(0.0 * unit(A[1]), length_AB, num_sub_segments+1)

            sub_segments = [sub_segments; [A .+ AB_unit .* range_AB[j] for j in eachindex(range_AB)]] 

            sub_segments = [round.(sub_segments[i], sigdigits=8) for i in eachindex(sub_segments)]  #round to help unique function work

    end 

    coords = unique(sub_segments)  
    
    return coords

end


function define_rotation_matrix(A, B)


    AB = B - A

    ΔX = AB[1]
    ΔY = AB[2]

    χ = atan(ΔY, ΔX)


    γ = [cos(χ) sin(χ) 0
         -sin(χ) cos(χ) 0
         0  0 1]

    Γ = zeros(Float64, (6, 6))

    Γ[1:3, 1:3] .= γ
    Γ[4:6, 4:6] .= γ

    return Γ

end

#from InstantFrame 
function beam_shape_function(q1, q2, q3, q4, L, x)
    
    a0 = q1
    a1 = q2
    a2 = 1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3 = 1/L^3*(2*q1+q2*L-2*q3+q4*L)
    w = a0 + a1*x + a2*x^2 + a3*x^3

    return w

end

#from InstantFrame
function discretized_element_global_coords(node_i_coords, Γ, x)

    local_element_discretized_coords = [zeros(Float64, 3) for i in eachindex(x)]

    [local_element_discretized_coords[i][1] = x[i] for i in eachindex(x)]

    global_element_discretized_coords = [Γ'[1:3, 1:3] * (local_element_discretized_coords[i]) .+  node_i_coords for i in eachindex(x)]

    return global_element_discretized_coords

end

function define_circular_segment(center, radius, θ_start, θ_end, n)

    x(θ) = radius * cos(θ)
    y(θ) = radius * sin(θ)

    θc = range(θ_start, θ_end, length = n-1)
 

    xc = center[1] .+ [x(θc[i]) for i in eachindex(θc)]
    yc = center[2] .+ [y(θc[i]) for i in eachindex(θc)]

    coords = [[xc[i], yc[i]] for i in eachindex(xc)]

    return coords

end



function generate_xy_coordinates_from_dxf_JSON(data, n, n_r, reverse_arc_angles)

    #read LINES and ARCS
    cross_section = Vector{Union{Line, Arc}}(undef, size(data)[1])

    for i in eachindex(data)

        if data[i].type == "LINE"

            cross_section[i] = Line([data[i].start.x, data[i].start.y], [data[i].end.x, data[i].end.y])

        elseif data[i].type == "ARC"

            cross_section[i] = Arc([data[i].center.x, data[i].center.y], data[i].radius, data[i].start_angle, data[i].end_angle)

        end

    end


    #reverse start and end angles

    # if reverse_arc_angles == true
        # j = 1
        # for i in eachindex(data)

        #     if data[i].type == "ARC"

        #         if reverse_arc_angles[j] == true

        #                 cross_section[i].angle_start = data[i].end_angle 
        #                 cross_section[i].angle_end = data[i].start_angle 

        #         end
            
        #         j += 1
                
        #     end

        # end

        #     

        # end

    # end

    # end

    xy_all = Vector{Vector{Float64}}(undef, 0)

    θ_start_all = Vector{Float64}(undef, 0)
    θ_end_all = Vector{Float64}(undef, 0)

    xy_segments = Vector{Vector{Vector{Float64}}}(undef, 0)

    for i in eachindex(cross_section)

        if typeof(cross_section[i]) == Line

            if i == 1
            
                xy_segment = LinesCurvesNodes.subdivide_line_segments([cross_section[i].point_start, cross_section[i].point_end], n)

            else

                xy_segment = LinesCurvesNodes.subdivide_line_segments([cross_section[i].point_start, cross_section[i].point_end], n)[2:end]

            end


            xy_all = vcat(xy_all, xy_segment)

            xy_segments = [xy_segments; [xy_segment]]

        elseif typeof(cross_section[i]) == Arc

            center = cross_section[i].center
            radius = cross_section[i].radius

            θ_start = deg2rad(cross_section[i].angle_start) 

            θ_end = deg2rad(cross_section[i].angle_end) 

            # if θ_start >= 3π/2

            #     if (θ_end > 0) & (θ_end < π/2)

            #         θ_end =  θ_end + 2π

            #     elseif isapprox(θ_start, 2π) 

            #         θ_start =  0.0

            #     end

            # end

            if θ_end < θ_start

                if θ_end >= 0.0

                    θ_end = θ_end + 2π

                end

            end

            push!(θ_start_all, θ_start)
            push!(θ_end_all, θ_end)

            xy_arc = LinesCurvesNodes.define_circular_segment(center, radius, θ_start, θ_end, n_r)

            if reverse_arc_angles[i] == true

                xy_arc = reverse(xy_arc)[2:end]

            else

                xy_arc = xy_arc[2:end]

            end

            xy_all = vcat(xy_all, xy_arc)

            xy_segments = [xy_segments; [xy_arc]]

        end

    end

    return xy_all, xy_segments,  θ_start_all,  θ_end_all 

end

#for cross-sections with varying thickness, find the array range for each line segment that has the same thickness, used in CUFSM.Show
function find_linesegments(linewidths)

    index_start = 1

    linesegment_ranges = Vector{Vector{Int}}(undef, 0)
    linewidth_segments = Vector{Float64}(undef, 0)
    for i in eachindex(linewidths[1:end-1])

        if linewidths[i] != linewidths[i+1]

            index_end = i
            push!(linesegment_ranges, [index_start, index_end])
            push!(linewidth_segments, linewidths[i])
            index_start = i+1

        end

        if i == (length(linewidths) - 1)
            index_end = i + 1
            push!(linesegment_ranges, [index_start, index_end])
            push!(linewidth_segments, linewidths[i+1])
        end

    end

    return linesegment_ranges, linewidth_segments

end

#for cross-sections with varying thickness, used in CUFSM.Show
function combine_points_into_linesegments(linesegment_ranges, x, y)

    linesegments = Vector{Vector{Vector{Float64}}}(undef, size(linesegment_ranges)[1])
    # linesegments = Vector{Vector{Vector{Float64}}}(undef, size(linesegment_ranges)[1])

    for i in eachindex(linesegment_ranges)

        # linesegments[i] = cross_section_coords[linesegment_ranges[i][1]:linesegment_ranges[i][2]]
      
        nodes = Vector{Float64}[]
        for j = linesegment_ranges[i][1]:linesegment_ranges[i][2]

            push!(nodes, [x[j], y[j]])

            if (j+1) < length(x)  #stops 1 element short for a closed section, taken care of in Show.section to plot last closed section element
                push!(nodes, [x[j+1], y[j+1]])
            end

            # if (j+1) < length(x)  #stops 1 element short for a closed section, taken care of in Show.section to plot last closed section element
            # push!(nodes, [x[j+1], y[j+1]])
            # end

        end

   
        if (i == size(linesegment_ranges)[1]) & (length(x) > size(linesegment_ranges)[1])   #get last node in an open cross-section, for section plots
            push!(nodes, [x[end], y[end]])
        end
        

        nodes = unique(nodes)

        linesegments[i] = nodes


    end

    return linesegments

end


function find_elements_containing_nodes(node_set, elements)

    element_set = Int[]
    for i in eachindex(node_set)

        element_index = findfirst(node_number -> node_number==node_set[i], elements[:, 2:end])
        if !isnothing(element_index)
            push!(element_set, elements[element_index[1], 1])
        end

    end

    return unique(element_set)

end



end # module

