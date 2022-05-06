using LinesCurvesNodes

coords = [[0.0, 0.0], [3.6, 4.2], [7.3, 2.5]]

num_nodes = size(coords)[1]
num_segments = num_nodes - 1

sub_segments = Array{Vector{Vector{Float64}}}(undef, 0)

using LinearAlgebra

for i = 1:num_segments

        A = coords[i]
        B = coords[i+1]

        AB = B-A

        AB_unit = AB/norm(AB)

        length_AB = norm(AB)

        range_AB = range(0.0, length_AB, num_sub_segments+1)

        sub_segments = [sub_segments; [A .+ AB_unit .* range_AB[j] for j in eachindex(range_AB)]] 

        

end 

coords = unique(sub_segments)

# coords = unique(round.(sub_segments, sigdigits=8)) 



num_sub_segments = 2
coords = LinesCurvesNodes.subdivide_line_segments(coords, num_sub_segments)