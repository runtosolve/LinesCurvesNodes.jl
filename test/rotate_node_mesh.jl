using LinesCurvesNodes

Z = 0.0:1.0:12.0
X = zeros(Float64, length(Z))
Y = zeros(Float64, length(Z))

nodes = [X Y Z]

nodes = LinesCurvesNodes.rotate_nodes(nodes, rotation_axis="y", rotation_center=[0.0, 0.0, 0.0], θ=π/2)