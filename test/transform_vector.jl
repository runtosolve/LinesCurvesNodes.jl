using LinesCurvesNodes, Unitful 

#with units
L = 23.0u"inch"
start_node = [0.0u"inch", 0.0u"inch"]
θ = π/4

LinesCurvesNodes.transform_vector(L, start_node, θ)


#without units 
L = 23.0
start_node = [0.0, 0.0]
θ = π/4

LinesCurvesNodes.transform_vector(L, start_node, θ)
