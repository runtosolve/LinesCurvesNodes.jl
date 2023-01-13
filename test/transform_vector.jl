using LinesCurvesNodes, Unitful 


L = 23.0u"inch"
start_node = [0.0u"inch", 0.0u"inch"]
θ = π/4

LinesCurvesNodes.transform_vector(L, start_node, Θ)
