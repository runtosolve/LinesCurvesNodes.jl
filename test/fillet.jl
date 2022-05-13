
using LinesCurvesNodes

A = [0.0, 0.0]
B = [0.504, -0.720]
C = [2.914, -0.720]


r=0.118
n = 3

fillet = LinesCurvesNodes.generate_fillet(A, B, C, r, n)

convert(Vector{Vector{Float64}}, fillet)

answer =  [[0.4687679896444109, -0.6696685566348728]
[0.49496385452130676, -0.6966439937443945]
[0.528315875030654, -0.7140089087122952]
[0.5654373562656577, -0.7200000000000001]]


r=0.177
A = [2.9146, -0.72]
B = [2.9146, 7.22]
C = [5.4146, 7.22]

fillet = LinesCurvesNodes.generate_fillet(A, B, C, r, n)

answer = [ 2.9146   7.043
2.93831  7.1315
3.0031   7.19629
3.0916   7.22]