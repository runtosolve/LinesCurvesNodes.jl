
using LinesCurvesNodes

A = [0.0, 0.0]
B = [0.504, -0.720]
C = [2.914, -0.720]


r=0.118
n = 3

fillet = LinesCurvesNodes.generate_fillet(A, B, C, r, n)


answer = [ 0.468768  -0.669669
0.494964  -0.696644
0.528316  -0.714009
0.565437  -0.72]


r=0.177
A = [2.9146, -0.72]
B = [2.9146, 7.22]
C = [5.4146, 7.22]

fillet = LinesCurvesNodes.generate_fillet(A, B, C, r, n)