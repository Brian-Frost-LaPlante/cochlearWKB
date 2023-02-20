using Random, LinearAlgebra

M = 1000; N = 1000; # width by height of full image, 2-D case

P = 0.5 # prob of keeping

K = Int(floor(N*P)) # number of A-Scans to pick out, downsampled in KxM
A = Matrix(1.0I,M,M)

rowsToKeep = sort(randperm(N)[1:K])
A = A[rowsToKeep,:]
A
