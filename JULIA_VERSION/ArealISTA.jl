using DelimitedFiles, Plots, Images, ImageView
using DSP, Wavelets, Random, LinearAlgebra

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/zdata.txt"

# for now I only care about one frequency and one amplitude

sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp4/freq19/mapReal.txt"
sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp4/freq19/mapImag.txt"
y = readdlm(yFile,',',Float32)
z = readdlm(zFile,',',Float32)

slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0

wt = wavelet(WT.db5)
wSlice = dwt(slice,wt)

P = 0.5

M = length(y)
N = length(z)

K = Int(floor(P*M))  # number of A-Scans to keep

A = Matrix(1.0I,M,M) # identity map for the MxN image
rowsToKeep = sort(randperm(N)[1:K])
A = A[rowsToKeep,:]  # sensing matrix, deletes some A-Scans

uSlice = A*slice     # this is the measured undersampled image, K x N
wGuess = zeros(M,N)  # initialized guess in the wavelet domain

lambda = .1 # l1 weight
L = 1 # 1/step size
numSteps = 100000
obj = zeros(numSteps,1)
MSE = zeros(numSteps,1)

for step = 1:numSteps
    global errors
    global wGuess

    grad1_space = transpose(A)*(A*idwt(wGuess,wt) - uSlice)
    grad1       = dwt(grad1_space,wt)/L

    wGuess      = wGuess - grad1

    for n = 1:N
        for m = 1:M
            if abs(wGuess[n,m]) < lambda/L
                wGuess[n,m] = 0
            end
        end
    end

    norm1 = 0.5*norm(A*idwt(wGuess,wt) - uSlice,2)^2
    norm2 = norm(wGuess,1)
    obj[step] = norm1 + lambda*norm2
    MSE[step] = norm(wGuess-wSlice,2)/norm(wSlice,2)

    if (step % 100) == 0 
        println("STEP " * string(step))
        println("MSE: " * string(MSE[step]))
        println("")
    end
end

p1 = plot(1:numSteps,obj)
title!("Objective Functions")

p2 = plot(1:numSteps,MSE)
title!("MSE")

plot(p1,p2,layout=(2,1))
