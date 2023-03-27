using DelimitedFiles, Plots, Images, ImageView
using DSP, Wavelets, Random, LinearAlgebra

lambda = 5e3 # l1 weight
L = 0.5 # 1/step size
numSteps = 10000

logfile = "ArealLogs_A3P20l5e3_54.txt"
fp = open(logfile,"w")

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 54/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 54/zdata.txt"

# for now I only care about one frequency and one amplitude

for ff = 5:5:25 

println("FREQ "*string(ff)*"\n")
write(fp,"FREQ "*string(ff)*"\n")

global sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 54/amp3/freq"*string(ff)*"/mapReal.txt"
global sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 54/amp3/freq"*string(ff)*"/mapImag.txt"
y = readdlm(yFile,',',Float32)
z = readdlm(zFile,',',Float32)

global slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0

#y = y[1:2:end]; z = z[1:2:end]; slice = slice[1:2:end,1:2:end];

wt = wavelet(WT.db5)
wSlice = dwt(slice,wt,3)

P = 0.2

M = length(y)
N = length(z)

K = Int(floor(P*M))  # number of A-Scans to keep

A = Matrix(1.0I,M,M) # identity map for the MxN image
rr = randperm(N)
rowsToKeep = sort(rr[1:K])
A = A[rowsToKeep,:]  # sensing matrix, deletes some A-Scans


uSlice = A*slice     # this is the measured undersampled image, K x N
global wGuess = zeros(M,N)  # initialized guess in the wavelet domain

obj = zeros(numSteps,1)
MSE = zeros(numSteps,1)

stopEarly = false

for step = 1:numSteps
    if ~stopEarly
        global errors

        grad1_space = transpose(A)*(A*idwt(wGuess,wt,) - uSlice)
        grad1       = dwt(grad1_space,wt,3)/L

        wGuess      = wGuess - grad1

        for n = 1:N
            for m = 1:M
            if ~(m<125 && n<125)
                    if abs(wGuess[n,m]) < lambda/L
                        wGuess[n,m] = 0
                    else
                        wGuess[n,m] = sign(wGuess[n,m])*(abs(wGuess[n,m]) - lambda/L)
                    end
            end
            end
        end

        #norm1 = 0.5*norm(A*idwt(wGuess,wt,3) - uSlice,2)^2
        #norm2 = norm(wGuess,1)
        #obj[step] = norm1 + lambda*norm2
        MSE[step] = 100*norm(wGuess-wSlice,2).^2/norm(wSlice,2)^2
        
        if (step % 10) == 0 
            println("STEP " * string(step))
            println("MSE: " * string(MSE[step]))
            println("")
            if MSE[step] > MSE[step-5]
                stopEarly = true
            end
        write(fp,"STEP " * string(step) * "\nMSE: " * string(MSE[step])*"\n")
        end
    end
end

write(fp,"\n")
end
close(fp)
#p1 = plot(1:numSteps,obj)
#title!("Objective Functions")

#p2 = plot(1:numSteps,MSE)
#title!("MSE")

#plot(p1,p2,layout=(2,1))
imshow(abs.(wGuess))
