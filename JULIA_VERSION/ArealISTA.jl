using DelimitedFiles, Plots#, Images, ImageView
using DSP, Wavelets, Random, LinearAlgebra

lambda = 5e4 # l1 weight
L = .5 # 1/step size
numSteps = 10000
numAvg = 10

angle = "transverse"

logfile = "ArealLogs_db5_"*angle*".txt"
fp = open(logfile,"w")

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/"*angle*"/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/"*angle*"/zdata.txt"
y = readdlm(yFile,',',Float32)
z = readdlm(zFile,',',Float32)
wt = wavelet(WT.db5)
wSlice = dwt(slice,wt,3)

P = 0.25

M = length(y)
N = length(z)

K = Int(floor(P*M))  # number of A-Scans to keep



# for now I only care about one frequency and one amplitude

NMSE = 0

for aa = 1:4
    write(fp,"AMP "*string(aa)*"\n")
    println(fp,"AMP "*string(aa)*"\n")
for ff = 1:25 #5:5:25 
    global sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/"*angle*"/amp"*string(aa)*"/freq"*string(ff)*"/mapReal.txt"
    global sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/"*angle*"/amp"*string(aa)*"/freq"*string(ff)*"/mapImag.txt"
    global slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
    slice[isnan.(slice)] .= 0
    println("FREQ "*string(ff)*"\n")
    write(fp,"FREQ "*string(ff)*"\n")

    for nn = 1:numAvg


        #y = y[1:2:end]; z = z[1:2:end]; slice = slice[1:2:end,1:2:end];

        A = Matrix(1.0I,M,M) # identity map for the MxN image
        rr = randperm(N)
        rowsToKeep = sort(rr[1:K])
        A = A[rowsToKeep,:]  # sensing matrix, deletes some A-Scans

        nZero = norm(slice[:,sort(rr[K+1:end])],2)^2/norm(slice,2)^2
        println("ZEROED NMSE: "*string(100*nZero))
        write(fp,"ZEROED NMSE: "*string(100*nZero)*"\n")

        uSlice = A*slice     # this is the measured undersampled image, K x N
        global wGuess = zeros(M,N)  # initialized guess in the wavelet domain

        obj = zeros(numSteps,1)
        MSE = zeros(numSteps,1)

        stopEarly = false

        for step = 1:numSteps
            if ~stopEarly
                global errors
                global NMSE

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
                        NMSE = NMSE + MSE[step-5]
                        stopEarly = true
                    end
                elseif step == numSteps
                    NMSE = NMSE + MSE[step]
                    write(fp,"STEP " * string(step) * "\nMSE: " * string(MSE[step])*"\n")
                end
            end
        end
        println("CURRENT AVERAGE: "*string(NMSE/nn))
        write(fp,"\n")
    end
    write(fp,"AVERAGED NMSE " * string(NMSE/numAvg)*"\n")
    write(fp,"\n")
end
    write(fp,"\n")
end
close(fp)
#p1 = plot(1:numSteps,obj)
#title!("Objective Functions")

#p2 = plot(1:numSteps,MSE)
#title!("MSE")

#plot(p1,p2,layout=(2,1))
#imshow(abs.(wGuess))
