using DelimitedFiles, Plots, Images, ImageView
using Wavelets

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 65/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 65/zdata.txt"
for a in 1:4
    open("waveletErrors_65.txt","a") do io
        println(io,"AMPLITUDE "*string(a))
        close(io)
    end
    for f in 1:25
        open("waveletErrors_65.txt","a") do io
            println(io,"FREQUENCY "*string(f))
            close(io)
        end
    
        sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 65/amp" * string(a) * "/freq" * string(f) * "/mapReal.txt"
        sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/angle 65/amp" * string(a) * "/freq" * string(f) * "/mapImag.txt"

        y = readdlm(yFile,',',Float32)
        z = readdlm(zFile,',',Float32)


        slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
        slice[isnan.(slice)] .= 0

        wavt = wavelet(WT.db5)

        ww = dwt(slice,wavt,3)

        for P in [0.01,0.05,0.1]

            open("waveletErrors_65.txt","a") do io
                println(io,"KEEP PROBABILITY "*string(P))
                println(io,"")
                close(io)
            end

            N = Int(floor(P*length(y)*length(z)))

            local wT = threshold(ww[:], BiggestTH(),N)
            wT = reshape(wT,(1000,1000))
            local sliceR = idwt(wT,wavt)
            nError = 100 * norm(sliceR-slice,2)./norm(slice,2) # normalized error
            
            open("waveletErrors_65.txt","a") do io
                println(io,"ERROR: "*string(nError))
                println(io,"")
                println(io,"")
                close(io)
            end
        end
    end
end
