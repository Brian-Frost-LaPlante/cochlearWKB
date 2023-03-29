using DelimitedFiles, Plots, Images, ImageView
using Wavelets

j = im

global saveFile = "logFiles/sparsityErrors/transverse.txt"

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/zdata.txt"
for a in 1:4
    open(saveFile,"w") do io
        println(io,"AMPLITUDE "*string(a))
        close(io)
    end
        
    println("AMPLITUDE "*string(a))
    
    for f in 1:25
        open(saveFile,"a") do io
            println(io,"FREQUENCY "*string(f))
            close(io)
        end
        
        println("FREQUENCY "*string(f))
    
        sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp" * string(a) * "/freq" * string(f) * "/mapReal.txt"
        sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp" * string(a) * "/freq" * string(f) * "/mapImag.txt"

        y = readdlm(yFile,',',Float32)
        z = readdlm(zFile,',',Float32)
        slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
        slice[isnan.(slice)] .= 0

        for kSize = 1:20

            if kSize == 1
                wavt = wavelet(WT.db1)
            elseif kSize == 2
                wavt = wavelet(WT.db2)
            elseif kSize == 3
                wavt = wavelet(WT.db3)
            elseif kSize == 4
                wavt = wavelet(WT.db4)
            elseif kSize == 5
                wavt = wavelet(WT.db5)
            elseif kSize == 6
                wavt = wavelet(WT.db6)
            elseif kSize == 7
                wavt = wavelet(WT.db7)
            elseif kSize == 8
                wavt = wavelet(WT.db8)
            elseif kSize == 9
                wavt = wavelet(WT.db9)
            elseif kSize == 10
                wavt = wavelet(WT.db10)
            elseif kSize == 11
                wavt = wavelet(WT.Daubechies{11}())
            elseif kSize == 12
                wavt = wavelet(WT.Daubechies{12}())
            elseif kSize == 13
                wavt = wavelet(WT.Daubechies{13}())
            elseif kSize == 14
                wavt = wavelet(WT.Daubechies{14}())
            elseif kSize == 15
                wavt = wavelet(WT.Daubechies{15}())
            elseif kSize == 16
                wavt = wavelet(WT.Daubechies{16}())
            elseif kSize == 17
                wavt = wavelet(WT.Daubechies{17}())
            elseif kSize == 18
                wavt = wavelet(WT.Daubechies{18}())
            elseif kSize == 19
                wavt = wavelet(WT.Daubechies{19}())
            elseif kSize == 20
                wavt = wavelet(WT.Daubechies{20}())
            end

            println("KERNEL SIZE "*string(kSize))
            
            ww = dwt(slice,wavt,3)

            keepSearching = true

            P = 0.0001
            while keepSearching

                open(saveFile,"a") do io
                    println(io,"KEEP PROBABILITY "*string(P))
                    println(io,"")
                    close(io)
                end
                
                println("KEEP PROBABILITY "*string(P))
                
                N = Int(floor(P*length(y)*length(z)))

                local wT = threshold(ww[:], BiggestTH(),N)
                wT = reshape(wT,(1000,1000))
                local sliceR = idwt(wT,wavt)
                nError = norm(sliceR-slice,2).^2 /norm(slice,2).^2 # normalized error
                
                open(saveFile,"a") do io
                    println(io,"ERROR: "*string(nError))
                    println(io,"")
                    println(io,"")
                    close(io)
                end
                   
                println("ERROR: "*string(nError))
                if nError < 1
                    keepSearching = false
                else
                    P = P + 0.0001
                end
            end
        end
    end
end
