using DelimitedFiles, Plots, Images, ImageView
using Wavelets

j = im

sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp4/freq15/mapReal.txt"
sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp4/freq15/mapImag.txt"

slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0
        
wavt = wavelet(WT.haar)

ww = dwt(slice,wavt,3)
ww = ww[1:125,1:125]
#wT = threshold(ww[:], BiggestTH(),10000)
plot(sort(abs.(ww[:])))
