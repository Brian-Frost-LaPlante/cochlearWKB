using DelimitedFiles, Plots, Images, ImageView
using Wavelets

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/ydata.txt"
sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp4/freq21/mapReal.txt"
sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp4/freq21/mapImag.txt"

y = readdlm(yFile,',',Float32)
z = readdlm(zFile,',',Float32)


slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0

wavt = wavelet(WT.db5)

ww = dwt(slice,wavt)

P = 0.5
N = Int(floor(P*length(y)*length(z)))

wT = threshold(ww[:], BiggestTH(),N)
wT = reshape(wT,(1000,1000))
sliceR = idwt(wT,wavt)
sum(sum(abs.(sliceR-slice)))

# P = 0.5 -> SE 0.50
# P = 0.1 -> SE 16.9
# P= 0.01 -> SE 9180.5
