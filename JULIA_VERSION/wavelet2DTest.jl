using DelimitedFiles, Plots, Images, ImageView
using Wavelets

j = im

yFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/ydata.txt"
zFile = "/home/brian/cochlearWKB/interpolatedAreaCSVs/zdata.txt"

a = 4; f = 15;

sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp" * string(a) * "/freq" * string(f) * "/mapReal.txt"
sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/amp" * string(a) * "/freq" * string(f) * "/mapImag.txt"

y = readdlm(yFile,',',Float32)
z = readdlm(zFile,',',Float32)


slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0

wavt = wavelet(WT.db5)

ww = dwt(slice,wavt)
N = 1000
M = 1000

#wT = reshape(wT,(N,M))
#imshow(abs.(ww[1:Int(floor(N/4)),1:Int(floor(M/4))]))
imshow(abs.(ww))
