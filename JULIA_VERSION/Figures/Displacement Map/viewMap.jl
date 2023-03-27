using DelimitedFiles, Images, ImageView

j = im

sliceFileR = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp2/freq20/mapReal.txt"
sliceFileI = "/home/brian/cochlearWKB/interpolatedAreaCSVs/transverse/amp2/freq20/mapImag.txt"


slice = readdlm(sliceFileR,',',Float32) + j*readdlm(sliceFileI,',',Float32)
slice[isnan.(slice)] .= 0

imshow(abs.(slice))

