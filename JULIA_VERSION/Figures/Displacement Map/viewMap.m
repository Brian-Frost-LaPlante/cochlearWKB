fR1 = "~/cochlearWKB/interpolatedAreaCSVs/transverse/amp3/freq10/mapReal.txt";
fI1 = "~/cochlearWKB/interpolatedAreaCSVs/transverse/amp3/freq10/mapImag.txt";

R1 = load(fR1); I1 = load(fI1);

fR2 = "~/cochlearWKB/interpolatedAreaCSVs/transverse/amp4/freq10/mapReal.txt";
fI2 = "~/cochlearWKB/interpolatedAreaCSVs/transverse/amp4/freq10/mapImag.txt";

R2 = load(fR2); I2 = load(fI2);
%%
cplxImg1 = R1 + 1j*I1;
cplxImg2 = R2 + 1j*I2;

%max(max(abs(cplxImg2)))
cplxImg2(1,1) = max(max(abs(cplxImg1)));
cplxImg2(1:151,1:491) = 0;
cplxImg2(1:346, 1:176) = 0;
cplxImg1(737:end, 1:49) = 0;
%%
figure;
subplot(2,2,1)
imagesc(abs(cplxImg1))
colormap("parula")
colorbar
subplot(2,2,2)
imagesc(abs(cplxImg2))
colormap("parula")
colorbar

subplot(2,2,3)
imagesc(angle(cplxImg1)/(2*pi))
colormap("hot")
colorbar
subplot(2,2,4)
imagesc(angle(cplxImg2)/(2*pi))
colormap("hot")
colorbar