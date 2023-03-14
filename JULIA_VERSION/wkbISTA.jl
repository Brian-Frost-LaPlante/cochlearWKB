# WKB 2-D model
# Parameters from Steele and Taber, 1979

using Plots
using DSP
using Wavelets
using Random
using LinearAlgebra

j = im # I must use j!

rho = 1e-6; h = 1; L = 35; # density, height of scalae, length of cochlea
M0 = 1.5e-6; R0 = 2e-3; S0 = 10e3; # mass, resistance and stiffness at base
M1 = 0; R1 = 0; S1 = -0.2; # space constants of exponential parameter decay

X = 1000 # number of points in longitudinal linspace
x = LinRange(0,L,X) # longitudinal linspace
M = M0*exp.(M1*x); R = R0*exp.(R1*x); S = S0*exp.(S1*x); # parameters across x

# For this test I decided to use the Alessandro method that I do not particularly like. It's easy!
function wkbStablePoint(x,omega,Z,rho,h,numSteps,ampOrder)
    j = im
    k_LW = sqrt.(-j*2*omega'*rho ./ Z)
    X = length(x); L = x[X];


    k = k_LW

        #for step = 1:numSteps
    #    alpha = k .*h ./ tanh.(k.*h)
    #    k = k_LW .* sqrt.(alpha)
    #    k = abs.(real(k)) - j*abs.(imag(k))
    #end

    k0 = k[1,:]'

    int_k = cumsum(k,dims = 1)*L/X

    if ampOrder == 1
        A = (k0*h ./ tanh.(k0*h)) .* sqrt.((k0*h .* sech.(k0*h).^2 + tanh.(k0*h)) ./ (k*h .* sech.(k*h).^2 + tanh.(k*h)))
    elseif ampOrder == 0
        A = sqrt.(k0 ./ k)
    end

    P = exp.(-j*int_k) .* A
    V = P ./ Z

    return int_k, P, V

end
# V is the fully sampled vector in space, VV is the undersampled vector.
# W is the fully sampled vector in wavelet domain, WW is the undersampled vector in wavelet domain
# Vp is the guessed vector
# Wp is guessed vector in wavelet domain
# The objective function is (1/2)*||M*Wp - WW||_2^2 + lambda*||Wp||_2
# Subgradient index i of 2-norm is v_i/||v||_2 if ||v_2 > epsilon, 0 else


lambda = 10
L = 1
numSteps  = 100000
errors = zeros(numSteps,1)

wt = wavelet(WT.db5)


for freq = [300,400,500,600,700,800]
	println("FREQ: "*string(freq))
	omega = 2*pi*freq # radian frequencies

	Z = (S ./ (j*omega')) .+ R .+ M .* (j*omega') # impedance
	Y = 1 ./ Z # admittance

	k, P, V = wkbStablePoint(x,omega,Z,rho,h,10,1)

	global V

	P = .33 # prob of keeping

	XX = Int(floor(X*P)) # number of A-Scans to pick out, downsampled in KxM

	indsToKeep = sort(randperm(X)[1:XX])

	global A = Matrix(1.0I,X,X)
	A = A[indsToKeep,:]


	global xx = A*x
	global VV = A*V

	global Wp = rand(Complex{Float64}, X) 

	global VW = dwt(V[:,1],wt,3)

	for step = 1:numSteps
	    grad1_space = transpose(A)*(A*idwt(Wp[:,1],wt,3) - VV)
	    grad1 = dwt(grad1_space[:,1],wt,3)/L

	    Wp = Wp-grad1

	    for i = 125:X
		if abs(Wp[i]) < lambda/L
		    Wp[i] = 0
		end
	    end

	    norm1 = 0.5 * norm(A*idwt(Wp[:,1],wt,3) - VV,2)^2
	    norm2 = norm(Wp,1)
	    
	    
	    #errors[step] = norm1 + lambda*norm2

	    errors[step] = norm(VW-Wp,2)/norm(VW,2)

	    if step%10000==0
		println(string(errors[step]))
	    end

	end
end

Vp = idwt(Wp[:,1],wt,3)

p1 = plot(1:numSteps,errors)
p2 = plot(x,20*log10.(abs.(Vp)),linewidth=5)
plot!(x,20*log10.(abs.(V)),linewidth = 1)

ylims!(-20,60)
title!("Thresholded Responses")

plot(p1,p2,layout=(2,1))
