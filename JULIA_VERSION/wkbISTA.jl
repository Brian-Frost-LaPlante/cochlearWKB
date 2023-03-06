# WKB 2-D model
# Parameters from Steele and Taber, 1979

using Plots
using DSP
using Wavelets
using Random
using LinearAlgebra

j = im # I must use j!

rho = 1e-6; h = 1; L = 25; # density, height of scalae, length of cochlea
M0 = 1.5e-6; R0 = 2e-3; S0 = 10e3; # mass, resistance and stiffness at base
M1 = 0; R1 = 0; S1 = -0.2; # space constants of exponential parameter decay

X = 5000 # number of points in longitudinal linspace
x = LinRange(0,L,X) # longitudinal linspace
M = M0*exp.(M1*x); R = R0*exp.(R1*x); S = S0*exp.(S1*x); # parameters across x

f = 2000 # frequencies in Hz
omega = 2*pi*f # radian frequencies

Z = (S ./ (j*omega')) .+ R .+ M .* (j*omega') # impedance
Y = 1 ./ Z # admittance

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

k, P, V = wkbStablePoint(x,omega,Z,rho,h,10,1)

P = 0.5 # prob of keeping

XX = Int(floor(X*P)) # number of A-Scans to pick out, downsampled in KxM

indsToKeep = sort(randperm(X)[1:XX])

A = Matrix(1.0I,X,X)
A = A[indsToKeep,:]

xx = A*x
VV = A*V

# V is the fully sampled vector, VV is the undersampled vector.
# Vp is the guessed vector
# The undersampling operation A is just indexing by indsToKeep
# The objective function is (1/2)*||A*V - VV||_2^2 + lambda*||W[VV]||_2
# Subgradient index i of 2-norm is v_i/||v||_2 if ||v_2 > epsilon, 0 else

Vp = rand(Complex{Float64}, X) 

lambda = 1000.0 
eta = 0.01
epsilon = 10
numSteps  = 1000
errors = zeros(numSteps,1)

wt = wavelet(WT.db5)

for step = 1:numSteps
    
    global errors

    grad1 = transpose(A)*(A*Vp - VV)
    norm1 = 0.5*norm(grad1,2)^2 

    WW = dwt(Vp,wt)
    grad2 = zeros(X,1)
    norm2 = norm(WW,1)
    for i = 1:X
        if abs(WW[i]) >epsilon
            grad2[i] = angle(WW[i])
        end
    end
    #norm2 = norm(WW,2)
    #if norm2 < epsilon
    #    grad2 = 0
    #else
    #    grad2 = WW./norm2
    #end

    errors[step] = norm1 + lambda*norm2
    grad = grad1 .+ lambda.*grad2
    global Vp = Vp - eta*grad
end

p1 = plot(1:numSteps,errors)
p2 = plot(x,20*log10.(abs.(Vp)),linewidth=5)
plot!(x,20*log10.(abs.(V)),linewidth = 1)

ylims!(-20,40)
title!("Thresholded Responses")

plot(p1,p2,layout=(2,1))
