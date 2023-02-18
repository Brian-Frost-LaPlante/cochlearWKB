# WKB 2-D model
# Parameters from Steele and Taber, 1979

using Plots
using DSP
using Wavelets
using Random

j = im # I must use j!

rho = 1e-6; h = 1; L = 35; # density, height of scalae, length of cochlea
M0 = 1.5e-6; R0 = 2e-3; S0 = 10e3; # mass, resistance and stiffness at base
M1 = 0; R1 = 0; S1 = -0.2; # space constants of exponential parameter decay

X = 5000 # number of points in longitudinal linspace
x = LinRange(0,L,X) # longitudinal linspace
M = M0*exp.(M1*x); R = R0*exp.(R1*x); S = S0*exp.(S1*x); # parameters across x

f = [800, 1600, 3200,6400] # frequencies in Hz
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

#p1 = plot(x,abs.(V),label = ["800 Hz" "1.6 kHz" "3.2 kHz" "6.4 kHz"])
#plot!(yscale=:log10)
#ylims!(1e-1,3e2)
#p2 = plot(x,unwrap(angle.(V),dims=1)/(2*pi),label = ["800 Hz" "1.6 kHz" "3.2 kHz" "6.4 kHz"])
#plot(p1,p2,layout=(2,1))


# using db5 wavelets
wt = wavelet(WT.db5)

#v1 = dwt(V[:,1],wt); v2 = dwt(V[:,2],wt)
#v3 = dwt(V[:,3],wt); v4 = dwt(V[:,4],wt)

#P = 0.10 # percentage of kept terms
#N = Int(floor(P*X))

#vT1 = threshold(v1, BiggestTH(),N); vT2 = threshold(v2, BiggestTH(),N); 
#vT3 = threshold(v3, BiggestTH(),N); vT4 = threshold(v4, BiggestTH(),N);

#VT = hcat(idwt(vT1,wt),idwt(vT2,wt),idwt(vT3,wt),idwt(vT4,wt))

#plot(x,abs.(VT))
#plot!(yscale=:log10)
#ylims!(1e-1,3e2)
#title!("Thresholded Responses")


P = 0.5 # probability of keeping a point
mask = rand(Float64,(10,))

mask = mask.<P

mask
