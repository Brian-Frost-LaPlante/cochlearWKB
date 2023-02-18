function [k, P, V] = WKB_stablePoint(x,omega,Z,rho,h,numSteps,ampOrder)
    % WKB_stablePoint gives the wavenumber, pressure and velocity at the BM using the
    % 2-D WKB stable-point solution method for k. 
    % INPUTS
    % x:        a 1 x X vector containing longitudinal coordinates on the BM
    % omega:    a 1 x F vector containing radian frequency coordinates 
    % Z:        a 1 x X vector containing the complex impedance of the BM
    % rho:      a real number, the density of fluid
    % h:        a real number, the height of the chamber
    % numSteps: an integer, number of iterations used in stable point
    % ampOrder: either 0 for sqrt(k/k0) or 1 for hyperbolic garbage
    % root-finding algorithm
    % OUTPUTS
    % k:        a X x F vector containing the complex wavenumber as a
    % function of both space and frequency
    % P:        a X x F vector containing the complex pressure at the BM as a
    % function of both space and frequency
    % V:        a X x F vector containing the complex velocity at the BM as a
    % function of both space and frequency

    k_LW = sqrt(-1j*2*omega.'*rho./Z);
    k = k_LW; % start with long-wave approx
    X = length(x); L = x(X);


    for step = 1:numSteps
        alpha_WKB=k.*h./tanh(k.*h);
        k=(k_LW).*sqrt(alpha_WKB);
        k=abs(real(k))-1j.*abs(imag(k));%force to pick solution of wave going from base to apex;
    end
    
    k0 = k(:,1);
    
    int_k = cumtrapz(k,2)*L/X;
    
    if ampOrder == 1
        A = (k0*h./tanh(k0*h)).* ...
            sqrt((k0*h.*sech(k0*h).^2 + tanh(k0*h))./...
            (k*h.*sech(k*h).^2 + tanh(k*h)));
    elseif ampOrder == 0
        A = sqrt(k0./k);
    end
    
    P =  exp(-1j*int_k).*A;
    V = P./Z;

end

