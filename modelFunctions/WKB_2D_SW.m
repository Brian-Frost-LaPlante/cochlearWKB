function [k, P, V] = WKB_2D_SW(x,omega,Z,rho,h)
    % WKB_2D_SW gives the wavenumber, pressure and velocity at the BM using the
    % 2-D WKB short-wave approximation. 
    % INPUTS
    % x:      a 1 x X vector containing longitudinal coordinates on the BM
    % omega:  a 1 x F vector containing radian frequency coordinates 
    % Z:      a 1 x X vector containing the complex impedance of the BM
    % rho:    a real number, the density of fluid
    % h:      a real number, the height of the chamber
    % OUTPUTS
    % k:      a X x F vector containing the complex wavenumber as a
    % function of both space and frequency
    % P:      a X x F vector containing the complex pressure at the BM as a
    % function of both space and frequency
    % V:      a X x F vector containing the complex velocity at the BM as a
    % function of both space and frequency

    k = -2*1j*omega.'*rho./(h*Z);
    k0 = k(:,1);
    X = length(x); L = x(X);

    int_k = cumtrapz(k,2)*L/X;

    A = (k0*h./tanh(k0*h)).* ...
        sqrt((k0*h.*sech(k0*h).^2 + tanh(k0*h))./...
        (k*h.*sech(k*h).^2 + tanh(k0*h)));

    P =  exp(-1j*int_k).*A;
    V = P./Z;
    
end

