function [k, P, V] = WKB_2D_LW(x,omega,Z,rho)
    % WKB_2D_LW gives the wavenumber, pressure and velocity at the BM using the
    % 2-D WKB long-wave approximation. 
    % INPUTS
    % x:      a 1 x X vector containing longitudinal coordinates on the BM
    % omega:  a 1 x F vector containing radian frequency coordinates 
    % Z:      a 1 x X vector containing the complex impedance of the BM
    % rho:    a real number, the density of fluid
    % OUTPUTS
    % k:      a X x F vector containing the complex wavenumber as a
    % function of both space and frequency
    % P:      a X x F vector containing the complex pressure at the BM as a
    % function of both space and frequency
    % V:      a X x F vector containing the complex velocity at the BM as a
    % function of both space and frequency

    k = sqrt(-1j*2*omega.'*rho./Z);
    k0 = k(:,1);
    X = length(x); L = x(X);

    int_k = cumtrapz(k,2)*L/X;
    
    P =  exp(-1j*int_k).*sqrt(k0./k);
    V = P./Z;
end

