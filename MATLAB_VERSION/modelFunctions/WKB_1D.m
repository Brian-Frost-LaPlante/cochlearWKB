function [k, P, V] = WKB_1D(x,omega,Z,rho,h, order)
    % WKB_1D gives the wavenumber, pressure and velocity at the BM using the
    % 1-D WKB approximation. 
    % INPUTS
    % x:      a 1 x X vector containing longitudinal coordinates on the BM
    % omega:  a 1 x F vector containing radian frequency coordinates 
    % Z:      a 1 x X vector containing the complex impedance of the BM
    % rho:    a real number, the density of fluid
    % h:      a real number, the height of the chamber
    % order:  either 0, 1 or 2. Determines the order of WKB approximation
    % used to determine the pressure and velocity
    % OUTPUTS
    % k:      a X x F vector containing the complex wavenumber as a
    % function of both space and frequency
    % P:      a X x F vector containing the complex pressure at the BM as a
    % function of both space and frequency
    % V:      a X x F vector containing the complex velocity at the BM as a
    % function of both space and frequency

    k = sqrt(-1j*2*omega.'*rho./(h*Z));
    k0 = k(:,1);
    X = length(x); L = x(X);
    int_k = cumtrapz(k,2)*L/X;

    if order == 0
        P =  2*rho*omega.'.*exp(-1j*x.*k)./k0;
        V = P./Z;
    elseif order == 1
        P =  2*rho*omega.'.*exp(-1j*int_k)./k0(:,1);
        V = P./Z;
    elseif order == 2
        P =  2*rho*omega.'.*exp(-1j*int_k)./sqrt(k0.*k);
        V = P./Z;
    end
end

