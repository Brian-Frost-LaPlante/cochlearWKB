function [k, P, V] = WKB_walkingRF(x,omega,Z,rho,h,numSteps,epsilon)
    % WKB_walkingRF gives the wavenumber, pressure and velocity at the BM using the
    % 2-D WKB stable-point solution method for k. 
    % INPUTS
    % x:        a 1 x X vector containing longitudinal coordinates on the BM
    % omega:    a 1 x F vector containing radian frequency coordinates 
    % Z:        a 1 x X vector containing the complex impedance of the BM
    % rho:      a real number, the density of fluid
    % h:        a real number, the height of the chamber
    % numSteps: an integer, number of iterations used in Newton-Raphson
    % epsilon:  The error threshold above which we switch our initial value
    % OUTPUTS
    % k:        a X x F vector containing the complex wavenumber as a
    % function of both space and frequency
    % P:        a X x F vector containing the complex pressure at the BM as a
    % function of both space and frequency
    % V:        a X x F vector containing the complex velocity at the BM as a
    % function of both space and frequency


    k_LW = sqrt(-1j*2*omega.'*rho./Z); % start with long-wave approx at base
    X = length(x); F = length(omega); L = x(X);


    C = -2*1j*omega.'*rho./Z; % RHS of k relation, so ktanh(kh) - C = 0
    
    k = nan(F,X);
    
    kx = k_LW(:,1); % this is k at x = 0, where long-wave holds well;
    
    for i = 1:X
        startingSW = false;
        for ff = 1:F
            if ~startingSW
                for n = 1:numSteps
                    f = kx(ff).*tanh(kx(ff)*h) - C(ff,i);
                    df = tanh(kx(ff)*h) + kx(ff)*h.*sech(kx(ff)*h).^2; % the derivative in k
                    kx(ff) = kx(ff) - f./df; % the Newton-Raphson update
                    kx(ff) = abs(real(kx(ff))) -1j*abs(imag(kx(ff)));
                end
                if abs(f) <= epsilon || isnan(f)
                    k(ff,i) = kx(ff);
                else
    
                    startingSW = true;
                    Ri = real(Z(1,i)); % resistance is real part of Z
                    gamma = Ri./(2*rho*h*omega(ff).');
                    ky = (pi/2)*gamma/h - 1j*(pi/2)*(1-gamma^2)/h;
                end
            end
            
            if startingSW   
                f = ky.*tanh(ky*h) - C(ff,i);
                df = tanh(ky*h) + ky*h.*sech(ky*h).^2; % the derivative in k
                ky = ky - f./df; % the Newton-Raphson update
                ky = abs(real(ky)) -1j*abs(imag(ky));
                if abs(f) <= epsilon || isnan(f)
                    k(ff,i) = ky;
                else
                    Ri = real(Z(1,i)); % resistance is real part of Z
                    gamma = Ri./(2*rho*h*omega(ff).');
                    ky = (pi/2)*gamma/h - 1j*(pi/2)*(1-gamma^2)/h;
                end
            end
    
        end
    end
       
    k = abs(real(k))-1j*abs(imag(k));
    k_to_int = k;
    k_to_int(isnan(k))=0;
    
    int_k = cumtrapz(k_to_int,2)*L/X;
    
    
    k0 = k(:,1);
    
    A = (k0*h./tanh(k0*h)).* ...
        sqrt((k0*h.*sech(k0*h).^2 + tanh(k0*h))./...
        (k*h.*sech(k*h).^2 + tanh(k*h)));
    
    P =  exp(-1j*int_k).*A;
    V = P./Z;




end

