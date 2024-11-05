function phi = phase_fun(k,kr,n,Rsamp)

    kz = abs(sqrt(k^2 - kr^2)); % longitudinal wavevector 
    
    theta = atan(kr/kz); % cone angle of the beam
    
    alpha = atan(sin(theta) / (n - cos(theta))); %opening angle of axicon
    
    phi = -1 * k * Rsamp * sin(theta); % axicon phase function
    
    phi = mod(phi,2*pi); 
        
end