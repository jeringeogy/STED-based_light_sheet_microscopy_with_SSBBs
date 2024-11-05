function F1 = BLAS(F, lambda, z, FX, FY, L)
    
    u_lim = 1 / (sqrt((2 / L * z)^2 + 1) .* lambda); % band-limiting condition

    ANG=exp(-1i*pi*lambda*z*(FX.^2+FY.^2)); % transfer function
    % The following command implements the band-limiting condition
    % When evaluating the bleaching value, the band-limiting conditon is avoided to ensure that the total power remain conserved
    ANG = ANG .* (sqrt(FX.^2 + FY.^2) <= u_lim);
    ANG=fftshift(ANG);
    
    % Angular spectrum propagation
    F1=fft2(fftshift(F)); % Fourier transform of the input signal
    F1=F1.*ANG; % Multiplying fft(input signal) with transfer function
    F1=ifftshift(ifft2(F1)); % Taking inverse Fourier transform 

end