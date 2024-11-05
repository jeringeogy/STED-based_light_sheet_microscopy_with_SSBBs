function LS = STED_atto647(Exc, Dep, lambda, L, N)

    k_fl = 1 / (3.9 * 10^-9); % spontaneous decay rate of Atto 647N
    k_vib = 1 / (5 * 10^-12); % vibrational relaxation rate of Atto 647N
    tau_STED = 250 * 10^-12; % pulse width of STED beam
    hc = 1.99 * 10^-25; % Planck's constant (h) and speed of light (c)
    sigma = 1 * 10^-16; % cross-section of stimulated emission for the dye
    I_S = k_fl * hc / (sigma * lambda); % saturation intensity (in per cm2)
    THRES = I_S * ((L * 10^2) / N)^2 ; % calculating the pixel saturation intensity 

    seta = Dep / (THRES); % Saturation factor
    gamma = (seta .* k_vib) ./ (seta .* k_fl + k_vib); % effective saturation factor
    eeta_ps  = (1 + gamma .* exp(-1 .* k_fl * tau_STED .* (1 + gamma))) ./ (1 + gamma); % probability of spontaneous decay
    LS = Exc .* eeta_ps; % effective PSF 
    
end
