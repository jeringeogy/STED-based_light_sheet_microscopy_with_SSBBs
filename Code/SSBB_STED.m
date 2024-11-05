%%% SIMULATION OF STED-BASED LIGHT SHEET USING SIDELOBE-SUPPRESSED BESSEL BEAMS%%%

close all; clear; clc;

N = 12000; % No of pixels
L= 400 * 10^-6; % Length of the element 
n = 1.5; % Refractive index
lambda1 = 647 * 10^-9; % Excitation wavelength
lambda2 = 750 * 10^-9; % Depletion Wavelength
k1 = 2*pi/lambda1; % Propagation constant 1
k2 = 2*pi/lambda2; % Propagation constant 2
R1 = zeros(N,N);
R2 = zeros(N,N);
z0 = 2.4 * 10^-6; % Rayleigh distance of a Gaussian beam with the same core radius
p_w = 250 * 10^-12; % pulse width of STED beam
p_p = 12.5 * 10^-9; % pulse period of STED beam
nl = 1.4; % Non-linearity rate of atto647N dye

%%%% SAMPLING SPACE %%%%
dx = L/N; 
x=-L/2:dx:L/2-dx; %x coordinate
y=x; %y coordinate
[X,Y]=meshgrid(x,y);
[XX,YY]=meshgrid(x,x);
fx=(-1/(2*dx):1/(L):1/(2*dx)-1/(L))*1; %Nyquist criteria
[FX,FY]=meshgrid(fx,fx);
Rsamp = sqrt((X).^2 + (Y).^2);
psi = atan2(YY,XX);

%%%% INPUT GAUSSIAN BEAM %%%%
w0 = 50 *10^-6; % 1/e2 radius of input beam
power = 1; % average power of STED beam
E0 = exp(-1 * (Rsamp).^2 / (w0/2).^2);
IP = sum(sum(abs(E0).^2)) / (power * p_p / p_w);
E0 = E0 / sqrt(IP);
pow1 = sum(sum(abs(E0).^2));


%%% EXCITATION BEAM %%%%
m = 0.7; % Excitation wave-vector ratio of zeroth-order SSBB
% Zeroth-order Bessel beam can be selected by giving m = 1
opt = 1; % 1 for excitation; 2 for depletion
[kr1, kr2, R1, R2, z] = opt_data(m,opt,N); % obtaining radial wavevectors and amplitudes

phi1 = phase_fun(k1,kr1,n,Rsamp); % Axicon 1 phase
phi2 = phase_fun(k1,kr2,n,Rsamp); % Axicon 2 phase

SSBB = (R1.*phi1) + (R2.*phi2); %Spatially multiplexing two axicon phase functions
SSBB = SSBB / max(max(SSBB));
SSBB1 = E0 .* exp(1i * 2 * pi * SSBB); % Gaussian beam illuminating SSBB_phase

for j = 1 : 3
    z1 = [z, z+z0, z+(2*z0)]; % different propagation planes (z = 0, z_R, and 2z_R)
    sty = ["b","g:","r-."];
    SSBB2 = BLAS(SSBB1, lambda1, z1(j), FX, FY, L); % Propagation simulation using band-limited angular spectrum
    I1 = abs(SSBB2(N/2,:)).^2; % line-profile of the output beam 
    I1 = I1 / max(I1);
    LS_E(:,j) = I1; % This variable holds the zeroth-order SSBB Y-Z profile. Propagation distance 'z' can be incremented in small steps to obtain fig. 5 in manuscript
    lim = 8;
    figure(1); plot_line(x,I1,sty(j),lim);title (sprintf('m = %.2f', m),'FontSize',16); % x-profile plot of excitation beam
    hold on;
end


%%%% DEPLETION BEAM %%%%
m = 0.8; % Depletion wave-vector ratio of first-order SSBB 
% First-order Bessel beam can be selected by giving m = 1
opt = 2; % 1 for excitation 2 for depletion
[kr_1, kr_2, R_1, R_2, z] = opt_data(m,opt,N); % obtaining radial wavevectors and amplitudes

phi_1 = phase_fun(k2,kr_1,n,Rsamp); % Axicon 1 phase
phi_2 = phase_fun(k2,kr_2,n,Rsamp); % Axicon 2 phase

SSB_B = (R_1.*phi_1) + (R_2.*phi_2); %Spatially multiplexing two axicon phase functions
SSB_B = SSB_B / max(max(SSB_B));
SSBB_1 = E0 .* exp(1i * 2 * pi * SSB_B) .* exp(1i * 1 * psi); % Gaussian beam illuminating SSBB_phase + SPP_phase

for j = 1 : 3
    z2 = [z, z+z0, z+(2*z0)]; % different propagation planes (z = 0, z_R, and 2z_R)
    sty1 = ["b","g:","r-."];
    SSBB_2 = BLAS(SSBB_1, lambda2, z2(j), FX, FY, L); % Propagation simulation using band-limited angular spectrum
    power = sum(sum(abs(SSBB_2).^2))
    I_1 = abs(SSBB_2(N/2,:)).^2; % line-profile of the output beam
    LS_D(:,j) = I_1; % This variable holds the first-order SSBB Y-Z profile. Propagation distance 'z' can be incremented in small steps to obtain fig. 5 in manuscript
    I_1 = I_1 / max(I_1);
    
    BL = abs(SSBB_2).^(2*nl); %cross-sectional bleaching for atto 647N
    BLEC(j) = sum(sum(BL));% calculating bleaching value

    lim = 8;
    figure(2); plot_line(x,I_1,sty1(j),lim);title (sprintf('m = %.2f', m),'FontSize',16); % x-profile plot of depletion beam
    hold on;
end

%%%% EFFECTIVE PSF %%%%
LS = STED_atto647(LS_E,LS_D,lambda2,L,N); % Calculating the effective PSF with excitation and depletion beams

for j = 1 : 3
    INT = LS(:,j);
    INT = INT / max(INT);
    sty1 = ["b","g:","r-."];
    lim = 4;
    figure(3); plot_line(x,INT,sty1(j),lim);title (sprintf('m = %.2f', m),'FontSize',16); % x-profile plot of effective PSF
    hold on;
end



