close all; clear; clc;

N=12000; % No of pixels
lambda1 = 647 * 10^-9; % Excitation wavelength
lambda2 = 750 * 10^-9; % Depletion wavelength
L= 400 * 10^-6; % Length of the element 
fwhm1 = 830 * 10^-9; % excitation beam FWHM
fwhm2 = 1254 * 10^-9; % this gives 720 nm null width for depletion beam
w01 = fwhm1 / 1.18; % radius
w02 = fwhm2 / 1.18; % radius
k1 = 2*pi/lambda1;
k2 = 2*pi/lambda2;
n = 1.5;
z01 = pi * w01^2 / lambda1;
z02 = pi * w02^2 / lambda2;
z = 0 * 10^-6;


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
avg_pow = 200 * 10^-3;
p_w = 250 * 10^-12; % pulse width of STED beam
p_p = 12.5 * 10^-9; % pulse perion of STED beam
power = avg_pow * p_p / p_w;

%%% EXCITATION BEAM %%%%

for j = 1 : 3
    z1 = [z, z+z01, z+(2*z01)];
    sty = ["b","g:","r-."];
    wz = w01 * sqrt(1 + (z1(j)/z01)^2);
    Rz = z1(j) * (1 + (z01/z1(j))^2);
    GB = (w01/wz) * exp (-(Rsamp).^2 / wz^2) * exp(-1i*k1*z1(j));
    I1 = abs(GB(N/2,:)).^2;
    I1 = I1 / max(I1);
    LS_E(:,j) = I1;
    lim = 8;
    figure(1); plot_line(x,I1,sty(j),lim);title ('Conventional-STED','FontSize',16);
    hold on;
end

%%%% DEPLETION BEAM %%%%

for j = 1 : 3
    p = 0;
    l = 1;
    z2 = [z, z+z01, z+(2*z01)];
    U00 = 1/(1 + 1i*z2(j)/z02) .* exp(-Rsamp.^2/w02^2./(1 + 1i*z2(j)/z02));
    w = w02 * sqrt(1 + z2(j).^2/z02^2);
    R = sqrt(2)*Rsamp./w;

    % Lpl from OT toolbox (Nieminen et al., 2004)
    Lpl = nchoosek(p+l,p) * ones(size(R));   % x = R(r, z).^2
    for m = 1:p
        Lpl = Lpl + (-1)^m/factorial(m) * nchoosek(p+l,p-m) * R.^(2*m);
    end
    LG = U00.*R.^l.*Lpl.*exp(1i*l*psi).*exp(-1i*(2*p + l + 1)*atan(z/z02));
    pow = sum(sum(abs(LG).^2)) / power;
    LG = LG / sqrt(pow);
    I_1 = abs(LG(N/2,:)).^2;
    LS_D(:,j) = I_1;
    I_1 = I_1 / max(I_1);
    lim = 8;
    figure(2); plot_line(x,I_1,sty(j),lim);title ('Conventional-STED','FontSize',16);
    hold on;
end

% 
LS = STED_atto647(LS_E,LS_D,lambda2,L,N);


for j = 1 : 3
    INT = LS(:,j);
    INT = INT / max(INT);
    sty1 = ["b","g:","r-."];
    lim = 4;
    figure(3); plot_line(x,INT,sty1(j),lim);title ('Conventional-STED','FontSize',16);
    hold on;
end





