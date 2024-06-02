% SCENARIO : 
% 1) There is an AP, RIS and UE in the nearby vicinity. 

% 2) No direct path between AP and UE. 

% 3) RIS has 16x16 sub-atoms

% 4) We assume AP as the Tx and UE as the Rx(Downlink scenario)

% AIM: 
% 1) Generate channels-  AP2RIS and RIS2UE

% 2) Calculate the optimum phase-shifts required for the RIS sub-atoms to
% get maximum SNR at the UE.

% 3) Calculate the Direction of Arrival (DoA) at the UE

close all; clearvars; clc;
rng(2024);
global azimGrid; %#ok<GVMIS>
global elevGrid; %#ok<GVMIS>

% Beam Steering direction
refl_azim = 60*pi/180;
refl_elev = 30*pi/180;

% Frequency of operation and Wavelength
fc = ((5.15 +5.875)/2)*1e9 ;
lambda = physconst('LightSpeed')/fc; 


% The width and height of an RIS element
d = lambda/4;


% Number of elements in the horizontal and vertical dimensions(We assume a square RIS nH = nV)
nH = 16; 


% Area of an element
A = d^2;


% Set the average intensity attenuations
mu = db2pow(-55);


% Set transmit power in dBm
PdBm = 30;


% Set the noise power in dBm
sigma2dBm = -174 + 10*log10(10e6) + 10;


% Compute the transmit power over the noise power in linear scale
Psigma2 = db2pow(PdBm - sigma2dBm);


% Defining h_AP2RIS
impinging_azim = pi/4;
impinging_elev = pi/4;
arv1 = exp(-1i*pi*(0:(nH-1))*sin(impinging_azim)*cos(impinging_elev)).';
arv2 = exp(-1i*pi*(0:(nH-1))*sin(impinging_elev)).';
arv  = kron(arv1,arv2);
h_AP2RIS  = arv *exp(1i*2*pi*rand);


% Defining h_RIS2UE

arv1 = exp(-1i*pi*(0:(nH-1))*sin(refl_azim)*cos(refl_elev)).';
arv2 = exp(-1i*pi*(0:(nH-1))*sin(refl_elev)).';
arv  = kron(arv1,arv2);
% h_RIS2UE  = sqrt(A*mu) * arv ;%* sqrt(0.5)*(randn + 1i*randn);
h_RIS2UE  = arv *exp(1i*2*pi*rand);


% Compute optimim Array response of RIS
Psi_optimum = -angle(h_AP2RIS.*h_RIS2UE);
Psi_suboptimum_1bit = (pi/2) *sign(Psi_optimum);
Psi_suboptimum_1bit_greedy = (pi/2) *sign(Psi_optimum);

%Compute the SNR with an optimized RIS
SNR_RIS = Psigma2*(sum(abs(h_AP2RIS.*h_RIS2UE),1)).^2;


% Compute the SNR with an sub-optimized RIS
SNR_suboptRIS = Psigma2*abs(sum(h_AP2RIS.*exp(1i*Psi_optimum).*h_RIS2UE,1)).^2;


%Compute the SNR with a random RIS configuration (see Footnote 3)
SNR_noOpt = Psigma2*abs(sum(h_AP2RIS.*h_RIS2UE,1)).^2;


% Compute the SNR at different azim, elev angles for suboptimal_psi value.
BF = BF_pattern(nH,Psi_suboptimum_1bit,h_AP2RIS);

% Plotting the direction of beamforming 
set(groot,'defaultAxesTickLabelInterpreter','latex');
surf(azimGrid*180/pi, elevGrid*180/pi, BF,'FaceAlpha',1);
shading interp;
xlim([-90,90]);ylim([-90,90]);
xlabel('X: Azim angle(deg) ($\varphi$)','Interpreter','latex','rotation', 20);
ylabel('Y: Elev angle(deg) ($\theta$)','Interpreter','latex','rotation', -35);
zlabel('Z: Normalised BF gain(lin scale)','Interpreter','latex');
axis square

function BF = BF_pattern(nH,Psi,h_AP2RIS)

    global angleGrid; %#ok<GVMIS>
    global azimGrid; %#ok<GVMIS>
    global elevGrid; %#ok<GVMIS>
    gradations = 60;

    angleGrid= linspace(-pi/2,pi/2,gradations);
    [azimGrid,elevGrid] = meshgrid(angleGrid, angleGrid);
    BF = zeros(size(azimGrid));


    h_subopt = conj(h_AP2RIS).*exp(-1i*Psi);
    % R = h_subopt*h_subopt';


    for i = 1:length(angleGrid)
        for j = 1:length(angleGrid)

            arrayResponseVector1 = exp(-1i*pi*(0:(nH-1))*sin(azimGrid(i,j))*cos(elevGrid(i,j))).';
            arrayResponseVector2 = exp(-1i*pi*(0:(nH-1))*sin(elevGrid(i,j))).';
            arrayResponseVector = kron(arrayResponseVector1,arrayResponseVector2);
            % BF(i,j) = 1/real((arrayResponseVector'/R)*arrayResponseVector);
            BF(i,j) = abs(sum(arrayResponseVector.* h_subopt ));
            
        end
        disp([num2str(i) ' out of ' num2str(length(angleGrid)) ]);
    end

end
