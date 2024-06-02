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

%Frequency of operation and Wavelength
fc = ((5.15 +5.875)/2)*1e9 ;
lambda = physconst('LightSpeed')/fc; 


%The width and height of an RIS element
d = lambda/4;


%Number of elements in the horizontal and vertical dimensions(We assume a square RIS nH = nV)
nH = 16; 


%Area of an element
A = d^2;


%Set the average intensity attenuations
mu = db2pow(-55);


%Set transmit power in dBm
PdBm = 30;


%Set the noise power in dBm
sigma2dBm = -174 + 10*log10(10e6) + 10;


%Compute the transmit power over the noise power in linear scale
Psigma2 = db2pow(PdBm - sigma2dBm);


%Generate a grid for the elements
gridPoints = (0:nH-1)*d;
[X,Y] = meshgrid(gridPoints,gridPoints);
locations = X(:)+1i*Y(:);


%Total number of elements
N = length(locations);


%Compute the spatial correlation matrix
R = zeros(N,N);

for m = 1:N
    for l = 1:N
        
        R(m,l) = sinc(2*abs(locations(m)-locations(l))/lambda);
        
    end
end


%Generate channel realizations
Rsqrtm = sqrtm(R);
h1 = sqrt(A*mu) * Rsqrtm * (randn(N,1) + 1i*randn(N,1))/sqrt(2);
h2 = sqrt(A*mu) * Rsqrtm * (randn(N,1) + 1i*randn(N,1))/sqrt(2);


%Compute the SNR with an optimized RIS
SNR_RIS = Psigma2*(sum(abs(h1.*h2),1)).^2;


% sub-optimum RIS phase shift array with 1-bit granularity 
% subopt_psi = -sign(angle(h1.*h2))*pi/2;
subopt_psi = -angle(h1.*h2);

% Compute the SNR with an sub-optimized RIS
SNR_suboptRIS = Psigma2*abs(sum(h1.*exp(1i*subopt_psi).*h2,1)).^2;


%Compute the SNR with a random RIS configuration (see Footnote 3)
SNR_noOpt = Psigma2*abs(sum(h1.*h2,1)).^2;


% Compute the SNR at different azim, elev angles for suboptimal_psi value.
angleGrid = linspace(-pi/2,pi/2,361);
[azimGrid,elevGrid] = meshgrid(angleGrid, angleGrid);

BF = BF_pattern(azimGrid,elevGrid,nH,subopt_psi);


% Plotting the direction of beamforming 

surfl(azimGrid, elevGrid, BF);
colormap(pink);    % change color map
shading interp;    % interpolate colors across lines and faces


function BF = BF_pattern(azimGrid,elevGrid,nH,Psi)
    
    Psi = exp(1i*Psi);
    angleGrid = linspace(-pi/2,pi/2,361);
    BF = zeros(size(azimGrid));

    % Azimuth, Elevation  angle of the source in radians
    azimAngles = pi/4; elevAngles = -pi/4;
    

    % Array response matrix
    arrayResponseVector1 = exp(-1i*pi*(0:(nH-1))*sin(azimAngles)*cos(elevAngles)).';
    arrayResponseVector2 = exp(-1i*pi*(0:(nH-1))*sin(elevAngles)).';
    
    A = kron(arrayResponseVector1,arrayResponseVector2);
    % A = A.*Psi;
    Y2 = A*A';


    for i = 1:length(angleGrid)
        for j = 1:length(angleGrid)

            arrayResponseVector1 = exp(-1i*pi*(0:(nH-1))*sin(azimGrid(i,j))*cos(elevGrid(i,j))).';
            arrayResponseVector2 = exp(-1i*pi*(0:(nH-1))*sin(elevGrid(i,j))).';
            arrayResponseVector = kron(arrayResponseVector1,arrayResponseVector2);
            % arrayResponseVector = arrayResponseVector.* Psi;
            BF(i,j) = 1/real((arrayResponseVector'/Y2)*arrayResponseVector);
            
        end
    end

end
