% This Matlab script generates Figure 4 in the paper:
%
% Emil Björnson, Luca Sanguinetti, “Rayleigh Fading Modeling and Channel
% Hardening for Reconfigurable Intelligent Surfaces,” IEEE Wireless
% Communications Letters, To appear.
%
% Download article: https://arxiv.org/pdf/2009.04723.pdf
%
% This is version 1.0 (Last edited: 2021-01-01)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.


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
values_RIS = Psigma2*(sum(abs(h1.*h2),1)).^2;

%Compute the SNR with a random RIS configuration (see Footnote 3)
values_noOpt = Psigma2*abs(sum(h1.*h2,1)).^2;

