close all; clearvars; clc;

theta1 = linspace(-pi/2,pi/2,20);
phi1 = linspace(-pi/2,pi/2,20);

[theta,phi] = meshgrid(theta1, phi1);
BF = cos(theta).*phi; 
MagE        = reshape(BF,length(phi1),length(theta1));
r           = reshape(BF,length(phi1),length(theta1));

[X, Y, Z]   = sph2cart(phi, theta, r./max(max(r)));

surfHdl = surf(X,Y,Z,MagE, 'FaceColor','interp');
set(surfHdl,'LineStyle','none','FaceAlpha',1.0,'Tag','3D polar plot');


function [X, Y, Z]= sph2cart(phi, theta, r)


    Z  = r.*cosd(theta);
    X  = r.*sind(theta).*cosd(phi);
    Y  = r.*sind(theta).*sind(phi);
end