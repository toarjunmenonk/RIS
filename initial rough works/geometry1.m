close all; clearvars; clc;

% cd 'C:\Users\IIITB\OneDrive - iiit-b\Desktop\IRS';
%% VARIABLE PARAMETERS 

relected_azim = -(pi/180)*50;%pi/8; % It should vary from 0 to -pi/2. 0 = partition wall, -pi/2 is the surface where RIS is kept.
relected_elev = -(pi/180)*0; % It should vary from -pi/2 to +pi/2. -pi/2 = vertically down along the roomwall, pi/2 = vertically up along the roomwall.

impinging_azim = (pi/180)*45; % It should vary from pi/2 to 0. 0 = partition wall, pi/2 is the surface where RIS is kept.
impinging_elev = (pi/180)*0; % It should vary from -pi/2 to +pi/2. -pi/2 = vertically down along the roomwall, pi/2 = vertically up along the roomwall.
height_of_ap   = 1.5;
height_of_ue   = 1.5;
print_logs = 1;
choice_RIS_res = 3;% 1: 1-bit, 2: 2-bit , 3: optimum

if ~(relected_azim<=0 && relected_azim >=-pi/2)
   error('Error. Reflected Azim angle must be less than 0 and greater than -pi/2')
end
if ~(relected_elev<=pi/2 && relected_elev>=-pi/2)
   error('Error. Reflected Elev angle is not valid')
end


if ~(impinging_azim<=pi/2 && impinging_azim >=0)
   error('Error. Reflected Azim angle must be less than 0 and greater than -pi/2')
end
if ~(impinging_elev<=pi/2 && impinging_elev>=-pi/2)
   error('Error. Reflected Elev angle is not valid')
end


%% Draw RIS Plane

nH = 16;
nV = 16;
fc = ((5.150 +5.875)/2)*1e9 ;
lambda = physconst('LightSpeed')/fc; 
global gradations; %#ok<*GVMIS>
gradations = 600;

width_RIS = nH*lambda/2;
length_RIS = nV*lambda/2;
y = linspace(-width_RIS/2,width_RIS/2,nH);
a=0;b=1;c=0;d=0;
x = linspace(-length_RIS/2,length_RIS/2,nV);
[z2,x2]=meshgrid(y ,x );
y2 = (1/b).*(d+(a.*x2)+(c.*z2));
height_of_RIS = 1.5;
z2 = z2 + height_of_RIS ;
surf( x2,y2,z2,'FaceColor','none');xlabel("<-- +ve azim [x axis (in meter)] --> -ve azim ");ylabel("y axis (in meter)");zlabel("z axis (in meter)");
hold on;

%% Draw x,y,z axes lines 
view(3);
k = 3.14;
line([-k,k],[0,0],[0,0],'Color', 'k', 'LineWidth', 2);
line([0,0],[0,k],[0,0],'Color', 'k', 'LineWidth', 2);
line([0,0],[0,0],[0,k],'Color', 'k', 'LineWidth', 2);

%% Draw obstacle (wall) between UE and AP 
door_width = 0.5; % All distances in meters
room_height = 3;
wall_width = 0.1;
room_length = 5;
h = fill3([0 0 wall_width wall_width], [ 0  room_length  room_length  0]+door_width,         [0 0 0  0]  ,'b', ...
          [0 0 wall_width wall_width], [ 0  room_length  room_length  0]+door_width,         [room_height room_height room_height room_height],'b', ...
          [0 0 0 0],                   [ 0  room_length  room_length  0]+door_width,         [ 0 0 room_height room_height],'b',...
          [wall_width wall_width wall_width wall_width],[ 0  room_length  room_length  0]+door_width,[ 0 0 room_height room_height],'b',...
          [0 0 wall_width wall_width], [0 0 0 0]+door_width,[0 room_height room_height 0],'b',...
          [0 0 wall_width wall_width], [room_length room_length room_length room_length]+door_width,[0 room_height room_height 0],'b');

set(h,'FaceAlpha',0.3) ;

%% Draw LoS ray - RIS to UE(Negative azim angle)


[x_source,y_source,z_source]=getline(relected_azim,relected_elev);
quiver3(0,0,height_of_RIS,x_source,y_source,z_source,1,'LineWidth',3,'Color','r');
hold on;

text(x_source,y_source,z_source+height_of_ue,'UE');hold on;
stem3(x_source,y_source,height_of_ue,"filled",'LineWidth',2,'Color','g');

%% Draw LoS ray - AP to RIS(Positive azim angle)


[x_source,y_source,z_source]=getline(impinging_azim,impinging_elev);
quiver3(x_source,y_source,height_of_ap,-x_source,-y_source,0,1,'LineWidth',3,'Color','r');
hold on;

text(x_source,y_source,z_source+height_of_ap,'AP');hold on;

stem3(x_source,y_source,height_of_ap,"filled",'LineWidth',2,'Color','g');


%% Generating Channels
h_AP2RIS = getchannel(impinging_azim,impinging_elev,nH,nV);
h_RIS2UE = getchannel(relected_azim,relected_elev,nH,nV);

% Compute optimim Array response of RIS
Psi_optimum = -angle(h_AP2RIS.*h_RIS2UE);
Psi_suboptimum_1bit = (pi/2) *sign(Psi_optimum);

Psi_suboptimum_2bit = floor(2*(1+Psi_optimum/pi));
Psi_suboptimum_2bit(Psi_suboptimum_2bit == 0) = -3*pi/4;
Psi_suboptimum_2bit(Psi_suboptimum_2bit == 1) = -pi/4;
Psi_suboptimum_2bit(Psi_suboptimum_2bit == 2) = pi/4;
Psi_suboptimum_2bit(Psi_suboptimum_2bit == 3) = 3*pi/4;



switch(choice_RIS_res)
    case 1
        BF = BF_pattern(nH,Psi_suboptimum_1bit,h_AP2RIS);
        plot_BF_pattern(BF,height_of_RIS);
        res = " | 1- bit res";
    case 2
        BF = BF_pattern(nH,Psi_suboptimum_2bit,h_AP2RIS);
        plot_BF_pattern(BF,height_of_RIS);
        res = " | 2- bit res";
    case 3
        BF = BF_pattern(nH,Psi_optimum,h_AP2RIS);
        plot_BF_pattern(BF,height_of_RIS);
        res = " | Optimum";
    otherwise
        disp("Wrong choice of RIS phase resolution");
end

txt = strcat("Impinging(AP2RIS): Azimuth(deg) = ",num2str(180*impinging_azim/pi),", Elev(deg) = ",num2str(180*impinging_elev/pi),res);
title(txt);

if(print_logs)
    foldername = strcat("ImpingingAz",num2str(180*impinging_azim/pi),"Elev",num2str(180*impinging_elev/pi),"ReflAz",num2str(180*relected_azim/pi),"Elev",num2str(180*relected_elev/pi));
    mkdir(foldername);
    cd (foldername) ;
    writematrix(Psi_optimum,"Psi_opt.txt");
    writematrix(Psi_suboptimum_1bit,"Psi_1bit.txt");
    writematrix(Psi_suboptimum_2bit,"Psi_2bit.txt");
    writematrix(h_AP2RIS,"h_ap2ris.txt");
    writematrix(h_RIS2UE,"h_ris2ue.txt");
    cd ..\ ;
end



%% LOCAL FUNCTIONS
function [x,y,z] = getline(impinging_azim,impinging_elev)

    [x,y,z] = sph2cart(pi/2-impinging_azim,impinging_elev,3.14);

end

function h = getchannel(azim, elev,nH,nV)

    arv1 = exp(-1i*pi*(0:(nH-1))*sin(azim)*cos(elev)).'; % RIS Reflector Array Phase Response
    arv2 = exp(-1i*pi*(0:(nV-1))*sin(elev)).';
    arv  = kron(arv1,arv2);
    h  = arv *exp(1i*2*pi*rand);
    
end
function BF = BF_pattern(nH,Psi,h_AP2RIS)

    global angleGrid;
    global azimGrid; 
    global elevGrid; 
    global gradations; %#ok<*GVMIS>

    angleGrid= linspace(-pi/2,pi/2,gradations);
    [azimGrid,elevGrid] = meshgrid(angleGrid, angleGrid);
    BF = zeros(size(azimGrid));


    h_subopt = conj(h_AP2RIS).*exp(-1i*Psi);

    for i = 1:length(angleGrid)
        for j = 1:length(angleGrid)

            arrayResponseVector1 = exp(-1i*pi*(0:(nH-1))*sin(azimGrid(i,j))*cos(elevGrid(i,j))).';
            arrayResponseVector2 = exp(-1i*pi*(0:(nH-1))*sin(elevGrid(i,j))).';
            arrayResponseVector = kron(arrayResponseVector2,arrayResponseVector1);
            BF(i,j) = abs(sum(arrayResponseVector.* h_subopt ));
        end
        disp([num2str(i) ' out of ' num2str(length(angleGrid)) ]);
    end
    BF =BF./256;
    BF = BF.^2;
end
function plot_BF_pattern(BF,h)

    global azimGrid; 
    global elevGrid; 
    norm_BF = BF/max(max(BF));
    magnification_factor = 8;
    [X, Y, Z]   = sph2cart(azimGrid,elevGrid, magnification_factor*norm_BF);
    colormap turbo
    BFsurface = surf(X,Y,Z+h, norm_BF,'EdgeColor','none');

    view([125 20])
    set(BFsurface,'LineStyle','none','FaceAlpha',0.5,'Tag','3D polar plot');
    
end

