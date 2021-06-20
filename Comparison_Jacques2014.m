%% Based on inhomogeneous.m
clear all
close all

% Geometry - blood vessel in skin tissue
width_epidermis = 0.1; %100um
width_dermis = 2-width_epidermis; %2mm
diameter_bv = 0.25; %500um
width_laser = 0.8; %800um

xsize =  4;	% width of the region [mm]
ysize =  4;	% height of the region [mm]
dh = 0.005;         % discretisation size [mm]
vmcmesh = createRectangularMesh(xsize, ysize, dh);
%% Give varying optical parameters
% Water
% Set constant background coefficients (WATER)
vmcmedium.absorption_coefficient = 3.56*10^-5;     % absorption coefficient [1/mm]

% ** Changing scattering coefficient of water to 1 works **
vmcmedium.scattering_coefficient = 1;      % scattering coefficient [1/mm]


vmcmedium.scattering_anisotropy = 1;       % anisotropy parameter g of 
                                             % the Heneye-Greenstein scattering 
                                             % phase function [unitless]
vmcmedium.refractive_index = 1;            % refractive index [unitless]

% Resize the fields in vmcmedium so that they match the number of elements in the mesh
vmcmedium = createMedium(vmcmesh,vmcmedium);

% epidermis
% Select elements from the mesh
width = width_epidermis*2+width_dermis;                 % [mm]
height = ysize;
centercoord = [xsize/2-width_dermis/2,0];     % [mm]
elements_of_the_epidermis = findElements(vmcmesh, 'rectangle', centercoord, width, height);

% Assign a unique absorption coefficient to the selected elements
vmcmedium.absorption_coefficient(elements_of_the_epidermis) = 1.6;
vmcmedium.scattering_coefficient(elements_of_the_epidermis) = 37.6;
vmcmedium.scattering_anisotropy(elements_of_the_epidermis) = 0.9;

% Dermis
% Select elements from the mesh
width = width_dermis;                 % [mm]
height = ysize;
centercoord = [xsize/2-width_dermis/2,0];     % [mm]
elements_of_the_dermis = findElements(vmcmesh, 'rectangle', centercoord, width, height);

% Assign a unique absorption coefficient to the selected elements
vmcmedium.absorption_coefficient(elements_of_the_dermis) = 0.046;
vmcmedium.scattering_coefficient(elements_of_the_dermis) = 35.7;
vmcmedium.scattering_anisotropy(elements_of_the_dermis) = 0.9;

% Blood vessel
% Select elements from the mesh
radius = diameter_bv;                 % [mm]
centercoord = [1,0];     % [mm]
elements_of_the_bv = findElements(vmcmesh, 'circle', centercoord, radius);

% Assign a unique absorption coefficient to the selected elements
vmcmedium.absorption_coefficient(elements_of_the_bv) = 23.1;
vmcmedium.scattering_coefficient(elements_of_the_bv) = 9.4;
vmcmedium.scattering_anisotropy(elements_of_the_bv) = 0.9;


%% View geometry (scattering coefficient)
figure
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', vmcmedium.scattering_coefficient, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Scattering coeff';
hold off

%% Create a light source - laser beam 
lightsize = width_laser;
lightstart = floor((ysize/2-lightsize/2)/dh);
lightstop = ceil((ysize/2+lightsize/2)/dh);
vmcboundary.lightsource(lightstart:lightstop) = {'direct'};

%% Run the Monte Carlo simulation
options.photon_count = 1e6;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% Plot the solution - log10 of fluence
figure
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', log10(solution.element_fluence*100), 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'log( Fluence ) [J/cm^2]';
hold off

