%% Based on inhomogenous.m
clear all
close all

%% Geometry: a fluid filled space, with a layer of bone, then a layer of neural tissue
width_bone = 0.12;
width_neural = 0.1;
dist_neural = 0.5;

xsize =  width_neural/2+dist_neural;	% width of the region [mm]
ysize =  0.3;	% height of the region [mm]
dh = 0.002;         % discretisation size [mm]
vmcmesh = createRectangularMesh(xsize, ysize, dh);
%% Give varying optical parameters6
% For fluid filled space
% Set constant background coefficients
vmcmedium.absorption_coefficient = 0.96;     % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient = 0.0;      % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy = 0.85;       % anisotropy parameter g of 
                                             % the Heneye-Greenstein scattering 
                                             % phase function [unitless]
vmcmedium.refractive_index = 1.33;            % refractive index [unitless]

% Resize the fields in vmcmedium so that they match the number of elements in the mesh
vmcmedium = createMedium(vmcmesh,vmcmedium);

% Bone
% Select elements from the mesh (rectangle)
width = width_bone;                 % [mm]
height = ysize;
centercoord = [xsize/2-width_neural/2,0];% [xsize/2-width_neural-width_bone/2,0];     % [mm]
elements_of_the_bone = findElements(vmcmesh, 'rectangle', centercoord, width, height);

% Assign a unique absorption coefficient to the selected elements
vmcmedium.scattering_coefficient(elements_of_the_bone) = 2;

% Neural tissue
% Select elements from the mesh
width = width_neural;                 % [mm]
height = ysize;
centercoord = [xsize/2-width_neural/2,0];     % [mm]
elements_of_the_neural = findElements(vmcmesh, 'rectangle', centercoord, width, height);

% Assign a unique absorption coefficient to the selected elements
vmcmedium.scattering_coefficient(elements_of_the_neural) = 1;

%% Create a light source - Optical fibre with Gaussian distribution
lightsize = 0.2;
lightstart = floor((ysize/2-lightsize/2)/dh);
lightstop = floor((ysize/2+lightsize/2)/dh);
vmcboundary.lightsource(lightstart:lightstop) = {'gaussian'};
vmcboundary.lightsource_gaussian_sigma = 0.1;

%% View the geometry
figure
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', vmcmedium.scattering_coefficient, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Fluence [J/mm^2]';
hold off

%% Run the Monte Carlo simulation
options.photon_count = 1e6;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% Plot the solution
% Raw results
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Fluence [J/mm^2]';
hold off

% Results scaled to match the paper
figure
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence/40000000, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Fluence [J/mm^2]';
hold off

%% Absorption

% fluence = absorption/muA
% absorption = fluence * muA

myfluence = solution.element_fluence;
myfluence = myfluence.*vmcmedium.absorption_coefficient;
figure;
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', myfluence/40000000, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Absorption [J/mm^3]';
hold off
