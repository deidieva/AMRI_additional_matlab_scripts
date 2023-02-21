% Make a phantom consisting of N x N array of test tubes
% Tubes are along chosen direction, line 8

% Set dimensions
tube_radius = 0.0041;         % tube radius   (m)
% tube_inner_radius = 0.013;   % tube inner radius (m)
tube_length = 0.1;           % tube length   (m)
delta = 0.001;      % voxel size    (m)
direction = 'x';    % tube direction {'x', 'y', 'z'}
number_of_tubes = 1;             % Num,ber of tubes
Dx = 25e-2;         % Dimensions of the phantom     (m)
Dy = 25e-2;
Dz = 25e-2;

% Set MRI parameters (these are 1 x N*N vectors, the first element
% corresponds to the background material)
Rho = [0 1*ones(1, number_of_tubes*number_of_tubes)];
T1  = [0 0.145*ones(1, number_of_tubes*number_of_tubes)];
T2  = [0 0.058*ones(1, number_of_tubes*number_of_tubes)];
T2Star = [0 0.055*ones(1, number_of_tubes*number_of_tubes)];
ECon   = [0 0.3*ones(1, number_of_tubes*number_of_tubes)];
MassDen = [0 1045*ones(1, number_of_tubes*number_of_tubes)];

%% Voxelize the geometry
[X,Y,Z] = ndgrid(-Dx/2:delta:Dx/2, -Dy/2:delta:Dy/2, -Dz/2:delta:Dz/2); %create matrix 3d
A = zeros(size(X), 'uint8');

%y0 = linspace(-Dy/2+a, Dy/2-a, N);
%z0 = linspace(-Dz/2+a, Dz/2-a, N);
y0 = 0; % only meant for 1 element
z0 = 0; % only meant for 1 element
n = 1;
for y = y0
  for z = z0
    A( (Y - y).^2 + (Z - z).^2 <= tube_radius^2 ) = n;  % describes which element it is
    n = n+1;
  end
end
A(abs(X) > tube_length/2) = 0; % anything outside of phantom
%carve inner tube
% for y = y0
%   for z = z0
%     A( (Y - y).^2 + (Z - z).^2 <= tube_inner_radius^2 ) = 0;  
%     %n = n+1;
%   end
% end
A = A+1;


switch direction
  case 'y'
    A = permute(A, [3 1 2]);
  case 'z'
    A = permute(A, [2 3 1]);
end

%% Make the VObj struct
VObj.Gyro = 2.675380303797e+08;
VObj.Model = 'Normal';
VObj.Name = 'VObj_15x15_TubePhantom_';
VObj.Notes = ['Radius=' num2str(1e3*tube_radius) 'mm_Length=' num2str(1e3*tube_length) 'mm_Dir=' direction ];
VObj.Type = 'Water';
VObj.TypeNum = 1;
VObj.XDim = size(A,1);
VObj.XDimRes = delta;
VObj.YDim = size(A,2);
VObj.YDimRes = delta;
VObj.ZDim = size(A,3);
VObj.ZDimRes = delta;
VObj.Rho = Rho(A);
VObj.T1 = T1(A);
VObj.T2 = T2(A);
VObj.T2Star = T2Star(A);
VObj.ECon = repmat(ECon(A), [1 1 1 3]);
VObj.MassDen = MassDen(A);
VObj.ChemShift = 0;

%% Save VObj
% path = 'Z:\Desktop\MRiLab\MRiLab-master3\MRiLab-master3\MRiLab-master\Resources\VObj';
% filename = [VObj.Name '_' direction '.mat'];
% save([path filesep filename], 'VObj')
% [vol_handle]=VoxelPlotter(A); 
%visual effects (I recommend using the FigureRotator function from MATLAB
%Centeral
% view(3);
% daspect([1,1,1]);
% set(gca,'xlim',[0 gridesize], 'ylim',[0 gridesize], 'zlim',[0 gridesize]);
save new_tube_30mm_dia.mat VObj  % saved locally
