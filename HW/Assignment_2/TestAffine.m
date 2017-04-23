% Affine Transform Script
clear;
I = checkerboard;
%I =[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 0; 0 0 0 0 0 0; 0 0 0 00 0] % Define the image
figure;imagesc(I); axis('image'); colormap('gray');title('Original Image')
%Display the image
xdim = 256*2;
ydim = 256*2;
I_affine = zeros(xdim,ydim);
% Define the Transformation Matrix
T = [ 1 0 xdim/4; 0 1 ydim/4; 0 0 1] % Translation
R = [cos(20/180*pi) sin(20/180*pi) 0; -sin(20/180*pi) cos(20/180*pi) 0; 0 0 1] % Rotation
S = [1 0 0;0 1 0;0 0 1] % Scaleing
%H_x = [1 tan(10/180*pi) 0;0 1 0;0 0 1] % Shear in x only %tan(20/180*pi)
H_x = [1 0 0;0 1 0;0 0 1] % Shear in x only %tan(20/180*pi)
T_affine = T*S*R*H_x % Calculate the Affine Transformation Matrix

% Define indices of Coordinates to Transform
x = 1:length(I(:,1));
y = 1:length(I(1,:));
[X, Y] = meshgrid(x,y); % This automatically creates a matrix made up of the x and y coordinates
Z = ones(length(x),length(y));
clear Index;
Index(:,:,1) = X;
Index(:,:,2) = Y;
Index(:,:,3) = Z;
Index; % The columns of this 3D matrix are the coordinates to transform

% Perform Affine Transformation with Nearest Neighbor Interpolation
for j = 1:length(y),
for i = 1:length(x),
temp = Index(i,j,:); % temp is the column representing the coordinates
temp = squeeze(temp); % collapse to a 1X3 matrix
Test = T_affine*temp; % apply the transformation to the coordinates
I_affine(round(Test(2)),round(Test(1))) = I(temp(2),temp(1));
% round() performs the nearest neighbor interpolation
end
end

figure, imshow(I_affine);

I_affine = zeros(xdim,ydim);
% Perform Affine Tranformation with Bi-linear Interpolation
for j = 1:length(y),
for i = 1:length(x),
temp = Index(i,j,:); % temp is the column representing the coordinates
temp = squeeze(temp); % collapse to a 1 X 3 matrix
Test = T_affine*temp; % apply the tranformation to the coordinates
alpha = (Test(1)-floor(Test(1)));
beta = (Test(2)-floor(Test(2)));
w_11 = (1-alpha)*(1-beta);
w_12 = alpha*(1-beta);
w_22 = beta*alpha;
w_21 = (1-alpha)*(beta);
%if ceil(Test(1)) < length(x) && ceil(Test(2)) < length(y)
I_affine(floor(Test(2)),floor(Test(1))) = I_affine(floor(Test(2)),floor(Test(1)))+ w_11*I(i,j);
I_affine(floor(Test(2)),ceil(Test(1))) = I_affine(floor(Test(2)),ceil(Test(1)))+w_12*I(i,j);
I_affine(ceil(Test(2)),floor(Test(1))) = I_affine(ceil(Test(2)),floor(Test(1)))+w_21*I(i,j);
I_affine(ceil(Test(2)),ceil(Test(1))) = I_affine(ceil(Test(2)),ceil(Test(1)))+w_22*I(i,j); % load the values of the original matrix into the new coordinates
%end
% w_xx perform the bi-linear interpolation
end
end

figure, imshow(I_affine);