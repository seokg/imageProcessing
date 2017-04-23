%% 3. Geometric Transform

% translation in pixel
Tr_x = 100; Tr_y = 200;
% scale
Sc_x = 2; Sc_y = 1;
% rotation angle in radian
Ang = pi/18;

% translation matrix
Tr  = [1 0 Tr_x; 
       0 1 Tr_y; 
       0 0 1    ];

% scale matrix
S   = [Sc_x 0    0; 
       0    Sc_y 0; 
       0    0    1];
   
% rotation matrix
Rot = [cos(Ang) sin(Ang) 0 ; 
      -sin(Ang) cos(Ang) 0 ; 
       0        0        1  ];
   
% transfrom matrix
AffineT = Tr * S * Rot;
   
   
% initializing the geometry of the new image rectangle
width =  size(R,2); height = size(R,1);
scWidth = width*Sc_x; scHeight = height*Sc_y;

newHeight = floor(scWidth * sin(Ang) + scHeight * sin(pi/2-Ang))+Tr_y;
newWidth = floor(scWidth * cos(Ang) + scHeight * cos(pi/2-Ang))+Tr_x;

% find the mapping between original image and transformed image
% create homogeneous coordinate for final trasformed image
idxX = 1: newWidth; idxY = 1: newHeight;
[X,Y] = meshgrid(idxX, idxY);
Z = ones (1,newHeight*newWidth); 
idxTransImage = [X(1:newHeight*newWidth);Y(1:newHeight*newWidth);Z];

% initialized final transformed image
transImageR = zeros(newHeight, newWidth); 
transImageG = zeros(newHeight, newWidth); 
transImageB = zeros(newHeight, newWidth); 

% finding the corresponding index of transformed image
idxOriginImage = round(pinv(AffineT)*idxTransImage); 

% find index that is within the range
idxofidxOriginImage = find((idxOriginImage(1,:) > 0 & idxOriginImage(1,:) <= width)&(idxOriginImage(2,:) >0 & idxOriginImage(2,:) <= height));

% changing to row vector
rowIdx =sub2ind(size(R), idxOriginImage(2,idxofidxOriginImage), idxOriginImage(1,idxofidxOriginImage));
rowIdx2 = sub2ind(size(transImageR), idxTransImage(2,idxofidxOriginImage), idxTransImage(1,idxofidxOriginImage));

% give correct value for the projected area
transImageR(rowIdx2) = R(rowIdx);
transImageG(rowIdx2) = G(rowIdx);
transImageB(rowIdx2) = B(rowIdx);

transImage = zeros(newHeight, newWidth,3);
transImage(:,:,1) = transImageR;
transImage(:,:,2) = transImageG;
transImage(:,:,3) = transImageB;

figure, imshow(R);
figure, imshow(transImage);