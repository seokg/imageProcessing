%% Open Image files
clear;

file = fopen(fullfile('Test_images\Gray_baboon_256x256.raw'),'rb');
Gray_baboon = fread(file,fliplr([256, 256]),'*uint8')';
fclose(file);

file = fopen(fullfile('Test_images\Gray_barbara_720x576.raw'),'rb');
Gray_barbara = fread(file,fliplr([720, 567]),'*uint8')';
fclose(file);

file = fopen(fullfile('Test_images\Gray_lenna_256x256.raw'),'rb');
Gray_lenna = fread(file,fliplr([256, 256]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\Color_baboon_256x256.raw'),'rb');
Color_baboon=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\Color_barbara_720x576.raw'),'rb');
Color_barbara=fread(file,fliplr([720, 576*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\Color_lenna_256x256.raw'),'rb');
Color_lenna=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

% Lenna tr1,2,3
file=fopen(fullfile('Test_images\lenna_tr1_256x256.raw'),'rb');
lenna_tr1=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\lenna_tr2_256x256.raw'),'rb');
lenna_tr2=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\lenna_tr3_256x256.raw'),'rb');
lenna_tr3=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

clear ans file;
%% 1.Math preliminary

H1 = [0.0751 0.1238 0.0751; 0.1238 0.2044 0.1238; 0.0751 0.1238 0.0751];
H2 = [0.1667 0.6667 0.1667; 0.6667 -3.3333 0.6667; 0.1667 0.6667 0.1667];

Gray_baboon_H1_result = conv2(double(H1), double(Gray_baboon));
Gray_barbara_H1_result = conv2(H1, double(Gray_barbara));
Gray_lenna_H1_result = conv2(H1, double(Gray_lenna));

Gray_baboon_H2_result = conv2(H2, double(Gray_baboon));
Gray_barbara_H2_result = conv2(H2, double(Gray_barbara));
Gray_lenna_H2_result = conv2(H2, double(Gray_lenna));

figure, imshow(uint8([Gray_baboon_H1_result Gray_baboon_H2_result]) ), title('H1 H2 Kernel baboon');
figure, imshow(uint8([Gray_barbara_H1_result Gray_barbara_H2_result]) ), title('H1 H2 Kernel barbara');
figure, imshow(uint8([Gray_lenna_H1_result Gray_lenna_H2_result]) ), title('H1 H2 Kernel lenna');

clear H1 H2 Gray_baboon_H1_result Gray_barbara_H1_result Gray_lenna_H1_result Gray_baboon_H2_result Gray_barbara_H2_result Gray_lenna_H2_result
%% 2. Color Coordinate Transform

[R G B] = RGBchan( Color_baboon );
a = zeros(256,256);
R_baboon = cat(3, R, a, a);
G_baboon = cat(3, a, G, a);
B_baboon = cat(3, a, a, B);

RGB_baboon = cat(3,R,G,B);
result_baboon= [RGB_baboon R_baboon G_baboon B_baboon];
figure, imshow(uint8(result_baboon)); title('RGB image');

clear R_baboon G_baboon  B_baboon RGB_baboon result_baboon;

%% RGB to HSI
R1 = double(R/255.0);
G1 = double(G/255.0);
B1 = double(B/255.0);

I = 1/3 *(R1+G1+B1);
S=1- (3./(R1+G1+B1)+0.000001).*(min(min(R1,G1),B1));
H = acosd ((2*R1-G1-B1) ./(2*sqrt( ( (R1-G1).^2 + (R1-B1).*(G1-B1) ) )+0.000001)); % 1/(2*pi)
H(B1>G1)=360-H(B1>G1);
H=H/360;

z = zeros(size(H));
HSI=zeros(size(R,1),size(R,2),3);
HSI(:,:,1)=H;
HSI(:,:,2)=S;
HSI(:,:,3)=I;

dispH = cat(3,H,z,z);
dispS = cat(3,S,z,z);
dispI = cat(3,I,z,z);

figure,imshow([HSI dispH dispS dispI]);title('HSI Image');

clear H S I HSI dispH dispS dispI
clear B1 G1 R1
%% RGB to YCbCr
R_vec = reshape(R, [1,size(R,1)*size(R,2)]);
G_vec = reshape(G, [1,size(G,1)*size(G,2)]);
B_vec = reshape(B, [1,size(B,1)*size(B,2)]);

YCbCr_vec = [0.299 0.587 0.114; -0.169 -0.331 0.500; 0.500 -0.419 -0.081]*[R_vec;G_vec;B_vec];
Y = reshape(YCbCr_vec(1,:), [size(R,1),size(R,2)]);
Cb = reshape(YCbCr_vec(2,:), [size(G,1),size(G,2)]) + 127;
Cr = reshape(YCbCr_vec(3,:), [size(B,1),size(B,2)]) + 127;
%%
YCbCr=zeros(size(R,1),size(R,2),3);
YCbCr(:,:,1)=Y;
YCbCr(:,:,2)=Cb;
YCbCr(:,:,3)=Cr;
dispY = cat(3,Y,Y,Y);
dispCb = cat(3,Cb,z,z);
dispCr = cat(3,Cr,z,z);
figure, imshow(uint8([YCbCr dispY dispCb dispCr]));title('YCbCr Image');

clear R_vec G_vec B_vec YCbCr_vec Y Cb Cr YCbCr dispY dispCr dispCb
%% 2.2 
[R1, G1, B1]= RGBchan(lenna_tr1);
R1 = R1/255; G1 = G1/255; B1 = B1/255;
[H1, S1, I1] = rgb2HSI(R1,G1,B1);
HSI1=zeros(size(R1,1),size(R1,2),3);
dispHSI1 = cat(3,H1,S1,I1);
dispH1 = cat(3,H1,z,z);
dispS1 = cat(3,S1,z,z);
dispI1 = cat(3,I1,z,z);


[R2, G2, B2]= RGBchan(lenna_tr2);
R2 = R2/255; G2 = G2/255; B2 = B2/255;
[H2, S2, I2] = rgb2HSI(R2,G2,B2);
dispHSI2 = cat(3,H2,S2,I2);
dispH2 = cat(3,H2,z,z);
dispS2 = cat(3,S2,z,z);
dispI2 = cat(3,I2,z,z);

[R3, G3, B3]= RGBchan(lenna_tr3);
R3 = R3/255; G3 = G3/255; B3 = B3/255;
[H3, S3, I3] = rgb2HSI(R3,G3,B3);
dispHSI3 = cat(3,H3,S3,I3);
dispH3 = cat(3,H3,z,z);
dispS3 = cat(3,S3,z,z);
dispI3 = cat(3,I3,z,z);

figure, imshow(([dispHSI1 dispH1,dispS1,dispI1 ;dispHSI2 dispH2,dispS2,dispI2 ;dispHSI3 dispH3,dispS3,dispI3 ]))

clear dispHSI1 dispH dispS1 dispI1 dispHSI2 dispH2 dispS2 dispI2 dispHSI3 dispH3 dispS3 dispI3 
clear R1 G1 B1 H1 S1 I1 R2 G2 B2 H2 S2 I2 R3 G3 B3 H3 S3 I3
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

originImage = zeros(size(R,1),size(R,2),3);
originImage(:,:,1) = R;
originImage(:,:,2) = G;
originImage(:,:,3) = B;
figure, imshow(uint8(originImage)), title('Orignal Image');
figure, imshow(uint8(transImage)), title('Transformed Image');

% clear transImage transImageR transImageG transImageB
% clear originImage R G B
% clear rowIdx rowIdx2

% clear Tr_x Tr_y Sc_x Sc_y Ang 
% clear Tr S Rot AffineT 
% clear width height scWidth scHeight newHeight newWidth idxX idxY X Y Z
% clear idxTransImage idxOriginImage idxofidxOriginImage
%% 4. Image down and up Sampling
% use pyramid and bell kernels
% down sample a test image of 256*256 to an image of 128* 128
% up-sample the down smapled image to 384*384

[R, G, B] = RGBchan(Color_lenna);
R1 = R / 255; G1 = G / 255; B1 = B / 255;

preSize = size(R1);
downScale = 2;
downSize = (1/downScale).*preSize;

pyramid = 1/4 * [1 2 1; 2 4 2; 1 2 1];
bell = 1/16*[1 3 3 1;3 9 9 3; 3 9 9 3; 1 3 3 1];


% use lowpass filter to blur the image
lowpass = fir1(preSize(1), 1/downScale);
blurImageR = conv2(R1,lowpass);
blurImageG = conv2(G1,lowpass);
blurImageB = conv2(B1,lowpass);

blurImage = zeros(size(blurImageR,1), size(blurImageR,2), 3);
blurImage(:,:,1) = blurImageR;
blurImage(:,:,2) = blurImageG;
blurImage(:,:,3) = blurImageB;

downSampleImage = zeros(downSize(1), downSize(2), 3);

downSampleImage(:,:,1) = blurImage(1:downScale:end, 1+((size(lowpass,2)-1)/2):downScale:end-round((size(lowpass,2)-1)/2),1 );
downSampleImage(:,:,2) = blurImage(1:downScale:end, 1+((size(lowpass,2)-1)/2):downScale:end-round((size(lowpass,2)-1)/2),2 );
downSampleImage(:,:,3) = blurImage(1:downScale:end, 1+((size(lowpass,2)-1)/2):downScale:end-round((size(lowpass,2)-1)/2),3 );

figure(), imshow(downSampleImage), title('Down Sampled Image');

clear blurImage blurImageB blurImageR blurImageG lowpass
clear preSize downScale downSize
% figure(), imshow(imresize(R1,[128,128]));

%%

upSampleImage = zeros(384,384,3);
idxX = 1: 128; idxY = 1: 128;
[X,Y] = meshgrid(idxX, idxY);
idx1 = [X(1:128*128);Y(1:128*128)];
idxDown = sub2ind(size(downSampleImage(:,:,1)),idx1(1,:),idx1(2,:));
idx2 = [X(1:128*128)*3;Y(1:128*128)*3];
idxUp = sub2ind(size(upSampleImage(:,:,1)), idx2(1,:), idx2(2,:));
upSampleImage(idxUp) = downSampleImage(idxDown);
upSampleImage(idxUp+384*384) = downSampleImage(idxDown+128*128);
upSampleImage(idxUp+384*384*2) = downSampleImage(idxDown+128*128*2);


figure, imshow((upSampleImage))
%%
% pyramidUpSampleR = zeros(386,386); pyramidUpSampleB = zeros(386,386,1); pyramidUpSampleG = zeros(386,386,1);
pyramidUpSampleR= conv2(upSampleImage(:,:,1),pyramid);
pyramidUpSampleG= conv2(upSampleImage(:,:,2), pyramid);
pyramidUpSampleB= conv2(upSampleImage(:,:,3),pyramid);

pyramidUpSample = zeros(386, 386,3);
pyramidUpSample(:,:,1) = pyramidUpSampleR;
pyramidUpSample(:,:,2) = pyramidUpSampleG;
pyramidUpSample(:,:,3) = pyramidUpSampleB;


figure, imshow(pyramidUpSample), title('Upsample Pyramid');

% bellUpSample= conv2(upSampleImage,bell);
bellUpSampleR= conv2(upSampleImage(:,:,1),bell);
bellUpSampleG= conv2(upSampleImage(:,:,2), bell);
bellUpSampleB= conv2(upSampleImage(:,:,3),bell);

bellUpSample = zeros(387, 387,3);
bellUpSample(:,:,1) = bellUpSampleR;
bellUpSample(:,:,2) = bellUpSampleG;
bellUpSample(:,:,3) = bellUpSampleB;


figure, imshow(bellUpSample), title('Upsample Bell');

% %%
% convR = conv2(pyramid,R);
% convG = conv2(pyramid,G);
% convB = conv2(pyramid,B);
% 
% downSampPyramid = zeros(128,128,3);
% downSampPyramid(:,:,1) = convR(2:2:end-1,2:2:end-1);
% downSampPyramid(:,:,2) = convG(2:2:end-1,2:2:end-1);
% downSampPyramid(:,:,3) = convB(2:2:end-1,2:2:end-1);
% 
% imshow(uint8(downSampPyramid));
% % impyramid
% %# Initializations:
% 
% scale = [1/2 1/2];              %# The resolution scale factors: [rows columns]
% oldSize = size(R);                   %# Get the size of your image
% newSize = max(floor(scale.*oldSize(1:2)),1);  %# Compute the new image size
% 
% %# Compute an upsampled set of indices:
% 
% rowIndex = min(round(((1:newSize(1))-0.5)./scale(1)+0.5),oldSize(1));
% colIndex = min(round(((1:newSize(2))-0.5)./scale(2)+0.5),oldSize(2));
% 
% %# Index old image to get new image:
% 
% outputImage = R(rowIndex,colIndex);


%% 5. Quantization
%% 5.1 2 / 4 / 6 bit uniform quantization


[R1, G1, B1] = RGBchan(Color_lenna);
nbitQuant = 2; 
[ finalImage2 ] = unifromQuantization( nbitQuant, R1, G1, B1);
nbitQuant = 4; 
[ finalImage4 ] = unifromQuantization( nbitQuant, R1, G1, B1);
nbitQuant = 6; 
[ finalImage6 ] = unifromQuantization( nbitQuant, R1, G1, B1);

orgImage =  zeros(size(R1,1),size(R1,2),3);
orgImage(:,:,1) = R1;
orgImage(:,:,2) = G1;
orgImage(:,:,3) = B1;



figure, imshow(uint8([ finalImage2  finalImage4  finalImage6 orgImage])),title('2 3 6 bit Qunatization and Original Image')
[ PSNR1 ] = findPNSR( orgImage, finalImage2 );
[ PSNR2 ] = findPNSR( orgImage, finalImage4 );
[ PSNR3 ] = findPNSR( orgImage, finalImage6 );

%% 5.1 2 / 4 / 6 bit uniform quantization

% compress
input = 0:255;
fx = zeros(length(input),1);
for i = 1 : length(input)
    x = 0 : i-1;
    numeraPu = sum((1/255)^(1/3)*x);
    denomiPu= sum((1/255)^(1/3)*[0:255]);
    fx(i,1) = 255* numeraPu/denomiPu;
end



% plot(0:255,fx)

[R1, G1, B1] = RGBchan(Color_lenna);
orgImage =  zeros(size(R1,1),size(R1,2),3);
orgImage(:,:,1) = R1;
orgImage(:,:,2) = G1;
orgImage(:,:,3) = B1;

[ finalImage2 ] = companderQuantization( 2, R1, G1, B1 ,fx);
[ finalImage4 ] = companderQuantization( 4, R1, G1, B1 ,fx);
[ finalImage6 ] = companderQuantization( 6, R1, G1, B1 ,fx);
figure, imshow(uint8([finalImage2 finalImage4 finalImage6 orgImage])),title('compander')
[ comPSNR1 ] = findPNSR( orgImage, finalImage2 );
[ comPSNR2 ] = findPNSR( orgImage, finalImage4 );
[ comPSNR3 ] = findPNSR( orgImage, finalImage6 );

%% 5.2
randNoise = 32*rand(size(R1))-ones(size(R1))*16;

[R1, G1, B1] = RGBchan(Color_lenna);
orgImage =  zeros(size(R1,1),size(R1,2),3);
orgImage(:,:,1) = R1;
orgImage(:,:,2) = G1;
orgImage(:,:,3) = B1;

R1 = R1+randNoise;
G1 = G1+randNoise;
B1 = B1+randNoise;

noiseImage =  zeros(size(R1,1),size(R1,2),3);
noiseImage(:,:,1) = R1;
noiseImage(:,:,2) = G1;
noiseImage(:,:,3) = B1;

[ quantImg2 ] = unifromQuantization( 2, R1, G1, B1);
[ quantImg4 ] = unifromQuantization( 4, R1, G1, B1);
[ quantImg6 ] = unifromQuantization( 6, R1, G1, B1);
[ PSNR1 ] = findPNSR( orgImage, quantImg2 );
[ PSNR2 ] = findPNSR( orgImage, quantImg4 );
[ PSNR3 ] = findPNSR( orgImage, quantImg6 );
figure, imshow(uint8([quantImg2 quantImg4 quantImg6 noiseImage orgImage])), title('2 3 6 bit Qunatization, Noise Image, and Original Image');


rR = round(R1);
rR = rR.*~(rR>255) + (rR>255)*255;
rR = rR.*~(rR<=0)+ (rR<=0);

rG = round(G1);
rG = rG.*~(rG>255) + (rG>255)*255;
rG = rG.*~(rG<=0) + (rG<=0);

rB = round(B1);
rB = rB.*~(rB>255) + (rB>255)*255;
rB = rB.*~(rB<=0)+ (rB<=0);

[ finalImage2 ] = companderQuantization( 2, rR, rG, rB ,fx);
[ finalImage4 ] = companderQuantization( 4, rR, rG, rB ,fx);
[ finalImage6 ] = companderQuantization( 6, rR, rG, rB ,fx);
figure, imshow(uint8([finalImage2 finalImage4 finalImage6 noiseImage orgImage])),title('compander')
[ comPSNR1 ] = findPNSR( orgImage, finalImage2 );
[ comPSNR2 ] = findPNSR( orgImage, finalImage4 );
[ comPSNR3 ] = findPNSR( orgImage, finalImage6 );
