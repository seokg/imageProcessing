
%% 1.Math preliminary
clear;

H1 = [0.0751 0.1238 0.0751; 0.1238 0.2044 0.1238; 0.0751 0.1238 0.0751];
H2 = [0.1667 0.6667 0.1667; 0.6667 -3.3333 0.6667; 0.1667 0.6667 0.1667];

file = fopen(fullfile('Test_images\Gray_baboon_256x256.raw'),'rb');
Gray_baboon = fread(file,fliplr([256, 256]),'*uint8')';
fclose(file);

file = fopen(fullfile('Test_images\Gray_barbara_720x576.raw'),'rb');
Gray_barbara = fread(file,fliplr([720, 567]),'*uint8')';
fclose(file);

file = fopen(fullfile('Test_images\Gray_lenna_256x256.raw'),'rb');
Gray_lenna = fread(file,fliplr([256, 256]),'*uint8')';
fclose(file);

Gray_baboon_H1_result = conv2(double(H1), double(Gray_baboon));
Gray_barbara_H1_result = conv2(H1, double(Gray_barbara));
Gray_lenna_H1_result = conv2(H1, double(Gray_lenna));

Gray_baboon_H2_result = conv2(H2, double(Gray_baboon));
Gray_barbara_H2_result = conv2(H2, double(Gray_barbara));
Gray_lenna_H2_result = conv2(H2, double(Gray_lenna));


[Gray_baboon_row, Gray_baboon_col]= size(Gray_baboon_H2_result);
[Gray_barbara_row, Gray_barbara_col]= size(Gray_barbara_H2_result);
[Gray_lenna_row, Gray_lenna_col]= size(Gray_lenna_H2_result);

Problem_1_H1_result = zeros(Gray_barbara_row, Gray_baboon_col+Gray_barbara_col+Gray_lenna_col);

Problem_1_H1_result(1:Gray_baboon_row,1:Gray_baboon_col) = Gray_baboon_H1_result;
Problem_1_H1_result(1:Gray_barbara_row,1+Gray_baboon_col:Gray_baboon_col+Gray_barbara_col)= Gray_barbara_H1_result;
Problem_1_H1_result(1:Gray_lenna_row,1+Gray_baboon_col+Gray_barbara_col:Gray_baboon_col+Gray_barbara_col+Gray_lenna_col) = Gray_lenna_H1_result;

Problem_1_H2_result = zeros(Gray_barbara_row, Gray_baboon_col+Gray_barbara_col+Gray_lenna_col);

Problem_1_H2_result(1:Gray_baboon_row,1:Gray_baboon_col) = Gray_baboon_H2_result;
Problem_1_H2_result(1:Gray_barbara_row,1+Gray_baboon_col:Gray_baboon_col+Gray_barbara_col)= Gray_barbara_H2_result;
Problem_1_H2_result(1:Gray_lenna_row,1+Gray_baboon_col+Gray_barbara_col:Gray_baboon_col+Gray_barbara_col+Gray_lenna_col) = Gray_lenna_H2_result;

figure, title('Problem_1_H1_result');
imshow(uint8(Problem_1_H1_result));

figure, title('Problem_1_H2_result');
imshow(uint8(Problem_1_H2_result));

%% 2. Color Coordinate Transform

file=fopen(fullfile('Test_images\Color_baboon_256x256.raw'),'rb');
Color_baboon=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\Color_barbara_720x576.raw'),'rb');
Color_barbara=fread(file,fliplr([720, 576*3]),'*uint8')';
fclose(file);

file=fopen(fullfile('Test_images\Color_lenna_256x256.raw'),'rb');
Color_lenna=fread(file,fliplr([256, 256*3]),'*uint8')';
fclose(file);

a = zeros(256,256);

% 
% R = Color_baboon(:,1:3:end);
% R_baboon = cat(3, R, a, a);
% 
% G = Color_baboon(:,2:3:end);
% G_baboon = cat(3, a, G, a);
% 
% B = Color_baboon(:,3:3:end);
% B_baboon = cat(3, a, a, B);
% 
% RGB_baboon = cat(3,R,G,B);
% result_baboon= [RGB_baboon R_baboon G_baboon B_baboon];
% 
% figure, imshow(result_baboon)


R = double(Color_baboon(:,1:3:end))/256;
R_baboon = cat(3, R, a, a);

G = double(Color_baboon(:,2:3:end))/256;
G_baboon = cat(3, a, G, a);

B = double(Color_baboon(:,3:3:end))/256;
B_baboon = cat(3, a, a, B);

RGB_baboon = cat(3,R,G,B);
result_baboon= [RGB_baboon R_baboon G_baboon B_baboon];
figure, imshow(result_baboon)

% % RGB to HSI
% I = 1/3 (R+G+B)
% S = 1 - min(R,G,B)/I
% H = 1/2*pi *inv_cos( (2*R-G-B)/(2*sqrt((R-G)^2 + (R-B)*(G-B))) )

I = 1/3 *(R+G+B);
I_baboon = cat(3, I, I, I);

S = 1 - min(min(R,G),B)./I;

H = 1/(2*pi) * acos (double(2*R-G-B)./(2*sqrt( double((R-G).^2 + (R-B).*(G-B)) )));

figure, imshow([H S I])


% RGB to YCbCr
R_vec = reshape(R, [1,size(R,1)*size(R,2)]);
G_vec = reshape(G, [1,size(G,1)*size(G,2)]);
B_vec = reshape(B, [1,size(B,1)*size(B,2)]);

YCbCr_vec = [0.299 0.587 0.114; -0.169 -0.331 0.500; 0.500 -0.419 -0.081]*[R_vec;G_vec;B_vec];
Y = reshape(YCbCr_vec(1,:), [size(R,1),size(R,2)]);
Cb = reshape(YCbCr_vec(2,:), [size(G,1),size(G,2)]);
Cr = reshape(YCbCr_vec(3,:), [size(B,1),size(B,2)]);
figure, imshow([Y Cb Cr])

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

[R1, G1, B1]= RGBchan(lenna_tr1);
[H1, S1, I1] = rgb2HSI(R1,G1,B1);
[R2, G2, B2]= RGBchan(lenna_tr2);
[H2, S2, I2] = rgb2HSI(R2,G2,B2);
[R3, G3, B3]= RGBchan(lenna_tr3);
[H3, S3, I3] = rgb2HSI(R3,G3,B3);

figure, imshow([H1,S1,I1 ;H2,S2,I2 ;H3,S3,I3 ])
%% 3. Geometric Transform

Tr_x = 100;
Tr_y = 200;
Sc_x = 2;
Sc_y = 3;
Ang = pi;

Tr = [1 0 Tr_x; 0 1 Tr_y; 0 0 1];
S = [Sc_x 0 0; 0 Sc_y 0; 0 0 1];
Rot = [cos(Ang) sin(Ang) 0 ; -sin(Ang) cos(Ang) 0 ; 0 0 1];
%% 4. Image down and up Sampling
% use pyramid and bell kernels
% down sample a test image of 256*256 to an image of 128* 128
% up-sample the down smapled image to 384*384

pyramid = 1/4 * [1 2 1; 2 4 2; 1 2 1];
bell = 1/16*[1 3 3 1;3 9 9 3; 3 9 9 3; 1 3 3 1];

%% 5. Quantization

% % input image(gray)
% file=fopen(fullfile('Test_images\Color_baboon_256x256.raw'),'rb');
% lena=fread(file,fliplr([256, 256]),'*uint8')';
% fclose(file);
% % input image(color)
% file=fopen(fullfile('Test_images\Color_lenna_256x256'),'rb');
% lena=fread(file,fliplr([256, 256*3]),'*uint8')';
% fclose(file);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%%%%% source code %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % output image (at once)
% figure, title('Problem_1_H1_result');
% imshow(uint8(H1_result));
% figure, title('Problem_2_R');
% imshow(uint8(lena_R));
% figure, title('Problem_2_G');
% imshow(uint8(lena_G));
% figure, title('Problem_2_B');
% imshow(uint8(lena_B));