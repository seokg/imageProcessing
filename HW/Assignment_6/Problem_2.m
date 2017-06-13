function Problem_2()

    close all;
    clear all;
   
    
    imgdir = uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\bridge_gray_256x256.raw'),'rb');
    
    bridge = fread(file,fliplr([256,256]),'*uint8')';
    fclose(file);
    figure, imshow(bridge,[]); title('bridge image original');
%% Bit-Plane 
% change to graycode
graycodeImg = zeros(size(bridge,1),size(bridge,2),8);
graycodeImg = graycodeImage(bridge);
figure, imshow(graycodeImg(:,:,1),[]); title('2.1 Gray coded images bit 7 plane');
figure, imshow(graycodeImg(:,:,2),[]); title('2.1 Gray coded images bit 6 plane');
figure, imshow(graycodeImg(:,:,3),[]); title('2.1 Gray coded images bit 5 plane');
figure, imshow(graycodeImg(:,:,4),[]); title('2.1 Gray coded images bit 4 plane');
figure, imshow(graycodeImg(:,:,5),[]); title('2.1 Gray coded images bit 3 plane');
figure, imshow(graycodeImg(:,:,6),[]); title('2.1 Gray coded images bit 2 plane');
figure, imshow(graycodeImg(:,:,7),[]); title('2.1 Gray coded images bit 1 plane');
figure, imshow(graycodeImg(:,:,8),[]); title('2.1 Gray coded images bit 0 plane');

%% 2.1 apply bit-plane coding to Gray coded bit planes by using the white block skipping (WBS)
% with N = 2,4 and 8
% N = 2
[data2, numzero2 ]=encodeWBS(graycodeImg,1,2);
br2 = bitrate(data2,numzero2,2,1)
compratio2 = 8/br2;
% N = 4
[data4, numzero4 ]=encodeWBS(graycodeImg,2,2);
br4= bitrate(data4,numzero4,2,2)
compratio4 = 8/br4;

% N = 8
[data8, numzero8 ]=encodeWBS(graycodeImg,4,2);
br8 = bitrate(data8,numzero8,4,2)
compratio8 = 8/br8;

disp('compression ratio with N=2')
disp(compratio2)
disp('compression ratio with N=4')
disp(compratio4)
disp('compression ratio with N=8')
disp(compratio8)

% recover the orignal image
% N = 2
decodeImg2 =decodeWBS(data2, 1, 2, 8, 256, 256);
recoverImg2 = binaryImage(decodeImg2);
figure,imshow(recoverImg2, []);title('2.1recovered image N=2');

% N = 4
decodeImg4 =decodeWBS(data4, 2, 2, 8, 256, 256);
recoverImg4 = binaryImage(decodeImg4);
figure,imshow(recoverImg4, []);title('2.1recovered image N=4');

% N = 8
decodeImg8 =decodeWBS(data8, 4, 2, 8, 256, 256);
recoverImg8 = binaryImage(decodeImg8);
figure, imshow(recoverImg8, []);title('2.1recovered image N=8');

%% compare the compression ratio for N = 2,4, and 8
end

function output = de2bin(input)
%% deci to binary number
    output = [];

    while( input > 0)
        remainder = mod(input,2);
        output = [remainder output];
        input = floor(input/2);
%     remainder = mod(input,2);
%     if(input ~= 1 && input~=0)
%         output = de2bin(floor(input/2));
%     end
%     output = [output remainder];
    end
end

function output = bin2de(input)
%% binary to deci number
    len = length(input);
    output = 0;
    idx = 1;
    for i = len:-1:1
        output=output+input(idx)*2^(i-1);
        idx = idx + 1;
    end
end

function output = cvt2GrayCode(input)
%% convert binary to Graycode
%     zero padding for 8 bits input
    if length(input) ~= 8
       numzeropad = 8 - length(input);
       padding = zeros(1,numzeropad);
       input = [padding input];
    end
    
    output = size(input);    
    output(1)=input(1);
    
    for i = 2: length(input)
        output(i) = mod(input(i) +input (i-1),2);
    end
end

function output = cvt2Binary(input)
%% convert  Graycode to binary

    output = size(input);    
    output(1)=input(1);
    for i = 2: length(input)
        output(i) = mod(input(i) +output (i-1),2);
    end

end

function output = graycodeImage(input)
output = zeros(size(input,1),size(input,2),8);
for i = 1: size(input,1)
    for j = 1:size(input,2)
         currGraycode= cvt2GrayCode(de2bin(double(input(i,j))));
         for k = 1:8
              output(i,j,k)= currGraycode(k);
         end          
    end
end
end

function output = binaryImage(input)
output = zeros(size(input,1),size(input,2));
for i = 1: size(input,1)
    for j = 1:size(input,2)
        array = [];
        for k = 1:size(input,3)
            array = [ array input(i,j,k)];
        end
         currBinary= bin2de(cvt2Binary(array));

         output(i,j)= currBinary;
    end
end
end

function [output, numzero] = encodeWBS(inputimage, row, col)
    numzero=0;
    [height, width, channel] = size(inputimage);
    hStep = height / row;
    wStep = width / col;
    output = [];
    for l = 1: channel
    for i = 1: hStep
        for j = 1: wStep
            blockImg = inputimage( (i-1)*row+1 : i*row, (j-1)*col+1 : j*col, l);
            bool = checkZero(blockImg);
            if bool == 0
               output = [output 0];
               numzero = numzero+1;
            elseif bool ~=0
                appendBlock = [];
                for k = 1:row
                appendBlock = [appendBlock blockImg(k,:)];
                end
                output = [output 1 appendBlock];
            end
            
        end
    end
    end
    
end

function output =decodeWBS(inputdata, row, col, channel, width, height)
    dataLen = length(inputdata);
    dataLenPerChannel = dataLen / channel;
    hStep = height / row;
    wStep = width / col;
    output = zeros(height,width,channel);
    n = 1;
    for i = 1 : channel
        for j = 1 : hStep
            for k = 1 : wStep
                if inputdata(n) == 1
                    hetoroData = [];
                    for a = 1: row*col
                        hetoroData = [hetoroData inputdata(n+a)];
                    end
                     hetoroData2 = zeros(row,col);
                    for a = 1: row
                        for b = 1: col
                        hetoroData2(a,b) = hetoroData (b+(a-1)*col);
                        end
                    end
                    
                    output((j-1)*row+1 : j*row, (k-1)*col+1 : k*col, i) =hetoroData2;
                    n = n + row*col + 1; 
                elseif inputdata(n) == 0
                    output((j-1)*row+1 : j*row, (k-1)*col+1 : k*col, i) = zeros(row,col);
                    n = n + 1;
                end
                
            end
        end
    end
end

function output = checkZero(input)
    [x, y ]=size(input);
    
    for a = 1: x
        for b = 1: y
            if input(a,b) ~= 0
                output = 1; % something is none zero
                return
            end
        end
    end
    output = 0; %everything is zero
    
end

function output = bitrate(data, numzero,p,q)
    N = p*q;
    datalen = length(data);
    pzero = numzero/ datalen;
    pnonzero = (datalen- numzero) / datalen;
    
    output = ((1-pzero)*(N+1) + pzero*1) / N;    
end

function output = compRatio(origin, compress)
    output = origin / compress;
end