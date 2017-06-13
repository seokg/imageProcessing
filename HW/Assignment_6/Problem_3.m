function Problem_3()
    close all;
    clear all;
   
    
    imgdir = uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\cameraman_gray_256x256.raw'),'rb');
    
    camera = fread(file,fliplr([256,256]),'*uint8')';
    fclose(file);
    figure, imshow(camera,[]);title('bridge orignal image');
%% Differntial Pulse Coding Modulation
%% 3.1 DPCM with quantization parameter of 4, 8 , and 16
[DPCMImg4, nonquanterror4] = encodeDPCM(camera, 4 );
[DPCMImg8, nonquanterror8] = encodeDPCM(camera, 8);
[DPCMImg16, nonquanterror16] =encodeDPCM(camera,16);

figure, imshow(DPCMImg4(2:256,2:256),[]);title('DPCM Image 4');
figure, imshow(DPCMImg8(2:256,2:256),[]);title('DPCM Image 8');
figure, imshow(DPCMImg16(2:256,2:256),[]);title('DPCM Image 16');

reconImg4 = decodeDPCM( DPCMImg4, 4 );
reconImg8 = decodeDPCM( DPCMImg8, 8 );
reconImg16 = decodeDPCM( DPCMImg16, 16 );
figure, imshow(reconImg4(1:256,1:256),[]);title('reconstructed Image 4');
figure, imshow(reconImg8(1:256,1:256),[]);title('reconstructed Image 8');
figure, imshow(reconImg16(1:256,1:256),[]);title('reconstructed Image 16');

%% 3.2 PSNR of the reconstructed image

psnr_Q4 = PSNR(reconImg4, double(camera));
psnr_Q8 = PSNR(reconImg8, double(camera));
psnr_Q16 = PSNR(reconImg16, double(camera));
disp('psnr value for quantization parameter 4')
disp(psnr_Q4);
disp('psnr value for quantization parameter 8')
disp(psnr_Q8);
disp('psnr value for quantization parameter 16')
disp(psnr_Q16)

%% 3.3 Mean Square Value of coding error of u(n) and e(n)
MSE_u_Q4 = MSE(reconImg4, double(camera));
MSE_u_Q8 = MSE(reconImg8, double(camera));
MSE_u_Q16 = MSE(reconImg16, double(camera));
disp('u(n) MSE value for quantization parameter 4')
disp(MSE_u_Q4);
disp('u(n) MSE value for quantization parameter 8')
disp(MSE_u_Q8);
disp('u(n) MSE value for quantization parameter 16')
disp(MSE_u_Q16)
MSE_e_Q4 = MSE(DPCMImg4,nonquanterror4);
MSE_e_Q8 = MSE(DPCMImg8,nonquanterror8);
MSE_e_Q16 = MSE(DPCMImg16,nonquanterror16);
disp('e(n) MSE value for quantization parameter 4')
disp(MSE_e_Q4);
disp('e(n) MSE value for quantization parameter 8')
disp(MSE_e_Q8);
disp('e(n) MSE value for quantization parameter 16')
disp(MSE_e_Q16)

end

function output = quantization(image, bit)
%% quatization
% image: input image
% bit: bit number
    level = 2^bit;
    step = 256*2/level;
    start = -256 + step/2;
    startQuant = -256;
    [row,col]= size(image);
    
    for i = 1: level
         for j = 1: row
            for k = 1: col
                
                if image(j,k) < (startQuant+step) && image(j,k) >= startQuant
                    image(j,k) = start;
                end

            end
         end       
        start = start + step;
        startQuant = startQuant +step;
    end

    output = image;



end

function [output,outputnonquant] = encodeDPCM( input, quant )
    input = double(input);
    output = zeros(size(input,1),size(input,2));
    output(1,:) = input(1,:);
    output(:,1) = input(:,1);
  
    input = input;
    
    for i=2:size(input,1)
        for j=2:size(input,2)
                % predict
                ubarprime = 0.95*input(i,j-1)-0.95*0.95*input(i-1,j-1)+0.95*input(i-1,j);
            
                % error 
                error = (input(i,j)-(ubarprime));
                % Quantification
%                 errorprime = round(error/quant);
%                 output(i,j) = errorprime;
                output(i,j) = error;
        end
    end
    
    outputnonquant = output;
    output = quantization(output, quant);
end

function output = decodeDPCM( input, quant )
    input = double(input);
    output = zeros(size(input,1),size(input,2));
    output(1,:) = input(1,:);
    output(:,1) = input(:,1);
    
    for i=2:size(input,1)
        for j=2:size(input,2)
            sum = -0.95*0.95*output(i-1,j-1)+0.95*output(i-1,j)+0.95*output(i,j-1);
%             output(i,j)=input(i,j)*quant+sum;
            output(i,j)=input(i,j)+sum;

        end
    end

end

function output =PSNR(noise, ref)
    
    mse = MSE(noise,ref);
    peaksnr = 10*log10(255.^2/mse);
    output = peaksnr;
end

function output = MSE(a,b)
    output = sum(sum((a - b).^2)) / (size(a,1) * size(a,2));
end

