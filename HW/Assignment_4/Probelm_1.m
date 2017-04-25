function Probelm_1()
    % template for ee535 DIP
    % Insert the cod in the designated area below
    
    %% loading directory for image files
    imgdir= uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\Test_images\Gray_baboon_256x256.raw'),'rb');
    gray_image = fread(file, fliplr([512,512*3],'*uint8'));
    fclose(file);
    
    
 
    %% -------------- insert code below --------------%%
    output_image = InnerFunction(gray_image); % sample code - delete
    
    figure; imshow (output_image,[]); title('Problem_1.1');
    figure; imshow (output_image,[]); title('Problem_1.1');
    
    % ...



end

function OUTPUT = InnerFunction(INPUT)
    %% functions code

end