% Authors: Mariah Costa, Peter Niedbalski
% Last edited: 4/18/2021

% Version of code used in MC ISMRM 2021 project
% This code prompts user for list and raw files for RADIAL KOOSHBALL data,
% subject health status, acqusition style, and noise threshold factor.

% Computes key images via both linear and interleaved methods, extracts attenuation
% maps, corrected images, C1, C1 uncertainty, and angles

% Common abbreviations found in variables:
% "H" used interchangeably with "L" = Linear method to calculate C1 ie halved ratio
% "I" - interleaved method to calculate C1

%%
clear all;clc;close all; tic
je= jet; % set color scale for figures
je(1,:) = [0 0 0];
set(0,'DefaultFigureVisible','on')

% Set stop condition andf counter
Fin=1;
count = 0;


% While loop will run while user has provided new data. User will be
% prompted to select data set and it will be processed. Then the user will
% be asked if they want to process more data. If no, while loop closes
while Fin == 1
    count=count+1;
    % User selects files, paths
    [sin, sinpath]=uigetfile('*.sin', 'Select SIN file','C:\Users\CosG3M\Documents\MATLAB\Human data\3D Radial healthy and CF');
    [list, listpath]=uigetfile('*.list', 'Select LIST file','C:\Users\CosG3M\Documents\MATLAB\Human data\3D Radial healthy and CF');
    % Get locations of files
    sinloc=strcat(sinpath,sin);
    listloc=strcat(listpath,list);
         % Ask user for subject health status
    tag = questdlg('Subject status','Data tag','Healthy','CF','LAM','Other');
    
    % Ask user for Acquisition style
    AcqStyle = questdlg('Aquisition style?','Ordering','Golden Means','Other','Other');
    
    % Set reconstruction variables dependent on acquisition style
    switch AcqStyle
        case 'Golden Means'
            A=80;
            B=0.1;
            DataforExport=1;
        case 'Other'
            A=0;
            B=0;
            DataforExport=0;
    end
        
    % Default
    Dim1_L_factor = 12;
    Dim1_R_factor = 9; 
    Dim2_L_factor = 12;
    Dim2_R_factor = 14;
    Dim3_T_factor = 14;
    Dim3_B_factor = 10;
    
    prompt = {'Enter subject number:'};
    dlgtitle = 'Subject number';
    dims = [1 35];
    definput = {'307'};
    SubjNum = inputdlg(prompt,dlgtitle,dims,definput);
    SubjNum = string(SubjNum);
    SubjNum= str2double(SubjNum);  
    
    
    % Ask user if using both left and right or just left - (just right
    % can be added later but no need right now)
    % Ask user ROI style
    % FIXME - this can be better optimized, I'm just not sure how yet.
    % Ideally, the threshold for noise vs lungs should be
    % threshold = 2 * standard deviation of the noise
    % But in practicality this resulted in high variablity of the quality of
    % masks. So right now the user has to find the optimal threshold factor
    % (ie 2, 3, 10 etc) to make the best mask for each dataset
   
    
    % Recon data
    [Raw,Recon] = philips_xenon_kooshball_recon(A,B,DataforExport,listloc,sinloc);
    
    % Pull original reconstructed image, keyhole images, and number of
    % projections from philips_xenon_kooshball_recon script
    C1_H=Recon.Half_Decay;
    C1_I=Recon.Int_Decay;
    OG_image=abs(Recon.Image);
    npro=Raw.npro;
    % Display originally reconned image
    figure
    imslice(OG_image)
    title('OG image')
    Imsize=size(C1_H);
    % Save workspace
%     save('WorkSpace.mat')
%     save_all_figs('FileType','-djpeg')
    
    %%
    % Ask user for noise threshold factor - dependent on individual dataset
    prompt = {'Enter theshold 1:'};
    dlgtitle = 'Threshold Input';
    dims = [1 35];
    definput = {'5'};
    factor_thresh = inputdlg(prompt,dlgtitle,dims,definput);
    factor_thresh = string(factor_thresh);
    factor_thresh = str2double(factor_thresh);
    
    
    %%
    
    % Preset size of variables to be filled later to speed code
    C1_H_vals_i = zeros(Imsize);
    C1_I_vals_i = zeros(Imsize);
    C1_H_vals_j = zeros(Imsize);
    C1_I_vals_j = zeros(Imsize);
    C1_H_vals_k = zeros(Imsize);
    C1_I_vals_k = zeros(Imsize);
    
    AttMap_H_vals_i = zeros(Imsize);
    AttMap_I_vals_j = zeros(Imsize);
    AttMap_H_vals_j = zeros(Imsize);
    AttMap_I_vals_k = zeros(Imsize);
    AttMap_H_vals_k = zeros(Imsize);
    AttMap_I_vals_k = zeros(Imsize);
    
    noise_vals_i = zeros(Imsize);
    noise_vals_j = zeros(Imsize);
    noise_vals_k = zeros(Imsize);

    SigValsOG_i = zeros(Imsize);
    SigValsOG_j = zeros(Imsize);
    SigValsOG_k = zeros(Imsize);
    %% Calculate Attenuation map, Corrected image - PN
    AttMap_H = 1 - (1/npro).*(1-C1_H.^npro)./(1-C1_H);
    AttMap_I = 1 - (1/npro).*(1-C1_I.^npro)./(1-C1_I);
    AttMap_I(isnan(AttMap_I)) = 0;
    
    % Correct image via attenuation map
%     CorrIm_H = OG_image./(1-AttMap_H);
%     CorrIm_I = OG_image./(1-AttMap_I);
    
    %% Find C1 uncertainty for each slice
    
    for i = 1:Imsize(1) % for every slice
        
        Image_i=squeeze((OG_image(i,:,:)));
        C1_H_i=squeeze(abs(C1_H(i,:,:)));
        C1_I_i=squeeze(abs(C1_I(i,:,:)));
        AttMap_H_slice_i = squeeze(abs(AttMap_H(i,:,:)));
        AttMap_I_slice_i = squeeze(abs(AttMap_I(i,:,:)));       
        %% Create noise mask
        ImSize2D = size(Image_i);
        
        % Create regions of interest (ROIs) in which to sample background noise
        % Dependednt on data set - some images are skewed to left or right
        
        % If using both left and right ROI
            % PShape 1 and 2 are rectangle ROIs on left and right of image
            PShape1 = poly2mask([2 round(ImSize2D(1)/Dim1_L_factor) round(ImSize2D(1)/Dim1_L_factor) 2 2],[round(ImSize2D(1)/Dim1_R_factor) round(ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor)],ImSize2D(1),ImSize2D(1));
            PShape2 = poly2mask(ImSize2D(1)-[2 round(ImSize2D(1)/Dim1_L_factor) round(ImSize2D(1)/Dim1_L_factor) 2 2],[round(ImSize2D(1)/Dim1_R_factor) round(ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor)],ImSize2D(1),ImSize2D(1));
            Noise_Mask = PShape1 + PShape2;
            
            if SubjNum == 1002 || SubjNum == 155 || SubjNum == 167
            Dim1_L_factor = 9;
            PShape1 = poly2mask([2 round(ImSize2D(1)/Dim1_L_factor) round(ImSize2D(1)/Dim1_L_factor) 2 2],[round(ImSize2D(1)/Dim1_R_factor) round(ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor) round((Dim1_R_factor-1)*ImSize2D(1)/Dim1_R_factor)],ImSize2D(1),ImSize2D(1));
            Noise_Mask = PShape1;

            end
            
            % gather noise values in ROI in one variable
            ROI_Noise_i = abs(Image_i(Noise_Mask==1));
            % Calculate noise threshold
            thresh_i = factor_thresh*std(ROI_Noise_i);
        
        % Create lung mask for slice, eroode and dilate
        mask_i = Image_i;
        mask_i(mask_i<thresh_i) = 0;
        mask_i(mask_i>=thresh_i)=1;
        [mask_i, BNmask_i]= erode_dilate(mask_i,1,1);
        
        % for every pixel in the 2D slice        
        for m = 1:Imsize(1)
            for n = 1:Imsize(1)
                % if the pixel is in the lung
                if BNmask_i(m,n) == 1
                    % store C1, Attenuation map values
                    C1_H_vals_i(i,m,n)=C1_H_i(m,n);
                    C1_I_vals_i(i,m,n)=C1_I_i(m,n);                    
                    AttMap_H_vals_i(i,m,n)=AttMap_H_slice_i(m,n);
                    AttMap_I_vals_i(i,m,n)=AttMap_I_slice_i(m,n);
                    SigValsOG_i(i,m,n) = Image_i(m,n);
                    % if the pixel is background noise
                elseif BNmask_i(m,n) == 0
                    % store noise values in separate variable
                    noise_vals_i(i,m,n)=Image_i(m,n);
                    % FIXME - experiment with noise vals vs NoiseVals
                end    
            end
        end
       
%         sanity check - view one slice's mask, original image, noise mask
        if i == 35
            figure
            imagesc(BNmask_i);
            figure
            imagesc(abs(Image_i));
            figure
            imshow(Noise_Mask)
        end
        %% j 
        
        Image_j=squeeze((OG_image(:,i,:)));
        C1_H_j=squeeze(abs(C1_H(:,i,:)));
        C1_I_j=squeeze(abs(C1_I(:,i,:)));
        AttMap_H_slice_j = AttMap_H(:,i,:);
        AttMap_I_slice_j = AttMap_I(:,i,:);
        
        %% Create noise mask
        ImSize2D = size(Image_j);
        
        % Create regions of interest (ROIs) in which to sample background noise
        % Dependednt on data set - some images are skewed to left or right
        
        % If using both left and right ROI
            % PShape 1 and 2 are rectangle ROIs on left and right of image
            PShape1 = poly2mask([2 round(ImSize2D(1)/Dim2_L_factor) round(ImSize2D(1)/Dim2_L_factor) 2 2],[round(ImSize2D(1)/Dim2_R_factor) round(ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor)],ImSize2D(1),ImSize2D(1));
            PShape2 = poly2mask(ImSize2D(1)-[2 round(ImSize2D(1)/Dim2_L_factor) round(ImSize2D(1)/Dim2_L_factor) 2 2],[round(ImSize2D(1)/Dim2_R_factor) round(ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor)],ImSize2D(1),ImSize2D(1));
            Noise_Mask = PShape1 + PShape2;
            % if using only left ROI
            if SubjNum == 1002 || SubjNum == 155 || SubjNum == 167
                Dim2_L_factor = 8;
                PShape1 = poly2mask([2 round(ImSize2D(1)/Dim2_L_factor) round(ImSize2D(1)/Dim2_L_factor) 2 2],[round(ImSize2D(1)/Dim2_R_factor) round(ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor) round((Dim2_R_factor-1)*ImSize2D(1)/Dim2_R_factor)],ImSize2D(1),ImSize2D(1));
                Noise_Mask = PShape1;
            end
            
        % gather noise values in ROI in one variable
        ROI_Noise_j = abs(Image_j(Noise_Mask==1));
        % Calculate noise threshold
        thresh_j = factor_thresh*std(ROI_Noise_j);
        
        % Create lung mask for slice, eroode and dilate
        mask_j = Image_j;
        mask_j(mask_j<thresh_j) = 0;
        mask_j(mask_j>=thresh_j)=1;
        [mask_j, BNmask_j]= erode_dilate(mask_j,1,1);
        
        % for every pixel in the 2D slice        
        for m = 1:Imsize(1)
            for n = 1:Imsize(1)
                % if the pixel is in the lung
                if BNmask_j(m,n) == 1
                    % store C1, Attenuation map values
                    C1_H_vals_j(i,m,n)=C1_H_j(m,n);
                    C1_I_vals_j(i,m,n)=C1_I_j(m,n);                    
                    AttMap_H_vals_j(i,m,n)=AttMap_H_slice_j(m,n);
                    AttMap_I_vals_j(i,m,n)=AttMap_I_slice_j(m,n);
                    SigValsOG_j(i,m,n) = Image_j(m,n);
                    % if the pixel is background noise
                elseif BNmask_j(m,n) == 0
                    % store noise values in separate variable
                    noise_vals_j(i,m,n)=Image_j(m,n);
                    % FIXME - experiment with noise vals vs NoiseVals
                end    
            end
        end
        
        %         sanity check - view one slice's mask, original image, noise mask
        if i == 35
            figure
            imagesc(BNmask_j);
            figure
            imagesc(abs(Image_j));
            figure
            imshow(Noise_Mask)
        end
        
        %% k
        
        Image_k=squeeze((OG_image(:,:,i)));
        C1_H_k=squeeze(abs(C1_H(:,:,i)));
        C1_I_k=squeeze(abs(C1_I(:,:,i)));
        AttMap_H_slice_k = AttMap_H(:,:,i);
        AttMap_I_slice_k = AttMap_I(:,:,i);
        
        % Create noise mask
        ImSize2D = size(Image_k);
        
        % Create regions of interest (ROIs) in which to sample background noise
        % Dependednt on data set - some images are skewed to left or right
        
            PShape1 = poly2mask([round(ImSize2D(1)/Dim3_T_factor) round(ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor)],[2 round(ImSize2D(1)/Dim3_B_factor) round(ImSize2D(1)/Dim3_B_factor) 2 2],ImSize2D(1),ImSize2D(1));
            PShape2 = poly2mask([round(ImSize2D(1)/Dim3_T_factor) round(ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor)],ImSize2D(1)-[2 round(ImSize2D(1)/Dim3_B_factor) round(ImSize2D(1)/Dim3_B_factor) 2 2],ImSize2D(1),ImSize2D(1));
            Noise_Mask = PShape1 + PShape2;
            
            if SubjNum == 1018
                Dim3_L_factor = 12;
                Dim3_R_factor = 18;
                PShape1 = poly2mask([2 round(ImSize2D(1)/Dim3_L_factor) round(ImSize2D(1)/Dim3_L_factor) 2 2],[round(ImSize2D(1)/Dim3_R_factor) round(ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor)],ImSize2D(1),ImSize2D(1));
                PShape2 = poly2mask(ImSize2D(1)-[2 round(ImSize2D(1)/Dim3_L_factor) round(ImSize2D(1)/Dim3_L_factor) 2 2],[round(ImSize2D(1)/Dim3_R_factor) round(ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor) round((Dim3_R_factor-1)*ImSize2D(1)/Dim3_R_factor)],ImSize2D(1),ImSize2D(1));
                Noise_Mask = PShape1 + PShape2;

            elseif SubjNum == 1015 || SubjNum == 155
                Dim3_T_factor = 5;
                PShape2 = poly2mask([round(ImSize2D(1)/Dim3_T_factor) round(ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor)],ImSize2D(1)-[2 round(ImSize2D(1)/Dim3_B_factor) round(ImSize2D(1)/Dim3_B_factor) 2 2],ImSize2D(1),ImSize2D(1));
                Noise_Mask = PShape2;

            elseif SubjNum == 99
                Dim3_T_factor = 11;
                PShape1 = poly2mask([round(ImSize2D(1)/Dim3_T_factor) round(ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor)],[2 round(ImSize2D(1)/Dim3_B_factor) round(ImSize2D(1)/Dim3_B_factor) 2 2],ImSize2D(1),ImSize2D(1));
                Noise_Mask = PShape1 ;
            elseif SubjNum == 307
                Dim3_T_factor = 5;
                PShape2 = poly2mask([round(ImSize2D(1)/Dim3_T_factor) round(ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor) round((Dim3_T_factor-1)*ImSize2D(1)/Dim3_T_factor)],ImSize2D(1)-[2 round(ImSize2D(1)/Dim3_B_factor) round(ImSize2D(1)/Dim3_B_factor) 2 2],ImSize2D(1),ImSize2D(1));
                Noise_Mask = PShape2;
            end


        % gather noise values in ROI in one variable
        ROI_Noise_k = abs(Image_k(Noise_Mask==1));
        % Calculate noise threshold
        thresh_k = factor_thresh*std(ROI_Noise_k);
        
        % Create lung mask for slice, eroode and dilate
        mask_k = Image_k;
        mask_k(mask_k<thresh_k) = 0;
        mask_k(mask_k>=thresh_k)=1;
        [mask_k, BNmask_k]= erode_dilate(mask_k,1,1);
        
        % for every pixel in the 2D slice        
        for m = 1:Imsize(1)
            for n = 1:Imsize(1)
                % if the pixel is in the lung
                if BNmask_k(m,n) == 1
                    % store C1, Attenuation map values
                    C1_H_vals_k(i,m,n)=C1_H_k(m,n);
                    C1_I_vals_k(i,m,n)=C1_I_k(m,n);                    
                    AttMap_H_vals_k(i,m,n)=AttMap_H_slice_k(m,n);
                    AttMap_I_vals_k(i,m,n)=AttMap_I_slice_k(m,n);
                    SigValsOG_k(i,m,n) = Image_k(m,n);
                    % if the pixel is background noise
                elseif BNmask_k(m,n) == 0
                    % store noise values in separate variable
                    noise_vals_k(i,m,n)=Image_k(m,n);
                    % FIXME - experiment with noise vals vs NoiseVals
                end
            end
        end
        
        %         sanity check - view one slice's mask, original image, noise mask
        if i == 35
            figure
            imagesc(BNmask_k);
            figure
            imagesc(abs(Image_k));
            figure
            imshow(Noise_Mask)
        end
        
        
    end
    
    % Now that all values have been allocated for that slice, find average
    % C1, standard deviation of C1, C1 uncertainty, SNR, and average attenuation
    C1_H = cat(1,C1_H_vals_i,C1_H_vals_j,C1_H_vals_k);
    C1_I = cat(1,C1_I_vals_i,C1_I_vals_j,C1_I_vals_k);
    
    SizeC1_H = nonzeros(C1_H);
    SizeC1_I = nonzeros(C1_I);
    
    SigValsOG = cat(1, SigValsOG_i,SigValsOG_j,SigValsOG_k);
    noise_vals = cat(1,noise_vals_i,noise_vals_j,noise_vals_k);
    
    AvgC1_H = mean(nonzeros(C1_H));
    AvgC1_I = mean(nonzeros(C1_I));

    std_H =std(nonzeros(C1_H));
    std_I =std(nonzeros(C1_I));
 
    UncH = std_H/AvgC1_H*100;
    UncI = std_I/AvgC1_I*100;
 
    SNR_OG = (mean(nonzeros(SigValsOG)))/std(nonzeros(noise_vals));

    %% Prepare data and labels for export
    
    % Label data location
    sin_label=strcat('A',sinloc);
    % Prepare average c1, standard dev, and C1 uncertainty for export
    myarray = [SubjNum AvgC1_H, AvgC1_I, std_H, std_I, UncH, UncI, SNR_OG, factor_thresh];
    DataforExport = num2cell(myarray);
    DataLabels = ["Subject Number","AvgC1_H", "AvgC1_I", "std_H", "std_I", "UncH", "UncI", "SNR (OG)", "thresh factor"];
    DataLabels = cellstr(DataLabels);
    Header = ["Sin Location", "List location","Subject type","Acquisition Style"];
    Header = cellstr(Header);
    
    % Turn data into cells for Excel for each dataset (each count)
    DataLoc=strcat('B',num2str(count+1));
    Cell_sinLoc=strcat('B',num2str(count+1));
    Cell_listloc=strcat('C',num2str(count+1));
    TagLoc=strcat('D',num2str(count+1));
    AcqStyleLoc=strcat('E',num2str(count+1));
    
    Cell_sin_label=cellstr(sin_label);
    Cell_list_label=cellstr(listloc);
    Cell_tag = cellstr(tag);
    Cell_AcqStyle = cellstr(AcqStyle);
    
    % Name Excel file based on the name of the sin file and subject health
    % status - ie Uncertainty _ sin file_ healthy
%     filename = strcat('Unc','_',sin(1:8),'_',tag,'.xlsx');
%     % if multiple data sets
%     if count > 1
%         % Name excel file accordingly
%         filename = 'Unc_MultiFile.xlsx';
%     end
filename = 'Unc_MultiFile.xlsx';
    % Designate folder for Excel file to live
    % FIXME - update so others can easily use
    folder = 'C:\Users\CosG3M\Documents\MATLAB\Human data\C1_Unc Results';
    fullFileName = fullfile(folder, filename);
    
   % Export to excel
    xlswrite(fullFileName,DataforExport,2,DataLoc)
    xlswrite(fullFileName,DataLabels,2,'B1')
    xlswrite(fullFileName,Header,1,'B1')
    
    xlswrite(fullFileName,Cell_sin_label,1,Cell_sinLoc)
    xlswrite(fullFileName,Cell_list_label,1,Cell_listloc)
    xlswrite(fullFileName,Cell_tag,1,TagLoc)
    xlswrite(fullFileName,Cell_AcqStyle,1,AcqStyleLoc)
    
    % Ask user if they want to process more data
    Fin = questdlg('More data?','Continue?', 'Yes','No','Quit');
    % Handle response
    switch Fin
        case 'Yes'
            Fin = 1;
        case 'No'
            Fin = 2;
        case 'Quit'
            Fin=2;
    end
    
end

toc

%% Figures - kept commented as default, user can uncomment if they want
% imslice(C1_half_slices) %choose slice 51
% caxis manual; caxis([.999 1]);
% axis off
% axis square; colorbar
% set(gca,'YTickLabel',[],'XTickLabel',[]);

% imslice(C1_Int_slices) %choose 51
% caxis manual; caxis([.999 1]); %use .999 for exageration - same caxis as halved
% makes it more obvious interleaved has higher variance
% use .94 for more continous image standalone
% axis off
% axis square; colorbar
% set(gca,'YTickLabel',[],'XTickLabel',[]);


% figure
% imagesc(C1_H_slice.*mask1)
% title('Extracted C_1 values via Linear Keyhole')
% caxis manual; caxis([.998 1]);
% colormap(je)
% colorbar
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% 
% figure
% imagesc(C1_I_slice.*mask1)
% title('Extracted C_1 values via Interleaved Keyhole')
% %caxis manual; caxis([.997 1]);
% caxis manual; caxis([.89 1]);
% colormap(je)
% colorbar
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% %
% 
% figure('Name','og image')
% imagesc(Image_slice)
% colormap(gray)
% colorbar
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% 
% figure('Name','Corrected Image Lin')
% imagesc(squeeze(CorrIm_H_slice).*mask1)
% colormap(gray)
% colorbar
% % caxis manual
% % caxis([0 200]);
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% 
% figure('Name','Corrected Image Interleaved')
% imagesc(squeeze(CorrIm_I_slice).*mask1)
% colormap(gray)
% colorbar
% % caxis manual
% % caxis([Imagemin Imagemax]);
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% 
% %% Percent change
% 
% figure('Name','Percent change  Lin')
% imagesc((PC_H).*mask1)
% colormap(je)
% colorbar
% caxis manual; caxis([0 300]);
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% 
% figure('Name','percent change I ')
% imagesc((PC_I).*mask1)
% colormap(je)
% colorbar
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% % caxis manual; caxis([Imagemin Imagemax]);
% 
% %% Display attn map
% AttMap_H_slice=squeeze((AttMap_L(slice,:,:)));
% figure('Name','att map L ')
% imagesc((AttMap_H_slice).*mask1)
% colormap(je)
% colorbar
% axis off; axis square
% set(gca,'YTickLabel',[],'XTickLabel',[]);
% % caxis manual; caxis([Imagemin Imagemax]);


%% Show attenuation maps

% imslice(AttMap_L.*BNmask,'Linear Attenuation Map')
% colormap(je)
% 
% imslice(AttMap_Int.*BNmask,'Interleaved Attenuation Map')
% colormap(je)
% 
% imslice(CorrIm_L.*BNmask,'Linear Corrected Image')
% imslice(CorrIm_Int.*BNmask,'Interleaved Corrected Image')
        
        
        