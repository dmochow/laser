% 05/14/18
% check alignment between BOLD, T1, and segmentation

clear all; close all; clc

%%
TR=2.8;
subjIndx='S02';
basePath=['../data/' subjIndx '/NII/'];
boldFilename=fullfile(basePath,'dsbold_e1_tlrc_al+tlrc.BRIK');
resampledAnatFilename=fullfile(basePath,'resampled_anat+tlrc.BRIK');
resampledSegFilename=fullfile(basePath,'resampled_Classes+tlrc.BRIK');
brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');

NUMBER=3;

switch NUMBER
    
    
    case 1
        %%
        % compare skull stripped anatomical with bold
        bold = BrikLoad (boldFilename);
        anat = BrikLoad (resampledAnatFilename);
        muBold=mean(bold,4);
        slice=40;
        
        figure;
        while 1
            imagesc(muBold(:,:,slice)); colormap bone
            str=input('Hit me','s');
            if strcmp(str,'e')
                break
            else
                imagesc(anat(:,:,slice));
            end
            pause
        end
        
    case 2
        %%
        % compare skull stripped anatomical with segmentation
        seg = BrikLoad (resampledSegFilename);
        anat = BrikLoad (resampledAnatFilename);
        slice=40;
        
        figure;
        while 1
            imagesc(anat(:,:,slice)); colormap bone
            str=input('Hit me','s');
            if strcmp(str,'e')
                break
            else
                imagesc(seg(:,:,slice));
            end
            pause
        end
        
      case 3
        %%
        % compare bold with segmentation
        bold = BrikLoad (boldFilename);
        seg = BrikLoad (resampledSegFilename);
        slice=40;
        
        figure;
        while 1
            imagesc(bold(:,:,slice)); colormap bone
            str=input('Hit me','s');
            if strcmp(str,'e')
                break
            else
                imagesc(seg(:,:,slice));
            end
            pause
        end
        
end
