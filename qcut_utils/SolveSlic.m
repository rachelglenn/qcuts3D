function [GVMean, suppixel,boundaries,PixNum, LabelLine,width, height,recon]=SolveSlic(image_now,suppix_num,patient)
    % Get superpixels
    [suppixel,numlabels] = superpixels3(image_now,suppix_num);
    %imSize = size(image_now);
    %imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
    %for plane = 1:imSize(3)
    %  BW = boundarymask(suppixel(:, :, plane));
    %  % Create an RGB representation of this plane with boundary shown
    %  % in cyan.
    %  imPlusBoundaries(:, :, :, plane) = imoverlay(uint8(image_now(:, :, plane)), BW, 'cyan');
    %end

    %implay(imPlusBoundaries,5)
    %implay(suppixel, 10);
    %filename = sprintf('%s/video_%d__%s.png',outdir,suppix_num, patient);
%     objWrite = VideoWriter('test');
%     open(objWrite);
%     for k = 1:imSize(3)
%       
%        writeVideo(objWrite, imPlusBoundaries(:,:,:,k));
%     end
%     close(objWrite);
    pixelIdxList = label2idx(suppixel);
    %numlabels = max(suppixel(:));
    GVMean = zeros(1,numlabels);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     wavelengthMin = 4/sqrt(2);
%     wavelength = 2.^(0:4) * wavelengthMin;
%     deltaTheta = 45;
%     orientation = 0:deltaTheta:(180-deltaTheta);
%     g = gabor(wavelength,orientation);
%     
%     
%     gabormag = imgaborfilt(Agray,g);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for superpixel = 1:numlabels
        memberPixelIdx = pixelIdxList{superpixel};
        GVMean(superpixel) = mean(image_now(memberPixelIdx));
    end

    LabelLine = suppixel(:) - 1;
    PixNum = numel(image_now);
    width = size(image_now,1); height = size(image_now,2);
    recon = [];
    recon = sup2pixel( numel(suppixel), (suppixel(:)-1), GVMean );
    recon = reshape(recon,size(suppixel));
    
    %Find superpixels on boundaries
    boundaries = [];
    for i=1:size(suppixel,3)
       for j=1:size(suppixel,2)
           sp = suppixel(:,j,i);
           rc = recon(:,j,i);
           boundaries = [boundaries; unique(sp(rc==min(rc(:))))];
       end
    end
    
    close all
    boundaries=unique(boundaries);

