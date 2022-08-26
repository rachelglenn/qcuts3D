function [output, coupling_out, suppix_out] = applyQCUTv3(image_now, suppix_num_in,patient, coupling)

%image_now = imadjustn(image_now);

agg = [];


i = 1;
topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
patientRef = topLevelFolder+ "/1025825";
phase = "/Ven.raw.nii.gz";
info = niftiinfo(patientRef + phase');
historef = niftiread(info);
%coupling = uint32(1):10:uint32(length( historef(1,1,:)  ));
disp(length(historef(1,1,:)) );

coupling_out = zeros(size(coupling).*size(suppix_num_in));
suppix_out = zeros(size(coupling).*size(suppix_num_in));

for suppix_num=suppix_num_in
    for lambda=coupling
        % Adjust image
        for n = 1 : length(image_now(1,1,:))
            %grayImage = im2single(mat2gray(image_now(:,:,n)));
            %grayImRef = im2uint8(mat2gray(historef(:,:,lambda)));
            grayImage = imlocalbrighten(image_now(:,:,n),lambda);
          
            %grayImage = locallapfilt(grayImage, single(lambda), 0.5);
            %if lambda ~= 1
            %    grayImage = histeq(grayImage, lambda);
            %end
            %grayImage = imadjust(grayImage,[lambda,0.8]);
            %grayImage = imhistmatch(grayImage,grayImRef,'method','polynomial');

            image_now(:,:,n) = grayImage;
        end
        [GVMean, suppixel, boundaries,PixNum, LabelLine,width, height,recon]=SolveSlic(image_now,suppix_num,patient);
        [neighbourhood,LF,max_label]=FindNeighbours(suppixel);
        
        ALL_DIST=DistFind(GVMean,max_label);
        ALL_DIST=ALL_DIST/max(ALL_DIST(:));
        
        H=AffinityAssign(neighbourhood,LF,ALL_DIST,max_label, lambda);
        
        potentials = zeros(size(GVMean));
        potentials(boundaries) = 1;
        numa = 1000;
      
        H_new=UpdateDiagonal(GVMean,[],H,1:max_label,potentials*numa);
        [SalMap,binarized]=QCUT(H_new,1,PixNum, suppixel,image_now);
        %disp("Binarized");
        %disp(size(binarized));
        agg = cat(4,agg,binarized);
        coupling_out(i) =  lambda;
        suppix_out(i) =  suppix_num;
        i = i+1;
    end
end
%gg = image_now>graythresh(image_now);
%agg = agg>graythresh(image_now*1.05);
output = agg;
%output_mean=mean(agg,4);

end