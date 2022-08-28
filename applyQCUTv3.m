function [output, coupling_out, suppix_out] = applyQCUTv3(image_now, suppix_num_in,patient, coupling)

%image_now = imadjustn(image_now);

agg = [];
i = 1;
coupling_out = zeros(size(coupling).*size(suppix_num_in));
suppix_out = zeros(size(coupling).*size(suppix_num_in));

for suppix_num=suppix_num_in
    for lambda=coupling
        out_image = TransformImage(image_now, lambda);
        [GVMean, suppixel, boundaries,PixNum, LabelLine,width, height,recon]=SolveSlic(out_image,suppix_num,patient);
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
agg = agg>graythresh(image_now*0.85);
output = agg;
%output_mean=mean(agg,4);

end