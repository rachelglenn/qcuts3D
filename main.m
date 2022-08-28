clear; clc;
close all;


addpath utils;
addpath 'qcut_utils';


phase = "/Art.raw.nii.gz";
topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';

% Get patient folders 
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..
% Optional fun : Print folder names to command window.
patientDice = zeros(size(length(subFolderNames)));
disp(subFolderNames);

%Setup parameters
suppix_num =[100, 1000,10000];
%suppix_num =[12000];


otherparams = [0.5, 1.0, 1.5,2.0];
otherparams = [25, 50, 75, 100, 125, 150];
% Prep dice list for patients
filename = sprintf('results/liver_dices.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'paitentID\tJaccard\tDice\tcoupling\tsuperRes\tbits\n');
fclose(fid);
%'1025825'
%subFolderNames = {'1052473','1053514','1053863','1081121','1082819','1084441','1084939','1091400','1096371','1101543','1104289','1106125','1120249','1124373','1128451','1131799','1158513','1193512','1208287','1260522','2082714','2085757','2150717','2212201','372868','539789','653097','699783','846260','867447','876094','927498','930104'};
% Perform image segmentation on images
for k = 1 : length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    
    info = niftiinfo(patient + phase');
    tmp = niftiread(info);
    %art = im2single(mat2gray(tmp));
    art = im2single(tmp);

    disp("PixelDimen");
    disp(info.PixelDimensions);
    
    truth = niftiread(patient + '/Truth.raw.nii.gz');
   
    outdir = append('results/', subFolderNames{k});
    mkdir(outdir);
    %filename = sprintf('%s/test.nii',outdir);
    %niftiwrite(art,filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    diceList = zeros(size(length(art(1,1,1,:))));
    jaccardList = zeros(size(length(art(1,1,1,:))));
    if info.BitsPerPixel == 16

        refImage = getRefImage(7, phase);
        %refImage = truth;
        
        [seg_qcut_out, coupling_out, suppix_out] = applyQCUTv3(art,suppix_num, subFolderNames{k}, otherparams);
        truth = int16(truth);
        seg_qcut = int16(seg_qcut_out);
    else
       [seg_qcut_out, coupling_out, suppix_out] = applyQCUTv3(art, suppix_num, subFolderNames{k}, otherparams);
       truth = single(truth);
       seg_qcut = single(seg_qcut_out);
    end 

    % Find the maximum dice with the parameters
    for n = 1 : length(seg_qcut(1,1,1,:))
        diceList(n) = generalizedDice(truth,seg_qcut(:,:,:,n));
        jaccardList(n) = jaccard(double(truth),double(seg_qcut(:,:,:,n)));
    end

    %disp("DiceList");
    %disp(diceList);
    [dicemax,index] = max(diceList);
    seg_nifit = seg_qcut(:,:,:, index);
    coupling_neigh = coupling_out(index);
    super_pix_res = suppix_out(index);

    info = niftiinfo(patient + phase');
    tmp = niftiread(info);
    art = im2single(mat2gray(tmp));
    % Save off images to display in paper

    disp("Sizes:");
    disp(size(seg_nifit));
    disp(size(art));
  
    %filename = sprintf('%s/liver_dice.txt',outdir);
    filename = sprintf('results/liver_dices.txt');
    fid = fopen(filename,'a+'); 

    fprintf(fid,'%s\t%f\t%f\t%d\t%f\t%d\n', subFolderNames{k}, ...
        jaccardList(index), diceList(index), coupling_neigh, ...
        super_pix_res,info.BitsPerPixel );
    fclose(fid);

%     topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
%     patientRef = topLevelFolder+ "/1025825";
%     phase = "/Ven.raw.nii.gz";
%     info = niftiinfo(patientRef + phase');
%     historef = niftiread(info);
%     %coupling = uint32(1):10:uint32(length( historef(1,1,:)  ));
%     disp(length(historef(1,1,:)) );

    for n = 1 : length(art(1,1,:))
        img = (art(:,:,n));
        imgtruth = truth(:,:,n);
        %filename = sprintf('%s/Before_%ds.png',outdir, n);
        %imwrite(uint8(img),filename); 
        segimg = seg_nifit(:,:,n);
        %img = imlocalbrighten(img,coupling_neigh);
        %grayImRef = im2uint8(mat2gray(historef(:,:,super_pix_res)));
        %img = imhistmatch(img,grayImRef,'method','polynomial');    
        imga = double(segimg);
        imgb = double(imgtruth);
        filename = sprintf('%s/Truth_%d.png',outdir, n);
        imwrite(double(imgtruth),filename);


        rgbImage = imoverlay(mat2gray(img),mat2gray(segimg), [1 0 0]);

        filename = sprintf('%s/Comb_%d.png',outdir, n);
        imwrite(rgbImage,filename);
        
        %filename = sprintf('%s/Pred_%d.png',outdir, n);
        %imwrite(uint8(img),filename);

%         if super_pix_res ~= 1
%             img = histeq(img, super_pix_res);
%         end
        %alpha = 4.0;
        %grayImage = locallapfilt(img, super_pix_res, alpha);

        
        
                   
         %filename = sprintf('%s/Before_%d.png',outdir, n);
         %imwrite(uint8(img),filename); 

         % Save histogram
         %imhist(img);
         %filename = sprintf('%s/Before_hist_%d.png',outdir, n);
         %disp(filename)
         %saveas(gcf,filename);

    end

   

    filename = sprintf('%s/QIS.nii',outdir);
    disp(filename);
    %niftiwrite(segnifit,filename,info, 'Compressed',true);
    %info = niftiinfo(patient + '/Pre.raw.nii.gz');
    %[a, b, c ] =size(segnifit);
    %infoinfo.Filesize = (a*b*c)*31;

%     if info.BitsPerPixel == 16
%         segnifit = zeros(size(art),'int16');
%     else
%         segnifit = zeros(size(art),'single');
%     end 
    info = niftiinfo(patient + phase);
    niftiwrite(seg_nifit,filename, info, 'Version', 'NIfTI1',  'Compressed',true);
end


function refImage = getRefImage(index, phase)

    topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
    patientRef = topLevelFolder+ "/1158513";
    info = niftiinfo(patientRef + phase');
    historef = niftiread(info);
    refImage = im2uint16(historef(:,:,40));
    invImage =  imcomplement(refImage);
    BInv = imcomplement(imreducehaze(invImage));%, 'Method','approx','ContrastEnhancement','boost'));
    
    
    %     Kdouble = double(refImage);                  % cast uint16 to double
%     kmult = 65535/(max(max(Kdouble(:))));
%     nbins = kmult;
%     histtmp = histogram(refImage, nbins);
%     ybin = histtmp.Values;
%     mu = 70;    % data mean
%     sd = 20;     % data std    
%     x_bin = (hg.BinEdges(1:end-1)+hg.BinEdges(2:end))/2;
%     y_pdf = 0.25* 1/(2*pi*sd)*exp(-(x_bin-mu).^2/(2*sd^2));
%     y_bin = ybin.*y_pdf;

%     
%     a = 0.3;%  standard deviation
%     b = 0.5; %  mean
%     ny = 3000;
%     y = a.*randn(ny,1) + b;   % data
%     mu = 70;    % data mean
%     sd = 20;     % data std
%     nbins = round(ny/20);
%     hg = histogram(refImage, nbins);
%     hold on;
%     y_bin = hg.Values;
%     x_bin = (hg.BinEdges(1:end-1)+hg.BinEdges(2:end))/2;
%     y_pdf =0.25* 1/(2*pi*sd)*exp(-(x_bin-mu).^2/(2*sd^2));
%     area_hist = trapz(x_bin, y_bin);
%     area_pdf = trapz(x_bin, y_pdf);
%     h_fit = plot(x_bin,y_pdf*area_hist/area_pdf,'LineWidth',3);
%     legend(h_fit, sprintf('mu %.3f, std %.3f', mu, sd));
%     title(sprintf('Gaussian fit to histogram of %d observations with %d bins', length(y), nbins));
%     disp(pd);
end
