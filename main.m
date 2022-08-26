clear; clc;
close all;


addpath utils;
addpath 'qcut_utils';


% input_loc = 'Data_SkullStripped/SkullStripped/';
% output_loc='Output_SkullStripped/';
% for i=1:80
%     path = strcat(input_loc, strcat( int2str(i),'.png'));
%     respath = strcat(input_loc, strcat( int2str(i),'.png'));
%     
%     sprintf( path);
%     im=imread (path);
%     %imshow(im);
%     fuzzyimage(im,i, output_loc);
% end
clc; clear;


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
suppix_num =[10, 50, 2000,12000];
suppix_num =[12000];

otherparams = [0.1, 0.3, 0.4, 0.5];
otherparams = [0.1];
% Prep dice list for patients
filename = sprintf('results/liver_dices.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'paitentID\tDice\t\n');
fclose(fid);
%'1025825'
subFolderNames = {'1052473','1053514','1053863','1081121','1082819','1084441','1084939','1091400','1096371','1101543','1104289','1106125','1120249','1124373','1128451','1131799','1158513','1193512','1208287','1260522','2082714','2085757','2150717','2212201','372868','539789','653097','699783','846260','867447','876094','927498','930104'};
% Perform image segmentation on images
for k = 1 : length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    
    info = niftiinfo(patient + phase');
    tmp = niftiread(info);
    art = im2single(mat2gray(tmp));



    disp("PixelDimen");
    disp(info.PixelDimensions);
    
    truth = niftiread(patient + '/Truth.raw.nii.gz');
   
    outdir = append('results/', subFolderNames{k});
    mkdir(outdir);
    %filename = sprintf('%s/test.nii',outdir);
    %niftiwrite(art,filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    diceList = zeros(size(length(art(1,1,1,:))));
    
    if info.BitsPerPixel == 16
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
        %diceList(n) = generalizedDice(truth,seg_qcut(:,:,:,n));
        diceList(n) = jaccard(double(truth),double(seg_qcut(:,:,:,n)));
        
    end
    disp("DiceList");
    disp(diceList);
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
    fprintf(fid,'%s\t%f\t%d\t%f\n', subFolderNames{k}, dicemax, coupling_neigh, super_pix_res );
    fclose(fid);

    topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
    patientRef = topLevelFolder+ "/1025825";
    phase = "/Ven.raw.nii.gz";
    info = niftiinfo(patientRef + phase');
    historef = niftiread(info);
    %coupling = uint32(1):10:uint32(length( historef(1,1,:)  ));
    disp(length(historef(1,1,:)) );

    for n = 1 : length(art(1,1,:))
        img = art(:,:,n);
        imgtruth = truth(:,:,n);
        %filename = sprintf('%s/Before_%ds.png',outdir, n);
        %imwrite(uint8(img),filename); 
        segimg = seg_nifit(:,:,n);
        img = imlocalbrighten(img,coupling_neigh);
        %grayImRef = im2uint8(mat2gray(historef(:,:,super_pix_res)));
        %img = imhistmatch(img,grayImRef,'method','polynomial');    
        imga = double(segimg);
        imgb = double(imgtruth);
        filename = sprintf('%s/Truth_%d.png',outdir, n);
        imwrite(double(imgtruth),filename);

        A = segimg;
        %out_matrix = A;
        B = im2bw(double(img), graythresh(double(img)));
        %imshow(A, []);
        A=A-min(A(:)); % shift data such that the smallest element of A is 0
        A=A/max(A(:)); % normalize the shifted data to 1    
        B = B -min(B(:));
        B = B/max(B(:));
        B = imadjust(B,[0.2 0.8]);

        rgbImage = imoverlay(uint8(img),double(A), [1 0 0]);

        filename = sprintf('%s/Comb_%d.png',outdir, n);
        imwrite(rgbImage,filename);
        
        filename = sprintf('%s/Pred_%d.png',outdir, n);
        imwrite(uint8(img),filename);

%         if super_pix_res ~= 1
%             img = histeq(img, super_pix_res);
%         end
        %alpha = 4.0;
        %grayImage = locallapfilt(img, super_pix_res, alpha);

        
        
                   
         filename = sprintf('%s/Before_%d.png',outdir, n);
         imwrite(uint8(img),filename); 
         imhist(img);
         filename = sprintf('%s/Before_hist_%d.png',outdir, n);
         disp(filename)
         saveas(gcf,filename);

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


