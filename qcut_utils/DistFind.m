function ALL_SIM=DistFind(GVMean,max_label)
ss = 0.28; % anything above 0.22 doesn't work
GVMean = mat2gray(GVMean);
M = squareform(pdist(GVMean'));
ALL_SIM =exp(-M.^2 ./ (2*ss^2))*1/ss^2*1/2*1/sqrt(pi);

mu = 0;
ss = 0.28;%1025825
sigma = ss;

%ALL_SIM = normpdf(M, mu, sigma);

%ALL_SIM=zeros(max_label,max_label);
%dist_func = @(x) 1./(0.00001+x);
%dist_func = @(x) exp(-(x)/(2*(ss)));

%for label_counter=1:max_label   
%    distance_all=(GVMean(label_counter:end)-GVMean(label_counter)).^2;
%    ALL_SIM(label_counter:end,label_counter)=dist_func(distance_all);%1./(0.00001+distance_all);
%    ALL_SIM(label_counter,label_counter:end)=dist_func(distance_all);%1./(0.00001+distance_all);
%end

