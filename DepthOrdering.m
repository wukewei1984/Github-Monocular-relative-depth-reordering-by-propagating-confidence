
addpath('C:\Users\hfu\Desktop\2015-2017 Monocular Depth Ordering\TCSVT begin at 2016 May-5G\牛玉婷 4.3 资料全\大论文\第三章\代码\代码\code\iccv07Final\src')
addpath('C:\Users\hfu\Desktop\2015-2017 Monocular Depth Ordering\TCSVT begin at 2016 May-5G\牛玉婷 4.3 资料全\大论文\第三章\代码\代码\code\iccv07Final\src\bp')
bnd_dir='C:\Users\hfu\Desktop\2015-2017 Monocular Depth Ordering\TCSVT begin at 2016 May-5G\牛玉婷 4.3 资料全\大论文\第三章\代码\中间变量及结果\bnd_mat\';
%seg_dir='F:\360浜戠洏\鐮旂┒鐢?713\2015 nyt\region_merge\dataset\segmentation\';

seg_dir='C:\Users\hfu\Desktop\2015-2017 Monocular Depth Ordering\TCSVT begin at 2016 May-5G\牛玉婷 4.3 资料全\大论文\第三章\代码\中间变量及结果\wseg_mat\';
tmp=dir([seg_dir,'*.mat']);
filename={tmp(:).name};
tmp1=dir([bnd_dir,'*.mat']);
filename1={tmp1(:).name};

for ni=1:1%length(filename)  
    %load(im_mat);% im
%     bnd_dir='D:\360浜戠洏\瀛︿範\褰撳墠宸ヤ綔\dataset\mydataset\images\';
% %seg_dir='F:\360浜戠洏\鐮旂┒鐢?713\2015 nyt\region_merge\dataset\segmentation\';
% 
% seg_dir='D:\360浜戠洏\瀛︿範\褰撳墠宸ヤ綔\dataset\mydataset\segmentations\';
tmp=dir([seg_dir,'*.mat']);
filename={tmp(:).name};
tmp1=dir([bnd_dir,'*.mat']);
filename1={tmp1(:).name};
    
    load([seg_dir filename{1,ni}]);
    load([bnd_dir filename1{1,ni}]);


  
%fileName = ['bndinfo' num2str(i) '.mat'];
%load(fileName);
%fileName1 = ['wseg' num2str(i) '.mat'];
%load(fileName1);   
%load('bndinfo3');
%load('wseg3');
%load('/home/nyt/code/data_multi_object/bndinfo79073');%bndinfo_5000121;bndinfo368078;bndinfo79073
%load('/home/nyt/code/gobal_depth/wseg79073');%wseg;wseg367087;wseg79073


%% Set some parameters

MIN_HV_DIST = 0.01;  % minimum distance between horizon and object contact 
SKY_DEPTH = 2*(1 ./ MIN_HV_DIST);


%% Get pixel images for ground, vertical, and sky

[tmp, glabels] = max(bndinfo.result.geomProb, [], 2); 
glabels((glabels>=2) & (glabels<=4)) = 2;
glabels(glabels==5) = 3;


labim = glabels(bndinfo.wseg);
vim = labim==2;
gim = labim==1;
sim = labim==3;
if (sum(gim(:))==0)
  clc;
  clear all;
  continue;
end
[imh, imw] = size(vim);

%% Make sure that horizon is above ground and below sky
scol = sum(sim, 2)>0;
gcol = sum(gim, 2)>0;
%v0 =0.5;
%v0 = (1-v0) * imh;

%try
if sum(scol)>0
minv = find(scol);  minv = minv(end-1);
maxv = find(gcol);  maxv = maxv(2);
v0=minv+0.2*(maxv-minv);
else
maxv = find(gcol);  maxv = maxv(2);  
v0=0.9*(maxv);  
end

v0 = 1 - v0/imh;
%disp(num2str(v0))

%figure(1), imagesc(cat(3, vim, gim, sim)), axis image


%% Get possible ground-vertical boundary pixels

statswseg = regionprops(bndinfo.wseg, 'BoundingBox','Centroid', 'PixelIdxList');
bbox = vertcat(statswseg.BoundingBox);
idx = {statswseg.PixelIdxList};

% bbox = [x1 y1 x2 y2]
bbox = [bbox(:, 1) bbox(:, 2) bbox(:, 1)+bbox(:,3) bbox(:,2)+bbox(:,4)];

vfilt = [ones(3, 1) ; zeros(3, 1)];
gfilt = 1-vfilt; 
boundaryim = (imfilter(double(vim), vfilt)==sum(vfilt(:))) & ...
    (imfilter(double(gim), gfilt)==sum(gfilt(:)));
boundaryim(end, :) = vim(end, :);
boundaryim = imdilate(boundaryim, ones(3, 1));

%figure(3), imagesc(boundaryim), axis image


bpts = cell(bndinfo.nseg, 1);
edges = bndinfo.edges.indices;
spLR = bndinfo.edges.spLR;
eim = zeros(bndinfo.imsize);
for k = 1:numel(edges)
    sp1 = spLR(k, 1);
    sp2 = spLR(k, 2);
    bndind = edges{k}(boundaryim(edges{k}));
    if ~isempty(bndind)
        if glabels(sp1)==2 % && glabels(sp2)==1
            bpts{sp1} = [bpts{sp1} ; bndind];
        elseif glabels(sp2)==2 % && glabels(sp1)==1
            bpts{sp2} = [bpts{sp2} ; bndind];
        end
    end
    eim(edges{k}) = 1;
end

%%
meano=zeros(bndinfo.nseg,1);
for i=1:bndinfo.nseg
    if isempty(bpts{i,1})
         continue;
    end
    [yo,xo] = ind2sub(bndinfo.imsize, bpts{i,1});
 
    % yo=sort(yo);
     %yoo=yo([round(length(yo)*0.8):1:end]);
     for k=1:length(yo)
     bzo{i,1}(k) = (1 ./ max(v0 - (imh-double(yo(k)))/imh, MIN_HV_DIST));
     
     end
     %meano(i,1)=0.5*(mean(bzo{i,1})-min(bzo{i,1}))+min(bzo{i,1});
     bzo{i,1}=sort(bzo{i,1});
     varo(i,1)=var(bzo{i,1});
     bzoo{i,1}=bzo{i,1}([1:1:round(length(yo)*0.2)]);
     meano(i,1)=mean(bzoo{i,1});
end


%%
stats = regionprops(labim,  'PixelIdxList');
if sum(sim(:))==0
    imvg=labim;
else
    imvg=labim;
    imvg(stats(3, 1).PixelIdxList)=2;
end
dx = uint8(imvg ~= imvg(:,[2:end end]));
dy = uint8(imvg ~= imvg([2:end end],:));
boundaries=dx|dy;

index=find(boundaries);
[yv,xv] = ind2sub(size(labim), index);

bz = zeros(numel(yv), 1);
for k=1:length(yv)
 bzv(k) = (1 ./ max(v0 - (imh-yv(k))/imh, MIN_HV_DIST)); 
end
imv=labim==2;
A=zeros(size(labim));
for i=1:length(xv)
    j=xv(i);
A(:,j)=bzv(i);
end
depth=zeros(size(labim));
depth(stats(2, 1).PixelIdxList)=A(stats(2, 1).PixelIdxList);
[yg,xg] = ind2sub(size(labim), stats(1, 1).PixelIdxList);
ygmin=min(yg); ygmax=max(yg);
yg1=[ygmin:1:ygmax]';
bzg = zeros(numel(yg1), 1);
for k=1:length(yg1)
   bzg(k) = (1 ./ max(v0 - (imh-yg1(k))/imh, MIN_HV_DIST));  
end
A=zeros(size(labim));
for i=1:length(yg1)
    j=yg1(i);
A(j,:)=bzg(i);
end

depth(stats(1, 1).PixelIdxList)=A(stats(1, 1).PixelIdxList);
%depth(stats(3, 1).PixelIdxList)=SKY_DEPTH;
indnear=find(depth==0);
depth(indnear)=1;
%depth(stats(3, 1).PixelIdxList)=0;
figure;imagesc(depth);
saveas(gcf,[num2str(ni) '_hypdepth.jpg']);
close(gcf);

%%

for i=1:bndinfo.nseg
    ind=find(wseg==i);
    [yd,xd] = ind2sub(bndinfo.imsize, ind);
    indd=find(yd==ceil(statswseg(i, 1).Centroid(2))|yd==floor(statswseg(i, 1).Centroid(2)));
    inddd=ind(indd);
    cendepth=depth(inddd);
    segdepth=depth(ind);
    meandepth(i)=mean(segdepth);
    vardepth(i)=std(cendepth);
end
meandepth=meandepth';
vardepth=vardepth';
segindex=[1:1:numel(meandepth)]';

if sum(sim(:))>0
    clears=unique(wseg(stats(3, 1).PixelIdxList));
else
   clears=[];
end
clearg=unique(wseg(stats(1, 1).PixelIdxList));
segindex([clearg;clears])=[];%segment conrespondonse new ordering in meandepth
meandepth([clearg;clears])=[];
meano([clearg;clears])=[];
vardepth([clearg;clears])=[];
mindepth=floor(min(meano,[],1));
maxdepth=ceil(max(meano,[],1));
m=(maxdepth-mindepth)/9.0;

x=zeros([length(meandepth),10]);
for i=1:length(meandepth)
    if meano(i)~=0
        a=meano(i);
    else
        a=meandepth(i);
    end

sigma=vardepth(i); 
y=mindepth:m:maxdepth;
%y=1:1:10;
x(i,:)=(1/((sqrt(2*pi))*sigma))*exp(-((y-a).^2)/(2*sigma.^2));
%for j=1:10
total=sum(x(i,:),2);
tot=repmat(total,1,10);
x(i,:)=x(i,:)./(tot+eps);
x(i,:)=roundn(x(i,:),-4) ;
if sum(x(i,:))==0
    x(i,round((a-mindepth)/m)+1)=1;
end
%x(i,j)=normcdf(y(j+1),a,sigma)-normcdf(y(j),a,sigma);
%end

%end
end



    
%% 
bptsex=cellfun('length',bpts)>0;
bptsind=bptsex(segindex);
posind=sum(x~=0,2)==1;

% dep_dir={'C:\Users\hfu\Desktop\2015-2017 Monocular Depth Ordering\TCSVT begin at 2016 May-5G\牛玉婷 4.3 资料全\大论文\第三章\代码\代码\code\gobal_depth\'}
load(['mij00']);
load(['mij01']);
ctol = 0.001;  % convergence tolerance
T = 0.5;       % annealing temperature
maxiter = Inf; % max iter for maxBeliefPropBethe
nnodes=length(meandepth);
for i=1:nnodes
    factor2var{i,1}=[i];
    if (bptsind(i)~=1)
    %if ((posind(i)==1&&bptsind(i)~=1)||sum(x(i,:),2)==0)
     % if sum(x(i,:),2)==0
        factors{i,1}=repmat(0.1,10,1);
    else
%         factors{i,1}=0.8*-log(x(i,:)+eps)'+repmat(0.2/10,10,1);
     factors{i,1}=0.8*(x(i,:))'+repmat(0.2/10,10,1);
    end
    
    %factors{i,1}=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';
end


lab = bndinfo.edges.boundaryType;
lab = lab(1:end/2) + 2*lab(end/2+1:end);
n=nnodes;
for i=1:nnodes
    t=i+1;
    %if (i==5||i==6||i==9)
    %       continue;
     %  end
    for j=t:nnodes
      % if (j==5||j==6||j==9)
       %     continue;
       %end
        seg1=segindex(i);
        seg2=segindex(j);
        all_fg=[];
        for k = 1:numel(lab)
            if (bndinfo.edges.spLR(k, 1)==seg1&&bndinfo.edges.spLR(k, 2)==seg2...
                    ||bndinfo.edges.spLR(k, 1)==seg2&&bndinfo.edges.spLR(k, 2)==seg1)
                fgind=lab(k);
                all_fg(k)=bndinfo.edges.spLR(k, fgind);
            end
        end
        a_fg=all_fg(find(all_fg~=0));
        if (numel(unique(a_fg))==1&&unique(a_fg)==seg1)
             n=n+1;
            factor2var{n,1} = [i j]; 
%            if(bptsex(i)==1&&posind(i)==1)
%               m=find(x(i,:)');
%               b=1.0/(10-m);1
%              mo=zeros(10,10);
%               for k=m+1:10
 %                mo(m,k)=b;
  %            end
  %             factors{n,1}=mo;
   %        else
            factors{n,1} = mij00; 
   %         end
        elseif (numel(unique(a_fg))==1&&unique(a_fg)==seg2)
             n=n+1;
             factor2var{n,1} = [i j]; 
 %         if(bptsex(j)==1&&posindmatlab.mat(j)==1)    
 %                m=find(x(j,:)');
 %               b=1.0/(10-m);
 %               mo=zeros(10,10);
 %              for k=m+1:10
 %                   mo(k,m)=b;
 %               end
 %               factors{n,1}=mo;
 %          else
            factors{n,1} = mij01; 
 %         end
               
        end
        
                
    end
end
    

nvals = [10*ones(nnodes, 1)]; % number of states for each node (variable)

[vals_max, bel_max] = maxBeliefPropBethe(factors, factor2var, nvals, ctol, T, maxiter);
bel_max = reshape(cell2mat(bel_max), [nvals(1) numel(nvals)])';
%% 
for i=1:10
     in=vals_max==i; 
      c=find(in);
     if ((sum(in)>1)&&(~any(meano(c)==0))) 
         num=numel(c);
         Ind=find(vals_max>i);
         vals_max(Ind)=vals_max(Ind)+num-1;
         [C,I]=sort(meano(c));
         vals_max(c(I))=[i:1:i+num-1]';
     end
 end
 ordering=zeros(numel(vals_max),1);
 [V,Index]=sort(vals_max);    
 ordering(Index)=[1:1:numel(vals_max)];
B=zeros(size(wseg));
stats = regionprops(wseg,  'PixelIdxList');
for i=1:length(segindex)
    t=segindex(i);
    B(stats(t, 1).PixelIdxList)=ordering(i);
end
save(['F:\test_vp_data\',dataname],'B');
figure;imagesc(B);
saveas(gcf,[num2str(ni) '_ourorder.fig']);
saveas(gcf,[num2str(ni) '_ourorder.jpg']);
close(gcf);

clear;
close all;
end
