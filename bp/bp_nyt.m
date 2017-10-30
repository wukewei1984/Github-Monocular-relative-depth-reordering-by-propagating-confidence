load('/home/nyt/code/data_multi_object/bndinfo79073');%bndinfo_5000121;bndinfo368078;bndinfo79073
load('/home/nyt/code/gobal_depth/wseg79073');%wseg;wseg367087;wseg79073
%load('bndinfo6');
%load('wseg6');

MIN_HV_DIST = 0.01; 
SKY_DEPTH = 2*(1 ./ MIN_HV_DIST);

[tmp, glabels] = max(bndinfo.result.geomProb, [], 2); 
glabels((glabels>=2) & (glabels<=4)) = 2;
glabels(glabels==5) = 3;
labim = glabels(bndinfo.wseg);

[imh,imw]=size(labim);
stats = regionprops(labim,  'PixelIdxList');
imvg=labim;
imvg(stats(3, 1).PixelIdxList)=2;
dx = uint8(imvg ~= imvg(:,[2:end end]));
dy = uint8(imvg ~= imvg([2:end end],:));
boundaries=dx|dy;
%v0=vp{3,1}(2)/imh;
v0=0.5;
%vl=vp{3,1}(1);
index=find(boundaries);
[yv,xv] = ind2sub(size(labim), index);
%indleft=x<=vl;
%xleft=x(indleft);
%yleft=y(indleft);
%indright=x>vl;
%xright=x(indright);
%yright=y(indright);

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
depth(stats(3, 1).PixelIdxList)=SKY_DEPTH;
indnear=find(depth==0);
depth(indnear)=2;
depth(stats(3, 1).PixelIdxList)=0;
imagesc(depth);

for i=1:max(max(wseg))
    ind=find(wseg==i);
    segdepth=depth(ind);
    meandepth(i)=mean(segdepth);
    vardepth(i)=var(segdepth);
end
 meandepth=meandepth';
 vardepth=vardepth';
 segindex=[1:1:numel(meandepth)]';
clearg=unique(wseg(stats(1, 1).PixelIdxList));
clears=unique(wseg(stats(3, 1).PixelIdxList));
segindex([clearg;clears])=[];%segment conrespondonse new ordering in meandepth
meandepth([clearg;clears])=[];
vardepth([clearg;clears])=[];
%segindex([clearg])=[];%segment conrespondonse new ordering in meandepth
%meandepth([clearg])=[];
%vardepth([clearg])=[];
mindepth=floor(min(meandepth,[],1));
maxdepth=ceil(max(meandepth,[],1));
m=(maxdepth-mindepth)/9.0;

for i=1:length(meandepth)
a=meandepth(i);sigma=vardepth(i); % ��ֵa=-6
y=mindepth:m:maxdepth;

x(i,:)=(1/((sqrt(2*pi))*sigma))*exp(-((y-a).^2)/(2*sigma.^2));
%x(i,1)=normcdf(1,a,sigma)
%x(i,2)=normcdf(2,a,sigma)-normcdf(1,a,sigma);
%x(i,3)=normcdf(3,a,sigma)-normcdf(2,a,sigma);
%x(i,4)=normcdf(4,a,sigma)-normcdf(3,a,sigma);
%x(i,5)=normcdf(5,a,sigma)-normcdf(4,a,sigma);
%x(i,6)=normcdf(6,a,sigma)-normcdf(5,a,sigma);
%x(i,7)=normcdf(7,a,sigma)-normcdf(6,a,sigma);
%x(i,8)=normcdf(8,a,sigma)-normcdf(7,a,sigma);
%x(i,9)=normcdf(9,a,sigma)-normcdf(8,a,sigma);
%x(i,10)=1.0-normcdf(9,a,sigma);


end
total=sum(x,2);
tot=repmat(total,1,10);
x=x./(tot+eps);


load('/home/nyt/code/gobal_depth/mij00');
load('/home/nyt/code/gobal_depth/mij01');
ctol = 0.001;  % convergence tolerance
T = 0.5;       % annealing temperature
maxiter = Inf; % max iter for maxBeliefPropBethe
nnodes=length(meandepth);
for i=1:nnodes
    factor2var{i,1}=[i];
    if sum(x(i,:),2)==0
        factors{i,1}=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';
    else
        factors{i,1}=x(i,:)';
    end
    
    %factors{i,1}=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';
end


lab = bndinfo.edges.boundaryType;
lab = lab(1:end/2) + 2*lab(end/2+1:end);
n=nnodes;
for i=1:nnodes
    t=i+1;
    for j=t:nnodes
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
            factors{n,1} = mij00; 
        elseif (numel(unique(a_fg))==1&&unique(a_fg)==seg2)
             n=n+1;
             factor2var{n,1} = [i j]; 
             factors{n,1} = mij01;   
        end
        
                
    end
end
    

nvals = [10*ones(nnodes, 1)]; % number of states for each node (variable)

[vals_max, bel_max] = maxBeliefPropBethe(factors, factor2var, nvals, ctol, T, maxiter);
bel_max = reshape(cell2mat(bel_max), [nvals(1) numel(nvals)])';
B=zeros(size(wseg));
stats = regionprops(wseg,  'PixelIdxList');
for i=1:length(segindex)
    t=segindex(i);
    B(stats(t, 1).PixelIdxList)=vals_max(i);
end
figure;imagesc(B);