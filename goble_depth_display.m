MIN_HV_DIST = 0.01; 
SKY_DEPTH = 2*(1 ./ MIN_HV_DIST);
[imh,imw]=size(labim);
stats = regionprops(labim,  'PixelIdxList');
imvg=labim;
imvg(stats(3, 1).PixelIdxList)=2;
dx = uint8(imvg ~= imvg(:,[2:end end]));
dy = uint8(imvg ~= imvg([2:end end],:));
boundaries=dx|dy;
%v0=vp{3,1}(2)/imh;
v0=0.5 ;
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
clearg=unique(wseg(stats(1, 1).PixelIdxList));
clears=unique(wseg(stats(3, 1).PixelIdxList));
meandepth([clearg clears])=[];
vardepth([clearg clears])=[];


figure;
for i=1:length(meandepth)
a=meandepth(i);sigma=vardepth(i); % ��ֵa=-6
y=0:0.00001:10;

x=(1/((sqrt(2*pi))*sigma))*exp(-((y-a).^2)/(2*sigma.^2));
plot(x,y,'b','LineWidth',1.5);
hold on; 
end



