function [a, dfstack2, Q ,cmapidx,amask,translations] = AuBallsetup(Dx,Dy,flag,imageFlag,sizeImg,sizeProb,alpha)
%
% usage: [a1, dfstack2, Q ] = AuBallsetup(Dx,Dy);
%a contaminated image ?????
% dfstack2 stacks croped
% Q ??
% cmapidx - croped id
% amask phobe
% traslations 
%
% image stored in the matrix A
%
% start gendata
if ~exist('imageFlag','var')
    imageFlag=1;
end
scaleP=1;
if imageFlag==1
load GoldBalls2;a = A;clear A;
end
if imageFlag==2
    %a1=mat2gray(imread('cameraman.tif'));
    
    a0=mat2gray(imread('cameraman.tif'));
    a1=double(imresize(rgb2gray(imread('peppers.png')),[256,256]))/255;
    a1=a0.*exp(1j*2*pi*(a1)/2.*scaleP);
    
    a=min(1,abs(a1)).*sign(a1); %???camera man ??????????? ??????
    clear a0, a1;
end

if ~exist('sizeImg','var')
    sizeImg=size(a,1);
end
% sizeImg - to squarre
a = imresize(real(a), [sizeImg sizeImg])+1j*imresize(imag(a), [sizeImg sizeImg]);

%a = a(end/2-128:end/2+127,end/2-128:end/2+127);

%s=0; a=repmat(a,[2^s,2^s]);


if ~exist('flag','var')
    flag=1;
end
% frame size
%
if exist('sizeProb','var')
nfx=sizeProb;
else

nfx = 64;
end

 %nfx nfy sizephobe
nfy = nfx;

% 
% translation
%
%if ~exist('Dx','var')
%    Dx=16;
%end
%if ~exist('Dy','var')
%    Dy=16;
%end

[Nx,Ny]=size(a);
% a contaminate fig
% number of frames
%
nsx=floor(Nx/Dx);
nsy=floor(Ny/Dy);
nframes=nsx*nsy;
%
% generate the illumination
%
[~, ~, rr]=make_grid(480,2);

%[~, ~, rr]=make_grid(320,2);


scale_r=35/3*5/3*sizeProb/64.;



r1 = scale_r*50/180; % ? select a phobe
r2 = scale_r;
%scale_r=35/3*5/1.5;r1 = scale_r*100/180;r2 = scale_r;

w = fftshift(rr>=r1 & rr <= r2);
aa = fftshift(ifft2(w));
aa = aa./max(abs(aa(:)));
amask = cropmat(aa,[nfx nfy]); %approximation! get a phobe from time domain
%%
% probe with %finite number of photons: tflux
tflux = 4.0445e7;

%tflux = 2e7;
amask=amask/norm(amask,'fro')*sqrt(tflux); 

if ~exist('alpha','var')
alpha=0;
end;

[mx,my]=meshgrid(1:size(amask,1),1:size(amask,2));

factor=(mx-size(amask,1)/2).^2+(my-size(amask,2)/2).^2; %factor dist to center
amask=ifft2(fftshift(fft2(amask)).*exp(1j*factor*alpha)); %dist to center - an elapsing zero now

%% out of focus alpha






maskframe = amask;  
%
% set up index maps for overlapping frames 
%
[ix,iy]=ndgrid(0:Dx:Nx-Dx,0:Dy:Ny-Dy);
[xx,yy]=meshgrid(1:nfx,1:nfy); %meshgrid and ndgrid
wrapNx=@(x) mod(x-1,Nx)+1;
wrapNy=@(x) mod(x-1,Ny)+1;
ixt=reshape(ix,1,1,numel(ix)); %get number to (1,1,:)
iyt=reshape(iy,1,1,numel(iy));
cmapidx=wrapNx(bsxfun(@plus,xx,ixt))+(wrapNy(bsxfun(@plus,yy,iyt))-1)*Nx;
% bug: cmapidx=wrapNx(bsxfun(@plus,yy,ixt))+(wrapNy(bsxfun(@plus,xx,iyt))-1)*Nx;
% cropped
%
[xx1,yy1]=ndgrid(floor((Nx-nfx)/2)+1+(0:nfx-1),floor((Ny-nfy)/2)+1+(0:nfy-1));
%cmapidx=wrapNx(bsxfun(@plus,xx1,ixt))+(wrapNy(bsxfun(@plus,yy1,iyt))-1)*Nx;
%

%%
%% hexagonal lattice

if flag==2 %flag lattice method
ix1=floor(nfx/2+(1:Dx:Dx*nsx)-Dx*nsx/2);
iy1=floor(nfy/2+(1:Dy:Dy*nsy)-Dy*nsy/2);

%make 2D lattice of positions
[ix,iy]=meshgrid(ix1,iy1);
% shear by 1/2 step to x every other row to get hexagonal lattice
% ... more uniform density...
xshift=floor(Dx/2)*(mod(1:length(ix1),2));
yshift=floor(Dy/2)*(mod(1:length(iy1),2));
ix=bsxfun(@plus,ix',xshift)';iy=bsxfun(@plus,iy',yshift)';
ixt=reshape(ix,1,1,numel(ix));iyt=reshape(iy,1,1,numel(iy));
cmapidx=wrapNx(bsxfun(@plus,xx1,ixt))+(wrapNy(bsxfun(@plus,yy1,iyt))-1)*Nx;
end

%%random
if flag==3
ix1=floor(nfx/2+(1:Dx:Dx*nsx)-Dx*nsx/2);
iy1=floor(nfy/2+(1:Dy:Dy*nsy)-Dy*nsy/2);

%make 2D lattice of positions
[ix,iy]=meshgrid(ix1,iy1);

% shear by 1/2 step to x every other row to get hexagonal lattice — ixt,iyt
% ... more uniform density...
xshift=floor(Dx/2)*(mod(1:length(ix1),2));
yshift=floor(Dy/2)*(mod(1:length(iy1),2));
ix=bsxfun(@plus,ix',xshift)';iy=bsxfun(@plus,iy',yshift)';
rng('default');
%randx=round(rand(size(ix))*2)-1;randy=round(rand(size(iy))*2)-1;

randx=round(rand(size(ix))*8)-1;randy=round(rand(size(iy))*8)-1; %random disturb 8?

ixt=reshape(ix+randx,1,1,numel(ix));iyt=reshape(iy+randy,1,1,numel(iy));
cmapidx=wrapNx(bsxfun(@plus,xx1,ixt))+(wrapNy(bsxfun(@plus,yy1,iyt))-1)*Nx;
end
%save all

% generate image stack
%imgstack=bsxfun(@times,a1(mapidx),mask1);
%
% generate cropped image stack
%


cropstack=bsxfun(@times,a(cmapidx),maskframe);

%
dfstack2 = abs(myfft2(cropstack)).^2; %cdi stack figures
maskstack = maskframe(:,:,ones(nframes,1));
% save all;
%
% Q is the extraction matrix and Q' is the overlap matrix
%
ncols = Nx*Ny;
nrows = numel(cmapidx);
%Q=[];
Q=sparse(1:nrows, cmapidx(:),(maskstack(:)),nrows,ncols);

translations=[squeeze(ixt) squeeze(iyt)];
%save all;

