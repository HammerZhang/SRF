
%Im=imread('./ORL/s01/01.jpg');
%Im=imread('./lena.bmp’);
Im=imread('Pic\\1.bmp');
figure;imshow(Im,[]);%,'border','tight');
I=double(Im);
[m,n]=size(Im);
%对每一幅图像都进行前4*r个svd分解的截取；
r=20;
[UI SI VI]=svds(I,r);
IM=UI*SI*VI';
I=uint8(IM);
figure,imshow(I,[]);%'border','tight');
p=0.6*m*n;
Omega=randsample(m*n,p);
b=I(Omega);
b=double(b)./255;

idx_row=mod(Omega-1,m)+1;
idx_col=ceil(Omega/m);
W = zeros(m,n);
for ii = 1:p
    W(idx_row(ii),idx_col(ii)) = 1;
end

D=zeros(m,n);
D(Omega)=b;
%DU=uint8(D*255);
DU=D*255;
figure,imshow(DU,[]);%'border','tight');

Dp = srf(D,W);  Dp = Dp*255; figure,imshow(uint8(Dp),[]);
err = norm(double(IM)-Dp,'fro')/norm(double(IM),'fro') ;

