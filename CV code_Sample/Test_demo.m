%% Matlad code implementing Chan-Vese model in the paper 'Active Contours Without Edges'(IEEE TIP 2001)
%   This method works well for bimodal images, for example the image 'three.bmp'

clear;
close all;
clc;
%Img=imread('test_images/cups.jpg');
%Img=imread('test_images/three.bmp');      % example that CV model works well
Img=imread('test_images/circle.jpg');      % example that CV model works well
%Img=imread('test_images/vessel.bmp');    % Warning: example image that CV model does NOT work well
%Img=imread('test_images/twoCells.bmp');   % Warning: example image that CV model does NOT work well
% 彩色图转灰度图
if size(Img, 3) > 1
    U = rgb2gray(Img);
else
    U = Img;
end

[nrow,ncol] =size(U); % get the size

% circle
% r=20;
% three
% r=40;
% twoCells
% r=40;
% cups
% r=80;
% vessel
% r=75;
r = 20;
delta_t = 5;
lambda=5;
nu=3;
h = 1;
epsilon = 0.4;
mu = 0.04; %0.01*255*255;

ic=nrow/2;
jc=ncol/2;

phi_0 = sdf2circle(nrow,ncol,ic,jc,r); % construct the signed distamnce function
figure; mesh(phi_0); title('Signed Distance Function')

I=double(U);

% 求g
sigma = 0.5;
% 初始化一个高斯滤波器
G = fspecial('gaussian', 15, sigma);
filted_result = conv2(I, G, 'same');
%求x, y两个方向上的梯度
[grad_x, grad_y] = gradient(filted_result);
grad = grad_x.^2 + grad_y.^2;
g = 1 ./ (1 + grad);

% iteration should begin from here
phi=phi_0;
figure(2);
subplot(1,2,1); mesh(phi);
subplot(1,2,2); imagesc(uint8(I));colormap(gray)
hold on;
plotLevelSet(phi,0,'r');

numIter = 30;
seg_region_old = zeros(size(U));
seg_region_new = zeros(size(U));

for k=1:1000,
    phi = evolution_cv(phi, mu, nu, lambda, delta_t, epsilon, numIter, g);   % update level set function
    if mod(k,2)==0
        pause(.5);
        figure(2); clc; axis equal; 
        title(sprintf('Itertion times: %d', k));
        subplot(1,2,1); mesh(phi);
        subplot(1,2,2); imagesc(uint8(I));colormap(gray)
        hold on; plotLevelSet(phi,0,'r');
        
        if k == 2
            seg_region_old = (phi < 0);
        else
            seg_region_new = (phi < 0);
            dif_pixNum = sum(sum(abs(seg_region_old - seg_region_new)));
            if dif_pixNum < 1 % 零水平集包围的区域不再变化，则终止迭代
                fprintf('Level set evolution is converged.\n');
                break;
            else
                seg_region_old = seg_region_new;
            end
        end
    end    
end;
