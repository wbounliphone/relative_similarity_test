function [X,Y,Z] = generate_data_vWacha(m,n,r,elongate,alpha,shape,means);

% Author: Wacha Bounliphone - wacha.bounliphone@centralesupelec.fr
% Copyright (c) 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1
    m=3001;
    n=3002;
    r=3003;
    elongate = [3,3,3];
    alpha =[0,pi/12,pi/2]; 
    means=[0,0,3];
    shape = [1,1,3];
end

if(~exist('alpha','var'))
    alpha =[0,pi/12,pi/4];
end

if(~exist('means','var'))
    means=[0,0,3];
end

if(~exist('shape','var'))
    shape =[1,1,2];
end

X = randn(m,2);
Y = randn(n,2);
Z = randn(r,2);
    
X(:,2)=X(:,2).*elongate(1);
Y(:,2)=Y(:,2).*elongate(2);
Z(:,2)=Z(:,2).*elongate(3);
       
     
alphaX = alpha(1);
RX = rotz(alphaX);

%    X = X*RX;
 X = RX*X';
 X = X';

alphaY = alpha(2);
RY = rotz(alphaY);
%    Y = Y*RY;
 Y = RY*Y';
 Y = Y';

alphaZ = alpha(3);
RZ = rotz(alphaZ);
%Z = Z*RZ;
 Z = RZ*Z';
 Z = Z';

 X=X*shape(1)+means(1);
 Y=Y*shape(2)+means(2);
 Z=Z*shape(3)+means(3);
     
end


function RY = rotz(a);
    RY = [cos(a) -sin(a);sin(a),cos(a)];
end

% figure; hold on;
% scatter(X(:,1),X(:,2),'b')
% scatter(Y(:,1),Y(:,2),'m')
% scatter(Z(:,1),Z(:,2),'g')
% % 
% figure();
% scatter(X(:,1),X(:,2),'b')
% axis equal
% xlim([-10 10])
% ylim([-10 10])
% cropAxis
% print('Figures\synthetic_x','-dpdf')
% figure();
% scatter(Y(:,1),Y(:,2),'m')
% axis equal
% xlim([-10 10])
% ylim([-10 10])
% cropAxis
% print('Figures\synthetic_y','-dpdf')
% figure();
% scatter(Z(:,1),Z(:,2),'g')
% axis equal
% xlim([-10 10])
% ylim([-10 10])
% cropAxis
% print('Figures\synthetic_z','-dpdf')
