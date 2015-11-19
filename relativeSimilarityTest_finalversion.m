function [pvalue]=relativeSimilarityTest_finalversion(X,Y,Z)
% X : target
% Y : source 1
% Z : source 2

% Author: Wacha Bounliphone - wacha.bounliphone@centralesupelec.fr
% Copyright (c) 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% selection of theBandwidth;
myX = pdist2(X,Y);
myX = myX(:);
theBandwidth(1) = sqrt(median(myX(:))/2);
myX = pdist2(X,Z);
myX = myX(:);
theBandwidth(2) = sqrt(median(myX(:))/2);
theBandwidth=mean(theBandwidth);
params.sig=theBandwidth;
localSig=params.sig;

%rbf_dot kernel
Kxx = rbf_dot(X,X,localSig);
Kyy = rbf_dot(Y,Y,localSig);
Kzz = rbf_dot(Z,Z,localSig);
Kxy = rbf_dot(X,Y,localSig);
Kxz = rbf_dot(X,Z,localSig);


%mean vector \mu of the joint distribution of MMD
mu = [ MMD_unbiased(Kxx,Kyy,Kxy) ; MMD_unbiased(Kxx,Kzz,Kxz) ];

%covariance matrix \sigma of the joint distribution of MMD
Sigma(1,1) = thevariance(Kxx,Kyy,Kxy);
Sigma(2,2) = thevariance(Kxx,Kzz,Kxz);
Sigma(1,2) = thecovariance(Kxx,Kxy,Kxz);
Sigma(2,1) = Sigma(1,2);


theta=pi/4;
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

mu = R*mu;
Sigma = R*Sigma*R';

t = mu(1);
sig = sqrt(Sigma(1,1));
t = t/sig;
pvalue = normcdf(-t);

end



%Wacha covariance
function theCov = thecovariance(Kxx,Kxy,Kxz);

Kxxnd = Kxx-diag(diag(Kxx));

m = size(Kxx,1);
n = size(Kxy,2);
r = size(Kxz,2);

u_xx=sum(sum(Kxxnd))*( 1/(m*(m-1)) );
u_xy=sum(sum(Kxy))/(m*n);
u_xz=sum(sum(Kxz))/(m*r);


ct1 = (1/(m*(m-1)*(m-1)))   * sum(sum(Kxxnd*Kxxnd));
ct2 =  u_xx^2;
ct3 = (1/(m*(m-1)*r))       * sum(sum(Kxxnd*Kxz));
ct4 =  u_xx*u_xz;
ct5 = (1/(m*(m-1)*n))       * sum(sum(Kxxnd*Kxy));
ct6 = u_xx*u_xy;
ct7 = (1/(n*m*r))           * sum(sum(Kxy'*Kxz));
ct8 = u_xy*u_xz;

zeta_1 = (ct1-ct2)-(ct3-ct4)-(ct5-ct6)+(ct7-ct8);
theCov = (4*(m-2))/(m*(m-1)) * zeta_1;

end


%Wacha variance
function theVar = thevariance(Kxx,Kyy,Kxy);

Kxxnd = Kxx-diag(diag(Kxx));
Kyynd = Kyy-diag(diag(Kyy));


m = size(Kxx,1);
n = size(Kyy,1);

u_xx=sum(sum(Kxxnd))*( 1/(m*(m-1)) );
u_yy=sum(sum(Kyynd))*( 1/(n*(n-1)) );
u_xy=sum(sum(Kxy))/(m*n);

Kyx = Kxy';

vt1 = (1/(m*(m-1)*(m-1)))   * sum(sum(Kxxnd*Kxxnd)); %aa
vt2 =  u_xx^2;
vt3 = (1/(m*(m-1)*n))       * sum(sum(Kxxnd*Kxy)); %ac
vt4 =  u_xx*u_xy;
vt5 = (1/(n*(n-1)*(n-1)))   * sum(sum(Kyynd*Kyynd)); %bb
vt6 = u_yy^2;
vt7 = (1/(n*(n-1)*m))       * sum(sum(Kyynd*Kyx)); %bd
vt8 = u_yy*u_xy;
vt9 = (1/(m*n*n))           * sum(sum(Kyx*Kxy)); %cc
vt10 = u_xy^2;
vt11 = (1/(n*m*m))          * sum(sum(Kxy*Kyx)); %dd
vt12 = u_xy^2;

zeta_1 = (vt1-vt2)... %aa
    -2*(vt3-vt4) ... %ac
    +(vt5-vt6) ... % bb
    -2*(vt7-vt8) ... % bd
    + (vt9 - vt10) ... %cc
    + (vt11 - vt12); %dd

theVar = (4*(m-2))/(m*(m-1)) * zeta_1;

end

