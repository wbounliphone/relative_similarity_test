function [MMDsquared_unbiased] = MMD_unbiased(Kxx,Kyy,Kxy);

% Author: Wacha Bounliphone - wacha.bounliphone@centralesupelec.fr
% Copyright (c) 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%length
m = length(Kxx);
n = length(Kyy);

sumKxxnd = sum(sum(Kxx - diag(diag(Kxx))));
sumKyynd = sum(sum(Kyy - diag(diag(Kyy)) ));
sumKxy = sum(sum(Kxy));

sumKxxnd = sumKxxnd /(m*(m-1));
sumKyynd = sumKyynd / (n*(n-1));
sumKxy = -sumKxy * 2/m/n;

%MMD unbiased 
MMDsquared_unbiased = sumKxxnd +  sumKyynd + sumKxy;

end