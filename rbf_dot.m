function [H,deg]=rbf_dot(patterns1,patterns2,deg);

%Note : patterns are transposed for compatibility with C code.

if(~exist('patterns2','var'))
    patterns2=patterns1;
end

size1=size(patterns1);
size2=size(patterns2);

%new vectorised version

G = sum((patterns1.*patterns1),2);
H = sum((patterns2.*patterns2),2);

Q = repmat(G,1,size2(1));
R = repmat(H',size1(1),1);

H = Q + R - 2*patterns1*patterns2';

if(~exist('deg','var'))
    deg = sqrt(median(H(:)));
end

H=exp(-H/2/deg^2);

%old version

if 0
  
H=zeros(size1(1),size2(1));
for k=1:size1(1)
  for l=1:size2(1)
    H(k,l)=(patterns1(k,:)-patterns2(l,:))*(patterns1(k,:)-patterns2(l,:))';
  end
end

end