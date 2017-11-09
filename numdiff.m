function [ gnum, Hnum ] = numdiff(myfunc, x, A, h, tasknum)
% Calculates gradient and hessian numerically
[n,~,~]=size(x);
E=eye(n);
[~,l]=size(h);

% initalizing matrices
f1=zeros(n,1,l);f2=f1;
g1=zeros(n,n,l);g2=g1;
gnum=zeros(n,l);
Hnum=zeros(n,n,l);
    
% numerical gradient calculation
for i=1:l
    for p=1:n
    [f1(p,1,i),g1(:,p,i),~]=myfunc(x(:,i)+(h(1,i)*E(:,p)),A,tasknum);
    [f2(p,1,i),g2(:,p,i),~]=myfunc(x(:,i)-(h(1,i)*E(:,p)),A,tasknum);
    end
end

% gradient calculation
for i=1:l
   gnum(:,i)=(f1(:,:,i)-f2(:,:,i))./(2*h(1,i));
end 


% hessian calculaion is based on gradient's anlytical calculations
if (nargout>1)
    for i=1:l 
        Hnum(:,:,i)=(g1(:,:,i)-g2(:,:,i))./(2*(h(1,i)));
    end
end 


end