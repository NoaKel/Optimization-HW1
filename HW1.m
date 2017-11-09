% Compares analytical and numerical calculations 
close all;

% Compares results function 1 f(x)=h(Ax)
x=rand(10,1);
A=rand(10,10);
cmp(@myfunc,@numdiff,x,A,1);

% Compares results function 2 f(x)=phi(h(x)) 

x=rand(6,1);
cmp(@myfunc,@numdiff,x,[],2);