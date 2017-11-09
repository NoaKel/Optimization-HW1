function [  ] = cmp(myfunc,numdiff,x,par,tasknum)
% Compares analytical and numerical calculations 

epsilon_machine =2*10^-16;
experiments_number=10000;
h=linspace(epsilon_machine*10^9,epsilon_machine*10^12,experiments_number);

% plot infinity norm of gradient as a function of h

[~,g,H]= myfunc(x,par,tasknum);
x_rep=repmat(x,size(h));
[gnum,Hnum]=numdiff(myfunc,x_rep,par,h,tasknum);

g_error=max(abs(repmat(g,size(h))-gnum));
[a,b]=min(g_error);

figure;hold on;
loglog(log(h),log(g_error));
text(log(h(b)),log(a),['\fontsize{10}\rightarrow minimal error for h= '  num2str(h(b))]);
xlabel('log(h)');
ylabel('log(infinity norm of gradient)');
title(['\fontsize{14}function ' num2str(tasknum) '- log infinity norm of gradient as a function of h']);
hold off;

% plot infinity norm of h as a function of h

H_error=max(abs(repmat(H,1,size(h))-Hnum));H_error=squeeze(H_error);H_error=max(H_error);
[c,d]=min(H_error);

figure;hold on;
loglog(log(h),log(H_error));
text(log(h(d)),log(c),['\fontsize{10}\rightarrow minimal error for h= '  num2str(h(d))]);
xlabel('log(h)');
ylabel('log(infinity norm of hessian)');
title(['\fontsize{14}function ' num2str(tasknum) '- log infinity norm of hessian as a function of h']);
hold off;

% plot gradient error for optimal h

h=((epsilon_machine)^(1/3)*max(abs(x)));
[gnum,Hnum]=numdiff(myfunc,x,par,h,tasknum); 

g_error=abs(g-gnum);
avg1=mean(g_error);
max1=max(g_error);
min1=min(g_error);
figure;hold on;
plot(g_error);
title({['\fontsize{14} function' num2str(tasknum) '- gradient error for optimal h'];['  mean error: ' num2str(avg1) '  minimal error: ' num2str(min1) '  maximal error: ' num2str(max1)]} );
hold off;

% plot hessian error for optimal h

H_error=abs(H-Hnum);

figure;hold on;
imagesc(H_error);colorbar
title(['\fontsize{14} function ' num2str(tasknum) '-  hessian error for optimal h']);
hold off;

end

