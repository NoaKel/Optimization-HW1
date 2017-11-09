function [ f,g,H ] = myfunc(x,par,tasknum)


if(tasknum==1)
% Calculates based on analytical calculations f(x)=h(Ax)
% function value, gradient and hessian.
    A=par;
    u= A*x; 
    f=u(1)^3+u(1)*u(2)+u(2)^2+sin(u(3))+u(4)*u(5)+exp(u(6))+(u(7)*u(9))+cos(u(8))+10*u(10);
    % example gradphi=[2u1^2+u2,u1+2u2, cosu3, u5, u4, exp(u6), u9, -sinu8, u7, 10]

    if (nargout>1)
        g= A'*[3*u(1)^2+u(2);u(1)+2*u(2); cos(u(3)); u(5); u(4); exp(u(6)); u(9); -sin(u(8)); u(7); 10];

        if (nargout>2)
            H=A'*[6*u(1),1,0,0,0,0,0,0,0,0 ;...
                1,2,0,0,0,0,0,0,0,0; ...
                0,0,-sin(u(3)),0,0,0,0,0,0,0; ...
                0,0,0,0,1,0,0,0,0,0; ...
                0,0,0,1,0,0,0,0,0,0; ...
                0,0,0,0,0,exp(u(6)),0,0,0,0; ...
                0,0,0,0,0,0,0,0,1,0; ...
                0,0,0,0,0,0,0,-cos(u(8)),0,0; ...
                0,0,0,0,0,0,1,0,0,0; ...
                0,0,0,0,0,0,0,0,0,0; ...
                ]*A;
        end
    end
end

if (tasknum==2)
% Calculates based on analytical calculations f(x)=phi(h(x)) 
% function value, gradient and hessian.


    
    h=x(1)^3+x(1)*x(2)+x(2)^2+exp(x(3)+x(4))+cos(x(5)+x(6));
    f=sin(h); % example phi=sin(x)

    if (nargout>1)

        grad_h=[3*x(1)^2+x(2);x(1)+2*(x(2));exp(x(3)+x(4)); exp(x(3)+x(4)); -sin(x(5)+x(6)); -sin(x(5)+x(6))];
        g= cos(h)*grad_h; 
        
        if (nargout>2)

            Hessian_h=[6*x(1),1,0,0,0,0;...
            1,2,0,0,0,0;...
            0,0,exp(x(3)+x(4)),exp(x(3)+x(4)),0,0;...
            0,0,exp(x(3)+x(4)),exp(x(3)+x(4)),0,0;...
            0,0,0,0,-cos(x(5)+x(6)),-cos(x(5)+x(6));...
            0,0,0,0,-cos(x(5)+x(6)),-cos(x(5)+x(6))];
            H=-sin(h)*(grad_h)*(grad_h)'+cos(h)*Hessian_h;  
        end
    end
end

end
