%=====================================================================
% 
% Stocastic Gradient vs AdaGrad
%
% Optimal Control of Parabolic PDE using SGD with fixed mesh
%
% Solve the optimal problem 
%
% min ||y-0||^2 + ||u||^2,  u in L2(0,T;L2[0,1]),
%
% ||.|| = L2 norm (over [0,1]x[0,T] 
%    
% alpha=1
% l=2
% tau> 1/l
%
%=====================================================================
% Somak Das November 2019
%=====================================================================


np=length(p);
ne=length(e);
nt=length(t);

T=1;
tau=1/10;
N=ceil(T/tau);          % partition in time

z=@(x,y,t) 0; %desired state

Z=zeros(N+1,np);
for i=1:N+1
    time=(i-1)*tau;
    for j=1:np
        x=p(1,j);
        y=p(2,j);
        Z(i,j)=z(x,y,time);
    end
end

U=(-2)*ones(N+1,np);
%U=zeros(N+1,np);

n=1e3;
D=data(n); %creates 'n' uniform iid data in (0,1)

tol=0.0001;
alpha=2;

% AdaGrad
grad_J=1;
B=zeros(N+1,np)+0.001;
Risk=[0 0]; %to plot omega vs norm of gradient of risk fn
c=0;

tic
while norm(grad_J)>tol && c<1000
    c=c+1;
    step=1;
        
    %randomly pick a realisation from the data set
    i=rand; i=i*(n); i=floor(i);
    if i==0
        i=n;
    end
        
    
    [Y,M,A]=primal_par_2D(D(i),U,tau,N,p,e,t);
    P=dual_par_2D(Y,M,A,tau,N,z,p,e,t);
            
    
    grad_J=alpha*U+P;
    B=sqrt(B.^2 + grad_J.^2); 
    U=U-step*grad_J./B;       
       
            
    newrow=[c, norm(grad_J)];
    Risk=[Risk; newrow];
    
            
    newrow=[c, norm(grad_J)];
    Risk=[Risk; newrow];
end
display('Adagrad:')
toc


display(c);


R=Risk(2:(c+1),1);
S=Risk(2:(c+1),2);
figure(1)
hold on
plot(R,S,'DisplayName','AdaGrad');

%SGD
U=zeros(N+1,np);

grad_J=1;

Risk=[0 0]; %to plot omega vs norm of gradient of risk fn
k=c;
c=0;
E=0;
tic
while norm(grad_J)>tol && c<k
    c=c+1;
    step=1;    
    
    %randomly pick a realisation from the data set
    i=rand; i=i*(n); i=floor(i);
    if i==0
        i=n;
    end
    
    [Y,M,A]=primal_par_2D(D(i),U,tau,N,p,e,t);
    P=dual_par_2D(Y,M,A,tau,N,z,p,e,t);
            
    grad_J=alpha*U+P;
    U=U-step*grad_J; 
    
    
    newrow=[c, norm(grad_J)];
    Risk=[Risk; newrow];
    
            
    newrow=[c, norm(grad_J)];
    Risk=[Risk; newrow];
end
display('SGD')
toc

display(c);


R=Risk(2:(c+1),1);
S=Risk(2:(c+1),2);
figure(1)
plot(R,S,'--','DisplayName','SGD');
hold off


