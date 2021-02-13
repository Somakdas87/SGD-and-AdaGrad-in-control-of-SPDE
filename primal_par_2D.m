function [Y,M,A]=primal_par_2D(omega,U,tau,N,p,e,t)

%=================================================================
% solution to the primal function in two dimension
%
% primal: solution to primal equation
% par: parabolic differential equation
% 2: two dimension
%
%=================================================================
% Somak Das July 2019
%=================================================================

%-------------------------------------------------
% Exporting triangular mesh 
np=length(p);
ne=length(e);
nt=length(t);
%-------------------------------------------------

%-------------------------------------------------
% Initializing 
A=zeros(np,np);
M=zeros(np,np);
G_trans=zeros(N+1,np);

Y=zeros(N+1,np);
%-------------------------------------------------


%-------------------------------------------------
% Initial value y(x,y,0) 
for j=1:np
    x1=p(1,j); y1=p(2,j);
    %Y(1,j)=x1*y1*(1-x1)*(1-y1);
    Y(1,j)=1*(-x1*(x1-0.5)^(2)*(x1-1)-y1*(y1-0.5)^(2)*(y1-1));
end
%-------------------------------------------------

% a(x)=x+y+w

%-------------------------------------------------
% Solution to primal equation 
for i=2:N+1
    for j=1:nt
        
        x1=p(1,t(1,j));    x2=p(1,t(2,j));    x3=p(1,t(3,j));
        y1=p(2,t(1,j));    y2=p(2,t(2,j));    y3=p(2,t(3,j));
        
        a10=x2*y3-x3*y2; a11=y2-y3; a12=x3-x2;
        a20=x3*y1-x1*y3; a21=y3-y1; a22=x1-x3;
        a30=x1*y2-x2*y1; a31=y1-y2; a32=x2-x1;
        
        delta=(0.5)*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)); % area of triangle
        
        b=[[a11,a21,a31];[a12,a22,a32]];
        a=b'*b;
        a=(1/(12*delta))*(5*(x1+y1)+(x2+y2)+(x3+y3)+3*omega)*a;
        
        m=[[2,1,1];[1,2,1];[1,1,2]];
        m=(delta/12)*m;
        
        u=[U(i,t(1,j)); U(i,t(2,j)); U(i,t(3,j))];
        
        time=(i-1)*tau;
        
        f1=exp(-time)*(2*x1*(1-x1)+2*y1*(1-y1)-x1*y1*(1-x1)*(1-y1));
        f2=exp(-time)*(2*x2*(1-x2)+2*y2*(1-y2)-x2*y2*(1-x2)*(1-y2));
        f3=exp(-time)*(2*x3*(1-x3)+2*y3*(1-y3)-x3*y3*(1-x3)*(1-y3));
        f=[f1 f2 f3];
        G=f*m;
        G=G';
        
        for k=1:3
            for l=1:3
                A(t(k,j),t(l,j))=A(t(k,j),t(l,j))+a(k,l);
                M(t(k,j),t(l,j))=M(t(k,j),t(l,j))+m(k,l);
            end
            G_trans(i,t(k,j))=G_trans(i,t(k,j))+G(k)+u(k);
        end
        
    end
        
       
    B=M+tau*A;
    F=(Y(i-1,:)*M+tau*G_trans(i,:));
    
    for j=1:ne    
        B(e(1,j),:)=zeros(np,1);
        B(:,e(1,j))=zeros(1,np);
        B(e(1,j),e(1,j))=1;
        F(e(1,j))=0;
    end
    
    Y(i,:)=B\F';
         
end

end