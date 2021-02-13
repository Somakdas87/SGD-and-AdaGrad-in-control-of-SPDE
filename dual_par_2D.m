function P=dual_par_2D(Y,M,A,tau,N,z,p,e,t)

%=================================================================
% solution to the dual function in two dimension
%
% dual: solution to primal equation
% par: parabolic differential equation
% 2: two dimension
%
%=================================================================
% Somak Das July 2019
%=================================================================


np=length(p);
ne=length(e);
nt=length(t);

G=zeros(N+1,np);
g=zeros(3,1);

P=zeros(N+1,np);

%----------------------------------------------------------
%terminal condition 
%p(T,v)=0
for j=1:np    
    P(N+1,j)=0;
end
%----------------------------------------------------------

for i=N:-1:1
    
    time=(i-1)*tau;
    
    for j=1:nt
        x1=p(1,t(1,j));    x2=p(1,t(2,j));    x3=p(1,t(3,j));
        y1=p(2,t(1,j));    y2=p(2,t(2,j));    y3=p(2,t(3,j));        
        
        delta=(0.5)*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)); % area of triangle
        
        m=[[2,1,1];[1,2,1];[1,1,2]];
        m=(delta/12)*m;
    
        for k=1:3
            g(1)= Y(i,t(1,j))-z(x1,y1,time);
            g(2)= Y(i,t(2,j))-z(x2,y2,time);
            g(3)= Y(i,t(3,j))-z(x3,y3,time);
        end
        g=g'*m;
        g=g';
        
        for k=1:3
            G(i,t(k,j))=G(i,t(k,j))+g(k);            
        end
        
    end
        
       
    B=M+tau*A;
    F=(tau*G(i,:)-P(i+1,:)*M);
    
    for j=1:ne    
        B(e(1,j),:)=zeros(np,1);
        B(:,e(1,j))=zeros(1,np);
        B(e(1,j),e(1,j))=1;
        F(e(1,j))=0;
    end
    
    P(i,:)=B\F';
         
end

end
