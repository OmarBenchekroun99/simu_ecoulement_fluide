Nx=50;
Nt=50;
kappa=0.01;
exec("my_cholesky.sce")

dt=1/Nt
dx=1/Nt

function y=phi0(x)
    if (x>=0) & (x<0.25) then
        y=0;
    elseif (x>=0.25) & (x<0.375) then
        y=2*(x-0.25);
    elseif (x>=0.375) & (x<0.5) then
        y=2*(0.5-x);
    elseif (x>=0.5) & (x<=1) then
        y=0;
    end
endfunction

phi_0=[]
for i = 1:Nt
    phi_0($+1)=phi0((i-1)*dx);
end


function y=conv(x)
    y=0.4*(x-0.25)
endfunction


N=zeros(Nt, Nt)
for i=1:Nt
    N(i,i)=1+2*kappa*(dt/dx^2);
    if i+1<=Nt then
        N(i,i+1)=-kappa*(dt/dx^2);
    else //donc i=Nt
        N(i,1)=-kappa*(dt/dx^2);
    end;
    if i-1>=1 then
        N(i,i-1)=-kappa*(dt/dx^2);
    else //donc i=1
        N(i,Nt)=-kappa*(dt/dx^2);
    end;
end;


M=zeros(Nt, Nt)
for j=1:Nt
    c=conv((j-1)*dx);
    M(j,j) = 1-c^2*dt^2/dx^2
    if j>=2 then
        M(j-1,j) = -c*dt/(2*dx) + c^2*dt^2/(2*dx^2);
    else
        M(Nt,1) = -c*dt/(2*dx) + c^2*dt^2/(2*dx^2);
    end;
    if j<=Nt-1 then
        M(j+1,j) = c*dt/(2*dx) + c^2*dt^2/(2*dx^2);
    else
        M(1,Nt) = c*dt/(2*dx) + c^2*dt^2/(2*dx^2);
    end;
end;



phi = phi_0;
fin=Nt;
for i=1:fin
    phi=my_cholesky(N, M*phi);
end
scf;

maillage=linspace(0,1,Nx)';

plot(maillage,[phi_0, phi]);
