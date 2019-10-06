
function dt=calcul_dt(u,dx)
  dt=dx/max(abs(u));
endfunction


function vort_s=solveur_1D(vort,Ux,Nx,kappa,dt,dx)
        //Implémentation du solveur 1D
        
        N=zeros(Nx, Nx)
        M = zeros(Nx, Nx)
        cst = dt/(dx)^2
        cst2 = dt/dx
        for n = 1:Nx;
        if (n == 1) then
           N(n, Nx) = -kappa*cst;
           M(n, Nx) = Ux(n)*0.5*cst2+Ux(n)^2*0.5*cst2^2;
        else
           M(n, n - 1) = Ux(n)*0.5*cst2+Ux(n)^2*0.5*cst2^2;
           N(n, n - 1) = -kappa*cst;
        end 
        N(n, n) = 1 + 2*kappa*cst;
        M(n, n) = 1 - Ux(n)^2*cst2^2;
        if (n == Nx) then
           N(n, 1) = -kappa*cst;
           M(n, 1) = -Ux(n)*0.5*cst2+Ux(n)^2*0.5*cst2^2;
       else
           N(n, n + 1) = -kappa*cst;
           M(n, n + 1) = -Ux(n)*0.5*cst2+Ux(n)^2*0.5*cst2^2;
       end
    end
    vort_s = umfpack(sparse(N), "\", M*vort)
endfunction

function vort = solveur_2D(vort,Ux,Uy,Nx,Ny,kappa,dt,dx,dy)
    //Implémentation du solveur 2D
    for i=1:Ny
        vort(i,:) = solveur_1D(vort(i,:)', Ux(i,:), Nx, kappa, dt, dx)'
    end
    for j = 1:Nx
        vort(:,j) = solveur_1D(vort(:,j), Uy(:,j), Ny, kappa, dt, dy)
    end
endfunction
