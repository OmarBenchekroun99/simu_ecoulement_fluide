
dx=Lx/Nx;
dy=Ly/Ny;

exec("dif-conv-f.sce")

maillage_x=linspace(0,(Nx-1)/Nx*Lx,Nx)';
maillage_y=linspace(0,(Ny-1)/Ny*Ly,Ny)';

cx=zeros(Ny,Nx); //composante x de la vitesse de convection
cy=zeros(Ny,Nx); //composante y de la vitesse de convection

phi=zeros(Ny,Nx);   //fonction à calculer
phi_i=zeros(Ny,Nx); //condtion initiale

//------------------------------------------
//TODO remplir les tableaux cx cy phi phi_i
//------------------------------------------
for j=1:Nx
    x=j*dx;
    for i=1:Ny
        y = i*dy;
        c = conv(y,x);
        cx(i,j) = c(1);
        cy(i,j) = c(2);
        phi_i(i,j) = phi_0(y,x);
    end

end

phi=phi_i;
dt=min(calcul_dt(cx,dx),calcul_dt(cy,dy));
Nt=floor(Tf/dt);
for k=1:Nt
    disp(k)
    //-----------------
    // TODO
    //----------------
    phi=solveur_2D(phi, cx, cy, Nx, Ny, kappa, dt, dx, dy);
end

//Graphe initial
G_init = scf()
surf(maillage_x, maillage_y, phi_i)
G_init.color_map = jetcolormap(100)

//Graphe final
G_final = scf()
surf(maillage_x, maillage_y, phi)
G_final.color_map = jetcolormap(100)
