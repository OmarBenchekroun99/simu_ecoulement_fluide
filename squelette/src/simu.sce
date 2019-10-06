
// Domain size and discretization
Lx = 1.0
Ly = 1.0
Nx = 128
Ny = 128
dx = Lx/Nx
dy = Ly/Ny
X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
cmm = [0,64]

// Simulation parameters
T     = 1.50
nu    = 1e-4  //1e-4
rho   = 30.0 //30
delta = 0.05


// Initialize vorticity
function [W] = init_vorticity(y,x)
    if y <= 0.5 then
        W = 2*%pi*delta*cos(2*%pi*x) + rho*(1 - tanh(rho*(y - 0.25))^2)
    else
        W = 2*%pi*delta*cos(2*%pi*x) - rho*(1 - tanh(rho*(0.75 - y))^2)
    end
endfunction


// Plot and dump fields (at every steps of the simulation so that we can generate a video)
function plot_fields(W, Ux, Uy, iteration)
    if (t~=0.80) & (t~=1.20) then
        return
    end
    fig = scf(0)
    clf()
    fig.color_map = [jetcolormap(64); whitecolormap(1)]
    fig.background = 65;
    cmm = [0,64]

    subplot(131)
    title("W(x,y)")
    colorbar(min(W), max(W), cmm)
    Sgrayplot(X, Y, W', colminmax=cmm)
    xlabel("x")
    ylabel("y")
    
    subplot(132)
    title("Ux(x,y)")
    colorbar(min(Ux), max(Ux), cmm)
    Sgrayplot(X, Y, Ux', colminmax=cmm)
    xlabel("x")
    ylabel("y")
    
    subplot(133)
    title("Uy(x,y)")
    colorbar(min(Uy), max(Uy), cmm)
    Sgrayplot(X, Y, Uy', colminmax=cmm)
    xlabel("x")
    ylabel("y")
        
    figname = sprintf("ite_%04d.png", iteration)
    xs2png(fig, figname)
endfunction

// Plot isocontours (only at t=0.80 and t=1.20)
function plot_isocontours(W, figname)
    if (t~=0.80) & (t~=1.20) then
        return
    end
    
    fig = scf(1)
    clf()
    fig.color_map = [jetcolormap(64); whitecolormap(1)]
    fig.background = 65;
    

    // Affichage des isocontours
    
    title("W(x,y)")
    
    colorbar(min(W), max(W), cmm)
    Sgrayplot(X,Y,W', colminmax=cmm)
    contour2d(X,Y,W', -36:6:36) //-70:10:70 pour la question suivante
    xlabel("x")
    ylabel("y")

    figname = sprintf("isocontours_%f.png", t)
    xs2png(fig, figname)
endfunction

function plot_champs(nom_figure)
    if (t~=0.80) & (t~=1.20) then
        return
    end
    
    fig = scf(2)
    clf()
    
    fig.color_map = [jetcolormap(64); whitecolormap(1)]
    fig.background = 65;
    
    title("Champs des vitesses")
    // L'utilisation de fchamp présente un problème
    //La densité des flèches est trop grande et la figure devient illisible
    //De fait nous avons réduit le nombre de point (divisé par 4) pour un soucis de 
    //lisibilité
    X2 = zeros(Nx/4,1)
    Ux2 = zeros(Ny/4,Nx/4)
    Uy2 = zeros(Ny/4,Nx/4)
    Y2 = zeros(Ny/4,1)
    //On définit une matrice INt qui prend les intensités de ces nouvelles valeurs 
    //pour créer l'échelle de couleurs
    Int = zeros(Ny/4, Nx/4)
    
    for x = 1:Nx/4
        X2(x) = X((x-1)*4 + 1)
        Y2(x) = Y((x-1)*4 + 1)
        for y = 1:Ny/4
            Ux2(x,y) = Ux((x-1)*4+1, (y-1)*4+1)
            Uy2(x,y) = Uy((x-1)*4+1, (y-1)*4+1)
            Int(x,y) = sqrt(Ux2(x,y)^2 + Uy2(x,y)^2)
        end
    end
    
    //Affichage du champ de vitesses
    
    champ1(X2, Y2, Ux2', Uy2')
    colorbar(min(Int), max(Int))
    xlabel("x")
    ylabel("y")
    
    nom_figure = sprintf("Champs_%f.png", t)
    xs2png(fig, nom_figure)
endfunction
// Load the Poisson solver and the advection-diffusion solver
dir  = get_absolute_file_path("simu.sce")
file = dir+"poisson/poisson.sce" 
exec(file, -1)
file = dir+"diff/dif-conv-f.sce" 
exec(file, -1)

// Figure setup (fig0 = fields, fig1 = isocontours)
figure(0, "position", [0,0,1200,400])
figure(1, "position", [0,0,800,800])

// Initialize vorticity and loop untill final time using adaptive timestep
t = 0.0
ite = 0
W = feval(Y, X, init_vorticity)
while t<T
    // compute velocity from vorticity
    [Ux, Uy] = poisson_curl_2d(W, Nx, Ny, Lx, Ly)

    // compute new timestep from stability criteria
    dt = min(calcul_dt(Ux,dx), calcul_dt(Uy,dy))
    if (t<0.80) & (t+dt>0.80) then
        dt = 0.80-t
    elseif (t<1.20) & (t+dt>1.20) then
        dt = 1.20-t
    elseif (t+dt>T) then
        dt = T-t
    end
    
    printf("\niteration %i, from t=%f to t=%f", ite, t, t+dt)
    plot_fields(W,Ux,Uy,ite)
    plot_isocontours(W,t)
    plot_champs(t)
    
    
    // advection-diffusion on vorticity
    W = solveur_2D(W, Ux, Uy, Nx, Ny, nu, dt, dx, dy);
 
    // Mise à jour de t et ite
    t = t+dt;
    ite = ite +1
end
plot_fields(W,Ux,Uy,ite)
plot_isocontours(W,t)
plot_champs(t)

printf("\nDone in %i iterations!\n", ite)
exit(0)
