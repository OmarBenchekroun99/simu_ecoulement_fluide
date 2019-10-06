// Initialise f(x,y)
function [w]=init_field_curl(y,x)
    w = 8*%pi^2*cos(2*%pi*x)*cos(2*%pi*y)
endfunction


// Solution de référence du problème laplacien(psi(x,y)) = f(x,y)
function [refx]=solution_field_ux(y,x)
    refx = - 2* %pi* cos(2*%pi*x)*sin(2*%pi*y) 
endfunction

// Solution de référence du problème laplacien(psi(x,y)) = f(x,y)
function [refy]=solution_field_uy(y,x)
    refy = 2* %pi* sin(2*%pi*x)*cos(2*%pi*y) 
endfunction


// Affichage de la fonction f, de la solution de référence et 
// de la solution obtenue, ainsi que de l'erreur commise
// F   = f(x,y)         -- fonction testée
// Ref = Psi_alpha(x,y) -- solution analytique
// Psi = Psi_star(x,y)  -- solution du solveur
function plot_error_curl(Ref_Ux, Ref_Uy, Ux, Uy)  
    fig = gcf()
    Err_Ux = abs(Ref_Ux - Ux)
    Err_Uy = abs(Ref_Uy - Uy)
    fig.color_map = jetcolormap(100)
    subplot(321)
    plot3d1(Y, X, Ux)
    colorbar(min(Ux), max(Ux))
    xtitle(" Ux = ux(x,y) ")
    subplot(322)
    plot3d1(Y, X, Uy)
    colorbar(min(Uy), max(Uy))
    xtitle(" Uy = uy(x,y) ")
    subplot(323)
    plot3d1(Y, X, Ref_Ux)
    colorbar(min(Ref_Ux), max(Ref_Ux))
    xtitle(" Ref = ux_star(x,y) ")
    subplot(324)
    plot3d1(Y, X, Ref_Uy)
    colorbar(min(Ref_Uy), max(Ref_Uy))
    xtitle(" Ref = Uy_star(x,y) ")
    subplot(325)
    plot3d1(Y, X, Err_Ux)
    colorbar(min(Err_Ux), max(Err_Ux))
    xtitle(" Error_Ux ")
    subplot(326)
    plot3d1(Y, X, Err_Uy)
    colorbar(min(Err_Uy), max(Err_Uy))
    xtitle(" Error_Uy ")
    xs2png(fig, "poisson_error_curl.png")
endfunction


// Fonction de test pour le solveur de Poisson
function test_poisson_curl(Lx, Ly, Nx, Ny)
    printf("::Testing poisson operator::")
    printf("\n  Domain size:    [%0.2f, %0.2f]", Lx, Ly)
    printf("\n  Discretization: [%i, %i]", Nx, Ny)
    
    // X[i] = i*dx avec dx = Lx/Nx et i=0..Nx-1
    // Y[i] = j*dy avec dy = Ly/Ny et j=0..Ny-1
    X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
    Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
    disp(size(X))
    disp(size(Y))
    printf("\n\n  Initializing field F(x,y).")
    W   = feval(Y, X, init_field_curl)
    
    printf("\n  Initializing reference solution Ref_Ux(x,y).")
    Ref_Ux = feval(Y, X, solution_field_ux)
    
    printf("\n  Initializing reference solution Ref_Uy(x,y).")
    Ref_Uy = feval(Y, X, solution_field_uy)

    dir  = get_absolute_file_path("test_poisson_curl.sce")
    file = dir+"poisson.sce" 
    printf("\n\n  Loading poisson_2d function from file %s%s%s.", char(39), file, char(39))
    exec(file, -1)
   
    printf("\n\n  Computing Poisson solution Psi(x,y).")
    [Ux, Uy] = poisson_curl_2d(W, Nx, Ny, Lx, Ly)

    printf("\n  Computing error |Psi-Ref|(x,y).")
    Err_Ux = abs(Ux-Ref_Ux)
    Err_Uy = abs(Uy-Ref_Uy)
    
    file = pwd()+"/poisson_error.png"
    printf("\n\n  Plotting everything to %s%s%s.", char(39), file, char(39))
    plot_error_curl(Ref_Ux, Ref_Uy, Ux, Uy)
    
    printf("\n\n")
    mErr_Ux = max(Err_Ux)
    mErr_Uy = max(Err_Uy)
    max_error = 1e-12
      
    if (mErr_Ux > max_error) & (mErr_Uy > max_error) then
        printf("  Maximal error is %.10ef, TEST FAILURE (max_error=%.10ef).\n", mErr_Ux, max_error)
        exit(1)
    else
        printf("  Maximal error is only %.10ef, TEST SUCCESS.\n", mErr_Ux)
        exit(0)
    end
endfunction


// Taille du domaine
Lx = 1.0
Ly = 1.0

// Discretisation du domaine
Nx = 64
Ny = 32

test_poisson_curl(Lx, Ly, Nx, Ny)
