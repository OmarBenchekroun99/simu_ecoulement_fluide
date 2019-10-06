
// Initialise f(x,y)
function [f]=init_field(y,x)
    f = sin(2*%pi*x)*sin(2*%pi*y)
endfunction


// Solution de référence du problème laplacien(psi(x,y)) = f(x,y)
function [ref]=solution_field(y,x)
    alpha = -1/(8*%pi^2)
    ref = alpha * sin(2*%pi*x)* sin (2*%pi*y)
endfunction


// Affichage de la fonction f, de la solution de référence et 
// de la solution obtenue, ainsi que de l'erreur commise
// F   = f(x,y)         -- fonction testée
// Ref = Psi_alpha(x,y) -- solution analytique
// Psi = Psi_star(x,y)  -- solution du solveur
function plot_error(F, Ref, Psi)  
    fig = gcf()
    err = abs(Ref - Psi)
    fig.color_map = jetcolormap(100)
    subplot(2,2,1)
    plot3d1(Y, X, F)
    colorbar(min(F), max(F))
    xtitle(" z = f(x,y) ")
    subplot(2,2,2)
    plot3d1(Y, X, Ref)
    colorbar(min(Ref), max(Ref))
    xtitle(" z = psi_alpha(x,y) ")
    subplot(2,2,3)
    plot3d1(Y, X, Psi)
    colorbar(min(Psi), max(Psi))
    xtitle(" z = psi_star(x,y) ")
    subplot(2,2,4)
    plot3d1(Y, X, err)
    colorbar(min(err), max(err))
    xtitle(" Error ")
    xs2png(fig, "poisson_error.png")
endfunction


// Fonction de test pour le solveur de Poisson
function test_poisson(Lx, Ly, Nx, Ny)
    printf("::Testing poisson operator::")
    printf("\n  Domain size:    [%0.2f, %0.2f]", Lx, Ly)
    printf("\n  Discretization: [%i, %i]", Nx, Ny)
    
    // X[i] = i*dx avec dx = Lx/Nx et i=0..Nx-1
    // Y[i] = j*dy avec dy = Ly/Ny et j=0..Ny-1
    X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
    Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
    
    printf("\n\n  Initializing field F(x,y).")
    F   = feval(Y, X, init_field)
    
    printf("\n  Initializing reference solution Ref(x,y).")
    Ref = feval(Y, X, solution_field)

    dir  = get_absolute_file_path("test_poisson.sce")
    file = dir+"poisson.sce" 
    printf("\n\n  Loading poisson_2d function from file %s%s%s.", char(39), file, char(39))
    exec(file, -1)

    printf("\n\n  Computing Poisson solution Psi(x,y).")
    Psi = poisson_2d(F, Nx, Ny, Lx, Ly)

    printf("\n  Computing error |Psi-Ref|(x,y).")
    Err = abs(Psi-Ref)
    
    file = pwd()+"/poisson_error.png"
    printf("\n\n  Plotting everything to %s%s%s.", char(39), file, char(39))
    plot_error(F, Ref, Psi)
    
    printf("\n\n")
    mErr = max(Err)
    max_error = 1e-12
      
    if (mErr > max_error) then
        printf("  Maximal error is %.10ef, TEST FAILURE (max_error=%.10ef).\n", mErr, max_error)
        exit(1)
    else
        printf("  Maximal error is only %.10ef, TEST SUCCESS.\n", mErr)
        exit(0)
    end
endfunction


// Taille du domaine
Lx = 1.0
Ly = 1.0

// Discretisation du domaine
Nx = 64
Ny = 32

test_poisson(Lx, Ly, Nx, Ny)
