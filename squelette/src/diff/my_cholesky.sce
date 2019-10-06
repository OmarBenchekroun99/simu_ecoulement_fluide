
function [A]=cholesky_fact(A)
    // retourne la matrice triangulaire inférieure L tel que A = L.tr(L)
    n = size(A)(1)
    L = zeros(n,n)
    L(1,1) = sqrt(A(1,1))
    for j = 2:n
        L(j,1) = A(1,j)/L(1,1)
    end
    for i = 2:n
        somme1 = A(i,i);
        for k = 1:i-1
            somme1 = somme1 - L(i,k)^2;
        end;
        L(i,i) = sqrt(somme1);
        
        for j = i+1:n
            somme2 = A(i,j);
            for k = 1:i+1
                somme2 = somme2 - L(i,k)*L(j,k);
            end;
            L(j,i) = somme2/L(i,i);
        end;
    end
    A = L
endfunction


function [y]=up_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
      y=x;
      y(n)=y(n)/A(n,n);
      for i=n-1:-1:1,
          //y(i)=y(i)-(A(i,i+1:n))*y(i+1:n);
          //}y(i)=y(i)/A(i,i)
          y(i)=(y(i)-(A(i+1:n,i))'*y(i+1:n))/A(i,i);
      end;
  end;
endfunction

function [y]=down_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
      y=x;
      y(1)=y(1)/A(1,1);
      for i=2:n
          y(i)=y(i)-A(i,1:i-1)*y(1:i-1);
          y(i)=y(i)/A(i,i);
      end;
  end;
endfunction

function [U]=my_cholesky(N,S)
    T = cholesky_fact(N)
    U = down_sweep_cholesky(T,S)
    U = up_sweep_cholesky(T,U)
endfunction
