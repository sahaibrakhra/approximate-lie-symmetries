read "HybridInvolutiveFormLHPDE.mpl";
read "NumericDiffGeometryToolsComplex.mpl";
with(PDEtools):
with(DEtools):

testsample := proc(tpoint::list,upoint::list,m,n,tol,ProlongationTimes)
              local NLS, DNLS, SymDimList, ApproxSymDim, i, j;

              NLS := [diff(u(x, t), t)*I + diff(u(x, t), x, x) + u(x, t)^2*v(x, t)^2*u(x, t) = 0, (-1)*diff(v(x, t), t)*I + diff(v(x, t), x, x) + u(x, t)^2*v(x, t)^2*v(x, t) = 0];
              DNLS := convert(DeterminingPDE(NLS, [xi(x, t, u, v), tau(x, t, u, v), eta(x, t, u, v), zeta(x, t, u, v)], integrabilityconditions = false), list);
              SymDimList := [];

              for i from tpoint[1] by m to tpoint[2] do
                  for j from upoint[1] by n to upoint[2] do
                      ApproxSymDim := GeometricInvolutiveForm(DNLS,[[x, t, u, v], [xi, tau, eta, zeta]],tol,ProlongationTimes,[0.1,i,j,j]);
                      SymDimList := [op(SymDimList),ApproxSymDim];
                  end do;
              end do;

              return SymDimList;
              end proc;