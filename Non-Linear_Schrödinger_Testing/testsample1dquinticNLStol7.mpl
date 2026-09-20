read "1dquinticNLSsafety.mpl";
Digits := 15;
SymDimResult := testsample([0,4],[0,4],1,0.2,10^(-7),7);

save SymDimResult, "1dquinticNLSresult_tol7.m":