read "1dquinticNLSsafety.mpl";
Digits := 15;
SymDimResult := testsample([0,4],[0,4],1,0.2,10^(-10),7);

save SymDimResult, "1dquinticNLSresult_tol10.m":