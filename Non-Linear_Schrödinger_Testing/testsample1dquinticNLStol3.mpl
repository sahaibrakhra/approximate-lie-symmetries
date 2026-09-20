read "1dquinticNLSsafety.mpl";
Digits := 15;
SymDimResult := testsample([0,4],[0,4],1,0.2,10^(-3),7);

save SymDimResult, "1dquinticNLSresult_tol3.m":