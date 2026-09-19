### Approximate Symmetry Project: Chapter 4 of the Zahra's thesis,
### 	Hibrid (Involutive Form) for Linear Homogeneous Partial differential equations
###  Authors: Zahra Mohammadi & Greg Ried
###  Last modified Jun 2020 by Zahra

# HybridIFofLHPDE := proc(sys::list,vars::list,splitvar)
# # add third argument which is either float or a list of function, for exmaple G(x,y,z), which can be used to find system placed in ApSys
# ## Input:  sys is a list of polynomially linear homogeneous differential equations
# ## vars is a list of [[x], [u]] where [x] is a list of independent variable names and [u] is a list of dependent variable names.
# ## option encrypted;
# uses difforder=PDEtools:-difforder,initialdata=DEtools:-initialdata, rifsimp=DEtools:-rifsimp, dsubs=PDEtools:-dsubs; 
# local ExSys, ApSys, i, AppCoeff, rExSys, IDExSys,LHrExSys, PLHrExSys, d, MindegreeDep, InvIndepMulti, leastComM, Listmondep, 
# 		Remdep, j, k, n, leastComDepL, SolvedrifExact, Monomialdep, MaxdegreeDep, MaxDegreeEx, MinDegreeEx,
# 		HFExSys, SimpApSys, q2, mono; 
# 	# Dorder := difforder(sys);
# 	#Split the input to the exact and approximate parts (alter)
# 	ExSys := [];
# 	ApSys := sys;
# 	flag := true;
# 	# ApSys := [];
# 	# for i from 1 to nops(sys) do  
# 	# 	AppCoeff := indets(sys[i],float);
# 	# 	if AppCoeff <> {} then
# 	# 		ApSys := [op(ApSys), sys[i]];
# 	# 	else 
# 	# 		ExSys := [op(ExSys), sys[i]];
# 	# 	end if;
# 	# end do;
# 	userinfo(1, HybridIFofLHPDE, "The list of Exact equations =", ExSys);
# 	userinfo(1, HybridIFofLHPDE, "The list of Approximate equations =", ApSys);
# 	# rExSys := DEtools:-rifsimp(ExSys);
# 	userinfo(1, HybridIFofLHPDE, "rif form of Exact equations =", rExSys[Solved]);
# 	IDExSys := DEtools:-initialdata(rExSys[Solved]);
# 	userinfo(1, HybridIFofLHPDE, "Initial Data form of Exact equations =", [IDExSys[Finite], IDExSys[Infinite]]);
# 	#predict q1 (order of prolongation), order Involutivity exact system
# 	LHrExSys := map(lhs, rExSys[Solved]);
# 	PLHrExSys:= PolynomialTools:-PDEToPolynomial(LHrExSys, vars[1], vars[2]);
# 	lprint(`PLHrExSys=`,PLHrExSys);
# 	d := [];
# 	MindegreeDep := [];
# 	## this loop will be used for finding lower bound
# 	InvIndepMulti := 1;
# 	for i in vars[1] do
# 		InvIndepMulti := InvIndepMulti*(i)^(-1);
# 	end do;
# 	#lprint(`inversemultiindep =` , InvIndepMulti);
# 	## this loop provids list of monomals regards to each dependend variables and lcm
# 	for j from 1 to nops(vars[2]) do
# 		leastComM := [];
# 		if nops(vars[2]) > 1 then
# 			Listmondep := select(has, PLHrExSys, vars[2][j]);
# 			#lprint(`Listmondep =`, Listmondep);
# 			Remdep :=	map(z -> z/(vars[2][j]), Listmondep);
# 		else 
# 			Remdep := PLHrExSys; 
# 		end if; 
	
# 		if nops(Remdep) = 1 then 
# 			leastComM := [op(leastComM), op(Remdep)];
# 			MindegreeDep := degree(op(leastComM));
# 			#lprint(`One leading for dep var`, vars[2][j], `with order =`, degree(op(leastComM)));
# 		else 
# 			for k from 1 to nops(Remdep)-1 do
#  	 			for n from k+1 to nops(Remdep) do 
#   	 	 			leastComM := [op(leastComM), lcm(Remdep[n], Remdep[k])];
#    				end do;
# 			end do;
# 			leastComDepL := map(z->InvIndepMulti*z, leastComM);
# 			lprint(`leastComDepL =`, leastComDepL);
# 			mono :=map(z->z*vars[2][j],leastComDepL);
# 			#lprint(`mono =`, mono);
# 			SolvedrifExact := rExSys[Solved];
#   		 	Monomialdep := DetectBorder(mono,SolvedrifExact,vars);
# 			#lprint(`MonomialdepOnbou =`, Monomialdep[MonoOnbound]);
# 			MindegreeDep := [op(MindegreeDep), Monomialdep[MininumDegree]];
# 			#lprint(`MindegreeDep =`, MindegreeDep);
# 		end if;
# 		#lprint(`leastComM=`, leastComM);
# 		MaxdegreeDep := map(degree, leastComM);
# 		d := [op(d), max(MaxdegreeDep) - 1];
# 	end do;
# 	MaxDegreeEx := max(d);
# 	#lprint(`MaxdegreeEx =`, MaxDegreeEx);
# 	MinDegreeEx := max(MindegreeDep);
# 	#lprint(`MinDegreeEx =`, MinDegreeEx);
# 	HFExSys := DifferentialHilbertFunction(IDExSys, s);
# 	SimpApSys := dsubs(rExSys[Solved], ApSys);
# 	q2 := difforder(SimpApSys);
# 	return table([ExactSystem = ExSys, rifExact = rExSys[Solved], SimplifiedApSys = SimpApSys, 
# 						MaxDegreeExact = MaxDegreeEx, MinDegreeExact = MinDegreeEx ]);
# end proc:

SplitExApSys:=proc(ExSys, ApSys, inflag, splitvar)
    # Input: Disjoint systems exact system ExSys, approximate system ApSys and a flag.
    # Output: [rExSys, SimpApSys, flag] where rExSys is in rif-form, SimpApSys is an approximate system simplified with respect to rExSys.
    local rExSys, SimpApSys, ExSimpApSys, i, AppCoeff, newApSys, outflag, ExSysOut, ApSysOut;
	ExSimpApSys := [];
	newApSys := [];
	if ExSys <> [] then 
	    rExSys := DEtools:-rifsimp(ExSys);
	    SimpApSys := PDEtools:-dsubs(rExSys[Solved], ApSys);
		userinfo(2, NumericJetTools,`rExSys =`, rExSys[Solved]);
	else 
	    rExSys := []; 
		SimpApSys := ApSys;
		userinfo(2, NumericJetTools,`rExSys =`, rExSys);
	end if;
	
	if splitvar = float then
	    for i from 1 to nops(SimpApSys) do  
		    AppCoeff := indets(SimpApSys[i],splitvar);
			if AppCoeff = {} then
		        ExSimpApSys := [op(ExSimpApSys),SimpApSys[i]];	
		    else
			    newApSys := [op(newApSys), SimpApSys[i]];
			end if;
	    end do;
	else
	    for i from 1 to nops(SimpApSys) do  
		    AppCoeff := has(SimpApSys[i],splitvar);
			if AppCoeff = false then
		        ExSimpApSys := [op(ExSimpApSys),SimpApSys[i]];
			else
			    newApSys := [op(newApSys), SimpApSys[i]];
		    end if;
	    end do;
	end if;
	if ExSimpApSys = [] then
	    outflag := true;
	else
	    outflag :=false;
    end if;

	userinfo(3, NumericJetTools,`ExSimpApSys =`, ExSimpApSys);
    if rExSys <> [] then
	    ExSysOut := [op(rExSys[Solved]),op(ExSimpApSys)];
	else
	    ExSysOut := [op(ExSimpApSys)];
	end if;
	ApSysOut := newApSys;
	return table([ExactSystem = ExSysOut, ApproximateSystem = ApSysOut, flagname = outflag]);
end proc:

ComputeInvolutiveOrder := proc(ExSys, vars)
local LHrExSys,d, PLHrExSys, MindegreeDep, InvIndepMulti, leastComM, Listmondep, Remdep, leastComDepL, SolvedrifExact, 
mono, Monomialdep, MaxdegreeDep, MansfieldInvDeg, ReducedInvDeg, j, k, n,
s, InitialdataExSys, HilbertFuncExSys; # Aug 31, 2022
# rExSys...
# list of indep vars for Poisson eq [x,y,z,u])
# list of dep vars  [xi, eta, zeta, phi]
# Input: rif-form rExSys of the exact system output of the first stage of the HybridGIF algorithm
# Derivative order that rExSys becomes involutive when prolonged to this order
# 	#predict q1 (order of prolongation), order Involutivity exact system

LHrExSys := map(lhs, DEtools:-rifsimp(ExSys)[Solved]);
userinfo(2, NumericJetTools,`listexactsys`, LHrExSys);
# Aug 31, 2022
InitialdataExSys := DEtools:-initialdata(ExSys);
HilbertFuncExSys := DifferentialHilbertFunction(InitialdataExSys, s);
userinfo(2, NumericJetTools,`DHS`,HilbertFuncExSys); 

PLHrExSys:= PolynomialTools:-PDEToPolynomial(LHrExSys, vars[1], vars[2]);
	userinfo(2, NumericJetTools,`PLHrExSys=`,PLHrExSys);   
	d := [];
	MindegreeDep := [];
	## this loop will be used for finding lower bound
	# InvIndepMulti := 1;
	# for i in vars[1] do
	# 	InvIndepMulti := InvIndepMulti*(i)^(-1);
	# end do;
	InvIndepMulti := -nops(vars[1]);
	#lprint(`inversemultiindep =` , InvIndepMulti);
	## this loop provids list of monomals regards to each dependend variables and lcm
	for j from 1 to nops(vars[2]) do
		leastComM := [];
		if nops(vars[2]) > 1 then
			Listmondep := select(has, PLHrExSys, vars[2][j]);
			#lprint(`Listmondep =`, Listmondep);
			Remdep :=	map(z -> z/(vars[2][j]), Listmondep);
		else 
			Remdep := PLHrExSys; 
		end if; 
	
		if nops(Remdep) = 1 then 
			leastComM := [op(leastComM), op(Remdep)];
			MaxdegreeDep := map(degree, leastComM);
			d := [op(d),max(MaxdegreeDep)];
			MindegreeDep := [op(MindegreeDep), degree(op(leastComM))];
			#lprint(`One leading for dep var`, vars[2][j], `with order =`, degree(op(leastComM)));
		else 
			for k from 1 to nops(Remdep)-1 do
 	 			for n from k+1 to nops(Remdep) do 
  	 	 			leastComM := [op(leastComM), lcm(Remdep[n], Remdep[k])];
   				end do;
			end do;
			MaxdegreeDep := map(degree, leastComM);
			d := [op(d), max(MaxdegreeDep) - 1];
			leastComDepL := map(z->InvIndepMulti*z, leastComM);
			userinfo(2, NumericJetTools,`leastComDepL =`, leastComDepL); #userinfo
			mono := map(z->z*vars[2][j],leastComDepL);
			userinfo(2, NumericJetTools,`mono =`, mono);
			SolvedrifExact := DEtools:-rifsimp(ExSys)[Solved];
			userinfo(2, NumericJetTools,`rExSys =`, DEtools:-rifsimp(ExSys)[Solved]);
  		 	Monomialdep := DetectBorder(mono,SolvedrifExact,vars);
			#lprint(`MonomialdepOnbou =`, Monomialdep[MonoOnbound]);
			MindegreeDep := [op(MindegreeDep), Monomialdep[MininumDegree]];
			#lprint(`MindegreeDep =`, MindegreeDep);
		end if;
		userinfo(2, NumericJetTools,`leastComM=`, leastComM);
		# MaxdegreeDep := map(degree, leastComM);
		# if nops(Remdep) = 1 then
		#     d := [op(d),max(MaxdegreeDep)];
		# else
		#     d := [op(d), max(MaxdegreeDep) - 1];
		# end if;
		userinfo(0, NumericJetTools,`Maxdegreelist =`, d);
		userinfo(0, NumericJetTools,`MinDegreelist =`, MindegreeDep);
	end do;
	MansfieldInvDeg := max(d);
	#lprint(`MaxdegreeEx =`, MaxDegreeEx);

	# Aug 31, 2022
	if type(HilbertFuncExSys[3], polynom) then
	    ReducedInvDeg := degree(HilbertFuncExSys[3])+1;
	else
	    ReducedInvDeg := max(MindegreeDep);
	end if;


	# MinDegreeEx := max(MindegreeDep);
	#lprint(`MinDegreeEx =`, MinDegreeEx);
	# HFExSys := DifferentialHilbertFunction(IDExSys, s);
	# SimpApSys := dsubs(rExSys[Solved], ApSys);
	# q2 := difforder(SimpApSys);
	return table([MansfieldInvolutiveDegree = MansfieldInvDeg, ReducedInvolutiveDegree = ReducedInvDeg]);
end proc:

####

DetectBorder:= proc(mono::list, sys::list, vars::list) 
	## Input: mono is a list of monomial after times 1/vars[1] icluds related dep variable, sys is a rif form of the exact system, vars list [[indeps],[deps]]
	## Output :
	local Notonbound, Onbound, m, DleastComDep, DLCM, Mindegree, syszero, i, monoaftertrans,
	indepvardegreelist;
	Notonbound := [];
	Onbound := [];
	indepvardegreelist := [];    # Aug 18, 2022
	syszero := map(z->lhs(z)=0, sys);
	userinfo(0, NumericJetTools,`syszero =` , syszero);
	for m from 1 to nops(mono) do
	    
		monoaftertrans := mono[m];
		for i from 1 to nops(vars[1]) do
		    monoaftertrans := monoaftertrans / vars[1][i];
			indepvardegreelist := [degree(monoaftertrans,vars[1][i]),op(indepvardegreelist)];    # Aug 18, 2022
		end do;
		DleastComDep := PolynomialTools:-PolynomialToPDE([monoaftertrans], vars[1], vars[2]);
		userinfo(0, NumericJetTools,`monomial after transformation =` , monoaftertrans);
		userinfo(0, NumericJetTools,`degree for variables =` , indepvardegreelist);
		userinfo(0, NumericJetTools,`indets =` , indets(indepvardegreelist ,negative));
		if indets(indepvardegreelist ,negative) = {} then
			#mono := lcmdep[m]*vars[2];
			# DleastComDep := PolynomialTools:-PolynomialToPDE([mono[m]], vars[1], vars[2]);
			DLCM := PDEtools:-dsubs(syszero,DleastComDep );
			userinfo(0, NumericJetTools,`DLCM =` , DLCM);
			if DLCM = 0 then
				Notonbound := [op(Notonbound), PolynomialTools:-PolynomialToPDE([mono[m]], vars[1], vars[2])];
			else 
				Onbound := [op(Onbound), PolynomialTools:-PolynomialToPDE([mono[m]], vars[1], vars[2])];
			end if;
		else 
			Onbound := [op(Onbound), PolynomialTools:-PolynomialToPDE([mono[m]], vars[1], vars[2])];
		end if;
	end do;
	userinfo(0, NumericJetTools,`Onbound after subtracking one degree of indep variables =` , Onbound);	
	userinfo(0, NumericJetTools,`Notonbound after subtracking one degree of indep variables =` , Notonbound);
	Mindegree := PDEtools:-difforder(Onbound)-1;
	return table([MonoOnbound = Onbound, MonoNotOnbound = Notonbound , MininumDegree = Mindegree]);
end proc:


# DetectBorder:= proc(mono::list, sys::list, vars::list) 
# 	## Input: mono is a list of monomial after times 1/vars[1] icluds related dep variable, sys is a rif form of the exact system, vars list [[indeps],[deps]]
# 	## Output :
# 	local Inside, Onbound, m, DleastComDep, DLCM, Mindegree, i, OnboundDeg;
# 	Inside := [];
# 	Onbound := [];
# 	OnboundDeg := [];
# 	for m from 1 to nops(mono) do
# 		if indets(degree(mono[m])-1 ,negative) = {} then
# 			#mono := lcmdep[m]*vars[2];
# 			DleastComDep := PolynomialTools:-PolynomialToPDE([mono[m]], vars[1], vars[2]);
# 			userinfo(2, NumericJetTools, `DleastComDep = `, DleastComDep);
# 			DLCM := PDEtools:-dsubs(sys,DleastComDep );
# 			userinfo(2, NumericJetTools, `DLCM = `, DLCM);
# 			if DLCM <> DleastComDep then
# 				Onbound := [op(Onbound), DleastComDep ];
# 			else 
# 				Inside := [op(Inside), DleastComDep ];
# 			end if;
# 		else 
# 			Inside := [op(Inside), lcmdep[m] ];
# 		end if;
# 	end do;
# 	userinfo(2, NumericJetTools,`Onbound after subtracking one degree of indep variables =` , Onbound);	
# 	userinfo(2, NumericJetTools,`Inside after subtracking one degree of indep variables =` , Inside);
# 	userinfo(2, NumericJetTools,`Onbound order =` , PDEtools:-difforder(Onbound));	
# 	if Inside <> [] then
# 	    Mindegree := PDEtools:-difforder(Inside) + nops(vars[1]) - 1;
# 	else
# 	    for i from 1 to nops(Onbound) do
# 		    OnboundDeg := [op(OnboundDeg), PDEtools:-difforder(Onbound[i])];
# 	    end do;
# 	    Mindegree := min(OnboundDeg);
# 	end if;
# 	return table([MonoOnbound = Onbound, MonoInside = Inside , MininumDegree = Mindegree]);
# end proc:



DifferentialHilbertFunction := proc(ID, S)
local DHF, w, K, j, s, L, nfree, Dimens, DiffDim;	
#option encrypted;	
uses difforder=PDEtools:-difforder,initialdata=DEtools:-initialdata, rifsimp=DEtools:-rifsimp; 
	
# Note ID must be ID for a system R, in leading linear form, with no leading nonlinear DE 
# Finite part
	# DiffDim := 0; Dimens  := 0;
	s := S;	
	userinfo(5, MapDETools, `Enter Diff Hilbert Function`); 	
	if  ID[Finite] <> [] then
		K := map(w -> difforder((lhs)(w)), ID[Finite]);
		DHF[Finite] := add(s^K[j], j = 1 .. nops(K));
		userinfo(5, MapDETools, `Finite data part of Differential Hilbert Function =`, DHF[Finite]);
	else DHF[Finite] := 0;
	end if;
	if ID[Infinite] <> [] then
		L     := map(w -> difforder((lhs)(w)), ID[Infinite]);
		nfree := map(w -> nops([op((rhs)(w))]), ID[Infinite]);
		userinfo(5, MapDETools,`L=`, L, `nfree =`, nfree);
		DHF[Infinite] := sum(s^L[k]*diff((1-s)^(-1), [s$(nfree[k]-1)])/(factorial(nfree[k]-1)),k = 1 .. nops(L));
		DiffDim := max(nfree); Dimens  := infinity;
		userinfo(5, MapDETools,`Infinite data part of Differential Hilbert Function =`, DHF[Infinite]);
		userinfo(5, MapDETools, `DiffDim =`, DiffDim, `Dimens=`, Dimens):
	else DHF[Infinite] := 0; DiffDim := 0;
	end if;
		
	if 	ID[Finite] = [] and ID[Infinite] = [] then 
		Dimens := 0; DiffDim := 0;
		DHF[Finite]:= 0; DHF[Infinite]:= 0;
		return [Dimens, DiffDim, DHF[Finite] + DHF[Infinite]];
	elif ID[Finite] <> [] and ID[Infinite] = [] then 
		Dimens := nops(ID[Finite]); DiffDim := 0;
		DHF[Infinite]:= 0;
		return [Dimens, DiffDim, DHF[Finite] + DHF[Infinite]];
	elif ID[Finite] = [] and ID[Infinite] <> [] then 
		Dimens := infinity;
		DHF[Finite]:= 0;
		return [Dimens, DiffDim, DHF[Finite] + DHF[Infinite]];   
	else 
		return [Dimens, DiffDim, DHF[Finite] + DHF[Infinite]];
	end if;	
	
end proc:

HybridGeometricInvolutiveForm := proc(sys::list,vars::list,splitvar)
# add third argument which is either float or a list of function, for exmaple G(x,y,z), which can be used to find system placed in ApSys
## Input:  sys is a list of polynomially linear homogeneous differential equations
## vars is a list of [[x], [u]] where [x] is a list of independent variable names and [u] is a list of dependent variable names.
## option encrypted;
uses difforder=PDEtools:-difforder,initialdata=DEtools:-initialdata, rifsimp=DEtools:-rifsimp, dsubs=PDEtools:-dsubs;
local ExSys, ApSys, i, AppCoeff, rExSys, flag, SplitSys;
		## IDExSys,LHrExSys, PLHrExSys, d, MindegreeDep, InvIndepMulti, leastComM, Listmondep, 
		## Remdep, j, k, n, leastComDepL, SolvedrifExact, Monomialdep, MaxdegreeDep, MaxDegreeEx, MinDegreeEx,
		## HFExSys, SimpApSys, q2, mono; 

	ExSys := [];
	ApSys := sys;
	flag := false;
	while flag = false do
	    SplitSys := SplitExApSys(ExSys, ApSys, flag, splitvar);
		ExSys := SplitSys[ExactSystem];
		userinfo(2, NumericJetTools,`ExSys =`, ExSys);
		ApSys := SplitSys[ApproximateSystem];
		flag := SplitSys[flagname];
		userinfo(2, NumericJetTools,`flag =`, flag);
	end do;
	return [ExactSystem = ExSys, ApproximateSystem = ApSys];
end proc: