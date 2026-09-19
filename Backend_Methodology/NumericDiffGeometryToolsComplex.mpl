userinfo(0, 'NumericDiffGeometryTools', `NumericDiffGeometryTools`);     

Digits := 15:

with(PDEtools): with(DEtools):
with(LinearAlgebra):
with(PolynomialTools):

#  Short procedures

N  := (n, m, d)    -> m*binomial(n+d  , d):
Nd := (n, m, d)    -> m*binomial(n+d-1, d):
Nc := (n, m, d, k) -> m*binomial(n+d-k-1, d-1):
Ncd := proc(n, m, d, k)
     local s, ks;
     # s := sum(Nc(n,m,d,ks), ks = 1 .. (k-1) );
     # = # vars of degree d of class < k
     # = # vars class 1, 2, ... , k-1
     s := 0;
     for ks from 1 to (k-1) do
     	s := s + m*binomial(n+d-ks-1, d-1);
     end do;
     return s;
end proc:

################################# PROLONGATION ROUTINES ###################################
Deriv := proc(L,vars)
  local newL,l,v,nder;
     newL := []:
     for l in L do
        for v in vars do
           nder:=diff(l,v):
           if not member(nder,newL) then
              newL := [op(newL),nder]:
           end if;
        end do;
     end do;
     newL;
 end:

Prolong := proc (sys::list, vars, r)
# Examples:
# Prolongs the system of derivative order dordsys to an equivalent one of order (dordsys + r)
# Note: it also differentiates differential equations of order less than
# dordsys up to (dordsys + r)
# Prolong([xi(x,y,z), diff(eta(x,y,z),z)], [x,y,z], 1);
# Prolong([xi(x),diff(xi(x),x,x)+xi(x)], [x], 3);
   local i, j, T, P, dordsys, dord;
   dordsys := PDEtools[difforder](sys);
   for i to nops(sys) do
       dord := PDEtools[difforder](sys[i]);
       P[i] := sys[i]; T := sys[i];
       for j to r+dordsys-dord do
           T := op(Deriv([T], vars));
           P[i] := P[i], T
       end do
    end do;
    return [seq(P[i],i = 1 .. nops(sys))];
end proc:


ReorderMtx := proc(A::Matrix, trans)
	local nrows, ncols, i, j, ColIndexFunc, IndexFunc;
		
	nrows := LinearAlgebra[RowDimension](A);
	ncols := LinearAlgebra[ColumnDimension](A);
	
	if args[2]=RevColOrder then	
		ColIndexFunc := (i,j) -> A[i, ncols + 1 - j];
		return Matrix(nrows, ncols, ColIndexFunc);
	end if;
	
	if args[2]=ReflectInAntiDiagonal then
		if nrows<>ncols 
			then error `Matrix not square, can't do this transform on non-square mtx` ;
		end if;
		IndexFunc := (i, j) -> A[nrows + 1 - j, nrows + 1 - i];
		return Matrix(nrows, IndexFunc);
	end if;
end proc:

PommaretRankingMtx := proc(indeps, deps)
	local n, m, i, j, RM;
	n := nops(indeps);
	m := nops(deps);
	RM := Matrix(n+m, n+m, datatype = integer[8]);
	for j from 1 to n do RM[1,j] := 1; end do;
	for j from n+1 to n+m do RM[1,j] := 0; end do;
	for i from 2 to n do RM[i,i-1] := -1; end do;
	for i from n+1 to n+m do RM[i,i] := 1; end do;
return convert(RM, listlist);
end proc:


# Reduction procedure 
GeoInvBasis := proc(A::Matrix, vardeg, TolSpec, RandomPointSpec)
	local NumRank, d, UU, sigma, VVt, j, k, r, dimKer, C, Brsp, Bker,
	dimKerCheck, dimKerCheckL, 
	n, m, N, Nd, Nc, B, randInt, randPoint, randpt,
	indeps, deps, tol,
	dimSymb, rankSymb,
	dj, sumbeta, beta,
	numpr, maxpr, 
	BetaM, SpaceM, RankSymM, BrspM, BkerM, SigmaM,
	DimKerM, DimRspM, DataTypeSpec, 
	BetaMaxM, rowd, cold, Asymspan, Asym, RankAsym, Ncd, DimSymM, SymM,
	JetSp, NullSt, # Aug 18, 2015
	ClassCartanTest, # Aug 19, 2015
	degenerate, dgen; # Aug 25, 2015

	indeps := vardeg[1];
	deps	  := vardeg[2];
	d	  := vardeg[3];
	tol    := rhs(TolSpec);
	randpt := rhs(RandomPointSpec);
	maxpr := 0;
	numpr := 0;
	
     n := nops(indeps);
     m := nops(deps);
     
     N  := (n, m, d)    -> m*binomial(n+d  , d);
     Nd := (n, m, d)    -> m*binomial(n+d-1, d);
     Nc := (n, m, d, k) -> m*binomial(n+d-k-1, d-1);
     Ncd := proc(n, m, d, k)
     	local s, ks;
     	# s := sum(Nc(n,m,d,ks), ks = 1 .. (k-1) );
     	# = # vars of degree d of class < k
     	# = # vars class 1, 2, ... , k-1
     	s := 0;
     	for ks from 1 to (k-1) do
     		s := s + m*binomial(n+d-ks-1, d-1);
     	end do;
     	return s;
	end proc:

	if Digits < 15 then userinfo(1, 'NumericDiffGeometryTools',`Warning: to take advantage of hardware float computation speed, you should set Digits := 15 ` ); end if;

     if Digits <= 15 then DataTypeSpec := complex[8];  # changed float[8] to complex[8]
     	else 
     	userinfo(1, 'NumericDiffGeometryTools',`Since Digits = `, Digits, `> 15 the program will execute in software floats (very slowly!!)`);
     	userinfo(1, 'NumericDiffGeometryTools',`It will execute in the higher number of Digits.  Your tolerance should be compatible with this choice`);
     	userinfo(1, 'NumericDiffGeometryTools',`For example:  Tolerance = 1-e(Digits - 5) `);
     	DataTypeSpec := anything; 
     end if;

	# d = diff degree of the input system
	if ColumnDimension(A) <> N(n,m,d) 
		then error(`Your input matrix should have num columns =`, N(n,m,d)) 	
	end if;
	userinfo(1, 'NumericDiffGeometryTools',`You have entered a matrix for`, RowDimension(A), `eqns with system degree`, d, `independent vars`, indeps, `dependent vars`, deps );

	DimKerM  := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	DimRspM  := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	SigmaM   := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	RankSymM := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	DimSymM  := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);	
	SymM     := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);	
	BrspM    := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	BkerM    := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	BetaM    := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);
	BetaMaxM := Matrix(1 .. maxpr+1, 1 .. d+1, datatype=anything, fill = ` `);

	
	NumRank := proc(sigmav::Vector[column], tol) 
	# NumRank returns Numerical rank based on a column vector of sing vals
	# Alternative criterion sigma[j][k]/sigma[j][1] > tol
		local k, rk;
		rk := 0;
		for k from 1 to RowDimension(sigmav) while sigmav[k] > tol do  rk:= k;  end do;
		return rk;
	end proc;
	# Substitute an appropriate random point in B
	randInt := rand(90000 .. 100000);
	if rhs(RandomPointSpec)=false then 
     	randPoint := [seq(indeps[j] = evalf(1/100000*randInt()), j = 1 .. n)];
     	else 
     	randPoint :=  [seq(indeps[j] = evalf(randpt[j]), j = 1 .. n)];
     end if;
     userinfo(1, 'NumericDiffGeometryTools',`randPoint = `, randPoint );
     B := Matrix(convert(evalf(subs(randPoint, A)), listlist),  datatype=DataTypeSpec, storage=rectangular);
     userinfo(2,'NumericDiffGeometryTools', ` Matrix form evaluated at rand. point. `, randPoint);
	# userinfo(1, 'NumericDiffGeometryTools',`A with rand point substituted = `, B);
	# degenerate is a flag to distinguish between the degenerate(true) and non-degenerate(false) cases
	degenerate := false;
	# dgen = the order at which case becomes degenerate, so 0 <= dgen <= d
	# the defaut value is -1
	dgen := -1;
	userinfo(3,'NumericDiffGeometryTools',`Entering SVD at iteration for degree =`, d);
	UU[d], sigma[d], VVt[d] := SingularValues(B, output = ['U', 'S', 'Vt'], 
				outputoptions['Vt'] = [datatype = DataTypeSpec, storage = rectangular]);
	userinfo(1, 'NumericDiffGeometryTools',`Sigma = `, sigma[d] );
	# Determine the numerical rank of C = rowsp rank
   	r := NumRank(sigma[d], tol);  		
   	# userinfo(1, 'NumericDiffGeometryTools',`At degree = `, d, `dim ker = `, N(n,m,d) - r, `dim rowspace =`, r);
   	dimKer[d] := N(n,m,d) - r;
   	DimKerM[numpr+1, d+1] := dimKer[d];
   	DimRspM[numpr+1, d+1] := r;
   	SigmaM[numpr+1, d+1] := sigma[d];
   	# Generic case Bker[d] <> NullSt and Brsp[d] <> JetSp
   	# There is one further theoretical case:  r = 0 = N(n,m,d); not acccounted for yet
	# Bases for rowsp and ker are for order d
   	Bker[d] := VVt[d][(r+1) .. N(n,m,d)];   Brsp[d] := VVt[d][1 .. r];
   	BkerM[numpr+1, d+1] := Bker[d];         BrspM[numpr+1, d+1] := Brsp[d]; 
	# A span set for the proj of Bker by deleting components for order d - 1 is
   	C := VVt[d][(r+1) .. N(n,m,d), 1 .. N(n,m,d-1) ];
	if r = N(n,m,d) or r = 0 then
		degenerate := true;
		dgen := d;
		dimKer[d] := r;
   		DimKerM[numpr+1, d+1] := dimKer[d];
   		DimRspM[numpr+1, d+1] := N(n,m,d) - r;
   		SigmaM[numpr+1, d+1] := sigma[d];
		if r = N(n,m,d) then
			Bker[d] := JetSp(d);                 Brsp[d] := NullSt;
   			BkerM[numpr+1, d+1] := Bker[d];    BrspM[numpr+1, d+1] := Brsp[d];
   		else
   			Bker[d] := NullSt;               Brsp[d] := JetSp(d);
   			BkerM[numpr+1, d+1] := Bker[d];    BrspM[numpr+1, d+1] := Brsp[d];
   		end if;
   		userinfo(2,'NumericDiffGeometryTools',`first exceptional case at diff order = `, d);
	end if;

   	userinfo(3,'NumericDiffGeometryTools', ` Bker[d]= `, evalf(Bker[d],2)); 	     
	userinfo(3,'NumericDiffGeometryTools', ` C = `, C);
	
	for j from d-1 to 0 by -1 while degenerate = false do
		# userinfo(1, 'NumericDiffGeometryTools',`Starting iteration for degree = `, j-1);
		UU[j], sigma[j], VVt[j] := SingularValues(C, output = ['U', 'S', 'Vt'], 
				outputoptions['Vt'] = [datatype = DataTypeSpec, storage = rectangular]);
		# r = Numerical rank of C
   		r := NumRank(sigma[j], tol);
   		userinfo(3,'NumericDiffGeometryTools', ` numrank = `, r);
   		dimKer[j] := r;
   		DimKerM[numpr+1, j+1] := dimKer[j];
   		DimRspM[numpr+1, j+1] := N(n,m,j) - r;
   		dimSymb[j+1] := dimKer[j+1] - dimKer[j];
   		DimSymM[numpr+1, j+2] := dimSymb[j+1];
   		rankSymb[j+1]:= Nd(n,m,j+1) - dimSymb[j+1];
   		RankSymM[numpr+1, j+2] := rankSymb[j+1];
   		userinfo(4, 'NumericDiffGeometryTools',`At degree = `, j+1, `dim Sym = `, dimSymb[j+1] , `rank Symbol =`, rankSymb[j+1]);
   		#userinfo(1, 'NumericDiffGeometryTools',`At degree = `, j+1, `dim Sym = `, dimSymb[j+1] , `rank Symbol =`,  (N(n,m,j+1)-dimKer[j+1])-(N(n,m,j)-dimKer[j])); 
   		userinfo(4, 'NumericDiffGeometryTools',`At degree = `, j, `  dim ker = `, r, `dim rowspace =`, N(n,m,j) - r);  	
   		# if r = N(n,m,j) or r = 0 then degenerate := true; dgen := j; break; end if;	
   		if r = N(n,m,j) or r = 0 then
			degenerate := true;
			dgen := j;
			if r = N(n,m,j) then
				Bker[j] := JetSp(j);               Brsp[j] := NullSt;
   				BkerM[numpr+1, j+1] := Bker[j];    BrspM[numpr+1, j+1] := Brsp[j];
   			else
   				Bker[j] := NullSt;                 Brsp[j] := JetSp(j);
   				BkerM[numpr+1, j+1] := Bker[j];    BrspM[numpr+1, j+1] := Brsp[j];
   			end if;
   			userinfo(1, 'NumericDiffGeometryTools',`Just before break at diff order = `, j, `dgen =`, dgen);
   			break;
		end if;

		# Bases for rowsp and ker are:
   		Bker[j] := VVt[j][1 .. r];         Brsp[j] := VVt[j][(r+1) .. N(n,m,j)];
   		BkerM[numpr+1, j+1] := Bker[j];   	BrspM[numpr+1, j+1] := Brsp[j];  
   		SigmaM[numpr+1, j+1] := sigma[j];
   		# A span set for the proj of Brsp by deleting components is
		C := VVt[j][1 .. r, 1 .. N(n,m,j-1)];
		userinfo(3,'NumericDiffGeometryTools', ` N(n,m,j-1) `, N(n,m,j-1));
		userinfo(4,'NumericDiffGeometryTools', ` C = `, C);
	end do;

	# Note once if degenerate = true on entry to the following loop, it remains true.
	# once degenerate = true, then 0 <= dgen <= d 
	for j from dgen-1 to 0 by -1 do	
		if r = N(n,m,dgen) then
			dimKer[j] := N(n,m,j);
   			DimKerM[numpr+1, j+1] := dimKer[j];
   			DimRspM[numpr+1, j+1] := 0;
   			dimSymb[j+1] := dimKer[j+1] - dimKer[j];
   			DimSymM[numpr+1, j+2] := dimSymb[j+1];
   			rankSymb[j+1]:= Nd(n,m,j+1) - dimSymb[j+1];
   			RankSymM[numpr+1, j+2] := rankSymb[j+1];	
   			Bker[j] := JetSp(j);              Brsp[j] := NullSt;
   			BkerM[numpr+1, j+1] := Bker[j];   BrspM[numpr+1, j+1] := Brsp[j];
   	     else
   			if r<>0 then error `r should be zero, but is = `, r end if;
   			dimKer[j] := 0;
   			DimKerM[numpr+1, j+1] := dimKer[j];
   			DimRspM[numpr+1, j+1] := N(n,m,j);
   			dimSymb[j+1] := dimKer[j+1] - dimKer[j];
   			DimSymM[numpr+1, j+2] := dimSymb[j+1];
   			rankSymb[j+1]:= Nd(n,m,j+1) - dimSymb[j+1];
   			RankSymM[numpr+1, j+2] := rankSymb[j+1];
   			Bker[j] := NullSt;                Brsp[j] := JetSp(j); 
   			BkerM[numpr+1, j+1] := Bker[j];   BrspM[numpr+1, j+1] := Brsp[j];
   		end if; 
   	end do; 	
	
	# A special case, not covered by the above loop, is to compute rankSymbol[0] and dimSymb[0]
	dimSymb[0] := dimKer[0];
	DimSymM[numpr+1, 1] := dimSymb[0];
	rankSymb[0]:= Nd(n,m,0) - dimSymb[0];
	RankSymM[numpr+1, 1] := rankSymb[0];
		
	# Compute the maximal Cartan beta Characters for involutivity of the symbol
	for dj from d to 0 by -1 do
		userinfo(3, 'NumericDiffGeometryTools',`Entering max beta loop`, ` dj= `, dj);
		sumbeta := 0;
		beta := NULL ;
		for k from n to 1 by -1 do
			userinfo(3, 'NumericDiffGeometryTools',` [dj, k] = `, [dj,k], `Nc(n,m,dj,k) =`, Nc(n,m,dj,k), `rankSymb[dj] =`, rankSymb[dj] );
			if    Nc(n,m,dj,k) <= rankSymb[dj] - sumbeta then  beta := Nc(n,m,dj,k), beta; sumbeta := sumbeta + Nc(n,m,dj,k);
			elif  Nc(n,m,dj,k) >  rankSymb[dj] - sumbeta and rankSymb[dj] > sumbeta then beta := rankSymb[dj] - sumbeta, beta; sumbeta := rankSymb[dj];
			else  beta := 0, beta;
			end if;
		end do;
		j := 'j';
		userinfo(4, 'NumericDiffGeometryTools',` j = `, j, `beta =`, beta);
		BetaMaxM[numpr + 1, dj+1] :=  [[beta], sum(j*[beta][j], j = 1 .. n)];
		userinfo(4, 'NumericDiffGeometryTools',`BetaMaxM ` , [numpr+1, dj+1] , ` = `, BetaMaxM[numpr + 1, dj+1] );			
	end do;
	beta := 'beta':
	
	# Compute the Cartan characters using the classical test
	ClassCartanTest := true;
	for dj from d to dgen+1 by -1 while ClassCartanTest = true do
		userinfo(3, 'NumericDiffGeometryTools',`Entering loop to compute Cartan's betas at degree`, ` dj= `, dj);
		userinfo(3, 'NumericDiffGeometryTools',`We will compute a span set for rowspace of the symbol mtx denoted by SymM[numpr+1,dj+1]`);
		rowd, cold := Dimensions(BrspM[numpr+1,dj+1]);
		userinfo(3, 'NumericDiffGeometryTools',`row and col dimensions of BrspM[numpr+1,dj+1] are`, rowd, cold); 
		if rowd=0 then # betas are all zero
			userinfo(3, 'NumericDiffGeometryTools',`Wow easy ... betas all zero `);
			BetaM[numpr + 1, dj+1] :=  [[seq(0, j = 1 .. n)], 0];
			userinfo(3, 'NumericDiffGeometryTools',`Row space for the symbol mtx is degenerate 0 x Nd(n,m,dj) mtx`);
			SymM[numpr+1,dj+1] := BrspM[numpr+1,dj+1][1 .. rowd, (N(n,m,dj-1)+1) .. N(n,m,dj)];
		else
			userinfo(3, 'NumericDiffGeometryTools',`We will compute a span set at degree dj for rowspace of the symbol mtx denoted by SymM[numpr+1,dj+1]`);
			userinfo(3, 'NumericDiffGeometryTools',`By deleting components of degree < dj from BrspM[numpr+1,dj+1] `);
			userinfo(3, 'NumericDiffGeometryTools',` The resulting matrix Asymspan, represents a span set for the row space of the symbol mtx at degree dj `);
			Asymspan := BrspM[numpr+1,dj+1][1 .. rowd, (N(n,m,dj-1)+1) .. N(n,m,dj)];
			SymM[numpr+1,dj+1] := BrspM[numpr+1,dj+1][1 .. rowd, (N(n,m,dj-1)+1) .. N(n,m,dj)];
			if ColumnDimension(Asymspan)<>Nd(n,m,dj) then error(`Col dimension of Asymspan <> Nd(n,m,dj) `); end if;
			userinfo(2, 'NumericDiffGeometryTools',`row and col dimensions of Span set for Asymspan(BrspM) are`, Dimensions(Asymspan) );
			userinfo(3, 'NumericDiffGeometryTools',`Prepare to compute NumRank(Asymspan) .. by deleting relevant class var cols from Asymspan `);
			
			Asym[n]    := Asymspan[1 .. rowd, (Ncd(n,m,dj,n)+1) .. Nd(n,m,dj)];			
			sigma[n]   := SingularValues(Asym[n], output = ['S'], outputoptions['S'] = [datatype = DataTypeSpec, storage = rectangular]);				
			RankAsym[n]:= NumRank(sigma[n], tol);
			beta[n] := RankAsym[n];
			userinfo(3, 'NumericDiffGeometryTools',`Asym`,[n], ` = `,  Asym[n], `sigma[n] = `, sigma[n],`RankAsym[n] = beta[n] = `, RankAsym[n] );
				
			for k from n to 2 by -1 do
				userinfo(3, 'NumericDiffGeometryTools',` [dj, k] = `, [dj,k], `Start loop to calculate beta `, ` Ncd(n,m,dj,k-1)+1  = `, Ncd(n,m,dj,k-1)+1);
				# Delete cols from Asymspan corr to class k-1 or less
				Asym[k-1]    := Asymspan[1 .. rowd, (Ncd(n,m,dj,k-1)+1) .. Nd(n,m,dj)];
				# userinfo(1, 'NumericDiffGeometryTools',`Asym`,[k-1], ` = `,  Asym[k-1]);			
				sigma[k-1]   := SingularValues(Asym[k-1], output = ['S'], outputoptions['S'] = [datatype = DataTypeSpec, storage = rectangular]);				
				RankAsym[k-1]:= NumRank(sigma[k-1], tol);			
				userinfo(4, 'NumericDiffGeometryTools',`Asym`,[k-1], ` = `,  Asym[k-1], `sigma[k-1] = `, sigma[k-1],`RankAsym[k-1] =`, RankAsym[k-1] );					
				# Calculate sum(beta[ks], ks = k .. n) = rank SpanSym(cols cls > k-1) = rank[k]
				# Therefore beta[n] = rank[n], beta[n-1]+beta[n] = rank[n-1] so beta[n-1] = rank[n-1]-beta[n]
				# beta[n-2]+beta[n-1]+beta[n] = rank[n-2] so beta[n-2] = rank[n-2] - beta[n] - beta[n-1] = rank[n-1] - rank[n-2]
				beta[k-1] := RankAsym[k-1] - RankAsym[k];
			end do;
			# if n=1 then beta[1] := RankSymM[numpr+1,dj+1] else beta[1] := RankSymM[numpr+1,dj+1] - RankAsym[1]; end if;
			userinfo(3, 'NumericDiffGeometryTools',`beta[1] = `, beta[1], `dj = `, dj, `RankSymM[numpr+1,dj+1] =`, RankSymM[numpr+1,dj+1], ` RankAsym[1] =`, RankAsym[1] );
			j := 'j';
			BetaM[numpr + 1, dj+1] :=  [[seq(beta[j], j = 1 .. n)], sum(j*beta[j], j = 1 .. n)];
			userinfo(4, 'NumericDiffGeometryTools',`BetaM ` , [numpr+1, dj+1] , ` = `, BetaM[numpr + 1, dj+1] );
		end if;			
	end do;

	# dodge the classical cartan test for now
	
	if ClassCartanTest=true and dgen = -1 then
		# Finally the case dj = 0 obviously corresponds to beta[j] = 0 all j
		BetaM[numpr + 1, 1] :=  [[seq(0, j = 1 .. n)], 0];
		userinfo(4, 'NumericDiffGeometryTools',`BetaM = `, BetaM);
		userinfo(4, 'NumericDiffGeometryTools',`At dj = 0 we compute a span set for rowspace of the symbol mtx denoted by SymM[numpr+1,1]`);
		rowd, cold := Dimensions(BrspM[numpr+1,1]);
		userinfo(4, 'NumericDiffGeometryTools',`Span set for Row space for the symbol mtx at dj = 0 is rowd x Nd(n,m,dj) mtx`);
		SymM[numpr+1,1] := BrspM[numpr+1,1][1 .. rowd, 1 .. N(n,m,dj)];
	end if;

	SpaceM := Matrix(1 .. numpr+1, 1 .. d+1, (i,j) -> Pi^(d+1-j) * D^(i-1) * R);
	userinfo(1, 'NumericDiffGeometryTools', `SpaceM =`, SpaceM);
	userinfo(1, 'NumericDiffGeometryTools',`For a system with indep variables indeps = [x1,x2,x3, ...] and dependent vars deps = [x1,x2,x3, ...] `);
	userinfo(1, 'NumericDiffGeometryTools',`And random point RandomPoint = randpt = [x10, x20, ... , xn0] `);
	userinfo(1, 'NumericDiffGeometryTools',`Information about specific systems are accessed by`);
	userinfo(1, 'NumericDiffGeometryTools',`G := GeoInvBasis(A, [indeps,deps, deg], Tolerance = 1e-9, ranking = PommaretRankingMtx(indeps, deps), RandomPoint = randpt); `);
	userinfo(1, 'NumericDiffGeometryTools',`And the details of the output are accessed via: `);
	userinfo(1, 'NumericDiffGeometryTools',`G[DimKerMtx];  G[DimRspMtx];  G[RankSymMtx]; G[BetaMtx];... etc `);
	userinfo(1, 'NumericDiffGeometryTools',`Where: `); 
	userinfo(1, 'NumericDiffGeometryTools',`DimKerMtx =  Matrix for dimensions of the kernels of the projections`);
	userinfo(1, 'NumericDiffGeometryTools',`DimRspMtx =  Matrix for dimensions of the rowspaces of the projections`);
	userinfo(1, 'NumericDiffGeometryTools',`BkerMtx   =  Matrix for bases of the kernels of the projections`);
	userinfo(1, 'NumericDiffGeometryTools',`BrspMtx    =  Matrix for bases of the rowspaces of the projections`);
	userinfo(1, 'NumericDiffGeometryTools',`DimSymMtx =  Matrix for dimensions of the kernels of the symbols`);
	userinfo(1, 'NumericDiffGeometryTools',`RankSymMtx = Matrix for ranks of of the symbols`);
	userinfo(1, 'NumericDiffGeometryTools',`SymMtx    =  Matrix for span set of the rowspaces of the symbols`);
	userinfo(1, 'NumericDiffGeometryTools',`SigmaMtx  =  Matrix for the singular values in the main SVD calculations`);
	userinfo(1, 'NumericDiffGeometryTools',`BetaMtx   =  Matrix for Cartan's beta's ... may need generic coordinates for this to be correct`);
			     
	return table([ `DimKerMtx` =  DimKerM,
				`DimRspMtx` =  DimRspM,
				`SigmaMtx`  =  SigmaM,
				`BrspMtx`   =  BrspM,
				`BkerMtx`   =  BkerM,
				`DimSymMtx` = DimSymM,
				`SymMtx`	  = SymM,
				`RankSymMtx`= RankSymM,
			     `BetaMaxMtx`=  BetaMaxM,
			     `BetaMtx`   =  BetaM ]);
end proc:

MatrixFormPDE := proc(sys::list, vars::list)
# Input:  system of polys (then vars = [x1, x2, ... , xn] if so return error maybe implemented in future ... 
#         system of linear pdes (should check that it is linear) and add the type checking vars = [indeps, deps]
#         system of pol nonlinear pdes (should check that it is poly nonlinear) and add the type checking vars = [indeps, deps]
#         that was put already in GIFd
		local indeps, deps, diffvars, Sys, j, A, b, t1, cranking, UserDefinedRanking,
			dord0, WP, LinDiffSysFlag, PolyDiffSysFlag, jetsys, jetvars,
			i;
			#Aug 12, 2022 add local i

		LinDiffSysFlag  := false;
		PolyDiffSysFlag := false;

		if not( map(whattype, vars)=[seq(symbol, j = 1 .. nops(vars))] or map(whattype, vars)=[list, list] ) then
		error ` vars specification for the system is incorrect ` 
		end if; 
		if map(whattype, vars)=[seq(symbol, j = 1 .. nops(vars))] then 
			userinfo(1, 'NumericDiffGeometryTools',`You have entered a system with vars = `, vars, `The input syntax for a system of polynomials `);
			if has(map(type, sys, polynom(anything, vars)), false) then error `The given system is not a system of polynomials in vars` 
			end if;
			error ` MatrixFormPDE not yet implemented for input systems of polynomials `;
		end if;

		# polysysflag := false;			
		indeps := vars[1];
		deps   := vars[2];
		userinfo(1, 'NumericDiffGeometryTools',`You have specified vars for a DE system with independent variables`, indeps , `and dependent variables`, deps);
			if not(map(whattype, indeps)=[seq(symbol, j = 1 .. nops(indeps))]) then error `Your specification of indeps is incorrect` end if;
			if not(map(whattype,   deps)=[seq(symbol, j = 1 ..   nops(deps))]) then error `Your specification of   deps is incorrect` end if;
		# Now we check its differential order and that its a linear system of PDE
		dord0 := PDEtools[difforder](sys);
		# All possible dependent variables and derivatives in the PDE system
		
		UserDefinedRanking := false;
		
		if nargs=3 # and lhs(args[3]) = 'ranking'
			then UserDefinedRanking := true;   cranking := rhs(args[3]);    
			userinfo(1,'MatrixFormPDE', ` User defined ranking = `, cranking);
		end if;
		if nargs=4 then 
			UserDefinedRanking := true;   cranking := rhs(args[3]);
			WP := args(4);			    
			userinfo(1,'MatrixFormPDE', ` User defined ranking = `, cranking, ` Point(s) entered = `, WP);
		end if;
		
	   	if nargs>=5 then error `Only 2, 3 oe 4 args allowed to MatrixFormPDE, you entered more than 4 `  end if;
	   
		if UserDefinedRanking=true then 
			diffvars := DEtools[checkrank]([seq(deps[j](op(indeps)), j = 1 .. nops(deps))], deps, indep = indeps, degree = dord0, ranking = cranking);
			else 
			diffvars := DEtools[checkrank]([seq(deps[j](op(indeps)), j = 1 .. nops(deps))], deps, indep = indeps, degree = dord0, ranking = PommaretRankingMtx(indeps, deps));
		end if;
		## diffvars := DEtools[checkrank]([seq(deps[j](op(indeps)), j = 1 .. nops(deps))], deps, indep = indeps, degree = dord0);
		# userinfo(1, 'NumericDiffGeometryTools',` diffvars = `, diffvars );
		# Check that the system is linear in diffvars
		## if map(type, sys, linear(diffvars))<>[seq(true, i = 1 .. nops(sys))] then error `The system is not a linear PDE system` end if;
		if map(type, sys, linear(diffvars))=[seq(true, i = 1 .. nops(sys))] then 
			userinfo(1, 'NumericDiffGeometryTools',`The system is a linear PDE system`);
			LinDiffSysFlag := true;
			elif map(type, sys, polynom(anything, diffvars))=[seq(true, i = 1 .. nops(sys))] then 
				userinfo(1, 'NumericDiffGeometryTools',`The system is a differential polynomial system of PDE`, map(type, sys, polynom(anything,diffvars)), diffvars);
			PolyDiffSysFlag := true
			else
			error ` System not a linear or a polynomially nonlinear system of DE `;	
		end if;
		####
		### if map(type, sys, polynom(anything, diffvars))=[seq(true, i = 1 .. nops(sys))] then 
		### 	userinfo(1, 'NumericDiffGeometryTools',`The system is a differential polynomial system of PDE`, map(type, sys, polynom(anything,diffvars)), diffvars);
		### end if;		
		
        	# dOrder := PDEtools[difforder](sys);
        	# userinfo(1,'NumericDiffGeometryTools', `Input system has order`, dOrder);

		UserDefinedRanking := false;

		if nargs=3 # and lhs(args[3]) = 'ranking'
			then UserDefinedRanking := true;   cranking := rhs(args[3]);    
			userinfo(2,'NumericDiffGeometryTools', ` User defined ranking = `, cranking);
		end if;
		
		userinfo(3,'NumericDiffGeometryTools', ` About to check nargs ` );
		
	   	if nargs>=4 then error `Only 2 or 3 args allowed to MatrixFormPDE, you entered more than 3 `  end if;
        
        	# Sys := map(z -> z = 0, sys);  # temp change Aug 13
        	Sys := sys;
		
		userinfo(3,'NumericDiffGeometryTools', ` Jet form of system is `, Sys);

		jetvars := map(ToJet, diffvars, [seq(deps[j](op(indeps)), j = 1 .. nops(deps))]);
		jetsys := map(ToJet, Sys, [seq(deps[j](op(indeps)), j = 1 .. nops(deps))]);
		userinfo(3,'NumericDiffGeometryTools', `jet vars = `, jetvars);
		userinfo(3,'NumericDiffGeometryTools', `jet sys = `, jetsys);
		
		if LinDiffSysFlag then 
			userinfo(3,'NumericDiffGeometryTools', `About to compute matrix form of linear DE system`);  t1:= time();
        		A, b := LinearAlgebra[LinearAlgebra:-GenerateMatrix](Sys, diffvars);   # temp change Aug 13
        		userinfo(3,'NumericDiffGeometryTools', `Matrix form of system computed in time`, time()-t1);
        		elif PolyDiffSysFlag then
        		userinfo(3,'NumericDiffGeometryTools', `About to compute Jacobian matrix form of poly diff DE system`);  t1:= time();
			# compute the jacobian of the sys
			A := VectorCalculus[Jacobian](jetsys, jetvars);
			userinfo(2,'NumericDiffGeometryTools', `JacobianMtx =`, jacobianMtx );
        	end if;
        			        
return (A, jetvars, b);
end proc:


VSpaces := proc(A::Matrix, vardeg, TolSpec, DataTypeSpec)
	local NumRank, d, UU, sigma, VVt, j, k, r, dimKer, C, Brsp, Bker,
	# dimKerCheck, dimKerCheckL, 
	n, m, N, Nd, Nc, B, randInt, randPoint, randpt,
	indeps, deps, tol,
	# dimSymb, rankSymb,
	# dj, sumbeta, beta,
	# numpr, maxpr, 
	BetaM, SpaceM, RankSymM, BrspM, BkerM, SigmaM,
	DimKerM, DimRspM,
	BetaMaxM, rowd, cold, Asymspan, Asym, RankAsym, Ncd, DimSymM, SymM;

	indeps := vardeg[1];
	deps	  := vardeg[2];
	d	  := vardeg[3];
	tol    := rhs(TolSpec);
	
     n := nops(indeps);
     m := nops(deps);
     
     N  := (n, m, d)    -> m*binomial(n+d  , d);
     Nd := (n, m, d)    -> m*binomial(n+d-1, d);

     NumRank := proc(sigmav::Vector[column], tol) 
	# NumRank returns Numerical rank based on a column vector of sing vals
	# Alternative criterion sigma[j][k]/sigma[j][1] > tol
		local k, rk;
		rk := 0;
		for k from 1 to RowDimension(sigmav) while sigmav[k] > tol do  rk:= k;  end do;
		return rk;
	end proc;

	UU, sigma, VVt := SingularValues(A, output = ['U', 'S', 'Vt'], 
				outputoptions['Vt'] = [datatype = DataTypeSpec, storage = rectangular]);
	# Determine the numerical rank of C = rowsp rank
   	r := NumRank(sigma, tol);  		
   	userinfo(1, 'NumericDiffGeometryTools',`At degree = `, d, `dim ker = `, N(n,m,d) - r, `dim rowspace =`, r);

   	# Bases for rowsp and ker are:
   	Brsp := VVt[1 .. r];  Bker := VVt[(r+1) .. N(n,m,d)];
  	userinfo(3, 'NumericDiffGeometryTools',`Brsp = `, Brsp,`dim VVt = `, Dimensions(VVt),  `Bker = `, Bker);
	return table([ `DimKerMtx` =  N(n,m,d) - r,
				`DimRspMtx` =  r,
				`SigmaMtx`  =  sigma,
				`BrspMtx`   =  Brsp,
				`BkerMtx`   =  Bker ]);
   	
end proc:


# This is a combined routine for both differential and polynomial systems 

GeometricInvolutiveForm := proc(sys::list, vars::list, tol, Dmax, RandomPointSpec, hybridoption:=false )
	local j, dsys, u, G, dord0, dord, diffvars, A, vrs, b, B, k, InvolSysIndex,
	polysysflag, indeps, deps, RandomPointList,
	Involflag, Involdeg, # Aug 31, 2015
	DimAtPoint, j1, k1,  # Aug 11, 2017 Aug 13, 2017
	i,                   # Aug 12, 2022
	HybridSysFlag,       # Aug 22, 2022
	dApSysModExSys,      # Aug 23, 2022
	dordModExSys,        # Aug 23, 2022
	dordExSys,           # Sep 20, 2022
	prolongtimes,
	                     # Aug 28, 2022 delete dApSys
	DHF, seriesDHF, DimSymEx, RankSymEx, RankSymAp, DimSymAp, DimSymSys, RankSymSys,                 # Jun 08, 2023
	newsys,              # Jun 11, 2023
	rankApSys, rankExSys,# Jun 14, 2023
	DimSol,q,            # Jun 17, 2023
	DimEx;               # Jun 24, 2023 
	# sys is a system of differential equations with vars := [indeps, deps]
	# or sys is a system of polynomials with vars := [x1, x2, ... ]
	# or sys is a hybrid system of exact part and approximate part  # Aug 23, 2022
	# 

    userinfo(0,'NumericDiffGeometryTools', `Geometric involutive form enterd`);

	if not( map(whattype, vars)=[seq(symbol, j = 1 .. nops(vars))] or map(whattype, vars)=[list, list]  ) then
		error ` vars specification for the system is incorrect ` 
	end if; 
    


	# Aug 22, 2022
    if map(whattype, sys)=[list, list] then
	    HybridSysFlag := true;
		newsys := sys;
	elif type(sys,list) then
    # Jun 11, 2023
	    if hybridoption=false then
            HybridSysFlag := false;
			newsys :=sys;
		else
		    newsys := HybridGeometricInvolutiveForm(sys, vars, float);
			if rhs(newsys[1]) = [] then
			    HybridSysFlag := false;
				newsys := sys;
			else
			    HybridSysFlag := true;
				newsys := [rhs(newsys[1]),rhs(newsys[2])];
			end if;
		end if;
	else
		error ` You have entered neither a list of a system nor a list of two lists of a hybrid system `   # Aug 23, 2022
    end if;

    userinfo(1, 'NumericDiffGeometryTools',`HybridSystemFlag`, HybridSysFlag, `new sys`, newsys);

	if map(whattype, vars)=[seq(symbol, j = 1 .. nops(vars))] then 
		userinfo(1, 'NumericDiffGeometryTools',`You have entered a system with vars = `, vars, `The input syntax for a system of polynomials `);
		if has(map(type, newsys, polynom(anything, vars)), false) then error `The given system is not a system of polynomials in vars` 
		end if;
		# userinfo(1, 'NumericDiffGeometryTools',`about to calculate dord0`);
		dord0 := max(map(degree, collect(newsys, vars)));
		userinfo(1, 'NumericDiffGeometryTools',` You have entered a system of polynomials of degree `, dord0, ` in the variables `, vars);
		polysysflag := true;
		indeps := vars;
		deps   := [u];
		else
			polysysflag := false;			
			indeps := vars[1];
			deps   := vars[2];
			userinfo(1, 'NumericDiffGeometryTools',`You have specified vars for a DE system with independent variables`, indeps , `and dependent variables`, deps);
			# if not(map(whattype, indeps)=[seq(symbol, j = 1 .. nops(indeps))]) then error `Your specification of indeps is incorrect` end if;
			# if not(map(whattype,   deps)=[seq(symbol, j = 1 ..   nops(deps))]) then error `Your specification of   deps is incorrect` end if;
			# Now we check its differential order and that its a linear system of PDE

            # Aug 22, 2022 if hybridflag then dord0 = difforder of ApSys
            
			# Aug 23, 2022
			if HybridSysFlag = true then 
			    dord0 := PDEtools[difforder](newsys[2]);
		    else
                dord0 := PDEtools[difforder](newsys);
			end if;

            userinfo(2, 'NumericDiffGeometryTools', `dsys order with no prolongation`, dord0);

			# All possible dependent variables and derivatives in the PDE system
			diffvars := DEtools[checkrank]([seq(deps[j](op(indeps)), j = 1 .. nops(deps))], deps, indep = indeps, degree = dord0);
			# userinfo(1, 'NumericDiffGeometryTools',` diffvars = `, diffvars );
			# Check that the system is linear in diffvars
			## if map(type, sys, linear(diffvars))<>[seq(true, i = 1 .. nops(sys))] then error `The system is not a linear PDE system` end if;

			# Aug 23, 2022
			if HybridSysFlag then
                if map(type, newsys[2], linear(diffvars))=[seq(true, i = 1 .. nops(newsys[2]))] then 
				    userinfo(1, 'NumericDiffGeometryTools',`The system is a linear PDE system`);
			    end if;
			    if map(type, newsys[2], polynom(anything, diffvars))=[seq(true, i = 1 .. nops(newsys[2]))] then 
				    userinfo(1, 'NumericDiffGeometryTools',`The system is a differential polynomial system of PDE`, map(type, newsys[2], polynom(anything,diffvars)), diffvars);
			    end if;
			else
			    if map(type, newsys, linear(diffvars))=[seq(true, i = 1 .. nops(newsys))] then 
				    userinfo(1, 'NumericDiffGeometryTools',`The system is a linear PDE system`);
			    end if;
			    if map(type, newsys, polynom(anything, diffvars))=[seq(true, i = 1 .. nops(newsys))] then 
				    userinfo(1, 'NumericDiffGeometryTools',`The system is a differential polynomial system of PDE`, map(type, newsys, polynom(anything,diffvars)), diffvars);
			    end if;
			end if;
						
	end if;

    # DimSol := [ ];  # Jun 17, 2023

	InvolSysIndex := [ ];  # changed from NULL to  [ ] Aug 13, 2017
	Involflag := false;

	
	# Sep 20, 2022
    if HybridSysFlag then
	    userinfo(0, 'NumericDiffGeometryTools', 'sys1', newsys[1]);
	    dordExSys := ComputeInvolutiveOrder(rifsimp(newsys[1])[Solved],vars)[ReducedInvolutiveDegree];
		# Jun 08, 2023
		# rifnewsys := rifsimp(newsys[1]);
		# userinfo(1, 'NumericDiffGeometryTools', 'initialdatasys1', rifnewsys);
		DHF := DifferentialHilbertFunction(initialdata(rifsimp(newsys[1])),s);
		userinfo(0, 'NumericDiffGeometryTools', `DiffHF`, DHF);
	end if;

	# Aug 29, 2022
	prolongtimes := Dmax;
	#Involdeg := Dmax;      # probably at least a misleading assignment of a variable
	for j from 0 to Dmax do
	    userinfo(0, 'NumericDiffGeometryTools', `start prolongation`);
		if polysysflag then dsys[j] := Prolong(PolynomialToPDE(newsys, indeps, deps), indeps, j);
			else 
			if HybridSysFlag then
			# if Dmax = 0 then
			#     dsys[j] := Prolong(newsys[2], indeps, j + dordExSys - dord0);     # Jun 14, 2023
			# else
			    dsys[j] := Prolong(newsys[2], indeps, j + dordExSys - dord0);  # Jun 08, 2023
			# end if;

                # dsys[j] := Prolong(sys[2], indeps, j + dordExSys - dord0 + 1);  # Sep 20, 2022, j+dordExSys-dord0+1
				# Aug 30,2022, j+3-dord0+1
			else
			    dsys[j] := Prolong(newsys, indeps, j);
			end if;
		end if;
		
		dord[j] := PDEtools[difforder](dsys[j]);
		userinfo(2, 'NumericDiffGeometryTools',`dsys order`, dord[j]);

		# Aug 30, 2022
		# if dord[j]<>dord0+j then error(`Calculated order <> dord0 + j at prolongation j =`, j); end if;
		userinfo(2, 'NumericDiffGeometryTools',`Just before MatrixFormPDE`);

		# Aug 22, 2022 add dsubs, dsysModExSys, dordModExSys[j] for substituted one
        # Aug 22, 2022 MatrixFormPDE for both dsys and dsysModExSys for hybrid cases, compare dord[j] and dordModExSys[j]

        # Aug 23, 2022 
		if HybridSysFlag then
		    dApSysModExSys[j] := dsubs(DEtools:-rifsimp(newsys[1])[Solved],dsys[j]);
			userinfo(1, 'NumericDiffGeometryTools',`dapsys after substitution=`, dApSysModExSys[j]);
            dordModExSys[j] := PDEtools[difforder](dApSysModExSys[j]);
            userinfo(1, 'NumericDiffGeometryTools',`dsys order after substitution=`, dordModExSys[j]);
            A, vrs, b := MatrixFormPDE(dApSysModExSys[j], [indeps, deps], ranking = PommaretRankingMtx(indeps, deps) );
			userinfo(3, 'NumericDiffGeometryTools',`Just after MatrixFormPDE`);
        else
            A, vrs, b := MatrixFormPDE(dsys[j], [indeps, deps], ranking = PommaretRankingMtx(indeps, deps) );
		    userinfo(3, 'NumericDiffGeometryTools',`Just after MatrixFormPDE`);
		end if;

        

		# A, vrs, b := MatrixFormPDE(dsys[j], [indeps, deps], ranking = PommaretRankingMtx(indeps, deps) );
		# userinfo(3, 'NumericDiffGeometryTools',`Just after MatrixFormPDE`);

        

		userinfo(4, 'NumericDiffGeometryTools',`A = `, A, `vrs =`, vrs, `b = `, b, ` indets A =`, indets(A), convert(deps, set));
		# Error checking that the matrix A or vector b don't contain any dep variables
		if (indets(A) union indets(b)) intersect convert(deps, set) <> {} then 
			error `Some error in input, dependent variables occur in either the coeff mtx of coeff vector of the system`; 
		end if;
		B := ReorderMtx(A, RevColOrder);
		if polysysflag then RandomPointList := [seq(1, k = 1 .. nops(indeps))];
			else 
				RandomPointList := RandomPointSpec;
		end if;
        userinfo(1,'NumericDiffGeometryTools', `MatrixInput`, B);
		# Aug 23, 2022
		if HybridSysFlag then
		    G[j] := GeoInvBasis(B, [indeps, deps, dordModExSys[j]], Tolerance = tol, RandomPoint = RandomPointList):
			# userinfo(1, 'NumericDiffGeometryTools', G[j][DimKerMtx], G[j][BetaMaxMtx], G[j][RankSymMtx], G[j][BrspMtx]);
            # userinfo(1, 'NumericDiffGeometryTools', G[j][DimSymMtx]);
			DimEx := 0;
			# Jun 11, 2023
			for k from 0 to dordModExSys[j] do
            # Jun 08, 2023
			    seriesDHF := series(DHF[3],s, k+1);       # d+1 terms of series
		        # userinfo(1, 'NumericDiffGeometryTools', `Series of DiffHF`, seriesDHF);
		        # DimSymEx := coeff(seriesDHF, s, k);       # dimension of symbols of Exact system at order d
	            # RankSymEx := k + 1 - DimSymEx;         # rank = # jet vars - dim
		        # userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ex`, DimSymEx, `Rank of Symbol of Ex`, RankSymEx);
                # RankSymAp := G[j][RankSymMtx][1,k+1];    
		        # DimSymAp := k + 1 - RankSymAp;     
		        # userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ap`, DimSymAp, `Rank of Symbol of Ap`, RankSymAp);     
                # RankSymSys := RankSymEx + RankSymAp;         # rank of joint system = rank of ex + rank of ap
		        # DimSymSys := k + 1 - RankSymSys;
		        # userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of joint system`, DimSymSys);
				# userinfo(1, 'NumericDiffGeometryTools', `GjDimKer`, G[j][DimKerMtx]);
				userinfo(1, 'NumericDiffGeometryTools', `J`, j, `K`, k);
				
                # Jun 17, 2023
				# rankApSys := N(nops(indeps),nops(deps),k) - G[j][DimKerMtx][1,k+1];
				userinfo(1, 'NumericDiffGeometryTools', `DimSymbol`, G[j][DimSymMtx]);
				# userinfo(1, 'NumericDiffGeometryTools', `rank of ApSys`, rankApSys);

				# seriesDHF := series(DHF[3],s, k1+1);       # d+1 terms of series
		        # userinfo(1, 'NumericDiffGeometryTools', `Series of DiffHF`, seriesDHF);
		        # # DimSymEx := coeff(seriesDHF, s, k);       # dimension of symbols of Exact system at order d

				# rankExSys := add(Nd(nops(indeps),nops(deps),q)-coeff(seriesDHF,s, q), q=0..k);

                # Jun 24, 2023
				DimEx := DimEx + coeff(seriesDHF,s,k);
                userinfo(1, 'NumericDiffGeometryTools', `DimExSys`, DimEx);
				
				rankExSys := N(nops(indeps),nops(deps),k) - DimEx;

				userinfo(1, 'NumericDiffGeometryTools', `rank of ExSys`, rankExSys);
			    DimSol[j,k] := G[j][DimKerMtx][1,k+1] - rankExSys;
				userinfo(1, 'NumericDiffGeometryTools', `Dimension of solution`, DimSol[j, k]);

                # Jun 22, 2023
				if j > 0 then
				userinfo(1, 'NumericDiffGeometryTools', `Determine dimsys=0`);
			    if DimSol[j-1, k] = DimSol[j, k] and k >= dordExSys and k >= dord0 then
			        Involflag := true;
					# if Dmax = 0 then 
					# userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j, `and degree = `, k,` >= dord[0]? `, dord[0],  `k = mtx col =`, k);
					# else
					userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j-1, `and degree = `, k,` >= dord[0]? `, dord[0],  `k = mtx col =`, k);
					# end if;
					InvolSysIndex := [op(InvolSysIndex), [j-1,k]];
			    end if;
				end if;
			
			end do;
                if Involflag = true then
		            prolongtimes := j-1;
			        # userinfo(1, 'NumericDiffGeometryTools', `Dimensionofapsys`, G[j][DimKerMtx]);
			    break;
		        end if;

        # seriesDHF := series(DHF[3],s, dordModExSys[j]+1);       # d+1 terms of series
		# userinfo(1, 'NumericDiffGeometryTools', `Series of DiffHF`, seriesDHF);
		# DimSymEx := coeff(seriesDHF, s^dordModExSys[j]);       # dimension of symbols of Exact system at order d
	    # RankSymEx := dordModExSys[j] + 1 - DimSymEx;         # rank = # jet vars - dim
		# userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ex`, DimSymEx, `Rank of Symbol of Ex`, RankSymEx);
        # RankSymAp := G[j][RankSymMtx][1,dordModExSys[j]+1];    
		# DimSymAp := dordModExSys[j] + 1 - RankSymAp;     
		# userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ap`, DimSymAp, `Rank of Symbol of Ap`, RankSymAp);     
        # RankSymSys := RankSymEx + RankSymAp;         # rank of joint system = rank of ex + rank of ap
		# DimSymSys := dordModExSys[j] + 1 - RankSymSys;
		# userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of joint system`, DimSymSys);

		# 	if j > 0 then 
		# # Aug 27, 2022
        #     for k from 1 to dordModExSys[j-1]+1 do
		# 	    userinfo(1,'NumericDiffGeometryTools',`k=`, k);
		#         if G[j-1][DimKerMtx][1,k] = G[j][DimKerMtx][1,k] and G[j-1][BetaMaxMtx][1,k][2]=G[j][RankSymMtx][1,k+1] and k >= dordModExSys[0]+1 # and G[j][RankSymMtx][1,k+1]<>0
		# 	        then userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j-1, `and degree = `, k - 1,` >= dord[0]? `, dordModExSys[0],  `k = mtx col =`, k);
		# 		        InvolSysIndex := [op(InvolSysIndex), [j-1,k,evalb(k >= dordModExSys[0]+1)]];  # Changed to list of lists not expr sequence Aug 13, 2017
		# 		        Involflag := true;
		#         end if;
		#     end do;
		#     end if;
		else
		    G[j] := GeoInvBasis(B, [indeps, deps, dord[j]], Tolerance = tol, RandomPoint = RandomPointList):

			# Aug 30, 2022
			userinfo(5, 'NumericDiffGeometryTools', G[j][BrspMtx], G[j][BkerMtx]);
			if j > 0 then 
            for k from 1 to dord[j-1]+1 do
		        if G[j-1][DimKerMtx][1,k] = G[j][DimKerMtx][1,k] and G[j-1][BetaMaxMtx][1,k][2]=G[j][RankSymMtx][1,k+1] and k >= dord[0]+1 # and G[j][RankSymMtx][1,k+1]<>0
			        then userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j-1, `and degree = `, k - 1,` >= dord[0]? `, dord[0],  `k = mtx col =`, k);
				        InvolSysIndex := [op(InvolSysIndex), [j-1,k,evalb(k >= dord[0]+1)]];  # Changed to list of lists not expr sequence Aug 13, 2017
				        Involflag := true;
		        end if;
		    end do;
		    end if;
			if Involflag = true then
		    prolongtimes := j-1; # Aug 29, 2022
			# Involdeg := j-1;
			## userinfo(1, 'NumericDiffGeometryTools',`G[j-1][DimSymMtx][1,k] = 0`, G[j-1][DimSymMtx][1,k] = 0, `G[j-1][DimKerMtx][1,k] =`, G[j-1][DimKerMtx][1,k]);
			## if G[j-1][DimSymMtx][1,k] = 0 then DimAtPoint := G[j-1][DimKerMtx][1,k];  # Return the dimension Aug 11, 2017
			##	else DimAtPoint := infinity;
			## end if;
			break;
		    end if;
        end if;
		
		# if j > 0 then 
        #     for k from 1 to dord[j-1]+1 do
		#         if G[j-1][DimKerMtx][1,k] = G[j][DimKerMtx][1,k] and G[j-1][BetaMaxMtx][1,k][2]=G[j][RankSymMtx][1,k+1] and k >= dord[0]+1 # and G[j][RankSymMtx][1,k+1]<>0
		# 	        then userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j-1, `and degree = `, k - 1,` >= dord[0]? `, dord[0],  `k = mtx col =`, k);
		# 		        InvolSysIndex := [op(InvolSysIndex), [j-1,k,evalb(k >= dord[0]+1)]];  # Changed to list of lists not expr sequence Aug 13, 2017
		# 		        Involflag := true;
		#         end if;
		#     end do;
		# end if;

		# if Involflag = true then
		#     prolongtimes := j-1; # Aug 29, 2022
		# 	# Involdeg := j-1;
		# 	## userinfo(1, 'NumericDiffGeometryTools',`G[j-1][DimSymMtx][1,k] = 0`, G[j-1][DimSymMtx][1,k] = 0, `G[j-1][DimKerMtx][1,k] =`, G[j-1][DimKerMtx][1,k]);
		# 	## if G[j-1][DimSymMtx][1,k] = 0 then DimAtPoint := G[j-1][DimKerMtx][1,k];  # Return the dimension Aug 11, 2017
		# 	##	else DimAtPoint := infinity;
		# 	## end if;
		# 	break;
		# end if;
	end do:
    
	# Aug 29, 2022
	# for j from 0 to Involdeg+1 do
	for j from 0 to prolongtimes+1 do
	    if HybridSysFlag then
		    userinfo(4, 'NumericDiffGeometryTools',` j = `, j, `degree = `, dordModExSys[j] );
		else
	        userinfo(4, 'NumericDiffGeometryTools',` j = `, j, `degree = `, dord[j] );
		end if;
	    userinfo(4, 'NumericDiffGeometryTools',` G[j][DimKerMtx], G[j][BetaMtx], G[j][BetaMaxMtx], G[j][RankSymMtx],  G[j][DimSymMtx], G[j][DimRspMtx]  `);
	    userinfo(4, 'NumericDiffGeometryTools',  G[j][DimKerMtx], G[j][BetaMtx], G[j][BetaMaxMtx], G[j][RankSymMtx],  G[j][DimSymMtx], G[j][DimRspMtx] );
		# Aug 29, 2022
	    userinfo(4, 'NumericDiffGeometryTools',` G[j][DimKerMtx], G[j][BetaMtx], G[j][BetaMaxMtx], G[j][RankSymMtx], G[j][DimRspMtx],  G[j][DimSymMtx], G[j][BrspMtx]`); # G[j][BrspMtx], G[j][BkerMtx];
	    userinfo(4, 'NumericDiffGeometryTools',  G[j][DimKerMtx], G[j][BetaMtx], G[j][BetaMaxMtx], G[j][RankSymMtx], G[j][DimRspMtx],  G[j][DimSymMtx], G[j][BrspMtx]);  # G[j][BrspMtx], G[j][BkerMtx] );
	end do;

	#InvolSysIndex := NULL;
	#for j from 0 to Dmax-1 dof
		#for k from dord[0]+1 to dord[j]+1 do
	#	for k from 1 to dord[j]+1 do
	#	if G[j][DimKerMtx][1,k] = G[j+1][DimKerMtx][1,k] and G[j][BetaMaxMtx][1,k][2]=G[j+1][RankSymMtx][1,k+1] and G[j+1][RankSymMtx][1,k+1]<>0
	#		then userinfo(1, 'NumericDiffGeometryTools',`System involutive at prolongation j = `, j, `and degree = `, k - 1,` >= dord[0]? `, dord[0],  `k = mtx col =`, k);
	#	userinfo(1, 'NumericDiffGeometryTools',`G[j][DimKerMtx][1,k] = G[j+1][DimKerMtx][1,k],  G[j][BetaMaxMtx][1,k][2]=G[j+1][RankSymMtx][1,k+1]`);
	#	print( G[j][DimKerMtx][1,k]  = G[j+1][DimKerMtx][1,k],  G[j][BetaMaxMtx][1,k][2]=G[j+1][RankSymMtx][1,k+1]);
		# userinfo(1, 'NumericDiffGeometryTools',`G[j][BrspMtx][1,k] =`, G[j][BrspMtx][1,k], `G[j][BkerMtx][1,k] =`, G[j][BkerMtx][1,k]);
	#	InvolSysIndex := InvolSysIndex, [j,k,evalb(k >= dord[0])];		
		# return [G[j][BrspMtx][1,k], G[j][BkerMtx][1,k]];
	#	end if;
	#	end do;
	#end do;
	
	### if evalb(map(whattype, InvolSysIndex) = [integer, integer, symbol])
	###    then jj := InvolSysIndex[1][1]+1;  kk := InvolSysIndex[2];
	###		if InvolSysIndex[3] = true and G[jj-1][DimSymMtx][1,kk] = 0 then DimAtPoint := G[jj-1][DimKerMtx][1,kk];
	###			else DimAtPoint := infinity;
	###		end if;
	### end if;
	if  Involflag <> true then 
				userinfo(0, 'NumericDiffGeometryTools',`System did not test involutive perhaps increase the max # of prolongations`);
				DimAtPoint := -1;  # ie. FAIL = -1
	    else 
			if HybridSysFlag = false then
				j1 := InvolSysIndex[1][1]+1;  k1 := InvolSysIndex[1][2];
				userinfo(0, 'NumericDiffGeometryTools',`InvolSysIndex = `, InvolSysIndex, `j1 =`, j1, `k1 = `, k1);
				if G[j1-1][DimSymMtx][1,k1] = 0 then DimAtPoint := G[j1-1][DimKerMtx][1,k1];
				else DimAtPoint := infinity;
				end if;
			# Jun 11, 2023
			else
			    j1 := InvolSysIndex[1][1]; k1 := InvolSysIndex[1][2];

				# Jun 22, 2023
				seriesDHF := series(DHF[3],s, k1+1);       # d+1 terms of series
		        userinfo(1, 'NumericDiffGeometryTools', `Series of DiffHF`, seriesDHF);
		        DimSymEx := coeff(seriesDHF, s, k1);       # dimension of symbols of Exact system at order d
	            RankSymEx := k1 + 1 - DimSymEx;         # rank = # jet vars - dim
		        userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ex`, DimSymEx, `Rank of Symbol of Ex`, RankSymEx);
                RankSymAp := G[j1][RankSymMtx][1,k1+1];    
		        DimSymAp := k1 + 1 - RankSymAp;     
		        userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of Ap`, DimSymAp, `Rank of Symbol of Ap`, RankSymAp);     
                RankSymSys := RankSymEx + RankSymAp;         # rank of joint system = rank of ex + rank of ap
		        DimSymSys := k1 + 1 - RankSymSys;
		        userinfo(1, 'NumericDiffGeometryTools', `Dimension of Symbol of joint system`, DimSymSys);
				userinfo(1, 'NumericDiffGeometryTools', `GjDimKer`, G[j][DimKerMtx]);
				userinfo(1, 'NumericDiffGeometryTools', `K`, k);

				if DimSymSys = 0 then 
				# Jun 18, 2023
				    DimAtPoint := DimSol[j1,k1];
				else DimAtPoint := infinity;
				end if;
				# # rankApSys := sum(Nd(nops(indeps),nops(deps),q) - 'G[j1][DimSymMtx][1,q+1]', q = 0..k1);
				# rankApSys := N(nops(indeps),nops(deps),k1) - G[j1][DimKerMtx][1,k1+1];
				# userinfo(1, 'NumericDiffGeometryTools', `DimSymbol`, G[j1][DimSymMtx]);
				# userinfo(1, 'NumericDiffGeometryTools', `rank of ApSys`, rankApSys);

				# seriesDHF := series(DHF[3],s, k1+1);       # d+1 terms of series
		        # userinfo(1, 'NumericDiffGeometryTools', `Series of DiffHF`, seriesDHF);
		        # # DimSymEx := coeff(seriesDHF, s, k);       # dimension of symbols of Exact system at order d

				# rankExSys := sum(Nd(nops(indeps),nops(deps),q)-coeff(seriesDHF,s, q), q=0..k1);
				# userinfo(1, 'NumericDiffGeometryTools', `rank of ExSys`, rankExSys);
			    # DimAtPoint := N(nops(indeps),nops(deps),k1) - rankApSys - rankExSys;
			end if;
	end if;
	
	## return table([ [InvolSysIndex], seq(G[j], j = 0 .. Dmax) ]);
    userinfo(0, 'NumericDiffGeometryTools',`InvolSysIndex = `, InvolSysIndex, `j1 =`, j1, `k1 = `, k1);
	#Aug 16, 2022
	
	# Jun 08, 2023
	return [RandomPointSpec, DimAtPoint, InvolSysIndex]; #, G[j1-1][DimKerMtx], G[j1-1][DimRspMtx], G[j1-1][DimSymMtx], G[Dmax][BrspMtx], G[Dmax][BkerMtx]]; 
	# Aug 16, 2022, add G[Dmax][DimRspMtx] to the output to show the rank of the approximate system G[Dmax][DimRspMtx]
end proc:

# Sep 20, 2020

GeometricInvolutiveFormGrid := proc(sys::list, vars::list, tol, Dmax, PointSpec)
	local x, t, xi, tau, phi, a, b, c, d, nx, nt, dx, dt, TOL, i, j, k, M, GIF, MList;
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	if nops(vars[1])<>2 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	if nops(vars[2])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# lprint(`PointSpec=`, PointSpec);

	if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	x   := vars[1][1];
	t   := vars[1][2];
	xi  := vars[2][1];
	tau := vars[2][2];
	phi := vars[2][3];
	
	# a <= x <= b
	# c <= t <= d
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	d   := PointSpec[1][4];

	nx  := PointSpec[2][1];
	nt  := PointSpec[2][2];
	dx := (b - a)/nx:
	dt := (d - c)/nt:

	MList := [ ];

	TOL := tol;

	if not(type(TOL,list)) then TOL := [TOL]; end if;
	TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	for k to nops(TOL) do
		M  := Matrix(1 .. nx+1, 1 .. nt+1);
		for i from 1 to nx+1 do		
			for j from 1 to nt+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				GIF:= GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], TOL[k], Dmax, [a+(i-1/2)*dx, c+(j-1/2)*dt]):
				#end if;
				M[i,j] := GIF[2]:
				## lprint("M[",i,j,"]=", M[i,j]):
			end do:
		end do:
		MList := [op(MList), M ];
	end do:
	return table([PointSpecialization=PointSpec,ToleranceList = TOL, SymDimMtx = MList])
end proc:

GeometricInvolutiveFormGrid1 := proc(sys::list, vars::list, tol, Dmax, PointSpec)
	local x, t, u, xi, tau, phi, a, b, c, d, nx, nt, dx, dt, TOL, i, j, k, M, GIF, MList,
	Nlist;       # Aug 18, 2022
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	if nops(vars[1])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	if nops(vars[2])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# lprint(`PointSpec=`, PointSpec);

	if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	x   := vars[1][1];
	t   := vars[1][2];
	u 	:= vars[1][3];
	xi  := vars[2][1];
	tau := vars[2][2];
	phi := vars[2][3];
	
	# a <= x <= b
	# c <= t <= d
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	d   := PointSpec[1][4];

	nx  := PointSpec[2][1];
	nt  := PointSpec[2][2];
	dx := (b - a)/nx:
	dt := (d - c)/nt:

	MList := [ ];

	TOL := tol;

	if not(type(TOL,list)) then TOL := [TOL]; end if;
	TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	for k to nops(TOL) do
		M  := Matrix(1 .. nx+1, 1 .. nt+1);
		for i from 1 to nx+1 do		
			for j from 1 to nt+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				GIF:= GeometricInvolutiveForm(sys, [[x, t, u], [xi, tau, phi]], TOL[k], Dmax, [a+(i-1)*dx, c+(j-1)*dt, 1]):
				#end if;
				M[i,j] := GIF[2]:
				Nlist[i,j] := GIF[6]:
				# Aug 16, 2022 add Nlist[i,j] for rank matrix
				## lprint("M[",i,j,"]=", M[i,j]):
			end do:
		end do:
		MList := [op(MList), M ];
	end do:
	return table([PointSpecialization=PointSpec,ToleranceList = TOL, SymDimMtx = MList, RankMatrix = Nlist])
	# Aug 16, 2022 add RankMatrix = Nlist to the output
end proc:

GeometricInvolutiveFormGrid2 := proc(sys::list, vars::list, tol, Dmax, PointSpec)
	local x, y, z, u, xi, eta, zeta, phi, a, b, c, d, e, nx, ny, dx, dy, TOL, i, j, k, M, GIF, MList;
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	if nops(vars[1])<>4 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	if nops(vars[2])<>4 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# lprint(`PointSpec=`, PointSpec);

	if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	# if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	x   := vars[1][1];
	y   := vars[1][2];
	z   := vars[1][3];
	u 	:= vars[1][4];
	xi  := vars[2][1];
	eta := vars[2][2];
	zeta := vars[2][3];
	phi  :=vars[2][4];
	
	# a <= x <= b
	# c <= y <= d
	# z = e
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	d   := PointSpec[1][4];
	e   := PointSpec[1][5];

	nx  := PointSpec[2][1];
	ny  := PointSpec[2][2];
	dx := (b - a)/nx:
	dy := (d - c)/ny:

	MList := [ ];

	TOL := tol;

	if not(type(TOL,list)) then TOL := [TOL]; end if;
	TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	for k to nops(TOL) do
		M  := Matrix(1 .. nx+1, 1 .. ny+1);
		for i from 1 to nx+1 do		
			for j from 1 to ny+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				GIF:= GeometricInvolutiveForm(sys, [[x, y, z, u], [xi, eta, zeta, phi]], TOL[k], Dmax, [a+(i-1/2)*dx, c+(j-1/2)*dy, e, 1]):
				#end if;
				M[i,j] := GIF[2]:
				## lprint("M[",i,j,"]=", M[i,j]):
			end do:
		end do:
		MList := [op(MList), M ];
	end do:
	return table([PointSpecialization=PointSpec,ToleranceList = TOL, SymDimMtx = MList, ker = GIF[4],rsp = GIF[5]])
end proc:

# Aug 22, 2022, GIF for poisson line
GeometricInvolutiveFormGrid3 := proc(sys::list, vars::list, tol, Dmax, PointSpec)
	local x, y, z, u, xi, eta, zeta, phi, a, b, c, d, e, nx, ny, dx, dy, TOL, i, j, k, M, GIF, MList;
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	if nops(vars[1])<>4 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	if nops(vars[2])<>4 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# lprint(`PointSpec=`, PointSpec);

	if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	# if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	x   := vars[1][1];
	y   := vars[1][2];
	z   := vars[1][3];
	u 	:= vars[1][4];
	xi  := vars[2][1];
	eta := vars[2][2];
	zeta := vars[2][3];
	phi  :=vars[2][4];
	
	# a <= x <= b
	# y = x
	# z = c
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	# d   := PointSpec[1][4];
	# e   := PointSpec[1][5];

	nx  := PointSpec[2][1];
	ny  := nx;
	dx := (b - a)/nx:
	dy := dx:

	MList := [ ];

	TOL := tol;

	if not(type(TOL,list)) then TOL := [TOL]; end if;
	TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	for k to nops(TOL) do
		M  := Vector(1 .. nx+1);
		for i from 1 to nx+1 do		
			# for j from 1 to ny+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				GIF:= GeometricInvolutiveForm(sys, [[x, y, z, u], [xi, eta, zeta, phi]], TOL[k], Dmax, [a+(i-1/2)*dx, a+(i-1/2)*dx, c, 1]):
				#end if;
				M[i] := GIF[2]:
				## lprint("M[",i,j,"]=", M[i,j]):
			# end do:
		end do:
		MList := [op(MList), M ];
	end do:
	return table([PointSpecialization=PointSpec,ToleranceList = TOL, SymDimMtx = MList])
end proc:

# Aug 27, 2022
GeometricInvolutiveFormGrid4 := proc(sys::list, vars::list, tol, Dmax, PointSpec)
	local x, y, u, a, b, c, d, nx, ny, dx, dy, TOL, i, j, k, M, GIF, MList,
	Nlist;       # Aug 18, 2022
	userinfo(0, 'NumericDiffGeometryTools', "sys is a system of LHPDE (later extend to poly sys etc)");
	userinfo(0, 'NumericDiffGeometryTools', "vars:  [[indeps], [deps]] 2 indep vars and the corresponding infinitesimal names");
	userinfo(0, 'NumericDiffGeometryTools', "PointSpec is a rectange [[a, b, c, d],[nx,nt]]");
	userinfo(0, 'NumericDiffGeometryTools', "Dmax", "is the max number of prolongations to be used by GIF");

	# if nops(vars[1])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 2 independent variables` end if;

	# if nops(vars[2])<>3 then error `Currently GeometricInvolutiveFormGrid only allows 3 dependent variables` end if;

	# lprint(`PointSpec=`, PointSpec);

	if not(type(PointSpec, list)) then error `currently PointSpec must be a list` end if;

	if nops(PointSpec[1])<>4 or nops(PointSpec[2])<>2 then error `currently must have: nops(PointSpec[1])=4, nops(PointSpec[2])=2` end if;

	x   := vars[1][1];
	y   := vars[1][2];
	u 	:= vars[2][1];
	
	
	# a <= x <= b
	# c <= t <= d
	a   := PointSpec[1][1];
	b   := PointSpec[1][2];
	c   := PointSpec[1][3];
	d   := PointSpec[1][4];

	nx  := PointSpec[2][1];
	ny  := PointSpec[2][2];
	dx := (b - a)/nx:
	dy := (d - c)/ny:

	MList := [ ];

	TOL := tol;

	if not(type(TOL,list)) then TOL := [TOL]; end if;
	TOL := evalf(TOL);

	# if map(type, tol, float) <> [seq(true, i = 1 .. nops(tol))] 
	# 	then error(`You must enter a list of numerical tolerances with`, nops(vars[1]), `entries`);
	# end if;

	for k to nops(TOL) do
		M  := Matrix(1 .. nx+1, 1 .. ny+1);
		for i from 1 to nx+1 do		
			for j from 1 to ny+1 do
				## lprint("About to enter GIF with i = ", i, "j = ", j):
				# GIF :=     GeometricInvolutiveForm(sys, [[x, t], [xi, tau, phi]], tol, 3, [a+(i-1)*dx, c+(j-1)*dt]);
				#if GIF[2] = -1 then  # i.e. failed to get involutive form
				GIF:= GeometricInvolutiveForm(sys, [[x, y], [u]], TOL[k], Dmax, [a+(i-1)*dx, c+(j-1)*dy, 1]):
				#end if;
				M[i,j] := GIF[2]:
				Nlist[i,j] := GIF[6]:
				# Aug 16, 2022 add Nlist[i,j] for rank matrix
				## lprint("M[",i,j,"]=", M[i,j]):
			end do:
		end do:
		MList := [op(MList), M ];
	end do:
	return table([PointSpecialization=PointSpec,ToleranceList = TOL, SymDimMtx = MList, RankMatrix = Nlist])
	# Aug 16, 2022 add RankMatrix = Nlist to the output
end proc:
