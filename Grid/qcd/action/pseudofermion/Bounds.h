#pragma once

NAMESPACE_BEGIN(Grid);

    template<class Field>
    void HighBoundCheck(LinearOperatorBase<Field> &HermOp, 
			Field &Phi,
			RealD hi)
    {
      // Eigenvalue bound check at high end
      PowerMethod<Field> power_method;
      auto lambda_max = power_method(HermOp,Phi);
      std::cout << GridLogMessage << "Pseudofermion action lamda_max "<<lambda_max<<"( bound "<<hi<<")"<<std::endl;
      assert( (lambda_max < hi) && " High Bounds Check on operator failed" );
    }

     template<class Field> void ChebyBoundsCheck(LinearOperatorBase<Field> &HermOp,
						 Field &GaussNoise,
						 RealD lo,RealD hi) 
    {
      int orderfilter = 1000;
      Chebyshev<Field> Cheb(lo,hi,orderfilter);

      GridBase *FermionGrid = GaussNoise.Grid();

      Field X(FermionGrid);
      Field Z(FermionGrid);

      X=GaussNoise;
      RealD Nx = norm2(X);
      Cheb(HermOp,X,Z);
      RealD Nz = norm2(Z);

      std::cout << "************************* "<<std::endl;
      std::cout << " noise                    = "<<Nx<<std::endl;
      std::cout << " Cheb x noise             = "<<Nz<<std::endl;
      std::cout << " Ratio                    = "<<Nz/Nx<<std::endl;
      std::cout << "************************* "<<std::endl;
      assert( ((Nz/Nx)<1.0) && " ChebyBoundsCheck ");
    }
      
    template<class Field> void InverseSqrtBoundsCheck(int MaxIter,double tol,
						       LinearOperatorBase<Field> &HermOp,
						       Field &GaussNoise,
						       MultiShiftFunction &PowerNegHalf) 
    {
      GridBase *FermionGrid = GaussNoise.Grid();

      Field X(FermionGrid);
      Field Y(FermionGrid);
      Field Z(FermionGrid);

      X=GaussNoise;
      RealD Nx = norm2(X);

      ConjugateGradientMultiShift<Field> msCG(MaxIter,PowerNegHalf);
      msCG(HermOp,X,Y);
      msCG(HermOp,Y,Z);

      RealD Nz = norm2(Z);

      HermOp.HermOp(Z,Y);
      RealD Ny = norm2(Y);

      X=X-Y;
      RealD Nd = norm2(X);
      std::cout << "************************* "<<std::endl;
      std::cout << " noise                         = "<<Nx<<std::endl;
      std::cout << " (MdagM^-1/2)^2  noise         = "<<Nz<<std::endl;
      std::cout << " MdagM (MdagM^-1/2)^2  noise   = "<<Ny<<std::endl;
      std::cout << " noise - MdagM (MdagM^-1/2)^2  noise   = "<<Nd<<std::endl;
      std::cout << "************************* "<<std::endl;
      assert( (std::sqrt(Nd/Nx)<tol) && " InverseSqrtBoundsCheck ");
    }

NAMESPACE_END(Grid);

