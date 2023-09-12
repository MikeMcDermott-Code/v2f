/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "v2fVE.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(v2fVE, 0);
addToRunTimeSelectionTable(RASModel, v2fVE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> v2fVE::Ts() const
{
    return
		max(k_/epsilonTilda_, Cnabla_*sqrt(nu0()/epsilonTilda_));//k_/epsilonTilda_;
}

tmp<volScalarField> v2fVE::Ls() const
{
    return
        CL_*max(pow(k_, 1.5)/
       (epsilonTilda_ + epsilon0_), Ceta_*pow(pow(nu0(),3)/(epsilonTilda_ + epsilon0_),0.25));
}
// Dampening function for w'2
tmp<volScalarField> v2fVE::fd() const
{
    return
        min(max(pow((3*v2_/(2*k_)), 0.5), scalar(0.3)/Da()), scalar(1.0));
}
// Total viscosity
tmp<volScalarField> v2fVE::nu0() const
{
    return
        nu()/beta_;
}
// Polymer viscosity
tmp<volScalarField> v2fVE::nuP() const
{
    return
        (1-beta_)*nu0();
}
// Local eddy viscosity
tmp<volScalarField> v2fVE::fN() const
{
    return
        nut_/nu0();
}
// Implicit polymer damping
tmp<volScalarField> v2fVE::Da() const
{
    return
         (scalar(1)+CN3_*pow(fp_*sqrt(L2_), n1_));
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2fVE::v2fVE
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

	// Units
    Time_("Time", dimTime,1), 		
    Length_("Length", dimLength, 1),    
	// Netownoan constants
    Cmu_(dimensionedScalar::lookupOrAddToDict("Cmu", coeffDict_, 0.2)),
    Ceps2_(dimensionedScalar::lookupOrAddToDict("Ceps2", coeffDict_, 1.9)),
    sigmaEps_(dimensionedScalar::lookupOrAddToDict("sigmaEps", coeffDict_, 1.3)),
    sigmak_(dimensionedScalar::lookupOrAddToDict("sigmak", coeffDict_, 1.0)),
    C1_(dimensionedScalar::lookupOrAddToDict("C1", coeffDict_, 1.4)),
    C2_(dimensionedScalar::lookupOrAddToDict("C2", coeffDict_, 0.3)),
    CL_(dimensionedScalar::lookupOrAddToDict("CL", coeffDict_, 0.23)),
    Ceta_(dimensionedScalar::lookupOrAddToDict("Ceta", coeffDict_, 70)),
    Cnabla_(dimensionedScalar::lookupOrAddToDict("Cnabla", coeffDict_, 6.0)),    
	// Viscoelastic constants
    CN1_(dimensionedScalar::lookupOrAddToDict("CN1", coeffDict_, 0.3)),
    CN2_(dimensionedScalar::lookupOrAddToDict("CN2", coeffDict_, 0.2)),
    CN3_(dimensionedScalar::lookupOrAddToDict("CN3", coeffDict_, 1)),
    n1_(dimensionedScalar::lookupOrAddToDict("n1", coeffDict_, 1)),

	// Viscoelastic parameters
    L2_(dimensionedScalar::lookupOrAddToDict("L2", coeffDict_, 900)),    
    lambda_(dimensionedScalar::lookupOrAddToDict("lambda", coeffDict_, 0.063291139)),
    beta_(dimensionedScalar::lookupOrAddToDict("beta", coeffDict_, 0.9)),
    kappa_(dimensionedScalar::lookupOrAddToDict("kappa", coeffDict_, 0.0001)),    
    UNewt_(dimensionedScalar::lookupOrAddToDict("UNewt", coeffDict_, 19.2)),


    Q_ // shear-rate minimum
    (
			"Qmin",
            dimless/dimTime,
            VSMALL
    ), 

    normaln_
    (
        IOobject
        (
            "normaln",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm(fvc::grad(wallDist(mesh_).y()) * fvc::grad(wallDist(mesh_).y()))
        /max(magSqr(fvc::grad(wallDist(mesh_).y())), SMALL)
    ), 

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilonTilda_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
    epsilonPlus_
    (
        IOobject
        (
            "epsilonPlus",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu0()*epsilonTilda_
    ), 

    // Shear rate
    QR_
    (
        IOobject
        (
            "QR",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
		max(sqrt(2.0)*mag(symm(fvc::grad(U_))), Q_)
	),         

	Utau_ // friction velocity
	(
		IOobject
		(
			"U_tau",
			runTime_.timeName(),
			U_.db(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		sqrt(max(QR_)*nu0())
	),  

	yPlus_
	(
		IOobject
		(
			"yPlus",
			runTime_.timeName(),
			U_.db(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		wallDist(mesh_).y()/nu0() 
	), 

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ), 	   

    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ), 

    N_
    (
        IOobject
        (
            "N",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ), 

///////////////////////////////////////

    Ceps1_
    (
        IOobject
        (
            "Ceps1",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.4*(1.0 + 0.045*min(sqrt(k_/v2_), scalar(100)))
    ),   
    
    C_ // Conformation tensor
    (
        IOobject
        (
            "C",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_        
	),

    NLT_ // Non-linear-term
    (
        IOobject
        (
            "NLT",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),  

    NLTPlus_ 
    (
        IOobject
        (
            "NLTPlus",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu0()*NLT_
    ),  

    trC_  // limit on tr(C)
    (
        IOobject
        (
            "trC",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        min(max(tr(C_), 3. + SMALL), L2_ - SMALL) 
	),

    fp_	// Peterlin function
    (
        IOobject
        (
            "fp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
		min(max((L2_ - 3)/(L2_ - trC_), scalar(1.) + SMALL), scalar(100.))
    ), 

    epsilonV_  // Viscoelastic stress work	
    (
        IOobject
        (
            "epsilonV",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
		(1/Time_)*(nuP()/lambda_)*fp_*tr(NLT_)/2
	),

    EtauP_  // Viscoelastic destruction term
    (
        IOobject
        (
            "EtauP",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		Ceps1_*epsilonV_/Ts()
	),

    tauP_ // Polymer stress tensor
    (
        IOobject
        (
            "tauP",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ), 
		mesh_
	),	
    
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReNut("nut", mesh_)
    ),
    
    eVyy_  // Transverse component stress work
    (
        IOobject
        (
            "eVyy",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),		
		CN2_*epsilonV_*(v2_/k_)
	), 
	
    DiffC_ // Artificial diffusion
    (
        IOobject
        (
            "DiffC",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ), 
		mesh_
	),		   

	DR_ // Drag reduction percentage
    (
        IOobject
        (
            "DR",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ), 
		mesh_
	),

    RStress_
    (
        IOobject
        (
            "RStress",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        N_ - nut_*twoSymm(fvc::grad(U_))
    ),     
    
    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), SMALL)),    
    fMin_(dimensionedScalar("fMin", f_.dimensions(), 0.0))       
{		
	volScalarField Cmu1 = scalar(1.0) - scalar(1.5)*v2_/k_;
	volScalarField Cmu2 = (2-fd())/(2+fd()) - scalar(0.5)*v2_/k_;
	dimensionedSymmTensor Id = symmTensor(scalar(1),scalar(0),scalar(0),scalar(1),scalar(0),scalar(1));
	dimensionedSymmTensor T11 = symmTensor(scalar(1),scalar(0),scalar(0),scalar(0),scalar(0),scalar(0));		

	N_ = (2.0/3.0)*Id*k_ + Cmu1*(Id/3 - normaln_)*k_ + Cmu2*(2*T11 + normaln_ - Id)*k_;	

	tauP_ = (1/Time_)*(nuP()/lambda_)*(fp_*C_-I);
	
    nut_ = Cmu_*v2_*Ts();
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

	DiffC_ = (Length_*Length_/Time_)*fvc::laplacian(kappa_, C_);
	DR_ = 100*(scalar(1) - pow((Time_*Time_/Length_/Length_)*max(QR_)*nu0(), 0.2574)* pow(max(mag(U_))/(UNewt_*Length_/Time_), -1.743));
 
    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> v2fVE::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            N_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> v2fVE::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           N_- nuEff()*dev(twoSymm(fvc::grad(U_)))-tauP_
        )
    );
}

tmp<fvVectorMatrix> v2fVE::divDevReff() const
{
    return
    (
      - fvm::laplacian(nuEff(), U_)
      - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
      + fvc::div(N_)
      - fvc::div(tauP_)     
    );
}

bool v2fVE::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());        
        C1_.readIfPresent(coeffDict());    
        C2_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        Ceta_.readIfPresent(coeffDict());   
        Cnabla_.readIfPresent(coeffDict());   
        CN1_.readIfPresent(coeffDict());   
        CN2_.readIfPresent(coeffDict());   
        CN3_.readIfPresent(coeffDict());   
        n1_.readIfPresent(coeffDict());   
        L2_.readIfPresent(coeffDict());                
		lambda_.readIfPresent(coeffDict());                 
		beta_.readIfPresent(coeffDict());                
        kappa_.readIfPresent(coeffDict());                
        UNewt_.readIfPresent(coeffDict());                
        return true;
    }
    else
    {
        return false;
    }
}

void v2fVE::correct()
{
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilonTilda_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }	

	normaln_ == symm(fvc::grad(wallDist(mesh_).y()) * fvc::grad(wallDist(mesh_).y()))
				/max(magSqr(fvc::grad(wallDist(mesh_).y())), SMALL);

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));
    volScalarField G("RASModel::G", nut_*S2);
    volScalarField v2fAlpha = 1.0/Ts()*((C1_ - 6.0)*v2_ - 2.0/3.0*k_*(C1_ - 1.0));
	Ceps1_ = 1.4*(1.0 + 0.045*min(sqrt(k_/v2_), scalar(10000)));
	epsilonV_ = (1/Time_)*(nuP()/lambda_)*fp_*tr(NLT_)/2;
	EtauP_ =  Ceps1_*epsilonV_/Ts();
	eVyy_ = CN2_*epsilonV_*(v2_/k_);
			
    // Dissipation rate equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonTilda_)
      + fvm::div(phi_, epsilonTilda_)
      + fvm::SuSp(-fvc::div(phi_), epsilonTilda_)
      - fvm::laplacian(DepsilonEff(), epsilonTilda_)
     ==
        Ceps1_*G/Ts()        
      - fvm::Sp(Ceps2_/Ts(), epsilonTilda_)
      - EtauP_
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilonTilda_, epsilon0_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp((epsilonTilda_)/k_, k_) - epsilonV_
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);
	
	// Relxation function equation
    tmp<fvScalarMatrix> fEqn
    (     
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(Ls()), f_)
      - 1.0/sqr(Ls())/k_*(v2fAlpha - C2_*G/Da())
    );

    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);	
	
    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      + fvm::SuSp(-fvc::div(phi_), v2_)      
      - fvm::laplacian(DkEff(), v2_)
      ==
        min(k_*f_, C2_*G - v2fAlpha)
      - fvm::Sp(6.0*epsilonTilda_/k_, v2_)
      - eVyy_
    );

    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, v2Min_);
        
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U_);     
    // Shear rate
    QR_ = max(sqrt(2.0)*mag(symm(fvc::grad(U_))), Q_);    
	// Conformation velocity gradient correlation
    volTensorField A = C_ & L;  
    // Identity tensor
	dimensionedSymmTensor I("Identity", dimless, symmTensor::I);
	// Peterlin function bounded
    fp_ = min(max((L2_ - 3)/(L2_ - trC_), scalar(1.) + SMALL), scalar(100.));

	// Mean flow direction
	dimensionedSymmTensor T11 = symmTensor(scalar(1),scalar(0),scalar(0),scalar(0),scalar(0),scalar(0));	
	    
	// NLT
	NLT_ = CN1_*fN()*tr(A)*(T11+CN2_*(v2_/k_)*normaln_);
    NLT_.relax();
    
	// Conformation tensor transport
    tmp<fvSymmTensorMatrix> CEqn
    (
        fvm::ddt(C_)
	  + fvm::div(phi_, C_)
      - twoSymm(A) //Mij      
      - NLT_ //NLTij
      - (Length_*Length_/Time_)*fvm::laplacian(kappa_, C_) // artificial diffusion
	  ==
	  -(1/Time_)*(1.0/lambda_)*(fp_*C_-I)
    );
	     
    CEqn().relax();
	solve(CEqn);

///////////////////
	Info<< max(Utau_) << nl << endl;
	Info<< max(U_.component(0)) << nl << endl;
	Info<< max(DR_) << nl << endl;

	// Re-calculate terms 
	NLTPlus_ == nu0()*NLT_;
	trC_ == min(max(tr(C_), 3. + SMALL), L2_ - SMALL);   
    fp_ == min(max((L2_ - 3)/(L2_ - trC_), scalar(1.) + SMALL), scalar(100.));  
	epsilonPlus_ == nu0()*epsilonTilda_; 
	QR_ == max(sqrt(2.0)*mag(symm(fvc::grad(U_))), Q_);
	Utau_ == sqrt(max(QR_)*nu0());
	yPlus_ == wallDist(mesh_).y()/nu0();
	
	volScalarField Cmu1 = scalar(1.0) - scalar(1.5)*v2_/k_;
	volScalarField Cmu2 = (2-fd())/(2+fd()) - scalar(0.5)*v2_/k_;
	dimensionedSymmTensor Id = symmTensor(scalar(1),scalar(0),scalar(0),scalar(1),scalar(0),scalar(1));

	N_ == (2.0/3.0)*Id*k_ + Cmu1*(Id/3.0 - normaln_)*k_ + Cmu2*(2*T11 + normaln_ - Id)*k_;		

	tauP_ == (1/Time_)*(nuP()/lambda_)*(fp_*C_-I);

	nut_ == Cmu_*v2_*Ts();
    nut_ == min(nut_, nuRatio()*nu());

	DiffC_ == (Length_*Length_/Time_)*fvc::laplacian(kappa_, C_);    
	RStress_ == N_ - nut_*twoSymm(fvc::grad(U_));  
	
	
	DR_ == 100*(scalar(1) - pow((Time_*Time_/Length_/Length_)*max(QR_)*nu0(), 0.2574)* pow(max(mag(U_))/(UNewt_*Length_/Time_), -1.743));	  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
