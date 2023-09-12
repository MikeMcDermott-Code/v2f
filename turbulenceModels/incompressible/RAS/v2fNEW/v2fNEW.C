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

#include "v2fNEW.H"
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

defineTypeNameAndDebug(v2fNEW, 0);
addToRunTimeSelectionTable(RASModel, v2fNEW, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> v2fNEW::Ts() const
{
    return
		max(k_/epsilonTilda_, Cnabla_*sqrt(nu()/epsilonTilda_));//(k_/epsilonTilda_)*sqrt(1.0+36.0*epsilonTilda_*nu()/pow(k_, 2));
}

tmp<volScalarField> v2fNEW::Ls() const
{
    return
        CL_*max(pow(k_, 1.5)/
       (epsilonTilda_), Ceta_*pow(pow(nu(),3)/(epsilonTilda_),0.25));
}
// Function for w'2
tmp<volScalarField> v2fNEW::fd() const
{
    return
        min(max(pow((3*v2_/(2*k_)), 0.5), scalar(0.25)), scalar(1.0));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2fNEW::v2fNEW
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

	// Turbulent coefficients
    Cmu_(dimensionedScalar::lookupOrAddToDict("Cmu", coeffDict_, 0.19)),
    Ceps2_(dimensionedScalar::lookupOrAddToDict("Ceps2", coeffDict_, 1.92)),
    sigmaEps_(dimensionedScalar::lookupOrAddToDict("sigmaEps", coeffDict_, 1.3)),
    sigmak_(dimensionedScalar::lookupOrAddToDict("sigmak", coeffDict_, 1.0)),
    C1_(dimensionedScalar::lookupOrAddToDict("C1", coeffDict_, 1.4)),
    C2_(dimensionedScalar::lookupOrAddToDict("C2", coeffDict_, 0.3)),
    CL_(dimensionedScalar::lookupOrAddToDict("CL", coeffDict_, 0.23)),
    Ceta_(dimensionedScalar::lookupOrAddToDict("Ceta", coeffDict_, 80)),
    Cnabla_(dimensionedScalar::lookupOrAddToDict("Cnabla", coeffDict_, 6.0)),

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
        nu()*epsilonTilda_
    ), 
        
    QR_ // Shear rate
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
		sqrt( mag(nu() *mag(QR_)))
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
		wallDist(mesh_).y()/nu()
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

    LS_
    (
        IOobject
        (
            "LS",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        CL_*max(pow(k_, 1.5)
       /(epsilonTilda_ + epsilon0_), Ceta_*pow(pow(nu(),3)/(epsilonTilda_ + epsilon0_),0.25))
    ), 
    
    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), VSMALL)),    
    fMin_(dimensionedScalar("fMin", f_.dimensions(), 0.0))       
{	        	
	volScalarField Cmu1 = scalar(1.0) - scalar(1.5)*v2_/k_;
	volScalarField Cmu2 = (2-fd())/(2+fd()) - scalar(0.5)*v2_/k_;		
	dimensionedSymmTensor Id = symmTensor(scalar(1),scalar(0),scalar(0),scalar(1),scalar(0),scalar(1));
	dimensionedSymmTensor T11 = symmTensor(scalar(1),scalar(0),scalar(0),scalar(0),scalar(0),scalar(0));		

	N_ = (2.0/3.0)*Id*k_ + (Cmu1*(Id/3 - normaln_)*k_ + Cmu2*(2*T11 + normaln_ - Id)*k_);	
	
	nut_ == Cmu_*v2_*Ts();	
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> v2fNEW::R() const
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


tmp<volSymmTensorField> v2fNEW::devReff() const
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
           N_- nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> v2fNEW::divDevReff() const
{
		dimensionedSymmTensor I
    ( 
        "Identity", 
        dimless, 
        symmTensor::I
    );
	
    return
    (
      - fvm::laplacian(nuEff(), U_)
      - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
      + fvc::div(N_)
    );
}

bool v2fNEW::read()
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
        return true;
    }
    else
    {
        return false;
    }
}


void v2fNEW::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
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
    
	// normal field
	normaln_ == symm(fvc::grad(wallDist(mesh_).y()) * fvc::grad(wallDist(mesh_).y()))
				/max(magSqr(fvc::grad(wallDist(mesh_).y())), SMALL);

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));
    volScalarField G("RASModel::G", nut_*S2);
    volScalarField v2fAlpha = 1.0/Ts()*((C1_ - 6.0)*v2_ - 2.0/3.0*k_*(C1_ - 1.0));
    volScalarField Ceps1 = 1.4*(1.0 + 0.045*min(sqrt(k_/v2_), scalar(10000)));

    // Dissipation rate equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonTilda_)
      + fvm::div(phi_, epsilonTilda_)
      + fvm::SuSp(-fvc::div(phi_), epsilonTilda_)
      - fvm::laplacian(DepsilonEff(), epsilonTilda_)
     ==
        Ceps1*G/Ts()        
      - fvm::Sp(Ceps2_/Ts(), epsilonTilda_)
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
        G - fvm::Sp(epsilonTilda_/k_, k_)
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
      - 1.0/sqr(Ls())/k_*(v2fAlpha - C2_*G)
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
        k_*f_//min(k_*f_, C2_*G - v2fAlpha)
      - fvm::Sp(6.0*epsilonTilda_/k_, v2_)
    );

    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, v2Min_);

	// Calculate max velocity and friction
	Info<< max(Utau_) << nl << endl;
	Info<< max(U_.component(0)) << nl << endl;

	// Re-calculate terms   
	epsilonPlus_ == nu()*epsilonTilda_; 
	QR_ == max(sqrt(2.0)*mag(symm(fvc::grad(U_))), Q_);
	Utau_ == sqrt( mag(nu() *mag(QR_)));
	yPlus_ == wallDist(mesh_).y()/nu();
		
	volScalarField Cmu1 = scalar(1.0) - scalar(1.5)*v2_/k_;
	volScalarField Cmu2 = (2-fd())/(2+fd()) - scalar(0.5)*v2_/k_;		
	dimensionedSymmTensor Id = symmTensor(scalar(1),scalar(0),scalar(0),scalar(1),scalar(0),scalar(1));
	dimensionedSymmTensor T11 = symmTensor(scalar(1),scalar(0),scalar(0),scalar(0),scalar(0),scalar(0));		

	N_ == (2.0/3.0)*Id*k_ + (Cmu1*(Id/3.0 - normaln_)*k_ + Cmu2*(2*T11 + normaln_ - Id)*k_);		

	//volTensorField OMEGA = skew(fvc::grad(U_));
	//volVectorField VORTICITY = 2*operator*(OMEGA)/mag(OMEGA);
	//VORTICITYMAG_ == sqrt(scalar(2))*operator*(skew(fvc::grad(U_)))/max(mag(skew(fvc::grad(U_))), Q_);
	//titj == symm(operator^(fvc::grad(wallDist(mesh_).y()), VORTICITYMAG_) * operator^(fvc::grad(wallDist(mesh_).y()), VORTICITYMAG_));
	
	nut_ == Cmu_*v2_*Ts();
    nut_ == min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();
        
	RStress_ == N_ - nut_*twoSymm(fvc::grad(U_)); 	
	LS_ ==	CL_*max(pow(k_, 1.5)/(epsilonTilda_), Ceta_*pow(pow(nu(),3)/(epsilonTilda_),0.25));   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
