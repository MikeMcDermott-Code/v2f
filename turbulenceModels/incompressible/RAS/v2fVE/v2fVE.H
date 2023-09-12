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

Class
    Foam::incompressible::RASModels::v2fVE

Description
    v2f Visco
    The default model coefficients correspond to the following:
    @verbatim
        v2fVECoeffs
        {
            Cmu         0.09;
            C1          1.45;
            C2          1.9;
            sigmaEps    1.3;
 	    sigmak	1.1;
        }
    @endverbatim

SourceFiles
    v2fVE.C

\*---------------------------------------------------------------------------*/

#ifndef v2fVE_H
#define v2fVE_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class v2fVE Declaration
\*---------------------------------------------------------------------------*/

class v2fVE
:
    public RASModel
{
    // Private data

	// Units
            dimensionedScalar Time_;
            dimensionedScalar Length_;

        // Newtonian coefficients
       	    dimensionedScalar Cmu_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmak_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar Cnabla_;
	
	// viscoelastic coefficients
            dimensionedScalar CN1_;
            dimensionedScalar CN2_;
            dimensionedScalar CN3_;
            dimensionedScalar n1_;

	// rheology
            dimensionedScalar L2_;
            dimensionedScalar lambda_;
            dimensionedScalar beta_;
            dimensionedScalar kappa_;
            dimensionedScalar UNewt_;

	// bounds
            dimensionedScalar Q_;

        // Fields
	    volSymmTensorField normaln_;
            volScalarField k_;
            volScalarField epsilonTilda_;
            volScalarField epsilonPlus_;
            volScalarField QR_;
            volScalarField Utau_;
            volScalarField yPlus_;
            volScalarField f_;	
            volScalarField v2_;	
	    volSymmTensorField N_;
            volScalarField Ceps1_;	
	    volSymmTensorField C_;
	    volSymmTensorField NLT_;
	    volSymmTensorField NLTPlus_;
            volScalarField trC_;
            volScalarField fp_;
	    volScalarField epsilonV_;
	    volScalarField EtauP_;

            volSymmTensorField tauP_;	
            volScalarField nut_;

	    volScalarField eVyy_;
            volSymmTensorField DiffC_;
            volScalarField DR_;
	    volSymmTensorField RStress_;


            dimensionedScalar v2Min_;
            dimensionedScalar fMin_;

    // Private member functions

	tmp<volScalarField> Ts() const;
	tmp<volScalarField> Ls() const;
	tmp<volScalarField> fd() const;
	tmp<volScalarField> nu0() const;
	tmp<volScalarField> nuP() const;
	tmp<volScalarField> fN() const;
	tmp<volScalarField> Da() const;
public:

    //- Runtime type information
    TypeName("v2fVE");

    // Constructors

        //- Construct from components
        v2fVE
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~v2fVE()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }


        //- Return the trace of C
        virtual tmp<volScalarField> trC() const
        {
            return trC_;
        }

        //- Return NLT
        virtual tmp<volSymmTensorField> NLT() const
        {
            return NLT_;
        }

        //- Return fp
        virtual tmp<volScalarField> fp() const
        {
            return fp_;
        }

        //- Return epsilonV
        virtual tmp<volScalarField> epsilonV() const
        {
            return epsilonV_;
        }

        //- Return EtauP
        virtual tmp<volScalarField> EtauP() const
        {
            return EtauP_;
        }

        //- Return the conformation tensor
        virtual tmp<volSymmTensorField> C() const
        {
            return C_;
        }

        //- Return the polymer stress tensor
        virtual tmp<volSymmTensorField> tauP() const
        {
            return tauP_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilonTilda_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
