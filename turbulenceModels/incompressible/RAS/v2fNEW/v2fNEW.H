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
    Foam::incompressible::RASModels::v2fNEW

Description
    Nagano and Hishida low-Reynolds k-epsilon turbulence model for
    incompressible flows. Added ft() turbulent diffusion term added on k
    and epsilon transport for stability.

    The default model coefficients correspond to the following:
    @verbatim
        v2fNEWCoeffs
        {
            Cmu         0.09;
            C1          1.45;
            C2          1.9;
            sigmaEps    1.3;
 	    sigmak	1.1;
        }
    @endverbatim

SourceFiles
    v2fNEW.C

\*---------------------------------------------------------------------------*/

#ifndef v2fNEW_H
#define v2fNEW_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class v2fNEW Declaration
\*---------------------------------------------------------------------------*/

class v2fNEW
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmak_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar Cnabla_;
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
            volScalarField nut_;
	    volSymmTensorField RStress_;
	    volScalarField LS_;

            dimensionedScalar v2Min_;
            dimensionedScalar fMin_;

    // Private member functions

	tmp<volScalarField> Ts() const;
	tmp<volScalarField> Ls() const;
	tmp<volScalarField> fd() const;

public:

    //- Runtime type information
    TypeName("v2fNEW");

    // Constructors

        //- Construct from components
        v2fNEW
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~v2fNEW()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
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

        //- Return the wall normal Reynolds stress 
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

        //- Return the elliptic function
        virtual tmp<volScalarField> f() const
        {
            return f_;
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
