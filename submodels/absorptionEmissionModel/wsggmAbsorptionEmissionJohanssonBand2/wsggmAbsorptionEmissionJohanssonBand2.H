/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::wsggmAbsorptionEmissionJohanssonBand2

Description

    wsggmAbsorptionEmissionJohanssonBand2 radiation absorption and emission coefficients
    for continuous phase.

    All the bands should have the same number of species and have to be entered
    in the same order.

    The coefficients for the species in the lookup table have to be specified
    for use in 1/m x P [bar].

    The look Up table file should be in the constant directory.

    band dictionary:
    \verbatim

    band0
    {
        species
        {

            grey
            {

                Tcommon         298.;
                Tlow            298.;
                Thigh           2500.;
                invTemp         false;
                loTcoeffs (0 0 0 0 0 0 0 0 0 0 0);
                hiTcoeffs (0.055 0.012 0.358 0.0731 -0.0466 -0.165 -0.0554 0.093 0.0598 0.0028 -0.0256);
            }
        }
    }

    band1
    {
        species
        {

            grey
            {

                Tcommon         298.;
                Tlow            298.;
                Thigh           2500.;
                invTemp         false;
                loTcoeffs (0 0 0 0 0 0 0 0 0 0 0);
                hiTcoeffs (0.88	-0.021 0.392 -0.212 0.0191 -0.291 0.644 -0.209 0.0784 -0.197 0.0662);
            }
        }
    }

    band2
    {
        species
        {

            grey
            {

                Tcommon         298.;
                Tlow            298.;
                Thigh           2500.;
                invTemp         false;
                loTcoeffs (0 0 0 0 0 0 0 0 0 0 0);
                hiTcoeffs (10 -1.6 0.142 -0.0831 0.0148 0.348 -0.294 0.0662 -0.122 0.118 -0.0295);
            }
        }
    }

    band3
    {
        species
        {

            grey
            {

                Tcommon         298.;
                Tlow            298.;
                Thigh           2500.;
                invTemp         false;
                loTcoeffs (0 0 0 0 0 0 0 0 0 0 0);
                hiTcoeffs (135 -35 0.0798 -0.037 0.0023 0.0866 -0.106 0.0305 -0.0127 0.0169 -0.0051);
            }
        }
    }

    band4 // clear gas (absorption = 0)
    {
        species
        {

            grey
            {

                Tcommon         298.;
                Tlow            298.;
                Thigh           2500.;
                invTemp         false;
                loTcoeffs (0 0 0 0 0 0 0 0 0 0 0);
                hiTcoeffs (0 0 0.9718 -0.259 -0.0104 -0.0214 0.1886 -0.0193 0.0035 -0.0593 0.006);
            }
        }
    }
}
    \endverbatim


SourceFiles
    wsggmAbsorptionEmissionJohanssonBand2.C

\*---------------------------------------------------------------------------*/

#ifndef wsggmAbsorptionEmissionJohanssonBand2_H
#define wsggmAbsorptionEmissionJohanssonBand2_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModelBand.H"
#include "HashTable.H"
#include "absorptionCoeffsJohansson.H"
#include "fluidThermo.H"

#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggmAbsorptionEmissionJohanssonBand2 Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmissionJohanssonBand2
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of species considered for absorptivity
        static const int nSpecies_ = 5;

        //- Maximum number of bands
        static const int maxBands_ = 10;

        //-  Absorption coefficients
        FixedList<FixedList<absorptionCoeffsJohanssen, nSpecies_>, maxBands_> coeffs_;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Hash table with species names
        HashTable<label> speciesNames_;

        //- Indices of species in the look-up table
        FixedList<label, nSpecies_> specieIndex_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Bands
        label nBands_;

        //- Pointer list of species being solved involved in the absorption
        UPtrList<volScalarField> Yj_;

        //- Access to soot coefficient from case file
        const scalar Csoot_;

public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmissionJohanssonBand2");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmissionJohanssonBand2(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggmAbsorptionEmissionJohanssonBand2();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

        // weighting coefficient for wsggm (reimplemented in radiativeIntensityRay)
                tmp<volScalarField> ggCoeffCont(const label bandI = 0) const;

            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

        //- Add contribution of ILambda to the total radiative intensity in
        //  direction i
        tmp<volScalarField> addIntensity
        (
            const label i,
            const volScalarField& ILambda
        ) const;

        void correct
        (
            volScalarField& a_,
            PtrList<volScalarField>& aLambda

        ) const;

        void correctNew
        (
            volScalarField& ggCoeff_,
            PtrList<volScalarField>& ggCoeffLambda

        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
