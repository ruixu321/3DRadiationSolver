/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::radiation::radiativeIntensityRayBand

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    radiativeIntensityRayBand.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeIntensityRay_H
#define radiativeIntensityRay_H

#include "absorptionEmissionModelBand.H"
#include "blackBodyEmission.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// Forward declaration of classes
class fvDOMBand;

/*---------------------------------------------------------------------------*\
                    Class radiativeIntensityRayBand Declaration
\*---------------------------------------------------------------------------*/

class radiativeIntensityRayBand
{
public:

    static const word intensityPrefix;


private:

    // Private data

        //- Refence to the owner fvDOMBand object
        const fvDOMBand& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Absorption/emission model
        const absorptionEmissionModel& absorptionEmissionModel_;

        //- Black body
        const blackBodyEmission& blackBody_;

        //- Total radiative intensity / [W/m2]
        volScalarField I_;

// WSGGM (added by Ivan Sikic 26/03/2015)
        //- radiation source term component(weighted intensity for Ru()) [W/m3]
        volScalarField kI_;

        //- Total radiative heat flux on boundary
        volScalarField Qr_;

        //- Incident radiative heat flux on boundary
        volScalarField Qin_;

        //- Emitted radiative heat flux on boundary
        volScalarField Qem_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to radiative intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;

        //- Global ray id - incremented in constructor
        static label rayId;

        //- My ray Id
        label myRayId_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        radiativeIntensityRayBand(const radiativeIntensityRayBand&);

        //- Disallow default bitwise assignment
        void operator=(const radiativeIntensityRayBand&);


public:

    // Constructors

        //- Construct form components
        radiativeIntensityRayBand
        (
            const fvDOMBand& dom,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const absorptionEmissionModel& absEmmModel_,
            const blackBodyEmission& blackBody,
            const label rayId
        );


    //- Destructor
    ~radiativeIntensityRayBand();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add radiative intensities from all the bands
            void addIntensity();


        // Access

            //- Return intensity
            inline const volScalarField& I() const;

// WSGGM (added by Ivan Sikic 26/03/2015)
            //- return weighted intensity
            inline const volScalarField& kI() const;

            //- Return const access to the boundary heat flux
            inline const volScalarField& Qr() const;

            //- Return non-const access to the boundary heat flux
            inline volScalarField& Qr();

            //- Return non-const access to the boundary incident heat flux
            inline volScalarField& Qin();

            //- Return non-const access to the boundary emmited heat flux
            inline volScalarField& Qem();

            //- Return const access to the boundary incident heat flux
            inline const volScalarField& Qin() const;

            //- Return const access to the boundary emmited heat flux
            inline const volScalarField& Qem() const;

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the radiative intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiativeIntensityRayBandI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
