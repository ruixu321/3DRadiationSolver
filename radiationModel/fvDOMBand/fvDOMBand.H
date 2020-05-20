/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::fvDOMBand

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter. Class modified by
    Ivan Sikic to solve 3 or 4 RTEs per direction (WSGGM).

    Available banded absorption models:
        wsggmAbsorptionEmissionSmithBand
        wsggmAbsorptionEmissionJohanssonBand

    i.e. dictionary
    \verbatim
        fvDOMBandCoeffs
        {
            nPhi        4;          // azimuthal angles in PI/2 on X-Y.
                                    //(from Y to X)
            nTheta      0;          // polar angles in PI (from Z to X-Y plane)
            convergence 1e-3;       // convergence criteria for radiation
                                    //iteration
            maxIter     4;          // maximum number of iterations
            cacheDiv    true;       // cache the div of the RTE equation.
            //NOTE: Caching div is "only" accurate if the upwind scheme is used
            //in div(Ji,Ii_h)
        }

        solverFreq   1; // Number of flow iterations per radiation iteration
    \endverbatim

    The total number of solid angles is  4*nPhi*nTheta.

    In 1D the direction of the rays is X (nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

SourceFiles
    fvDOMBand.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelfvDOMBand_H
#define radiationModelfvDOMBand_H

#include "radiativeIntensityRayBand.H"
#include "radiationModel.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class fvDOMBand Declaration
\*---------------------------------------------------------------------------*/

class fvDOMBand
:
    public radiationModel
{
    // Private data


        //- Incident radiation  [W/m2]
        volScalarField G_;

// WSGGM
// added by Ivan Sikic 26/03/2015
        //- source term component (contribution to fvDOMBand::Ru()) [W/m3]
        volScalarField kG_;

        //- Total radiative heat flux [W/m2]
        volScalarField Qr_;

         //- Emmited radiative heat flux [W/m2]
        volScalarField Qem_;

        //- Incident radiative heat flux [W/m2]
        volScalarField Qin_;

        //- Incident radiative heat flux [W/m2]
        volScalarField Qn_;

        //- Total absorption coefficient [1/m]
        volScalarField a_;

// WSGGM
        // added by Ivan Sikic 10/10/2014 (weighting wsggm coefficient)
        volScalarField ggCoeff_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nLambda_;

        //- Wavelength total absorption coefficient [1/m]
        PtrList<volScalarField> aLambda_;

// WSGGM
        // added by Ivan Sikic 10/10/2014 (wavelength weighting wsggm coefficient)
        PtrList<volScalarField> ggCoeffLambda_;
/*
        // added by Ivan Sikic 24/03/2016
        PtrList<volScalarField> planckMeanAbsorptionCoeff_;
*/
        //- Black body
        blackBodyEmission blackBody_;

        //- List of pointers to radiative intensity rays
        PtrList<radiativeIntensityRayBand> IRay_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- List of cached fvMatrices for rays
        List<PtrList<fvScalarMatrix> >fvRayDiv_;

        //- Cache convection div matrix
        bool cacheDiv_;

        //- Maximum omega weight
        scalar omegaMax_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        fvDOMBand(const fvDOMBand&);

        //- Disallow default bitwise assignment
        void operator=(const fvDOMBand&);

        //- Update nlack body emission
        void updateBlackBodyEmission();


public:

    //- Runtime type information
    TypeName("fvDOMBand"); // modified runtime type name to avoid conflict with regular fvDOMBand class


    // Constructors

        //- Construct from components
        fvDOMBand(const volScalarField& T);

        //- Construct from components
        fvDOMBand(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~fvDOMBand();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayIdLambdaId
            (
                const word& name,
                label& rayId,
                label& lambdaId
            ) const;

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;


        // Access

            //- Ray intensity for rayI
            inline const radiativeIntensityRayBand& IRay(const label rayI) const;

            //- Ray intensity for rayI and lambda bandwidth
            inline const volScalarField& IRayLambda
            (
                const label rayI,
                const label lambdaI
            ) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;

            //- Const access to total absorption coefficient
            inline const volScalarField& a() const;

//WSGGM
            // added by Ivan Sikic 10/10/2014 ( const access to weighting coefficient)
            inline const volScalarField& ggCoeff() const;

            //- Const access to wavelength total absorption coefficient
            inline const volScalarField& aLambda(const label lambdaI) const;

// WSGGM
            // added by Ivan Sikic 10/10/2014 ( const access to wavelength weighting coefficient)
            inline const volScalarField& ggCoeffLambda(const label lambdaI) const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;
// WSGGM
// added by Ivan Sikic 26/03/2015
            //- Const access to incident radiation field
            inline const volScalarField& kG() const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& Qr() const;

            //- Const access to incident radiative heat flux field
            inline const volScalarField& Qin() const;

            //- Const access to emitted radiative heat flux field
            inline const volScalarField& Qem() const;

            //- Const access to black body
            inline const blackBodyEmission& blackBody() const;

            //- Const access to cached fvMatrix
            inline const fvScalarMatrix& fvRayDiv
            (
                const label lambdaI,
                const label rayId
            ) const;

            //- Caching div(Ji, Ilamda)
            inline bool cacheDiv() const;

            //- Return omegaMax
            inline scalar omegaMax() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDOMBandI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
