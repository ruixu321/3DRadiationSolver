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

\*---------------------------------------------------------------------------*/

#include "wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOMBand.H"
//#include "wideBandAbsorptionEmission.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    TName_("T")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, dict),
    TName_(dict.lookupOrDefault<word>("T", "T"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        const scalarField& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        refValue() =
            4.0*physicoChemical::sigma.value()*pow4(Tp)*emissivity()/pi;
        refGrad() = 0.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
(
    const wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOMBand& dom(refCast<const fvDOMBand>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchI = patch().index();

    if (dom.nLambda() == 0)
    {
        FatalErrorIn
        (
            "Foam::radiation::"
            "wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::updateCoeffs"
        )   << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;

    const vectorField n(patch().Sf()/patch().magSf());

    radiativeIntensityRayBand& ray =
        const_cast<radiativeIntensityRayBand&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());

    ray.Qr().boundaryFieldRef()[patchI] += Iw*nAve;

    const scalarField Eb
    (
        dom.blackBody().bLambda(lambdaId).boundaryField()[patchI]
    );

    scalarField temissivity = emissivity();

// WSGGM
// added by Ivan Sikic 13/11/2014
    scalarField wsggmWeightingCoeffBC = dom.ggCoeffLambda(lambdaId).boundaryField()[patchI];

    scalarField& Qem = ray.Qem().boundaryFieldRef()[patchI];
    scalarField& Qin = ray.Qin().boundaryFieldRef()[patchI];

    // Use updated Ir while iterating over rays
    // avoids to used lagged Qin
    scalarField Ir = dom.IRay(0).Qin().boundaryField()[patchI];

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).Qin().boundaryField()[patchI];

    }

    forAll(Iw, faceI)
    {
        const vector& d = dom.IRay(rayId).d();

        if ((-n[faceI] & d) > 0.0)
        {
            // direction out of the wall
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] =
                (
                    Ir[faceI]*(1.0 - temissivity[faceI])
                 // + temissivity[faceI]*Eb[faceI]     
                  + temissivity[faceI]*Eb[faceI]*wsggmWeightingCoeffBC[faceI] // added WSGGM coeff
                )/pi;

            // Emmited heat flux from this ray direction
         //   Qem[faceI] = refValue()[faceI]*nAve[faceI]; // for grey radiation models only


            if(lambdaId == 0)
            {
                Qem[faceI]=0; // reinitialisation needed after every iteration (Ivan Sikic 07/04/2015)
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 1)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 2)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 3)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }

            else if(lambdaId == 4)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 5)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 6)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 7)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            // additional bands for mixture box model 30/08/2016
            else if(lambdaId == 8)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 9)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 10)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 11)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else if(lambdaId == 12)
            {
                Qem[faceI] += refValue()[faceI]*nAve[faceI];
            }
            else
            {
                Qem[faceI] = refValue()[faceI]*nAve[faceI];
            }

        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

            // Incident heat flux on this ray direction
           // Qin[faceI] = Iw[faceI]*nAve[faceI]; // for grey radiation models only

            if (lambdaId == 0)
            {
                Qin[faceI] = 0;
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 1)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 2)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 3)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }

            else if (lambdaId == 4)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 5)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 6)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 7)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }

            // additional bands for mixture box model 30/08/2016
            else if (lambdaId == 8)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 9)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 10)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 11)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }
            else if (lambdaId == 12)
            {
                Qin[faceI] += Iw[faceI]*nAve[faceI];
            }

            else
            {
                Qin[faceI] = Iw[faceI]*nAve[faceI];
            }
        }

    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBase::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wideBandDiffusiveRadiationWSGGMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
