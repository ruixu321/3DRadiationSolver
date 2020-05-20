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

#include "fvDOMBand.H"
#include "absorptionEmissionModelBand.H"
#include "scatterModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOMBand, 0);
        addToRadiationRunTimeSelectionTables(fvDOMBand);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOMBand::initialise()
{
    if (mesh_.nSolutionD() == 3)    //3D
    {
        nRay_ = 4*nPhi_*nTheta_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi/(2.0*nPhi_);
        scalar deltaTheta = pi/nTheta_;
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                // update 21/09/2016

                // linear angular distribution
                 scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
/*
                // logarithmic angular distribution
                scalar thetai = pow(pi, (2.0*n - 1.0)*deltaTheta/2.0/pi );
*/
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRayBand
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    else
    {
        if (mesh_.nSolutionD() == 2)    //2D (X & Y)
        {

            // if (mesh_.solutionD()[vector::Z] != -1)
            // {
            //     FatalErrorInFunction
            //         << "Currently 2D solution is limited to the x-y plane"
            //         << exit(FatalError);
            // }

            scalar thetai = piByTwo;
            scalar deltaTheta = pi;
            nRay_ = 4*nPhi_;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi/(2.0*nPhi_);
            label i = 0;
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRayBand
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
        else    //1D (X)
        {

            // if (mesh_.solutionD()[vector::X] != 1)
            // {
            //     FatalErrorInFunction
            //         << "Currently 1D solution is limited to the x-direction"
            //         << exit(FatalError);
            // }
            
            scalar thetai = piByTwo;
            scalar deltaTheta = pi;
            nRay_ = 2;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi;
            label i = 0;
            for (label m = 1; m <= 2; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRayBand
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE // NO_WRITE // modified 25/04/2016
                ),
                a_
            )
        );
    }

    Info<< "fvDOMBand : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

// added by Ivan Sikic 01/10/2014
// construct field of WSGGM weighting coefficients

    forAll(ggCoeffLambda_, lambdaI)
    {
	ggCoeffLambda_.set
	(
	    lambdaI,
	    new volScalarField
	    (
	        IOobject
	        (
	            "ggCoeffLambda_" + Foam::name(lambdaI) ,
	            mesh_.time().timeName(),
	            mesh_,
	            IOobject::NO_READ,
	            IOobject::AUTO_WRITE // NO_WRITE // modified 25/04/2016
	        ),
	        ggCoeff_
	    )
	);
    }

    if (cacheDiv_)
    {
        Info<< "Caching div fvMatrix..."<< endl;
        for (label lambdaI = 0; lambdaI < nLambda_; lambdaI++)
        {
            fvRayDiv_[lambdaI].setSize(nRay_);

            forAll(IRay_, rayId)
            {
                const surfaceScalarField Ji(IRay_[rayId].dAve() & mesh_.Sf());
                const volScalarField& iRayLambdaI =
                    IRay_[rayId].ILambda(lambdaI);

                fvRayDiv_[lambdaI].set
                (
                    rayId,
                    new fvScalarMatrix
                    (
                        fvm::div(Ji, iRayLambdaI, "div(Ji,Ii_h)")
                    )
                );
            }
        }
    }

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << nl;
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOMBand::fvDOMBand(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    kG_
    (
        IOobject
        (
            "kG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kG", dimMass/pow3(dimTime)/dimLength, 0.0)
    ),

    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            "Qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
    Qn_
    (
        IOobject
        (
            "Qn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qn", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE // modified 25/04/2016
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),

// added by Ivan Sikic 10/10/2014
// WSGGM
    ggCoeff_
    (
        IOobject
        (
            "ggCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("ggCoeff", dimless, 0.0)
    ),

    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),

// WSGGM
    ggCoeffLambda_(nLambda_), // added by Ivan Sikic 10/10/2014
    // planckMeanAbsorptionCoeff_(nLambda_), // added by Ivan Sikic 24/03/2016
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    fvRayDiv_(nLambda_),
    cacheDiv_(coeffs_.lookupOrDefault<bool>("cacheDiv", false)),
    omegaMax_(0)
{
    initialise();
}


Foam::radiation::fvDOMBand::fvDOMBand
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),

    kG_
    (
        IOobject
        (
            "kG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kG", dimMass/pow3(dimTime)/dimLength, 0.0)
    ),

    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            "Qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
    Qn_
    (
        IOobject
        (
            "Qn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qn", dimMass/pow3(dimTime), 0.0)
    ),

    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),

// added by Ivan Sikic 10/10/2014
// WSGGM
    ggCoeff_
    (
        IOobject
        (
            "ggCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("ggCoeff", dimless, 0.0)
    ),

    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),

// WSGGM
    ggCoeffLambda_(nLambda_), // added by Ivan Sikic 10/10/2014
    // planckMeanAbsorptionCoeff_(nLambda_), // added by Ivan Sikic 24/03/2016
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    fvRayDiv_(nLambda_),
    cacheDiv_(coeffs_.lookupOrDefault<bool>("cacheDiv", false)),
    omegaMax_(0)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOMBand::~fvDOMBand()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOMBand::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::fvDOMBand::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

// WSGGM
    absorptionEmission_->correctNew(ggCoeff_, ggCoeffLambda_); // Ivan Sikic 14/10/2014

    updateBlackBodyEmission();

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0.0;
    label radIter = 0;
    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;
        radIter++;
        maxResidual = 0.0;
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();

                maxResidual = max(maxBandResidual, maxResidual);

                if (maxBandResidual < convergence_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

    } while (maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOMBand::Rp() const
{

   Info << " TEST ICI " << endl;

    if(nLambda_ == 4) // WSGGM Smith
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 5) // WSGGM Cassol or Johansson
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 6) // Box model, CO2
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 7) // Box model, H2O (overlapping bands seperately)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    // additional bands for mixture box model 30/08/2016
    if(nLambda_ == 8)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 9)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7] + aLambda_[8]*ggCoeffLambda_[8])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 10)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7] + aLambda_[8]*ggCoeffLambda_[8] + aLambda_[9]*ggCoeffLambda_[9])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 11)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7] + aLambda_[8]*ggCoeffLambda_[8] + aLambda_[9]*ggCoeffLambda_[9] + aLambda_[10]*ggCoeffLambda_[10])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 12)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7] + aLambda_[8]*ggCoeffLambda_[8] + aLambda_[9]*ggCoeffLambda_[9] + aLambda_[10]*ggCoeffLambda_[10] + aLambda_[11]*ggCoeffLambda_[11])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

    if(nLambda_ == 13)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),

            4.0*(aLambda_[0]*ggCoeffLambda_[0]+aLambda_[1]*ggCoeffLambda_[1]+aLambda_[2]*ggCoeffLambda_[2]+aLambda_[3]*ggCoeffLambda_[3] + aLambda_[4]*ggCoeffLambda_[4] + aLambda_[5]*ggCoeffLambda_[5] + aLambda_[6]*ggCoeffLambda_[6] + aLambda_[7]*ggCoeffLambda_[7] + aLambda_[8]*ggCoeffLambda_[8] + aLambda_[9]*ggCoeffLambda_[9] + aLambda_[10]*ggCoeffLambda_[10] + aLambda_[11]*ggCoeffLambda_[11] + aLambda_[12]*ggCoeffLambda_[12])*physicoChemical::sigma
           // 4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
    }

   if(nLambda_==25)
   {
     return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Rp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
		    false
		),

		4*(aLambda_[0]*ggCoeffLambda_[0]+
			aLambda_[1]*ggCoeffLambda_[1]+
			aLambda_[2]*ggCoeffLambda_[2]+
			aLambda_[3]*ggCoeffLambda_[3]+
			aLambda_[4]*ggCoeffLambda_[4]+
			aLambda_[5]*ggCoeffLambda_[5]+
			aLambda_[6]*ggCoeffLambda_[6]+
			aLambda_[7]*ggCoeffLambda_[7]+
			aLambda_[8]*ggCoeffLambda_[8]+
			aLambda_[9]*ggCoeffLambda_[9]+
			aLambda_[10]*ggCoeffLambda_[10]+
			aLambda_[11]*ggCoeffLambda_[11]+
			aLambda_[12]*ggCoeffLambda_[12]+
			aLambda_[13]*ggCoeffLambda_[13]+
			aLambda_[14]*ggCoeffLambda_[14]+
			aLambda_[15]*ggCoeffLambda_[15]+
			aLambda_[16]*ggCoeffLambda_[16]+
			aLambda_[17]*ggCoeffLambda_[17]+
			aLambda_[18]*ggCoeffLambda_[18]+
			aLambda_[19]*ggCoeffLambda_[19]+
			aLambda_[20]*ggCoeffLambda_[20]+
			aLambda_[21]*ggCoeffLambda_[21]+
			aLambda_[22]*ggCoeffLambda_[22]+
			aLambda_[23]*ggCoeffLambda_[23]+
			aLambda_[24]*ggCoeffLambda_[24]
			)*physicoChemical::sigma
			)
			);


   }

}

// WSGGM
// modified by Ivan Sikic 26/03/2015
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::fvDOMBand::Ru() const
{

	const DimensionedField<scalar, volMesh>& kG =
		kG_();
	const DimensionedField<scalar, volMesh> E =
		absorptionEmission_->ECont()()();
	//  const DimensionedField<scalar, volMesh> a =
	//     a_.dimensionedInternalField();

	return kG - E;
	//  return  a*G - E;
}


void Foam::radiation::fvDOMBand::updateBlackBodyEmission()
{
	for (label j=0; j < nLambda_; j++)
	{
		blackBody_.correct(j, absorptionEmission_->bands(j));
	}
}

// WSGGM
// modified by Ivan Sikic 26/03/2015
void Foam::radiation::fvDOMBand::updateG()
{
	G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	kG_ = dimensionedScalar("zero",dimMass/pow3(dimTime)/dimLength, 0.0); //WSGGM
	Qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
	Qem_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
	Qin_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
	Qn_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

	forAll(IRay_, rayI)
	{
		IRay_[rayI].addIntensity();
		G_ += IRay_[rayI].I()*IRay_[rayI].omega();
		kG_ += IRay_[rayI].kI()*IRay_[rayI].omega(); // WSGGM
		Qr_.boundaryFieldRef() += IRay_[rayI].Qr().boundaryField();
		Qem_.boundaryFieldRef() += IRay_[rayI].Qem().boundaryField();
		Qin_.boundaryFieldRef() += IRay_[rayI].Qin().boundaryField();


		// radiation probes 

		const vector rayd = IRay_[rayI].dAve();


		const vector direction90 (0,0,-1);

		double gaugeAngle = cos(7.*M_PI/180.);


		if((rayd & direction90)>0){

			Qn_ += IRay_[rayI].I()*(rayd & direction90);
		}
	}
}


void Foam::radiation::fvDOMBand::setRayIdLambdaId
(
 const word& name,
 label& rayId,
 label& lambdaId
 ) const
{
	// assuming name is in the form: CHARS_rayId_lambdaId
	size_type i1 = name.find_first_of("_");
	size_type i2 = name.find_last_of("_");

	rayId = readLabel(IStringStream(name.substr(i1+1, i2-1))());
	lambdaId = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());
}


// ************************************************************************* //
