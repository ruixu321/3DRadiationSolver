/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
// #define DEBUG(x) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
// #define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;

#include "ThermoSurfaceFilmMeredith.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Pstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::ThermoSurfaceFilmMeredith<CloudType>::interactionTypeNames_
(
    IStringStream
    (
        "(absorb bounce splashBai)"
    )()
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::ThermoSurfaceFilmMeredith<CloudType>::interactionType
Foam::ThermoSurfaceFilmMeredith<CloudType>::interactionTypeEnum(const word& it) const
{
    forAll(interactionTypeNames_, i)
    {
        if (interactionTypeNames_[i] == it)
        {
            return interactionType(i);
        }
    }

    FatalErrorIn
    (
        "ThermoSurfaceFilmMeredith<CloudType>::interactionType "
        "ThermoSurfaceFilmMeredith<CloudType>::interactionTypeEnum"
        "("
            "const word& it"
        ") const"
    )   << "Unknown interaction type " << it
        << ". Valid interaction types include: " << interactionTypeNames_
        << abort(FatalError);

    return interactionType(0);
}


template<class CloudType>
Foam::word Foam::ThermoSurfaceFilmMeredith<CloudType>::interactionTypeStr
(
    const interactionType& it
) const
{
    if (it >= interactionTypeNames_.size())
    {
        FatalErrorIn
        (
            "ThermoSurfaceFilmMeredith<CloudType>::interactionType "
            "ThermoSurfaceFilmMeredith<CloudType>::interactionTypeStr"
            "("
                "const interactionType& it"
            ") const"
        )   << "Unknown interaction type enumeration" << abort(FatalError);
    }

    return interactionTypeNames_[it];
}


template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilmMeredith<CloudType>::tangentVector
(
    const vector& v
) const
{
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        vector vTest = rndGen_.sample01<vector>();
        tangent = vTest - (vTest & v)*v;
        magTangent = mag(tangent);
    }

    return tangent/magTangent;
}


template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilmMeredith<CloudType>::splashDirection
(
    const vector& tanVec1,
    const vector& tanVec2,
    const vector& nf,
    scalar thetaSiMin_,
    scalar thetaSiMax_,
    scalar thetaSiMean_,
    scalar thetaSiVar_,
    const vector& dropTan,
    const scalar& beta
) const
{
    // azimuthal angle [rad]
    const scalar phiSi = twoPi*rndGen_.sample01<scalar>();

    // ejection angle [rad]
    scalar thetaSi=0.0;
    if(normalDistribution_){
        thetaSi = pi/180.0*(sampleNormal(thetaSiMin_,thetaSiMax_,thetaSiMean_,thetaSiVar_));
    }
    else{
        thetaSi = pi/180.0*(rndGen_.sample01<scalar>()*(thetaSiMax_ - thetaSiMin_) + thetaSiMin_);
    }
//    DEBUG(thetaSi*180.0/pi);

    // direction vector of new parcel
    const scalar alpha = sin(thetaSi);
    const scalar dcorr = cos(thetaSi);
    const vector normal = alpha*(tanVec1*cos(phiSi) + tanVec2*sin(phiSi));
    vector dirVec = dcorr*nf;
    dirVec += normal;
    /*account for tangential impingement*/
    vector newDirVec = beta*dropTan+(1.-beta)*dirVec;

    return newDirVec/mag(newDirVec);
}
template<class CloudType>
Foam::scalar Foam::ThermoSurfaceFilmMeredith<CloudType>::sampleNormal(
        scalar minValue, 
        scalar maxValue, 
        scalar expectation, 
        scalar variance) const
{

    scalar a = erf((minValue - expectation)/variance);
    scalar b = erf((maxValue - expectation)/variance);

    scalar y = rndGen_.sample01<scalar>();
    scalar x = erfInv(y*(b - a) + a)*variance + expectation;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = min(max(x, minValue), maxValue);

    return x;
}
template<class CloudType>
Foam::scalar Foam::ThermoSurfaceFilmMeredith<CloudType>::erfInv(const scalar y) const
{
    scalar a_=0.147;
    scalar k = 2.0/(constant::mathematical::pi*a_) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a_;
    scalar x = sqrt(-k + sqrt(k*k - h));
    if (y < 0.0)
    {
        x *= -1.0;
    }
    return x;
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::absorbInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mass,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " absorbInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Parcel normal velocity
    const vector Un = nf*(Urel & nf);

    // Parcel tangential velocity
    const vector Ut = Urel - Un;

    filmModel.addSources
    (
        pp.index(),
        faceI,
        mass,                           // mass
        mass*Ut,                        // tangential momentum
        mass*mag(Un),                   // impingement pressure
        mass*p.hs()                     // energy
    );

    this->nParcelsTransferred()++;

    keepParticle = false;
    if(debug){
        Info << "mass absorbed " << p.position()[0] << " " << p.position()[1]  << " " << p.position()[2] << " " << mass << endl;

        static const vector center(0,0,0);
        static scalarList binRadius(7,0.0);
        binRadius[0]=0.022;
        binRadius[1]=0.040;
        binRadius[2]=0.060;
        binRadius[3]=0.089;
        binRadius[4]=0.134;
        binRadius[5]=0.190;
        binRadius[6]=0.490;
        //const scalarList binRadius(.0,.018,.040,.060,.089,.134,.190);
        static scalarList cummulatedMass(binRadius.size()+1,0.0);//last bin for adhered mass from initial impingement

        scalar x=p.position()[0];
        scalar y=p.position()[1];
        scalar radius=sqrt(pow(x-center.x(),2)+pow(y-center.y(),2));
        if(p.typeId()==3){
            for(label j=0;j<binRadius.size();j++){
                if(j==0){
                    if(radius<binRadius[j]){
                        cummulatedMass[j]+=mass;
                    }
                }
                else{
                    if(binRadius[j-1]<=radius&&radius<binRadius[j]){
                        cummulatedMass[j]+=mass;
                    }
                }

            }
        }
        else{
            cummulatedMass[7]+=mass;
            //DEBUG(mass);
        }
        static char buffer[256];
        sprintf(buffer,"SplashedMassLag %10.5f ",this->owner().mesh().time().value());
        Info << buffer;
        sprintf(buffer," % 10.5g ",cummulatedMass[7]);
        Info << buffer;
        for(label j=0;j<cummulatedMass.size()-1;j++){
            sprintf(buffer," % 10.5g ",cummulatedMass[j]);
            Info << buffer;
        }
        scalar totalMass=sum(cummulatedMass);
        sprintf(buffer," % 10.5g ",sum(cummulatedMass));
        Info << buffer;
        Info << endl;

        sprintf(buffer,"SplashedMassNorm %10.5f ",this->owner().mesh().time().value());
        Info << buffer;
        sprintf(buffer," % 10.5g ",cummulatedMass[7]/totalMass*100.0);
        Info << buffer;
        for(label j=0;j<cummulatedMass.size()-1;j++){
            sprintf(buffer," % 10.5g ",cummulatedMass[j]/totalMass*100.0);
            Info << buffer;
        }
        sprintf(buffer," % 10.5g ",sum(cummulatedMass)/totalMass*100.0);
        Info << buffer;
        Info << endl;
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::bounceInteraction
(
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
) const
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " bounceInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Flip parcel normal velocity component
    p.U() -= 2.0*nf*(Urel & nf);

    keepParticle = true;
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::drySplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " drySplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
    const scalar Wec = Adry_*pow(La, -0.183);

    if (debug)
    {
        Info<< "Velocity " << p.U() << " Type " << p.typeId() << " We " << We << endl;
    }

    //DEBUG(We);
    //DEBUG(Wec);
    if (We < Wec) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else // splash
    {
        // ratio of incident mass to splashing mass
        const scalar mRatio = drySplashRatioMin_ + (drySplashRatioMax_-drySplashRatioMin_)*rndGen_.sample01<scalar>();
        //const scalar mRatio = drySplashRatioMin_ + (drySplashRatioMax_-drySplashRatioMin_)*0.5;
        bool wet=false;
        splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle,wet);
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::wetSplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " wetSplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    vector& U = p.U();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    if(debug)
    {
        Info<< "Velocity " << p.U() << " Type " << p.typeId() << " We " << We << endl;
    }

    // Critical Weber number
    const scalar Wec = Awet_*pow(La, -0.183);

    //DEBUG(We);
    //DEBUG(Wec);
    if (We < 1) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else if ((We >= 1) && (We < 20)) // bounce
    {
        // incident angle of impingement
        const scalar theta = pi/2 - acos(U/mag(U) & nf);

        // restitution coefficient
        const scalar epsilon = 0.993 - theta*(1.76 - theta*(1.56 - theta*0.49));

        // update parcel velocity
        U = -epsilon*(Un) + 5/7*(Ut);

        keepParticle = true;
        return;
    }
    else if ((We >= 20) && (We < Wec)) // spread - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else    // splash
    {
        // ratio of incident mass to splashing mass
        // splash mass can be > incident mass due to film entrainment
        const scalar mRatio = wetSplashRatioMin_ + (wetSplashRatioMax_-wetSplashRatioMin_)*rndGen_.sample01<scalar>();
        //const scalar mRatio = wetSplashRatioMin_ + (wetSplashRatioMax_-wetSplashRatioMin_)*0.5;
        bool wet=true;
        splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle, wet);
    }
}

template<class CloudType>
Foam::scalar Foam::ThermoSurfaceFilmMeredith<CloudType>::a0(Foam::scalar We,Foam::scalar m, Foam::scalar b){
        return m*We+b;
    }

template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::splashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mRatio,
    const scalar We,
    const scalar Wec,
    const scalar sigma,
    bool& keepParticle,
    bool wet
)
{
    // Patch face velocity and normal
    const fvMesh& mesh = this->owner().mesh();
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // Retrieve parcel properties
    const scalar np = p.nParticle();
    const scalar m = p.mass()*np;
    const scalar d = p.d();
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;
    const vector& posC = mesh.C()[p.cell()];
    const vector& posCf = mesh.Cf().boundaryField()[pp.index()][faceI];

    // Determine direction vectors tangential to patch normal
    const vector dropTan = Ut/mag(Ut);
    const scalar beta = mag(Ut)/mag(Urel);
    const vector tanVec1 = tangentVector(nf);
    const vector tanVec2 = nf^tanVec1;

    // total mass of (all) splashed parcels
    const scalar mSplash = m*mRatio;

    // number of splashed particles per incoming particle
    //make a0 = f(We)

    const scalar Ns = a0_*(We/Wec - 1.0);
    //const scalar Ns = a0(We,a0Slope_,a0Intercept_)*(We/Wec - 1.0);
    //DEBUG(a0_);
    //DEBUG(Ns);

    // average diameter of splashed particles
    const scalar dBarSplash = 1/cbrt(6.0)*cbrt(mRatio/Ns)*d + ROOTVSMALL;
    //DEBUG(dBarSplash);

    // cumulative diameter splash distribution
    //const scalar dMax = 0.25*cbrt(mRatio)*d;
    //const scalar dMin = 0.01*dMax;
    const scalar dMax = 0.0016;
    const scalar dMin = 0.00008;
    const scalar K = exp(-dMin/dBarSplash) - exp(-dMax/dBarSplash);

    // surface energy of secondary parcels [J]
    scalar ESigmaSec = 0;

    // sample splash distribution to determine secondary parcel diameters
    scalarList dNew(parcelsPerSplash_);
    scalarList npNew(parcelsPerSplash_);
    forAll(dNew, i)
    {
        const scalar y = rndGen_.sample01<scalar>();
        //scalar y = 0.5;
        //if(i==1){
        //    y = rndGen_.sample01<scalar>();
        //}

        dNew[i] = max(dMin,-dBarSplash*log(exp(-dMin/dBarSplash) - y*K+SMALL));//have andy update to this
        npNew[i] = mRatio*np*pow3(d)/pow3(dNew[i])/parcelsPerSplash_;
        ESigmaSec += npNew[i]*sigma*p.areaS(dNew[i]);//kvm, we don't need to multiply by ptr->np here?
    }
    //DEBUG(ESigmaSec);

    // incident kinetic energy [J]
    const scalar EKIn = 0.5*m*magSqr(Urel); //this depends on number of particles in parcel--that's not right
    //DEBUG(EKIn);

    // incident surface energy [J]
    const scalar ESigmaIn = np*sigma*p.areaS(d);//kvm, we don't need to multiply by np here?
    //DEBUG(ESigmaIn);

    // dissipative energy
    const scalar Ed = max(EKInCoeff_*EKIn, 0.75*Wec/12*pi*sigma*sqr(d)*np); //second term controls energy available for splashed vel
    //DEBUG(Ed);

    // total energy [J]
    const scalar EKs = EKIn + ESigmaIn - ESigmaSec - Ed;
    //DEBUG(EKs);

    // switch to absorb if insufficient energy for splash
    if (EKs <= 0)
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
        return;
    }

    // helper variables to calculate magUns0
    const scalar logD = log(d);
    const scalar coeff2 = log(dNew[0]) - logD + ROOTVSMALL;
    scalar coeff1 = 0.0;
    forAll(dNew, i)
    {
        coeff1 += sqr(log(dNew[i]) - logD);
    }

    // magnitude of the normal velocity of the first splashed parcel
    const scalar magUns0 =
        sqrt(2.0*parcelsPerSplash_*EKs/mSplash/(1.0 + coeff1/sqr(coeff2)));

    // Set splashed parcel properties
    forAll(dNew, i)
    {
        vector dirVec;
        if(wet){
            dirVec = splashDirection(tanVec1, tanVec2, -nf, wetThetaSiMin_, wetThetaSiMax_,wetThetaSiMean_,wetThetaSiVar_,dropTan,beta);
        }
        else{
            dirVec = splashDirection(tanVec1, tanVec2, -nf, dryThetaSiMin_, dryThetaSiMax_,dryThetaSiMean_,dryThetaSiVar_,dropTan,beta);
        }

        // Create a new parcel by copying source parcel
        parcelType* pPtr = new parcelType(p);

        pPtr->origId() = pPtr->getNewParticleID();

        pPtr->origProc() = Pstream::myProcNo();

        if (splashParcelType_ >= 0)
        {
            pPtr->typeId() = splashParcelType_;
        }

        // perturb new parcels towards the owner cell centre
        pPtr->position() += 0.5*rndGen_.sample01<scalar>()*(posC - posCf);

        pPtr->nParticle() = mRatio*np*pow3(d)/pow3(dNew[i])/parcelsPerSplash_;

        pPtr->d() = dNew[i];

        pPtr->U() = dirVec*(mag(Cf_*Ut) + magUns0*(log(dNew[i]) - logD)/coeff2);

        // Apply correction to velocity for 2-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());

        // Add the new parcel
        this->owner().addParticle(pPtr);

        nParcelsSplashed_++;
    }

    // transfer remaining part of parcel to film 0 - splashMass can be -ve
    // if entraining from the film
    const scalar mDash = m - mSplash;
    static scalar cM=0;
    static scalar cMSplash=0;
    static scalar cMDash=0;
    if(p.typeId()==1){
        cM+=m;
        cMSplash+=mSplash;
        cMDash+=mDash;
        Info << " cM " << cM; 
        Info << " cMSplash " << cMSplash; 
        Info << " cMDash " << cMDash; 
        Info << endl;
        Info << "origId " <<p.origId() << " m " << m << " mSplash " << mSplash << " mDash " << mDash;
        Info << " nParticle " << p.nParticle(); 
        Info << endl;
    }
    absorbInteraction(filmModel, p, pp, faceI, mDash, keepParticle);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilmMeredith<CloudType>::ThermoSurfaceFilmMeredith
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceFilmModel<CloudType>(dict, owner, typeName),
    rndGen_(owner.rndGen()),
    thermo_
    (
        owner.db().objectRegistry::template lookupObject<SLGThermo>("SLGThermo")
    ),
    TFilmPatch_(0),
    CpFilmPatch_(0),
    interactionType_
    (
        interactionTypeEnum(this->coeffDict().lookup("interactionType"))
    ),
    deltaWet_(0.0),
    splashParcelType_(0),
    parcelsPerSplash_(0),
    Adry_(0.0),
    Awet_(0.0),
    Cf_(0.0),
    wetSplashRatioMin_(0.2),
    wetSplashRatioMax_(1.1),
    drySplashRatioMin_(0.2),
    drySplashRatioMax_(0.8),
    wetThetaSiMin_(5),
    wetThetaSiMax_(50),
    dryThetaSiMin_(5),
    dryThetaSiMax_(50),
    EKInCoeff_(0.8),
    a0_(5),
    a0Slope_(0),
    a0Intercept_(5),
    normalDistribution_(false),
    nParcelsSplashed_(0)
{
    Info<< "    Applying " << interactionTypeStr(interactionType_)
        << " interaction model" << endl;

    if (interactionType_ == itSplashBai)
    {
        this->coeffDict().lookup("deltaWet") >> deltaWet_;
        splashParcelType_ =
            this->coeffDict().lookupOrDefault("splashParcelType", -1);
        parcelsPerSplash_ =
            this->coeffDict().lookupOrDefault("parcelsPerSplash", 2);
        this->coeffDict().lookup("Adry") >> Adry_;
        this->coeffDict().lookup("Awet") >> Awet_;
        this->coeffDict().lookup("Cf") >> Cf_;
        wetSplashRatioMin_ =
            this->coeffDict().lookupOrDefault("wetSplashRatioMin", 0.2);
        wetSplashRatioMax_ =
            this->coeffDict().lookupOrDefault("wetSplashRatioMax", 1.1);
        drySplashRatioMin_ =
            this->coeffDict().lookupOrDefault("drySplashRatioMin", 0.2);
        drySplashRatioMax_ =
            this->coeffDict().lookupOrDefault("drySplashRatioMax", 0.8);
        wetThetaSiMin_ =
            this->coeffDict().lookupOrDefault("wetThetaSiMin", 40.0);
        wetThetaSiMax_ =
            this->coeffDict().lookupOrDefault("wetThetaSiMax", 85.0);
        wetThetaSiMean_ =
            this->coeffDict().lookupOrDefault("wetThetaSiMean", 40.0);
        wetThetaSiVar_ =
            this->coeffDict().lookupOrDefault("wetThetaSiVar", 10.0);
        dryThetaSiMin_ =
            this->coeffDict().lookupOrDefault("dryThetaSiMin", 80.0);
        dryThetaSiMax_ =
            this->coeffDict().lookupOrDefault("dryThetaSiMax", 89.0);
        dryThetaSiMean_ =
            this->coeffDict().lookupOrDefault("dryThetaSiMean", 85.0);
        dryThetaSiVar_ =
            this->coeffDict().lookupOrDefault("dryThetaSiVar", 2.0);
        EKInCoeff_ =
            this->coeffDict().lookupOrDefault("EKInCoeff", 0.8);
        a0_ =
            this->coeffDict().lookupOrDefault("a0", 5.0);
        a0Slope_ =
            this->coeffDict().lookupOrDefault("a0Slope", 0.0);
        a0Intercept_ =
            this->coeffDict().lookupOrDefault("a0Intercept", 5.0);
        normalDistribution_ = 
            this->coeffDict().template lookupOrDefault<Switch>("normalDistribution",false);
    }
}


template<class CloudType>
Foam::ThermoSurfaceFilmMeredith<CloudType>::ThermoSurfaceFilmMeredith
(
    const ThermoSurfaceFilmMeredith<CloudType>& sfm
)
:
    SurfaceFilmModel<CloudType>(sfm),
    rndGen_(sfm.rndGen_),
    thermo_(sfm.thermo_),
    TFilmPatch_(sfm.TFilmPatch_),
    CpFilmPatch_(sfm.CpFilmPatch_),
    interactionType_(sfm.interactionType_),
    deltaWet_(sfm.deltaWet_),
    splashParcelType_(sfm.splashParcelType_),
    parcelsPerSplash_(sfm.parcelsPerSplash_),
    Adry_(sfm.Adry_),
    Awet_(sfm.Awet_),
    Cf_(sfm.Cf_),
    nParcelsSplashed_(sfm.nParcelsSplashed_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilmMeredith<CloudType>::~ThermoSurfaceFilmMeredith()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoSurfaceFilmMeredith<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Retrieve the film model from the owner database
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel =
        const_cast<regionModels::surfaceFilmModels::surfaceFilmModel&>
        (
            this->owner().db().time().objectRegistry::template
                lookupObject<regionModels::surfaceFilmModels::surfaceFilmModel>
                (
                    "surfaceFilmProperties"
                )
        );

    const label patchI = pp.index();

    if (filmModel.isRegionPatch(patchI))
    {
        const label faceI = pp.whichFace(p.face());

        switch (interactionType_)
        {
            case itBounce:
            {
                bounceInteraction(p, pp, faceI, keepParticle);

                break;
            }
            case itAbsorb:
            {
                const scalar m = p.nParticle()*p.mass();
                absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);

                break;
            }
            case itSplashBai:
            {
                bool dry = this->deltaFilmPatch_[patchI][faceI] < deltaWet_;
//ilass                dry=false;

                if (dry)
                {
                    drySplashInteraction(filmModel, p, pp, faceI, keepParticle);
                }
                else
                {
                    wetSplashInteraction(filmModel, p, pp, faceI, keepParticle);
                }

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool ThermoSurfaceFilmMeredith<CloudType>::transferParcel"
                    "("
                        "parcelType&, "
                        "const polyPatch&, "
                        "bool&"
                    ")"
                )   << "Unknown interaction type enumeration"
                    << abort(FatalError);
            }
        }

        // transfer parcel/parcel interactions complete
        return true;
    }

    // parcel not interacting with film
    return false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const label primaryPatchI,
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchI,
        primaryPatchI,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, TFilmPatch_);

    CpFilmPatch_ = filmModel.Cp().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    SurfaceFilmModel<CloudType>::setParcelProperties(p, filmFaceI);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFaceI];
    p.Cp() = CpFilmPatch_[filmFaceI];
}


template<class CloudType>
void Foam::ThermoSurfaceFilmMeredith<CloudType>::info(Ostream& os)
{
    SurfaceFilmModel<CloudType>::info(os);

    label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");
    label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());

    os  << "    New film splash parcels         = " << nSplashTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsSplashed", nSplashTotal);
        nParcelsSplashed_ = 0;
    }
}


// ************************************************************************* //
