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

\*---------------------------------------------------------------------------*/

inline const Foam::radiation::radiativeIntensityRayBand&
Foam::radiation::fvDOMBand::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::fvDOMBand::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::fvDOMBand::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::fvDOMBand::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::fvDOMBand::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::fvDOMBand::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMBand::a() const
{
    return a_;
}

// WSGGM
// added by Ivan Sikic 10/10/2014
inline const Foam::volScalarField& Foam::radiation::fvDOMBand::ggCoeff() const
{

    return ggCoeff_;
}

inline const Foam::volScalarField& Foam::radiation::fvDOMBand::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}

// WSGGM
// added by Ivan Sikic 10/10/2014
inline const Foam::volScalarField& Foam::radiation::fvDOMBand::ggCoeffLambda
(
    const label lambdaI
) const
{
    return ggCoeffLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::fvDOMBand::G() const
{
    return G_;
}

// WSGGM
// added by Ivan Sikic 26/03/2015
inline const Foam::volScalarField& Foam::radiation::fvDOMBand::kG() const
{
    return kG_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMBand::Qr() const
{
    return Qr_;
}

inline const Foam::volScalarField& Foam::radiation::fvDOMBand::Qin() const
{
    return Qin_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOMBand::Qem() const
{
    return Qem_;
}


inline const Foam::radiation::blackBodyEmission&
Foam::radiation::fvDOMBand::blackBody() const
{
    return blackBody_;
}


inline const Foam::fvScalarMatrix& Foam::radiation::fvDOMBand::fvRayDiv
(
    const label rayId,
    const label lambdaI
) const
{
    return fvRayDiv_[lambdaI][rayId];
}


inline bool Foam::radiation::fvDOMBand::cacheDiv() const
{
    return cacheDiv_;
}


inline Foam::scalar Foam::radiation::fvDOMBand::omegaMax() const
{
    return omegaMax_;
}


// ************************************************************************* //
