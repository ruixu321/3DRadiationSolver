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

inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayBand::I() const
{
    return I_;
}

// WSGGM (added by Ivan Sikic 26/03/2015)
inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayBand::kI() const
{
    return kI_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayBand::Qr() const
{
    return Qr_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayBand::Qr()
{
    return Qr_;
}


inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRayBand::Qin() const
{
    return Qin_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayBand::Qin()
{
    return Qin_;
}


inline const Foam::volScalarField& Foam::radiation::
radiativeIntensityRayBand::Qem() const
{
    return Qem_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRayBand::Qem()
{
    return Qem_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRayBand::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRayBand::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayBand::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayBand::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayBand::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRayBand::omega() const
{
    return omega_;
}


inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRayBand::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
