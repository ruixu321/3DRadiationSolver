#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "basicReactingMultiphaseCloud.H"

#include "ThermoSurfaceFilmMeredith.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceFilmModelType(ThermoSurfaceFilmMeredith, basicReactingMultiphaseCloud);
}
