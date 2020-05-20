#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "basicReactingCloud.H"

#include "ThermoSurfaceFilmMeredith.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceFilmModelType(ThermoSurfaceFilmMeredith, basicReactingCloud);
}
