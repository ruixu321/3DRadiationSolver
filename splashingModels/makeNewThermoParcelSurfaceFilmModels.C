#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "basicThermoCloud.H"

#include "ThermoSurfaceFilmMeredith.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceFilmModelType(ThermoSurfaceFilmMeredith, basicThermoCloud);
}
