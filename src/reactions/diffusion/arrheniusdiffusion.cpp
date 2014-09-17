#include "arrheniusdiffusion.h"
#include "../../soluteparticle.h"

using namespace kMC;

ArrheniusDiffusion::ArrheniusDiffusion(SoluteParticle *reactant, int dx, int dy, int dz) :
    DiffusionReaction(reactant, dx, dy, dz)
{

}

double ArrheniusDiffusion::calcRate()
{
    return linearRateScale()*exp(beta()*reactant()->energy());
}
