#include "concentrationwall.h"


using namespace kMC;

ConcentrationWall::ConcentrationWall(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation)
{

}

void ConcentrationWall::loadConfig(const Setting &setting)
{
    (void) setting;
}

void ConcentrationWall::update()
{

}
