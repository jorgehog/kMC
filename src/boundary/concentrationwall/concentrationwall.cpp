#include "concentrationwall.h"

#include "../../kmcsolver.h"

#include "../../soluteparticle.h"

using namespace kMC;

ConcentrationWall::ConcentrationWall(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::ConcentrationWall)
{

}

ConcentrationWall::~ConcentrationWall()
{

}


void ConcentrationWall::update()
{
    performConcentrationBoundaryConditionStep();
}
