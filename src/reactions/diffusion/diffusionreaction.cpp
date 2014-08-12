#include "diffusionreaction.h"

#include "../../kmcsolver.h"

#include "../../boundary/boundary.h"

#include "../../soluteparticle.h"

#include "../../debugger/debugger.h"

#include "../../potential/potential.h"


using namespace kMC;

DiffusionReaction::DiffusionReaction(SoluteParticle *reactant, int dx, int dy, int dz) :
    Reaction(reactant)
{

    m_path[0] = dx;
    m_path[1] = dy;
    m_path[2] = dz;

    m_pathIndices[0] = dx + 1;
    m_pathIndices[1] = dy + 1;
    m_pathIndices[2] = dz + 1;

}

DiffusionReaction::~DiffusionReaction()
{

}

void DiffusionReaction::loadConfig(const Setting &setting)
{

    const Setting & rPowerSettings = getSetting(setting, "rPowers");
    const Setting & strengthSettings = getSetting(setting, "strengths");

    if (rPowerSettings.getLength() != strengthSettings.getLength())
    {
        KMCSolver::exit("Mismatch in diffusionreaction species parameter lengths.");
    }

    SoluteParticle::nSpecies(rPowerSettings.getLength());

    vector<double> rPowers;
    vector<double> strengths;

    for (uint i = 0; i < SoluteParticle::nSpecies(); ++i)
    {
        rPowers.push_back(rPowerSettings[i]);
        strengths.push_back(strengthSettings[i]);
    }

    setPotentialParameters(rPowers, strengths);

}

uint DiffusionReaction::xD() const
{
    Boundary::setupCurrentBoundary(x(), 0);
    return Boundary::currentBoundaries(0)->transformCoordinate(x() + m_path[0]);
}

uint DiffusionReaction::yD() const
{
    Boundary::setupCurrentBoundary(y(), 1);
    return Boundary::currentBoundaries(1)->transformCoordinate(y() + m_path[1]);
}

uint DiffusionReaction::zD() const
{
    Boundary::setupCurrentBoundary(z(), 2);
    return Boundary::currentBoundaries(2)->transformCoordinate(z() + m_path[2]);
}

void DiffusionReaction::setPotentialParameters(const vector<double> &rPowers, const vector<double> &strenghts)
{

    if (rPowers.size() != strenghts.size())
    {
        KMCSolver::exit("mismatch in species mixing data.");
    }

    SoluteParticle::nSpecies(rPowers.size());

    m_rPowers.reset();
    m_rPowers.set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

    m_strengths.reset();
    m_strengths.set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

    for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
    {
        for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
        {
            //geometric mean for powers and arithmetic mean for strengths
            m_rPowers(speciesA, speciesB)   = std::sqrt(rPowers.at(speciesA)*rPowers.at(speciesB));
            m_strengths(speciesA, speciesB) = (strenghts.at(speciesA) + strenghts.at(speciesB))/2;
        }
    }

}

string DiffusionReaction::getFinalizingDebugMessage() const
{
#ifndef KMC_NO_DEBUG

    if (!Debugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    s << Reaction::getFinalizingDebugMessage();


    if (!Debugger::lastCurrentReaction->isType("DiffusionReaction"))
    {
        return s.str();
    }

    const DiffusionReaction *lastReaction = (DiffusionReaction*)Debugger::lastCurrentReaction;

    if (lastReaction != NULL)
    {
        X = lastReaction->path(0);
        Y = lastReaction->path(1);
        Z = lastReaction->path(2);
    }

    s << "\nDestination of last active reaction site marked on current site:\n\n";
    s << reactant()->info(X, Y, Z);

    return s.str();
#else
    return "";
#endif
}

void DiffusionReaction::setupPotential()
{

    using std::pow;

    uint i, j, k;
    double rPower, strength;

    m_potential.reset_objects();
    m_potential.reset();
    m_potential.set_size(Site::neighborhoodLength(),
                         Site::neighborhoodLength(),
                         Site::neighborhoodLength());


    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            for (int z = -1; z <= 1; ++z)
            {
                if ((x == 0) && (y == 0) && (z == 0))
                {
                    continue;
                }

                i = x + 1;
                j = y + 1;
                k = z + 1;

                m_pathLengths(i, j, k) = sqrt(x*x + y*y + z*z);

            }
        }
    }

    for (uint i = 0; i < Site::neighborhoodLength(); ++i)
    {
        for (uint j = 0; j < Site::neighborhoodLength(); ++j)
        {
            for (uint k = 0; k < Site::neighborhoodLength(); ++k)
            {

                m_potential(i, j, k).set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

                for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
                {
                    for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
                    {

                        if (i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit())
                        {
                            m_potential(i, j, k)(speciesA, speciesB) = numeric_limits<double>::max();
                            continue;
                        }

                        rPower = m_rPowers(speciesA, speciesB);
                        strength = m_strengths(speciesA, speciesB);

                        m_potential(i, j, k)(speciesA, speciesB) = strength/std::pow(Site::originTransformVector(i)*Site::originTransformVector(i)
                                                                                     + Site::originTransformVector(j)*Site::originTransformVector(j)
                                                                                     + Site::originTransformVector(k)*Site::originTransformVector(k),
                                                                                     rPower/2);
                    }
                }
            }

        }
    }

}


void DiffusionReaction::execute()
{
    reactant()->changePosition(Site::boundaries(0, 0)->transformCoordinate(x() + m_path[0]),
            Site::boundaries(1, 0)->transformCoordinate(y() + m_path[1]),
            Site::boundaries(2, 0)->transformCoordinate(z() + m_path[2]));
}

const string DiffusionReaction::info(int xr, int yr, int zr, string desc) const
{

    (void) xr;
    (void) yr;
    (void) zr;
    (void) desc;

    stringstream s;
    s << Reaction::info(m_path[0], m_path[1], m_path[2], "D");

    s << "Reaction initiates diffusion to\n\n";

    if (destinationSite() == NULL)
    {
        s << "BOUNDARY";
    }

    else
    {
        s << Site::info(xD(), yD(), zD(), -m_path[0], -m_path[1], -m_path[2], "O");
    }

    s << "\nPath: " << m_path[0] << " " << m_path[1] << " " << m_path[2] << endl;

    string full_string = s.str();

    return full_string;

}

bool DiffusionReaction::isAllowed() const
{

    if (destinationSite() == NULL)
    {
        return false;
    }

    return !destinationSite()->isActive();
}

void DiffusionReaction::reset()
{
    Reaction::reset();
}

void DiffusionReaction::registerBetaChange(const double newBeta)
{
    (void) newBeta;

    BADAssBool(rate() != UNSET_RATE, "Beta should not be changed untill rates have been calculated.");

    _setRate(rate() * m_betaChangeScaleFactor);

}


mat        DiffusionReaction::m_rPowers   = {1.0};
mat        DiffusionReaction::m_strengths = {1.0};

double     DiffusionReaction::m_betaChangeScaleFactor;

field<mat> DiffusionReaction::m_potential;

cube::fixed<3, 3, 3> DiffusionReaction::m_pathLengths;
