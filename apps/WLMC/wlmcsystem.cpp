#include "wlmcsystem.h"
#include "wlmcwindow.h"

using namespace WLMC;

WLMCSystem::WLMCSystem(const uint nParticles,
                       const uint NX,
                       const uint NY,
                       const uint NZ,
                       function<double()> URNG) :
    m_nParticles(nParticles),
    m_NX(NX),
    m_NY(NY),
    m_NZ(NZ),
    m_volume(m_NX*m_NY*m_NZ),
    m_freeVolume(m_volume - nParticles),
    m_URNG(URNG)
{

}

void WLMCSystem::performMove(WLMCWindow *window)
{
    const uint particle = m_URNG()*m_nParticles;
    const uint destination = m_URNG()*m_freeVolume;

    uint xd, yd, zd;
    findDestination(destination, xd, yd, zd);

    double oldEnergy = getTotalEnergy();
    double newEnergy = oldEnergy + getEnergyDifference(particle, xd, yd, zd);

    uint oldBin = window->getBin(oldEnergy);
    uint newBin = window->getBin(newEnergy);

    if (!window->isLegal(oldBin) || !window->isLegal(newBin))
    {
        changePosition(xd, yd, zd);
    }

    double oldDOS = window->DOS(oldBin);
    double newDOS = window->DOS(newBin);


    bool accepted = true;

    if (oldDOS < newDOS)
    {
        accepted = (m_URNG() < oldDOS/newDOS);
    }

    if (accepted)
    {
        changePosition(xd, yd, zd);

        window->registerVisit(newBin);
    }

    else
    {
        window->registerVisit(oldBin);
    }




}

void WLMCSystem::findDestination(const uint destination, uint &xd, uint &yd, uint &zd)
{
    uint search = 0;

    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (!isOccupiedLoction(x, y, z))
                {
                    if (search == destination)
                    {
                        xd = x;
                        yd = y;
                        zd = z;

                        return;
                    }

                    search++;
                }
            }
        }
    }

}

