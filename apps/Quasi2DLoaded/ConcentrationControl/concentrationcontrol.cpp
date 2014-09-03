#include "concentrationcontrol.h"

#include "../quasidiffusionevents.h"

using namespace kMC;

uint ConcentrationControl::nSolvants() const
{
    return concentration()*m_movingWall->freeVolume();
}


void ConcentrationControl1D::registerPopulationChange(int change)
{

    double dC = 1.0/m_movingWall->freeVolume();

    m_concentrations(0) += change*dC;

    BADAssBool(concentration() >= 0 && concentration() <= 1,
               "Illegal concentration values initiated.", [&] ()
    {
        BADAssSimpleDump(change, concentration(), m_movingWall->freeVolume());
    });

}

void ConcentrationControl1D::diffuse(const double dt)
{
    vec newConcentration(m_nCells);

    double fac = dt/(m_dxSquared);


    newConcentration(m_nCells - 1) = m_concentrations(m_nCells - 1) + fac*(m_boundaryConcentration - 2*m_concentrations(m_nCells - 1) + m_concentrations(m_nCells - 2));

    for (int cell = m_nCells - 2; cell > 0; --cell)
    {
        newConcentration(cell) = m_concentrations(cell) + fac*(m_concentrations(cell + 1) - 2*m_concentrations(cell) + m_concentrations(cell - 1));
    }

    double boundaryCondition = m_concentrations(0);
    newConcentration(0) = m_concentrations(0) + fac*(m_concentrations(1) - 2*m_concentrations(0) + boundaryCondition);

    m_concentrations = newConcentration;

    if (m_movingWall->nTimesExecuted() % MainLattice::saveFileSpacing() == 0)
    {
        m_concentrations.save(KMCSolver::instance()->filePath() + "conc0.arma");
    }
}




ConcentrationControl3D::ConcentrationControl3D(const double cBoundary, const double diffusivity, const uint nCells, const double dH) :
    ConcentrationControl(cBoundary, diffusivity, dH/nCells),
    m_dH(dH),
    m_nCells(nCells)
{

}


uint ConcentrationControl3D::r() const
{
    return m_movingWall->span()/2 + m_nCells;
}

void ConcentrationControl3D::initialize()
{
    m_concentrationField.set_size(m_movingWall->length(), 2*m_nCells + m_movingWall->height(), r());

    m_concentrationField.fill(m_boundaryConcentration);
}


void ConcentrationControl3D::diffuse(const double dt)
{
    cout << "hey" << endl;
    (void) dt;

    int yCentered;
    uint z;

    uint s = m_movingWall->span();
    uint H = 2*m_nCells + s;

    for (int y = 0; y < H; ++y)
    {
        yCentered = y - (int)(m_nCells + s/2);

        z = sqrt((H*H)/4 - yCentered*yCentered);

        for (uint x = 0; x < m_movingWall->length(); ++x)
        {
            m_concentrationField(x, y, z) = sqrt(pow(((int)y - (int)yCentered), 2) + z*z);
        }
    }

    m_concentrationField.save(KMCSolver::instance()->filePath() + "conc3D.arma", raw_ascii);
    exit(1);
}
