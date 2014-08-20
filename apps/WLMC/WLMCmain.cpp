#include <kMC>
#include <BADAss/badass.h>
#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

#include "kmcwlmcsystem.h"

using namespace H5Wrapper;
using namespace libconfig;
using namespace kMC;


void initializeWLMC(KMCSolver *solver, const Setting & root);

int main()
{

    Config cfg;

    cfg.readFile("infiles/WLMC.cfg");

    const Setting & root = cfg.getRoot();

    KMCDebugger_SetEnabledTo(false);

    KMCSolver::enableDumpLAMMPS(false);
    KMCSolver::enableDumpXYZ(false);
    MainLattice::enableEventFile(false);

    KMCSolver* solver = new KMCSolver(root);

    initializeWLMC(solver, root);

    return 0;

}


void initializeWLMC(KMCSolver *solver, const Setting &root)
{
    const Setting &initCFG = getSetting(root, "Initialization");


    const Setting &output = getSetting(initCFG, "output");

    const string &path = getSetting<string>(output, "path");
    const string &filename = getSetting<string>(output, "filename");
    const string &name = getSetting<string>(output, "name");

    const bool &overwrite = getSetting<int>(output, "overwrite") == 1;

    Root hdf5root(path + "/" + filename);

    const uint &adaptiveWindows = getSetting<uint>(initCFG, "adaptiveWindows");

    const uint &nParticlesStart = getSetting<uint>(initCFG, "nParticlesStart");
    const uint &nParticlesStop = getSetting<uint>(initCFG, "nParticlesStop");

    BADAss(nParticlesStart, <=, nParticlesStop);

    const uint &movesPerSampling = getSetting<uint>(initCFG, "movesPerSampling");

    const uint &nbins = getSetting<uint>(initCFG, "nbins");
//    const uint &nbins = solver->volume();

    const uint &nbinsOverOverlap = getSetting<uint>(initCFG, "nbinsOverOverlap");
    const uint overlap = nbins/nbinsOverOverlap;

    const uint &nbinsOverMinWindowSize = getSetting<uint>(initCFG, "nbinsOverMinWindowSize");
    const uint minWindowSize = nbins/nbinsOverMinWindowSize;


    const double &flatnessCriterion = getSetting<double>(initCFG, "flatnessCriterion");

    const double &deflationLimit = getSetting<double>(initCFG, "deflationLimit");

    const double &flatnessGradientTreshold = getSetting<double>(initCFG, "flatnessGradientTreshold");

    const double &logfStart = getSetting<double>(initCFG, "logfStart");

    const double &logfFinal = getSetting<double>(initCFG, "logfFinal");

    KMCWLMCSystem *system;

    WLMC::Window *mainWindow;

    Member &dataGroup = hdf5root.addMember(name);
    Member &potentialGroup = dataGroup.addMember(DiffusionReaction::potentialString());
    Member &systemGroup = potentialGroup.addMember(solver->volume());

    for (uint nParticles = nParticlesStart; nParticles <= nParticlesStop; ++nParticles)
    {
        solver->insertRandomParticles(nParticles);

        system = new KMCWLMCSystem(solver,
                                   movesPerSampling,
                                   flatnessCriterion,
                                   overlap,
                                   minWindowSize,
                                   flatnessGradientTreshold,
                                   deflationLimit);

        //memory will be freed by the wrapper
        Member &particles = systemGroup.addMember(nParticles, overwrite);

        mainWindow = system->execute(nbins, adaptiveWindows, logfStart, logfFinal);

        vec DOS = arma::exp(mainWindow->logDOS());
        particles.addData("logDOS", mainWindow->logDOS());
        particles.addData("DOS", DOS);
        particles.addData("EBins", mainWindow->energies());
        particles.addData("nbins", nbins);

        particles.file()->flush(H5F_SCOPE_GLOBAL);

        delete system;
        delete mainWindow;

        solver->clearParticles();
    }


}
