#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

#include "diamondSquare/src/diamondSquare/diamondSquare.h"

using namespace libconfig;
using namespace kMC;


void initialize_diamondSquareSurface(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/diamondSquareSurface.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("diamondSquareSurface");

    KMCDebugger_SetEnabledTo(getSurfaceSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_diamondSquareSurface(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    delete solver;


    return 0;

}

enum RNG
{
    AllZero,
    Uniform,
    Normal
};

void initialize_diamondSquareSurface(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG  = getSurfaceSetting(root, "Initialization");

    double H                 = getSurfaceSetting<double>(initCFG, "H");

    double sigma             = getSurfaceSetting<double>(initCFG, "sigma");

    bool addition            = static_cast<bool>(
                               getSurfaceSetting<int>   (initCFG, "addition"));


    uint clearing            = getSurfaceSetting<uint>  (initCFG, "clearing");

    double occupancyTreshold = getSurfaceSetting<double>(initCFG, "treshold");


    double NX = static_cast<double>(solver->NX());

    double NY = static_cast<double>(solver->NY());


    //Hack to support non quadratic lattices
    if (NX > NY)
    {
        NY = NX;
    }
    else
    {
        NX = NY;
    }


    int power2 = ceil(std::log2(NX - 1));

    //fsund Diamond Square code:
    DiamondSquare generator(power2, RNG::Uniform, Seed::initialSeed);

    auto Bottom = generator.generate(H,              //Hurst Exponent
                                     {0, NX, 0, NY}, //Corners
                                     sigma,          //Deviation of random displacements
                                     0.5,            //RandomFactor: What the fuck is this?
                                     addition,       //Some weird shit I don't understand!
                                     true);          //Periodic Boundary Conditions

    auto Top    = generator.generate(H,
                                     {0, NX, 0, NY},
                                     sigma,
                                     0.5,
                                     addition,
                                     true);


    uint maxHeight = 0;
    uint topMax    = 0;
    uint bottomMax = 0;

    uint currentHeight;

    //We should be guaranteed to be in range ..!
    for (uint i = 0; i < solver->NX(); ++i)
    {
        for (uint j = 0; j < solver->NY(); ++j)
        {

            uint bottomZ = Bottom.at(i).at(j);
            uint topZ    = Top.at(i).at(j);

            currentHeight = bottomZ + topZ;

            if (currentHeight > maxHeight)
            {
                maxHeight = currentHeight;
            }

            if (bottomZ > bottomMax)
            {
                bottomMax = bottomZ;
            }

            if (topZ > topMax)
            {
                topMax = topZ;
            }

        }
    }

    //calculate how many percent of the total z = z' surface is populated

    vec occupancyBottom(bottomMax, fill::zeros);
    vec occupancyTop   (topMax,    fill::zeros);

    for (uint i = 0; i < solver->NX(); ++i)
    {
        for (uint j = 0; j < solver->NY(); ++j)
        {

            uint bottomZ = Bottom.at(i).at(j);
            uint topZ    = Top.at(i).at(j);

            for (uint z = 0; z < bottomZ; ++z)
            {
                occupancyBottom(z)++;
            }

            for (uint z = 0; z < topZ; ++z)
            {
                occupancyTop(z)++;
            }

        }
    }

    uint area = solver->NX()*solver->NY();

    occupancyBottom /= area;
    occupancyTop    /= area;


    //Going from top to bottom, ending when the population is lower than
    //the given treshold.

    uint bottomCutoff = 0;
    while (occupancyTop(bottomCutoff) > occupancyTreshold)
    {
        bottomCutoff++;
    }

    uint topCutoff = 0;
    while (occupancyTop(topCutoff) > occupancyTreshold)
    {
        topCutoff++;
    }


    //calculate the new box height

    uint newNZ =  maxHeight - (topCutoff + bottomCutoff) + clearing;

    uvec newN  = solver->NVec();
    newN(2) = newNZ;

    solver->setBoxSize(newN);


    Site * currentSite;

    for (uint x = 0; x < solver->NX(); ++x)
    {
        for (uint y = 0; y < solver->NY(); ++y)
        {

            uint bottomEnd = static_cast<uint>(Bottom.at(x).at(y));
            uint topEnd    = static_cast<uint>(Top.at(x).at(y));

            //Points below the cutoff are filtered out. Avoiding uints going to negative values (overflow)
            uint bottomSurface = (bottomEnd > bottomCutoff) ? (     (bottomEnd - bottomCutoff)) : 0;
            uint topSurface    =    (topEnd > topCutoff)    ? (newNZ - (topEnd - topCutoff)   ) : newNZ;

            for (uint z = 1; z < newNZ - 1; ++z)
            {
                currentSite = solver->getSite(x, y, z);

                //Fill in crystals below the bottom and above the top surface
                if ((z < bottomSurface) || (z >= topSurface))
                {
                    currentSite->spawnAsCrystal();
                }

                //Fill the cavity with solution
                else if (currentSite->isLegalToSpawn())
                {
                    if (KMC_RNG_UNIFORM() < solver->targetSaturation())
                    {
                        currentSite->activate();
                    }
                }

            }
        }

    }
}

