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

    uint clearing            = getSurfaceSetting<uint>  (initCFG, "clearing");

    uint maxSpan             = getSurfaceSetting<uint>  (initCFG, "maxSpan");

    double occupancyTreshold = getSurfaceSetting<double>(initCFG, "treshold");


    uint NX = solver->NX();

    uint NY = solver->NY();


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

    double rf = 1/std::sqrt(2);

    //fsund Diamond Square code:
    DiamondSquare BottomGenerator(power2, RNG::Uniform, Seed::initialSeed);
    DiamondSquare TopGenerator   (power2, RNG::Uniform, Seed::initialSeed+power2);

    auto _Bottom = BottomGenerator.generate(H,              //Hurst Exponent
                                            {},             //Random Corners
                                            sigma,          //Deviation of random displacements
                                            rf,             //RandomFactor: What the fuck is this?
                                            false,          //Addition: Some weird shit I don't understand!
                                            true);          //Periodic Boundary Conditions

    auto _Top    = TopGenerator.generate(1-H,
                                        {},
                                        sigma,
                                        rf,
                                        false,
                                        true);

    //Don't want to learn how to do this with std::min(fucked up shit)

    double bottomMin = INFINITY;
    double topMin    = INFINITY;

    double bottomMax = 0;
    double topMax = 0;

    for (uint i = 0; i < solver->NX(); ++i)
    {
        for (uint j = 0; j < solver->NY(); ++j)
        {

            double bottomVal = _Bottom.at(i).at(j);
            double topVal    = _Top.at(i).at(j);

            if (topVal < topMin)
            {
                topMin = topVal;
            }

            if (topVal > topMax)
            {
                topMax = topVal;
            }

            if (bottomVal < bottomMin)
            {
                bottomMin = bottomVal;
            }

            if (bottomVal > bottomMax)
            {
                bottomMax = bottomVal;
            }

        }
    }


    //Normalize generated boundary to [0 - maxSpan] and rid the std::vectors
    //Translate from continuous to discrete surface points.

    double bottomSpan = bottomMax - bottomMin;
    double topSpan = topMax - topMin;

    umat Bottom(NX, NY);
    umat Top(NX, NY);

    for (uint i = 0; i < solver->NX(); ++i)
    {
        for (uint j = 0; j < solver->NY(); ++j)
        {

            Bottom(i, j) = std::round((_Bottom.at(i).at(j) - bottomMin)/bottomSpan*maxSpan);
            Top(i, j)    = std::round((_Top.at(i).at(j) - topMin)/topSpan*maxSpan);


        }
    }

    _Bottom.clear();
    _Top.clear();


    //calculate how many percent of the total z = z' surface is populated

    vec occupancyBottom(maxSpan, fill::zeros);
    vec occupancyTop   (maxSpan, fill::zeros);

    for (uint i = 0; i < solver->NX(); ++i)
    {
        for (uint j = 0; j < solver->NY(); ++j)
        {

            for (uint z = 0; z < Bottom(i, j); ++z)
            {
                occupancyBottom(z)++;
            }

            for (uint z = 0; z < Top(i, j); ++z)
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

    uint newNZ = 2*maxSpan - (topCutoff + bottomCutoff) + clearing;

    solver->setBoxSize({NX, NY, newNZ});


    Site * currentSite;

    for (uint x = 0; x < solver->NX(); ++x)
    {
        for (uint y = 0; y < solver->NY(); ++y)
        {

            uint bottomEnd = Bottom(x, y);
            uint topEnd    = Top(x, y);

            //Points below the cutoff are filtered out. Avoiding uints going to negative values (overflow)
            uint bottomSurface = (bottomEnd > bottomCutoff) ? (     (bottomEnd - bottomCutoff)) : 0;
            uint topSurface    =    (topEnd > topCutoff)    ? (newNZ - (topEnd - topCutoff)) : newNZ;

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
                    if (KMC_RNG_UNIFORM() < solver->targetConcentration())
                    {
                        currentSite->activate();
                    }
                }

            }
        }

    }
}

