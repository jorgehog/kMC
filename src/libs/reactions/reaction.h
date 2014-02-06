#ifndef REACTION_H
#define REACTION_H

#include "../site.h"
#include <sys/types.h>

#include <libconfig_utils/libconfig_utils.h>

class Site;
class KMCSolver;

class Reaction
{
public:

    Reaction();

    virtual ~Reaction();

    virtual void calcRate() = 0;
    virtual bool isActive() = 0;
    virtual void execute() = 0;

    void setSite(Site* site) {
        reactionSite = site;
    }

    static void resetAll() {
        IDcount = 0;
    }

    const uint & ID() {
        return m_ID;
    }

    const double &  rate() {
        return m_rate;
    }

    const static double & getScale() {
        return mu;
    }

    static void setSolverPtr(KMCSolver * solver);

    static void loadReactionSettings(const Setting & setting);


protected:

    static uint NX;
    static uint NY;
    static uint NZ;

    static double beta;
    static double mu;

    uint m_ID;
    static uint IDcount;

    Site* reactionSite;

    double m_rate;


    const uint & x()
    {
        return reactionSite->x();
    }

    const uint & y()
    {
        return reactionSite->y();
    }

    const uint & z()
    {
        return reactionSite->z();
    }

    static KMCSolver* mainSolver;

};

#endif // REACTION_H
