#ifndef REACTION_H
#define REACTION_H

#include "../site.h"
#include <sys/types.h>


class Site;
class KMCSolver;

class Reaction
{
public:
    Reaction();

    virtual void calcRate() = 0;
    virtual bool isActive() = 0;
    virtual void execute() = 0;

    void setSite(Site* site) {
        reactionSite = site;
    }

    void setMainsolver(KMCSolver* solver);

    uint ID() {
        return m_ID;
    }

    double rate() {
        return m_rate;
    }


protected:

    uint NX;
    uint NY;
    uint NZ;

    uint m_ID;
    static uint IDcount;

    Site* reactionSite;

    double m_rate;


    uint x()
    {
        return reactionSite->x();
    }

    uint y()
    {
        return reactionSite->y();
    }

    uint z()
    {
        return reactionSite->z();
    }

    KMCSolver* mainSolver;

};

#endif // REACTION_H
