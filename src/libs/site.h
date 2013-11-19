#ifndef SITE_H
#define SITE_H

#include <armadillo>
using namespace arma;

class Reaction;

class Site
{
public:
    Site(uint _x, uint _y, uint _z);

    void setOrigin(uint _x, uint _y, uint _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    void addReaction(Reaction* reaction);

    void activate() {
        m_active = true;
    }

    void deactivate() {
        m_active = false;
    }

    bool active() {
        return m_active;
    }

private:

    uint x;
    uint y;
    uint z;

    bool m_active = false;

    std::vector<Reaction*> siteReactions;

};

#endif // SITE_H
