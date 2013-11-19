#include "site.h"
#include "kmcsolver.h"

Site::Site(uint _x, uint _y, uint _z) :
    x(_x),
    y(_y),
    z(_z)
{

}

void Site::addReaction(Reaction *reaction)
{
    siteReactions.push_back(reaction);
}
