#include "reaction.h"
#include "../site.h"

Reaction::Reaction(Site* site)
{
    site->addReaction(this);
}

