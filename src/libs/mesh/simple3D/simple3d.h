#ifndef SIMPLE3D_H
#define SIMPLE3D_H

#include "../mesh.h"

class simple3D : public Mesh
{
public:
    simple3D();

    // Mesh interface
public:
    void initialize();
};

#endif // SIMPLE3D_H
