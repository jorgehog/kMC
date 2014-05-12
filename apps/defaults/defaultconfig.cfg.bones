buildTrace = 0;

System = {

    BoxSize = [100, 100, 100];

    nNeighborsLimit = 2;


    SaturationLevel = 0.001;


    #0 = Periodic
    #1 = Edge
    #2 = ConcentrationWall
    Boundaries = {
    #            #back #front
         types = ([0,    0],   #X
                  [0,    0],   #Y
                  [0,    0]);  #Z

    };

};

Reactions = {

    beta = 1.0;
    scale = 1.0;

    Diffusion = {

        rPowers =   [1.0];
        strengths = [1.0];

    };

};

Initialization = {

    someValue = 0;

};

Solver = {

    nCycles = 1000000;
    cyclesPerOutput = 100;

    #seedType:
    #0 = from time
    #1 = use specific seed

    seedType = 0;
    specificSeed = 1394447431;

};
