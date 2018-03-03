Info<< "Reading field h\n" << endl;

volScalarField c
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading velocity field on bottom h\n" << endl;

volVectorField vBottom
(
    IOobject
    (
        "vBottom",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// *** THIN FILM TRANSPORT PROPERTIES ***
Info<< "Reading thinFilmTransportProperties\n" << endl;

IOdictionary thinFilmTransportProperties
(
    IOobject
    (
        "thinFilmTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading film viscosity\n" << endl;
dimensionedScalar mu
(
    thinFilmTransportProperties.lookup("mu")
);

Info<< "Reading film density\n" << endl;
dimensionedScalar mu
(
    thinFilmTransportProperties.lookup("rho")
);

Info<< "Reading film surface tension\n" << endl;
dimensionedScalar sigma
(
    thinFilmTransportProperties.lookup("sigma")
);

#include "readGravitationalAcceleration.H"

// *** THIN FILM GROWTH PROPERTIES ***

Info<< "Reading thinFilmGrowthProperties\n" << endl;

IOdictionary thinFilmGrowthProperties
(
    IOobject
    (
        "thinFilmGrowthProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading growth distance from front\n" << endl;
dimensionedScalar d
(
    thinFilmGrowthProperties.lookup("d")
);

Info<< "Reading growth rate\n" << endl;
dimensionedScalar a
(
    thinFilmGrowthProperties.lookup("a")
);

Info<< "Reading maximum film height\n" << endl;
dimensionedScalar hmax
(
    thinFilmGrowthProperties.lookup("hmax")
);