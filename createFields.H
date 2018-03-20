Info<< "Reading field h\n" << endl;

volScalarField h
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

volScalarField lap_h
(
    IOobject
    (
        "lap_h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Reading imposed velocity field

Info<< "Reading velocity field on bottom h\n" << endl;

volVectorField vBottom
(
    IOobject
    (
        "vBottom",
        runTime.timeName(), // IN THE CONSTANT DIRECTORY
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(vBottom) & mesh.Sf()
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
dimensionedScalar rho
(
    thinFilmTransportProperties.lookup("rho")
);

Info<< "Reading film surface tension\n" << endl;
dimensionedScalar sigma
(
    thinFilmTransportProperties.lookup("sigma")
);

Info<< "Reading gravitational acceleration\n" << endl;
//Note that the code assumes that g is perpendicular to the field.
dimensionedScalar g
(
    thinFilmTransportProperties.lookup("g")
);

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

// Creating the growth field...exponential growth proportional to h is multiplied by growth_factor

volScalarField growth_factor
(
    IOobject
    (
        "growth_factor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    a
);

volScalarField growth
(
    IOobject
    (
        "growth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h*growth_factor // It is updated later
);




// Other miscellaneous stuff
Info << "Cutoff value for h is " << SMALL << "\n" << endl;