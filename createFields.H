Info<< "Reading field c\n" << endl;

volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating pressure field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusion constant D\n" << endl;

dimensionedScalar D
(
    transportProperties.lookup("D")
);

Info<< "Reading base density of the media rho_0\n" << endl;

dimensionedScalar rho_0
(
    transportProperties.lookup("rho_0")
);

Info<< "Reading solute expansion coefficient beta \n" << endl;

dimensionedScalar beta
(
    transportProperties.lookup("beta")
);

#include "readGravitationalAcceleration.H"

Info<< "Creating density field\n" << endl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho_0*(1 + beta*c)
);