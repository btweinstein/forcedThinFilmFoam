/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    // Determine the growth field BEFORE incrementing time.
    forAll(mesh.C(), cell_i){
        cur_h = h[cell_i];

        if(cur_h < SMALL) growth[cell_i] = 0; // Zero check
        else if(cur_h < d) growth[cell_i] = a*h;
        else if (h < hmax) growth[cell_i] = d*h;
        else growth[cell_i] = 0; // No growth past hmax
    };

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Get the laplacian of h for the surface tension
        volScalarField lap_h('lap_h', fvc::laplacian(h));

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvm::div(h*vBottom)
              - ((rho*g)/(3*mu))*fvm::laplacian(Foam::pow(h,3),h)
              + (sigma/(3*mu))*fvm::laplacian(Foam::pow(h,3), lap_h)
              ==
              growth
            );
            hEqn.solve();
        }

        // Determine the growth field BEFORE incrementing time.
        forAll(mesh.C(), cell_i){
            cur_h = h[cell_i];

            if(cur_h < SMALL) growth[cell_i] = 0; // Zero check
            else if(cur_h < d) growth[cell_i] = a*h;
            else if (h < hmax) growth[cell_i] = d*h;
            else growth[cell_i] = 0; // No growth past hmax
        };

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    
    return 0;
}


// ************************************************************************* //
