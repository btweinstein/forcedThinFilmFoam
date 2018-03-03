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
  
    while (runTime.run())
    {

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //*** Calculate the concentration field ***

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix cEqn
            (
                fvm::ddt(c)
              - fvm::laplacian(D, c)
            );
            cEqn.solve();
        }

        //*** Calculate the resulting pressure field
        rho = rho_0*(1 + beta*c);
        volVectorField rhog("rhog", rho*g);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(p)
              - fvc::div(rhog)
            );
            //pEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
            //pEqn.setReference(0, 0);
            // Assumes that pressure at one boundary is set.
            pEqn.solve();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    
    return 0;
}


// ************************************************************************* //
