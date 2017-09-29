/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaBuoyancy.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaBuoyancy, 0);
addToRunTimeSelectionTable(RASModel, kOmegaBuoyancy, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaBuoyancy::kOmegaBuoyancy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    Cmu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alpha",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    nut_ = k_/(omega_ + omegaSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaBuoyancy::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kOmegaBuoyancy::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kOmegaBuoyancy::divDevReff() const
{
    return
    (
      - fvm::laplacian(nuEff(), U_)
      - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
    );
}


bool kOmegaBuoyancy::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaBuoyancy::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(omega_, omega0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G("RASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

////////////////////////////////////////////////////////////////////////
// Buoyancy correction -start (Brecht DEVOLDER, 19 september 2017)
////////////////////////////////////////////////////////////////////////
        
    // Access to the density
    volScalarField& rho_ = const_cast<volScalarField&>
    (
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>("rho")
    );

    // Mass flux
    surfaceScalarField rhoPhi = fvc::interpolate(rho_)*this->phi_;

    // Gravitational acceleration
    dimensionedVector g
    (
        "g",
        dimensionSet(0, 1, -2, 0, 0, 0, 0),
        vector(0, 0, -9.81)
    );
    
    // Constant coefficients
    scalar sigmaT = 0.85;	//turbulent Prandtl number (dimensionless)
    dimensionedScalar kSmall
    (
        "kSmall",
        k_.dimensions(),
        SMALL
    );
    
    // Buoyancy correction term
    volScalarField Gb("Gb", -nut_/sigmaT*(g & fvc::grad(rho_)));

////////////////////////////////////////////////////////////////////////
// Buoyancy correction -end (Brecht DEVOLDER, 19 september 2017)
////////////////////////////////////////////////////////////////////////

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(rho_, omega_)
      + fvm::div(rhoPhi, omega_)
      + fvm::SuSp(-fvc::div(rhoPhi), omega_)
      - fvm::laplacian(rho_*DomegaEff(), omega_)
     ==
        rho_*alpha_*G*omega_/k_
      - fvm::Sp(rho_*beta_*omega_, omega_)
    );

    omegaEqn().relax();

    // No longer needed: matrix completes at the point of solution
    // HJ, 17/Apr/2012
//     omegaEqn().completeAssembly();

    solve(omegaEqn);
    bound(omega_, omega0_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(rhoPhi, k_)
      + fvm::SuSp(-fvc::div(rhoPhi), k_)
      - fvm::laplacian(rho_*DkEff(), k_)
     ==
        rho_*G
      + fvm::Sp(Gb/max(k_, kSmall), k_) //buoyancy correction in k-eqn (Brecht DEVOLDER, 19 september 2017)
      - fvm::Sp(rho_*Cmu_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = k_/(omega_ + omegaSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
