#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import ref, man


#----------------------------------------------------------------------------
def read_controls( args, runTime, pimple ):
   
    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )

    nAlphaCorr = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaCorr" ) ) )

    nAlphaSubCycles = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaSubCycles" ) ) )

    if nAlphaSubCycles > 1 and nCorrPIMPLE != 1:
        from Foam.OpenFOAM import ext_Info, nl
        ref.ext_Info() << args.executable() << "FATAL ERROR: Sub-cycling alpha is only allowed for PISO operations, i.e. when the number of outer-correctors = 1" << ref.nl;
        import os; os_exit( 1 ) 
        pass
    
    return adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles

    
    
#----------------------------------------------------------------------------
def _createFields( runTime, mesh, g ):
        
    ref.ext_Info() << "Reading field p_rgh\n" << ref.nl
    p_rgh = man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.MUST_READ,
                                              ref.IOobject.AUTO_WRITE ),
                              mesh )

    ref.ext_Info() << "Reading field alpha1\n" << ref.nl
    alpha1 = man.volScalarField( man.IOobject( ref.word( "alpha1" ),
                                               ref.fileName( runTime.timeName() ),
                                               mesh,
                                               ref.IOobject.MUST_READ,
                                               ref.IOobject.AUTO_WRITE ),
                                 mesh )

    ref.ext_Info() << "Calculating field alpha1\n" << ref.nl
    alpha2 = man.volScalarField( ref.word( "alpha2" ), 1.0 - alpha1 )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    phi = man.createPhi( runTime, mesh, U )

    ref.ext_Info() << "Reading transportProperties\n" << ref.nl
    
    twoPhaseProperties = man.twoPhaseMixture (U, phi)
    
    rho10 = ref.dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase1Name() ).lookup( ref.word( "rho0" ) ) )
    rho20 = ref.dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase2Name() ).lookup( ref.word( "rho0" ) ) )
   
    psi1 = ref.dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase1Name() ).lookup( ref.word( "psi" ) ) )
    psi2 = ref.dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase2Name() ).lookup( ref.word( "psi" ) ) )

    pMin = ref.dimensionedScalar( twoPhaseProperties.lookup( ref.word( "pMin" ) ) )
    
    ref.ext_Info() << "Calculating field g.h\n" << ref.nl
    gh = man.volScalarField( ref.word( "gh" ), g & man.volVectorField( mesh.C(), man.Deps( mesh ) ) )
    
    ghf = man.surfaceScalarField( ref.word( "ghf" ), g & man.surfaceVectorField( mesh.Cf(), man.Deps( mesh ) ) )

    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.NO_READ,
                                          ref.IOobject.AUTO_WRITE ),
                        ( ( p_rgh + gh * ( alpha1 * rho10 + alpha2 * rho20 ) ) / ( 1.0 - gh * ( alpha1 * psi1 + alpha2 * psi2 ) ) ).ext_max( pMin ) ) #

    rho1 = rho10 + psi1 * p
    rho2 = rho20 + psi2 * p

    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.READ_IF_PRESENT,
                                            ref.IOobject.AUTO_WRITE ),
                              alpha1 * rho1 + alpha2 * rho2 )

    # Mass flux
    # Initialisation does not matter because rhoPhi is reset after the
    # alpha1 solution before it is used in the U equation.
    rhoPhi = man.surfaceScalarField( man.IOobject( ref.word( "rho*phi" ),
                                                   ref.fileName( runTime.timeName() ),
                                                   mesh,
                                                   ref.IOobject.NO_READ,
                                                   ref.IOobject.NO_WRITE ),
                                     man.fvc.interpolate( rho ) * phi )

    dgdt = alpha2.pos() * man.fvc.div( phi ) / alpha2.ext_max( 0.0001 )

    # Construct interface from alpha1 distribution
    interface = man.interfaceProperties( alpha1, U, twoPhaseProperties )

    # Construct incompressible turbulence model
    turbulence = man.incompressible.turbulenceModel.New( U, phi, twoPhaseProperties )

    return p_rgh, alpha1, alpha2, U, phi, twoPhaseProperties, rho10, rho20, psi1, psi2, pMin, \
           gh, ghf, p, rho1, rho2, rho, rhoPhi, dgdt, interface, turbulence
    

#--------------------------------------------------------------------------------------
def alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr ):

    alphaScheme = ref.word( "div(phi,alpha)" )
    alpharScheme = ref.word( "div(phirb,alpha)" )

    phir = phic*interface.nHatf()
    
    for gCorr in range( nAlphaCorr ):
        Sp = ref.volScalarField.DimensionedInternalField( ref.IOobject( ref.word( "Sp" ),
                                                                        ref.fileName( runTime.timeName() ),
                                                                        mesh ),
                                                           mesh,
                                                           ref.dimensionedScalar( ref.word( "Sp" ), dgdt.dimensions(), 0.0 ) )

        Su = ref.volScalarField.DimensionedInternalField( ref.IOobject( ref.word( "Su" ),
                                                                        ref.fileName( runTime.timeName() ),
                                                                        mesh ),
                                                          # Divergence term is handled explicitly to be
                                                          # consistent with the explicit transport solution
                                                          divU * alpha1.ext_min( 1.0 ) )
        for celli in range( dgdt.size() ):
            if dgdt[ celli ] > 0.0 and alpha1[ celli ] > 0.0:
                Sp[ celli ] -= dgdt[ celli ] * alpha1[ celli ]
                Su[ celli ] += dgdt[ celli ] * alpha1[ celli ]
                pass
            elif dgdt[ celli ] < 0.0 and alpha1[ celli ] < 1.0:
                Sp[ celli ] += dgdt[ celli ] * ( 1.0 - alpha1[ celli ] )
                pass
            pass

        phiAlpha1 = ref.fvc.flux( phi, alpha1, alphaScheme ) + ref.fvc.flux( - ref.fvc.flux( -phir, alpha2, alpharScheme ), alpha1, alpharScheme )
        
        ref.MULES.explicitSolve( ref.geometricOneField(), alpha1, phi, phiAlpha1, Sp, Su, 1.0, 0.0 )
        
        rho1f = ref.fvc.interpolate( rho1 )
        rho2f = ref.fvc.interpolate( rho2 )
        rhoPhi << phiAlpha1 * ( rho1f - rho2f ) + phi * rho2f

        alpha2 << 1.0 - alpha1

        pass

    
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "Liquid phase volume fraction = " << alpha1.weightedAverage( mesh.V() ).value() \
               << "  Min(alpha1) = " << alpha1.ext_min().value() << "  Min(alpha2) = " << alpha2.ext_min().value() << nl
    pass


#----------------------------------------------------------------------------------------
def alphaEqnsSubCycle( runTime, pimple, mesh, phi, alpha1, alpha2, rho, rho1, rho2, rhoPhi, dgdt, interface ):
    
    nAlphaCorr = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaCorr" ) ) )
    nAlphaSubCycles = ref.readLabel( pimple.dict().lookup( ref.word( "nAlphaSubCycles" ) ) )
    
    phic = ( phi() / mesh.magSf() ).mag() # mixed calculations
    phic << ( interface.cAlpha() * phic ).ext_min( phic.ext_max() )
    
    divU = ref.fvc.div( phi )
    
    if nAlphaSubCycles > 1:

        totalDeltaT = runTime.deltaT()
        rhoPhiSum = 0.0 * rhoPhi

        alphaSubCycle = ref.subCycle_volScalarField( alpha1, nAlphaSubCycles )
        for item in alphaSubCycle: 
            alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr )
            rhoPhiSum += ( runTime.deltaT() / totalDeltaT ) * rhoPhi
            pass
       # To make sure that variable in the local scope will be destroyed
       # - during destruction of this variable it performs some important actions
       # - there is a difference between C++ and Python memory management, namely
       # if C++ automatically destroys stack variables when they exit the scope,
       # Python relay its memory management of some "garbage collection" algorithm
       # that do not provide predictable behavior on the "exit of scope"
        del alphaSubCycle

        rhoPhi << rhoPhiSum
        pass
    else:
        alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr )
        pass

    if pimple.corr() == 1:
       interface.correct()
       pass
    
    pass


#--------------------------------------------------------------------------------------
def fun_UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, pimple ):
    muEff = ref.surfaceScalarField( ref.word( "muEff" ),
                                    twoPhaseProperties.muf() + ref.fvc.interpolate( rho * turbulence.ext_nut() ) )

    UEqn = ref.fvm.ddt( rho, U ) + ref.fvm.div( rhoPhi, U ) - ref.fvm.laplacian( muEff, U ) - ( ref.fvc.grad( U ) & ref.fvc.grad( muEff ) )
    
    UEqn.relax()

    if pimple.momentumPredictor():
       ref.solve( UEqn == \
                   ref.fvc.reconstruct( ( ref.fvc.interpolate( interface.sigmaK() ) * ref.fvc.snGrad( alpha1 ) - ghf * ref.fvc.snGrad( rho ) \
                                                                                                 - ref.fvc.snGrad( p_rgh ) ) * mesh.magSf() ) )
       pass
    
    return UEqn


#--------------------------------------------------------------------------------------
def fun_pEqn( runTime, mesh, pimple, UEqn, p, p_rgh, phi, U, rho, rho1, rho2, rho10, rho20, gh, ghf, dgdt, pMin, \
              psi1, psi2, alpha1, alpha2, interface ):
    rAU = 1.0/UEqn.A()
    
    rAUf = ref.fvc.interpolate( rAU )

    p_rghEqnComp = None

    if pimple.transonic():
        p_rghEqnComp = ref.fvm.ddt( p_rgh ) + ref.fvm.div( phi, p_rgh ) - ref.fvm.Sp( ref.fvc.div( phi ), p_rgh )
        pass
    else:
        p_rghEqnComp = ref.fvm.ddt( p_rgh ) + ref.fvc.div( phi, p_rgh ) - ref.fvc.Sp( ref.fvc.div( phi ), p_rgh ) 
        pass

    U << rAU * UEqn.H()

    phiU = ref.surfaceScalarField( ref.word( "phiU" ),
                                   ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) )

    phi << (phiU + ( ref.fvc.interpolate( interface.sigmaK() ) * ref.fvc.snGrad( alpha1 ) - ghf * ref.fvc.snGrad( rho ) ) * rAUf * mesh.magSf() )

    while (pimple.correctNonOrthogonal()):
        p_rghEqnIncomp = ref.fvc.div( phi ) - ref.fvm.laplacian( rAUf, p_rgh ) 
        
        ref.solve( ( alpha1.ext_max( 0.0 ) * ( psi1 / rho1 ) + alpha2.ext_max( 0.0 ) * ( psi2 / rho2 ) ) *p_rghEqnComp() + p_rghEqnIncomp,   ##
               mesh.solver( p_rgh.select( pimple.finalInnerIter() ) ) )

        if pimple.finalNonOrthogonalIter():
            dgdt << ( alpha2.pos() * ( psi2 / rho2 ) - alpha1.pos() * ( psi1 / rho1 ) ) * ( p_rghEqnComp & p_rgh )
            phi += p_rghEqnIncomp.flux()
            pass

    U += rAU * ref.fvc.reconstruct( ( phi() - phiU ) / rAUf ) # mixed calculations
    U.correctBoundaryConditions()

    p << ( ( p_rgh() + gh * ( alpha1() * rho10 + alpha2 * rho20 ) ) /( 1.0 - gh * ( alpha1() * psi1 + alpha2 * psi2 ) ) ).ext_max( pMin ) #

    rho1 << rho10 + psi1 * p
    rho2 << rho20 + psi2 * p

    ref.ext_Info() << "max(U) " << U.mag().ext_max().value() << ref.nl
    ref.ext_Info() << "min(p_rgh) " << p_rgh.ext_min().value() << ref.nl
    pass
    
    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    g = man.readGravitationalAcceleration( runTime, mesh)
    
    pimple = ref.pimpleControl( mesh )
    
    adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles = read_controls( args, runTime, pimple )

    cumulativeContErr = ref.initContinuityErrs()

    p_rgh, alpha1, alpha2, U, phi, twoPhaseProperties, rho10, rho20, psi1, psi2, pMin, \
                      gh, ghf, p, rho1, rho2, rho, rhoPhi, dgdt, interface, turbulence = _createFields( runTime, mesh, g )
    
    CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
    
    runTime = ref.setInitialDeltaT( runTime, adjustTimeStep, maxCo, CoNum )

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl

    while runTime.run():
        adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles = read_controls( args, runTime, pimple )

        CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )

        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )

        runTime.increment()

        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        # --- Outer-corrector loop
        while pimple.loop():
            alphaEqnsSubCycle( runTime, pimple, mesh, phi, alpha1, alpha2, rho, rho1, rho2, rhoPhi, dgdt, interface )

            ref.solve( ref.fvm.ddt( rho ) + ref.fvc.div( rhoPhi ) )

            UEqn = fun_UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, pimple )

            # --- Pressure corrector loop
            while pimple.correct(): 
                fun_pEqn( runTime, mesh, pimple, UEqn, p, p_rgh, phi, U, rho, rho1, rho2, rho10, rho20, gh, ghf, dgdt, pMin, \
                          psi1, psi2, alpha1, alpha2, interface )
                pass
            if pimple.turbCorr():
                turbulence.correct()
                pass
            pass
        rho <<  alpha1 * rho1 + alpha2 * rho2

        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n\n" << ref.nl
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher \n "     
   pass


#--------------------------------------------------------------------------------------

