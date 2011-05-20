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
def read_controls( args, runTime, mesh ):
    from Foam.finiteVolume.cfdTools.general.include import readPISOControls
    piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr = readPISOControls( mesh )
    
    from Foam.finiteVolume.cfdTools.general.include import readTimeControls
    adjustTimeStep, maxCo, maxDeltaT = readTimeControls( runTime )

    from Foam.OpenFOAM import word, readLabel
    nAlphaCorr = readLabel( piso.lookup( word( "nAlphaCorr" ) ) )

    nAlphaSubCycles = readLabel( piso.lookup( word( "nAlphaSubCycles" ) ) )

    if nAlphaSubCycles > 1 and nOuterCorr != 1:
        from Foam.OpenFOAM import ext_Info, nl
        ext_Info() << args.executable() << "FATAL ERROR: Sub-cycling alpha is only allowed for PISO, i.e. when the number of outer-correctors = 1" << nl;
        import os; os_exit( 1 ) 
        pass
    
    return piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr, adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles

    
    
#----------------------------------------------------------------------------
def _createFields( runTime, mesh, g ):
    from Foam.OpenFOAM import ext_Info, nl
    from Foam.OpenFOAM import IOdictionary, IOobject, word, fileName
    from Foam.finiteVolume import volScalarField
        
    ext_Info() << "Reading field p_rgh\n" << nl
    p_rgh = volScalarField( IOobject( word( "p_rgh" ),
                                        fileName( runTime.timeName() ),
                                        mesh,
                                        IOobject.MUST_READ,
                                        IOobject.AUTO_WRITE ),
                              mesh )

    ext_Info() << "Reading field alpha1\n" << nl
    alpha1 = volScalarField( IOobject( word( "alpha1" ),
                                       fileName( runTime.timeName() ),
                                       mesh,
                                       IOobject.MUST_READ,
                                       IOobject.AUTO_WRITE ),
                             mesh )

    ext_Info() << "Calculating field alpha1\n" << nl
    from Foam.OpenFOAM import scalar
    alpha2 = volScalarField( word( "alpha2" ), scalar( 1 ) - alpha1 )

    ext_Info() << "Reading field U\n" << nl
    from Foam.finiteVolume import volVectorField
    U = volVectorField( IOobject( word( "U" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    from Foam.finiteVolume.cfdTools.incompressible import createPhi
    phi = createPhi( runTime, mesh, U )

    ext_Info() << "Reading transportProperties\n" << nl
    
    from Foam.transportModels import twoPhaseMixture
    twoPhaseProperties = twoPhaseMixture (U, phi)
    
    from Foam.OpenFOAM import dimensionedScalar
    rho10 = dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase1Name() ).lookup( word( "rho0" ) ) )
    rho20 = dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase2Name() ).lookup( word( "rho0" ) ) )
   
    psi1 = dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase1Name() ).lookup( word( "psi" ) ) )
    psi2 = dimensionedScalar( twoPhaseProperties.subDict( twoPhaseProperties.phase2Name() ).lookup( word( "psi" ) ) )

    pMin = dimensionedScalar( twoPhaseProperties.lookup( word( "pMin" ) ) )
    
    ext_Info() << "Calculating field g.h\n" << nl
    gh = volScalarField( word( "gh" ), g & mesh.C() )
    
    from Foam.finiteVolume import surfaceScalarField
    ghf = surfaceScalarField( word( "ghf" ), g & mesh.Cf() )

    p = volScalarField( IOobject( word( "p" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.NO_READ,
                                  IOobject.AUTO_WRITE ),
                        ( ( p_rgh + gh * ( alpha1 * rho10 + alpha2 * rho20 ) ) / ( 1.0 - gh * ( alpha1 * psi1 + alpha2 * psi2 ) ) ).ext_max( pMin ) )

    rho1 = rho10 + psi1 * p
    rho2 = rho20 + psi2 * p

    rho = volScalarField( IOobject( word( "rho" ),
                                    fileName( runTime.timeName() ),
                                    mesh,
                                    IOobject.READ_IF_PRESENT,
                                    IOobject.AUTO_WRITE ),
                          alpha1 * rho1 + alpha2 * rho2 )

    # Mass flux
    # Initialisation does not matter because rhoPhi is reset after the
    # alpha1 solution before it is used in the U equation.
    from Foam import fvc
    rhoPhi = surfaceScalarField( IOobject( word( "rho*phi" ),
                                           fileName( runTime.timeName() ),
                                           mesh,
                                           IOobject.NO_READ,
                                           IOobject.NO_WRITE ),
                                 fvc.interpolate( rho ) * phi )

    dgdt = alpha2.pos() * fvc.div( phi ) / alpha2.ext_max( scalar( 0.0001 ) )

    # Construct interface from alpha1 distribution
    from Foam.transportModels import interfaceProperties
    interface = interfaceProperties( alpha1, U, twoPhaseProperties )

    # Construct incompressible turbulence model
    from Foam import incompressible
    turbulence = incompressible.turbulenceModel.New( U, phi, twoPhaseProperties )

    return p_rgh, alpha1, alpha2, U, phi, twoPhaseProperties, rho10, rho20, psi1, psi2, pMin, \
           gh, ghf, p, rho1, rho2, rho, rhoPhi, dgdt, interface, turbulence
    

#--------------------------------------------------------------------------------------
def alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr ):

    from Foam.OpenFOAM import word 
    alphaScheme = word( "div(phi,alpha)" )
    alpharScheme = word( "div(phirb,alpha)" )

    phir = phic*interface.nHatf()
    
    from Foam.finiteVolume import volScalarField
    from Foam.OpenFOAM import IOobject, word, fileName, dimensionedScalar, scalar, ext_Info, nl
    for gCorr in range( nAlphaCorr ):
        Sp = volScalarField.DimensionedInternalField( IOobject( word( "Sp" ),
                                                                fileName( runTime.timeName() ),
                                                                mesh ),
                                                      mesh,
                                                      dimensionedScalar( word( "Sp" ), dgdt.dimensions(), 0.0 ) )

        Su = volScalarField.DimensionedInternalField( IOobject( word( "Su" ),
                                                                fileName( runTime.timeName() ),
                                                                mesh ),
                                                      # Divergence term is handled explicitly to be
                                                      # consistent with the explicit transport solution
                                                      divU * alpha1.ext_min( scalar( 1 ) ) )
        for celli in range( dgdt.size() ):
            if dgdt[ celli ] > 0.0 and alpha1[ celli ] > 0.0:
                Sp[ celli ] -= dgdt[ celli ] * alpha1[ celli ]
                Su[ celli ] += dgdt[ celli ] * alpha1[ celli ]
                pass
            elif dgdt[ celli ] < 0.0 and alpha1[ celli ] < 1.0:
                Sp[ celli ] += dgdt[ celli ] * ( 1.0 - alpha1[ celli ] )
                pass
            pass

        from Foam import fvc
        phiAlpha1 = fvc.flux( phi, alpha1, alphaScheme ) + fvc.flux( - fvc.flux( -phir, alpha2, alpharScheme ), alpha1, alpharScheme )
        
        from Foam import MULES
        from Foam.OpenFOAM import geometricOneField
        MULES.explicitSolve( geometricOneField(), alpha1, phi, phiAlpha1, Sp, Su, 1.0, 0.0 )
        
        rho1f = fvc.interpolate( rho1 )
        rho2f = fvc.interpolate( rho2 )
        rhoPhi.ext_assign( phiAlpha1 * ( rho1f - rho2f ) + phi * rho2f )

        alpha2.ext_assign( scalar( 1 ) - alpha1 )

        pass

    
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "Liquid phase volume fraction = " << alpha1.weightedAverage( mesh.V() ).value() \
               << "  Min(alpha1) = " << alpha1.ext_min().value() << "  Min(alpha2) = " << alpha2.ext_min().value() << nl
    pass


#----------------------------------------------------------------------------------------
def alphaEqnsSubCycle( runTime, piso, mesh, phi, alpha1, alpha2, rho, rho1, rho2, rhoPhi, dgdt, interface, oCorr ):
    
    from Foam.OpenFOAM import word,readLabel
    nAlphaCorr = readLabel( piso.lookup( word( "nAlphaCorr" ) ) )
    nAlphaSubCycles = readLabel( piso.lookup( word( "nAlphaSubCycles" ) ) )
    
    from Foam.finiteVolume import surfaceScalarField
    phic = ( phi / mesh.magSf() ).mag()
    phic.ext_assign( ( interface.cAlpha() * phic ).ext_min( phic.ext_max() ) )
    
    from Foam import fvc
    divU = fvc.div( phi )
    
    if nAlphaSubCycles > 1:

        totalDeltaT = runTime.deltaT()
        rhoPhiSum = 0.0 * rhoPhi

        from Foam.finiteVolume import subCycle_volScalarField
        alphaSubCycle = subCycle_volScalarField( alpha1, nAlphaSubCycles )
        for item in alphaSubCycle: 
            alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr )
            rhoPhiSum.ext_assign( rhoPhiSum + ( runTime.deltaT() / totalDeltaT ) * rhoPhi )
            pass
       # To make sure that variable in the local scope will be destroyed
       # - during destruction of this variable it performs some important actions
       # - there is a difference between C++ and Python memory management, namely
       # if C++ automatically destroys stack variables when they exit the scope,
       # Python relay its memory management of some "garbage collection" algorithm
       # that do not provide predictable behavior on the "exit of scope"
        del alphaSubCycle

        rhoPhi.ext_assign( rhoPhiSum )
        pass
    else:
        alphaEqns( runTime, mesh, rho1, rho2, rhoPhi, phic, dgdt, divU, alpha1, alpha2, phi, interface, nAlphaCorr )
        pass

    if oCorr == 0:
       interface.correct()
       pass
    
    pass


#--------------------------------------------------------------------------------------
def fun_UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, momentumPredictor, oCorr, nOuterCorr ):
    from Foam.OpenFOAM import word
    from Foam.finiteVolume import surfaceScalarField
    from Foam import fvc
    muEff = surfaceScalarField( word( "muEff" ),
                                twoPhaseProperties.muf() + fvc.interpolate( rho * turbulence.ext_nut() ) )
    from Foam import fvm

    UEqn = fvm.ddt( rho, U ) + fvm.div( rhoPhi, U ) - fvm.laplacian( muEff, U ) - ( fvc.grad( U ) & fvc.grad( muEff ) )
    
    UEqn.relax()

    if momentumPredictor:
       from Foam.finiteVolume import solve
       solve( UEqn == \
                   fvc.reconstruct( ( fvc.interpolate( interface.sigmaK() ) * fvc.snGrad( alpha1 ) - ghf * fvc.snGrad( rho ) \
                                                                                                 - fvc.snGrad( p_rgh ) ) * mesh.magSf(),
                                     mesh.solver( U.select( oCorr == nOuterCorr-1 ) ) ) )
       pass
    
    return UEqn


#--------------------------------------------------------------------------------------
def fun_pEqn( runTime, mesh, UEqn, p, p_rgh, phi, U, rho, rho1, rho2, rho10, rho20, gh, ghf, dgdt, pMin, \
              psi1, psi2, alpha1, alpha2, interface, transonic, oCorr, nOuterCorr, corr, nCorr, nNonOrthCorr ):
    rUA = 1.0/UEqn.A()
    
    from Foam import fvc
    rUAf = fvc.interpolate( rUA )

    p_rghEqnComp = None

    from Foam import fvm
    if transonic:
        p_rghEqnComp = fvm.ddt( p_rgh ) + fvm.div( phi, p_rgh ) - fvm.Sp( fvc.div( phi ), p_rgh )
        pass
    else:
        p_rghEqnComp = fvm.ddt( p_rgh ) + fvc.div( phi, p_rgh ) - fvc.Sp( fvc.div( phi ), p_rgh ) 
        pass

    U.ext_assign( rUA * UEqn.H() )

    from Foam.finiteVolume import surfaceScalarField
    from Foam.OpenFOAM import word
    phiU = surfaceScalarField( word( "phiU" ),
                               ( fvc.interpolate( U ) & mesh.Sf() ) + fvc.ddtPhiCorr( rUA, rho, U, phi ) )

    phi.ext_assign(phiU + ( fvc.interpolate( interface.sigmaK() ) * fvc.snGrad( alpha1 ) - ghf * fvc.snGrad( rho ) ) * rUAf * mesh.magSf() )

    from Foam.finiteVolume import solve
    from Foam.OpenFOAM import scalar
    for nonOrth in range( nNonOrthCorr +1 ):
        p_rghEqnIncomp = fvc.div( phi ) - fvm.laplacian( rUAf, p_rgh ) 
        
        solve( ( alpha1.ext_max( scalar( 0 ) ) * ( psi1 / rho1 ) + alpha2.ext_max( scalar( 0 ) ) * ( psi2 / rho2 ) ) *p_rghEqnComp() + p_rghEqnIncomp,
               mesh.solver( p_rgh.select( oCorr == ( nOuterCorr - 1 ) and corr == ( nCorr-1 ) and nonOrth == nNonOrthCorr )  ) )

        if nonOrth == nNonOrthCorr:
            dgdt.ext_assign( ( alpha2.pos() * ( psi2 / rho2 ) - alpha1.pos() * ( psi1 / rho1 ) ) * ( p_rghEqnComp & p_rgh ) )
            phi.ext_assign( phi + p_rghEqnIncomp.flux() )
            pass

    U.ext_assign( U + rUA * fvc.reconstruct( ( phi - phiU ) / rUAf ) )
    U.correctBoundaryConditions()

    p.ext_assign( ( ( p_rgh + gh * ( alpha1 * rho10 + alpha2 * rho20 ) ) /( 1.0 - gh * ( alpha1 * psi1 + alpha2 * psi2 ) ) ).ext_max( pMin ) )

    rho1.ext_assign( rho10 + psi1 * p )
    rho2.ext_assign( rho20 + psi2 * p )

    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "max(U) " << U.mag().ext_max().value() << nl
    ext_Info() << "min(p_rgh) " << p_rgh.ext_min().value() << nl
    pass
    
    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    from Foam.OpenFOAM.include import setRootCase
    args = setRootCase( argc, argv )

    from Foam.OpenFOAM.include import createTime
    runTime = createTime( args )

    from Foam.OpenFOAM.include import createMesh
    mesh = createMesh( runTime )
    
    from Foam.finiteVolume.cfdTools.general.include import readGravitationalAcceleration
    g = readGravitationalAcceleration( runTime, mesh)
    
    piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr, \
           adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles = read_controls( args, runTime, mesh )

    from Foam.finiteVolume.cfdTools.general.include import initContinuityErrs
    cumulativeContErr = initContinuityErrs()

    p_rgh, alpha1, alpha2, U, phi, twoPhaseProperties, rho10, rho20, psi1, psi2, pMin, \
                      gh, ghf, p, rho1, rho2, rho, rhoPhi, dgdt, interface, turbulence = _createFields( runTime, mesh, g )
    
    from Foam.finiteVolume.cfdTools.incompressible import CourantNo
    CoNum, meanCoNum = CourantNo( mesh, phi, runTime )
    
    from Foam.finiteVolume.cfdTools.general.include import setInitialDeltaT
    runTime = setInitialDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "\nStarting time loop\n" << nl

    while runTime.run():
        piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr, \
               adjustTimeStep, maxCo, maxDeltaT, nAlphaCorr, nAlphaSubCycles = read_controls( args, runTime, mesh )

        from Foam.finiteVolume.cfdTools.incompressible import CourantNo
        CoNum, meanCoNum = CourantNo( mesh, phi, runTime )

        from Foam.finiteVolume.cfdTools.general.include import setDeltaT
        runTime = setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )

        runTime.increment()

        ext_Info() << "Time = " << runTime.timeName() << nl << nl
        
        from Foam import fvm, fvc
        from Foam.finiteVolume import solve
        # --- Outer-corrector loop
        for oCorr in range( nOuterCorr ):
            alphaEqnsSubCycle( runTime, piso, mesh, phi, alpha1, alpha2, rho, rho1, rho2, rhoPhi, dgdt, interface, oCorr )


            solve( fvm.ddt( rho ) + fvc.div( rhoPhi ) )

            UEqn = fun_UEqn( mesh, alpha1, U, p, p_rgh, ghf, rho, rhoPhi, turbulence, g, twoPhaseProperties, interface, momentumPredictor, oCorr, nOuterCorr )

            # --- PISO loop
            for corr in range( nCorr ): 
                fun_pEqn( runTime, mesh, UEqn, p, p_rgh, phi, U, rho, rho1, rho2, rho10, rho20, gh, ghf, dgdt, pMin, \
                          psi1, psi2, alpha1, alpha2, interface, transonic, oCorr, nOuterCorr, corr, nCorr, nNonOrthCorr )
                pass

        rho.ext_assign( alpha1 * rho1 + alpha2 * rho2 )

        turbulence.correct()

        runTime.write()

        ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n\n" << nl
        pass

    ext_Info() << "End\n" << nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "010700" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam1.7.0 or higher \n "     
   pass


#--------------------------------------------------------------------------------------

