package danfulea.phys.egsOutput;

import java.io.FileWriter;
import java.util.Calendar;
import java.util.Date;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EGS4Geom;
import danfulea.phys.egs.EGS4Grid;
import danfulea.phys.egs.EGS4Macro;
import danfulea.phys.egs.EGS4SrcEns;
import danfulea.phys.egs.EgsQuestion;

/**
 * DEMO CLASS<br>
 * THIS CODE SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY. 
 * IT IS INTENDED FOR USE IN CALCULATING QUANTITIES OF INTEREST FOR THICK-WALLED ION CHAMBERS EXPOSED TO PHOTON BEAMS ALTHOUGH IT MAY 
 * BE USED SIMPLY TO SCORE DOSE IN A CYLINDRICAL GEOMETRY. IN ADDITION TO SCORING DOSE, THE FOLLOWING QUANTITIES MAY ALSO BE SCORED: <p>
 * THE Aatt CORRECTION FACTOR (attenuation of primary beam) <br>
 * THE Ascat CORRECTION FACTOR (photon scatter)<br>
 * THE Apn CORRECTION FACTOR (for point-source beams, obtained by correlated sampling) <br>
 * THE STOPPING POWER RATIO (obtained by correlated sampling)<br>
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada). 
 * @author Dan Fulea, 07 NOV. 2005
 */
public class cavrz implements EgsQuestion{
	

	//"******************************************************************************
	//"
	//"
	//"                               ********************
	//"                               *                  *
	//"                               * cavrznrc.mortran *
	//"                               *                  *
	//"                               ********************
	//"
	//"
	//"       INTRODUCTION:
	//"       THIS CODE SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A
	//"       FINITE, RIGHT CYLINDRICAL GEOMETRY.
	//"       IT IS INTENDED FOR USE IN CALCULATING QUANTITIES OF INTEREST FOR
	//"       THICK-WALLED ION CHAMBERS EXPOSED TO PHOTON BEAMS ALTHOUGH IT MAY
	//"       BE USED SIMPLY TO SCORE DOSE IN A CYLINDRICAL GEOMETRY. IN ADDITION
	//"       TO SCORING DOSE, THE FOLLOWING QUANTITIES MAY ALSO BE SCORED:
	//"
	//"          1)THE Aatt  CORRECTION FACTOR    (attenuation of primary beam),
	//"          2)THE Ascat CORRECTION FACTOR    (photon scatter),
	//"          3)THE Apn   CORRECTION FACTOR    (for point-source beams, obtained
	//"                                             by correlated sampling),
	//"          4)THE Afl   CORRECTION FACTOR    (obtained by correlated sampling),
	//"          5)<s>g,w STOPPING POWER RATIO    (obtained by correlated sampling)

	//"*******************************************************************************
	//"                       INPUT/OUTPUT CONTROL INPUT
	//"                       **************************
	//"*******************************************************************************
	//"
	//"CARD IO1
	//"
	//"  IWATCH= off         (0)  FOR NORMAL OUTPUT
	//"        = interactions(1)  OUTPUT ON EVERY DISCRETE INTERACTION
	//"        = steps       (2)  OUTPUT ON EVERY ELECTRON/PHOTON STEP AS WELL
	//"        = deposited   (3)  PRINTS OUT ONLY WHEN ENERGY IS DEPOSITED
	//"        = graph       (4)  PRINTS OUT FILE FOR GRAPHICS
	//"
	//"  STORE INITIAL RANDOM NUMBERS
	//"        = no          (0)  DO NOT STORE THE INITIAL RANDOM NUMBERS
	//"        = last        (1)  STORE THE INITIAL RANDOM NUMBER FOR THE LAST HISTORY
	//"        = all deposited (2)STORE THE INITIAL RANDOM NUMBER FOR ALL
	//"                           THAT DEPOSIT ENERGY IN THE CAVITY
	//"        = all         (3)  STORE ALL THE INITIAL RANDOM NUMBERS
	//"
	//"  IRESTART
	//"        = first       (0)  first run for this data set
	//"        = restart     (1)  restart of a previous run
	//"        = analyze     (3)  just read in the raw data and do the statistical
	//"                           analysis
	//"        = start-RNS   (4)  read starting random numbers from a file
	//"        = parallel    (5)  Combine results from previous parallel runs
	//"
	//"  OUTPUT OPTIONS
	//"        = short              (0)  SHORT OUTPUT -JUST THE CAVITY SUMMARY
	//"                                  AND THE DOSE GRID.
	//"        = cavity details     (1)  ABOVE PLUS DETAILS FOR EACH CAVITY ZONE
	//"
	//"  STORE DATA ARRAYS
	//"        = yes             (0) Store data arrays for re-use
	//"        = no              (1) don't store them
	//;
	//"*******************************************************************************
	//"
	//" MONTE CARLO DELIMETERS:    :start Monte Carlo inputs:
	//"                            :stop Monte Carlo inputs:
	//"
	//"*******************************************************************************
	//"                       MONTE CARLO CONTROL INPUT
	//"                       *************************
	//"*******************************************************************************
	//"
	//"CARD MC1
	//"
	//"  NUMBER OF HISTORIES         (I)  # HISTORIES TO RUN
	//"                                   (MIN:100, DEFAULTS TO 20 000)
	//"
	//"  INITIAL RANDOM NO. SEEDS=  INTGER1, INTEGER2
	//"                  User-code can use RANLUX or RANMAR, depending on selection
	//"                  Default is RANLUX
	//"            RANLUX
	//"                  INTEGER1 is the luxury level, use 1 to 4, 4 taking longest
	//"                           default is 1 (set by $DEFAULT-LL in ranlux.macros)
	//"                  INTEGER2 selects the independent sequence to use, it
	//"                           can be from 1 to 1073741824 (2**30)
	//"            RANMAR
	//"                  INTEGER1 is a seed between 1 and 31328 (0 =>default 1802)
	//"                  INTEGER2 is a seed between 1 and 30081 (0 =>default 9937)
	//"                    Selection of unique INTEGER2 values guarantees independent
	//"                    sequences.
	//"            Note:
	//"                After the seeds are first input and used for initialization,
	//"                the variables jrng1 and jrng2 are just pointers used by the RNG
	//"     Note that the parameter sought does not change if the generator used
	//"     is changed.  This is so we don't need to change the input files
	//"     although it makes the meaning less clear!
	//"
	//"  MAX CPU HOURS ALLOWED       (I)  MAX CPU TIME ALLOWED IN HOURS, DEFAULT=999hr
	//"
	//"  IFULL
	//"         = dose and stoppers         (0) just calculate total dose and that due
	//"                                         to stoppers and discards.
	//"         = Aatt and Ascat            (1) Above plus Aatt, Ascat
	//"
	//"  STATISTICAL ACCURACY SOUGHT        (R) % statistical accuracy of the total
	//"                                         dose in the peak region that is sought
	//"                                         The program executes until this
	//"                                         accuracy is obtained ot the CPU time
	//"                                         runs out.
	//"
	//"  PHOTON REGENERATION
	//"         = yes (ifano = 1) the calculation is performed with regeneration
	//"                           of the parent photon after they have interacted. A
	//"                           typical setting when FANO conditions are examined.
	//"         = no (ifano = 0)  a normal calculation.
	//"         = no electrons from wall (ifano = 2) secondary electrons from
	//"                           interactions in the cavity wall are immediately
	//"                           eliminated.  Photons are not regenerated.
	//"
	//"*******************************************************************************
	//"*******************************************************************************
	//"
	//" GEOMRZ DELIMETERS:    :start geometrical inputs:
	//"                       :stop geometrical inputs:
	//"
	//"*******************************************************************************
	//"                  CYLINDRICAL GEOMETRY & MATERIAL INPUT
	//"                  *************************************
	//"
	//" EXTRA INPUT FOR CAVRZnrc:    METHOD OF INPUT= Cavity information
	//"
	//"                                JUST INPUT THE WALL THICKNESS, CAVITY AND
	//"                                ELECTRODE DIMENSIONS AND MATERIALS. THE
	//"                                GEOMETRY IS THEN SET UP AUTOMATICALLY.
	//"                                (USEFUL FOR SIMPLE,SYMMETRIC GEOMETRIES).
	//"                                3 PLANAR ZONES ARE ASSUMED AND THE CAVITY
	//"                                MATERIAL IS ASSUMED TO BE GAS.
	//"
	//"CARD CGM1
	//"
	//"  METHOD OF INPUT
	//"        = Groups              (0)   INPUT GROUPS OF SLABS OF EQUAL THICKNESS
	//"        = Individual          (1)   VERBOSE INPUT OF THE GEOMETRY AND MEDIA.
	//"        = Cavity information  (2)   SEE ABOVE
	//"
	//"*******************************************************************************
	//"
	//"CARD CGM2A (ONLY IF METHOD OF INPUT= Groups)
	//"
	//"  Z OF FRONT FACE        (R)   START OF FIRST SLAB (REAL)
	//"  NSLAB                  (M)   # PLANAR SLABS IN A GROUP (INTEGERS)
	//"  SLAB THICKNESS         (M)   THICKNESS OF EACH SLAB IN THE GROUP (REALS)
	//"
	//"*******************************************************************************
	//"
	//"CARD CGM2B (ONLY IF METHOD OF INPUT= Individual)
	//"
	//"  Z OF FRONT FACE        (R)   START OF FIRST SLAB (REAL)
	//"  DEPTH BOUDARIES        (M)   GEOMETRICAL Z-PLANE COORDINATES (REALS)
	//"
	//"*******************************************************************************
	//"
	//"CARD CGM3
	//"
	//"  RADII                  (M)   RADII OF CYLINDERS DEFINING THE GEOMETRY (REALS)
	//"
	//"*******************************************************************************
	//"                                MATERIAL INPUT
	//"                                **************
	//"*******************************************************************************
	//"
	//"CARD M1
	//"
	//"  MEDIA                  (M)   TYPE OF MATERIAL (FMT='24A1')
	//"                               BY CONVENTION THE PRIMARY WALL
	//"                               MATERIAL IS FIRST, THE CAVITY MATERIAL
	//"                               IS SECOND AND THE REST FOLLOWS
	//"
	//"*******************************************************************************
	//"
	//"CARD M2
	//"
	//"  DESCRIPTION BY= Regions(0)   USING THE IRL REGION NUMBER
	//"                = Planes (1)   USING THE IX, IZ PLANES
	//"
	//"*******************************************************************************
	//"
	//"CARD M3-A  (CHOICE OF CARD M3-A AND CARD M3-B)
	//"
	//"  MEDNUM                 (M)   THE MATERIAL NUMBER (INEGERS)
	//"                               (MEDNUM=0 TO SKIP)
	//"  START REGION           (M)   INITIAL GEOMETRICAL ZONE(IRL) (INTEGERS)
	//"  STOP REGION            (M)   FINAL GEOMETRICAL ZONE(IRL) (INTEGERS)
	//"                               ( >NREGLO TO INPUT MORE THAN ONE ZONE)
	//"                               DEFAULTS:   MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	//"                                           MEDNUM=1 FOR REGION=2,NREG
	//"
	//"*******************************************************************************
	//"
	//"CARD M3-B  (CHOICE OF CARD M3-A AND CARD M3-B)
	//"
	//"  MEDNUM                 (M)   THE MATERIAL NUMBER (INEGERS)
	//"                               (MEDNUM=0 TO SKIP)
	//"  START ZSLAB            (M)   INITIAL ZSLAB (IZ) (INTEGERS)
	//"  STOP ZSLAB             (M)   FINAL ZSLAB (IZ) (INTEGERS)
	//"  START RING             (M)   INITIAL RADIAL RING (IX) (INTEGERS)
	//"  STOP RING              (M)   FINAL RADIAL RING (IX) (INTEGERS)
	//"                               ( >NREGLO TO INPUT MORE THAN ONE ZONE)
	//"                               DEFAULTS:   MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	//"                                           MEDNUM=1 FOR REGION=2,NREG
	//"
	//"*******************************************************************************
	//"*******************************************************************************
	//"
	//"  CAVITY INPUT DELIMETERS:
	//"                      :start cavity inputs:
	//"                      :stop cavity inputs:
	//"
	//"*******************************************************************************
	//"                           CAVITY INFORMATION INPUTS
	//"                           *************************
	///"*******************************************************************************
	//"
	//"CARD CAV1-A     (ONLY IF   METHOD OF INPUT= Cavity information)
	//"
	//"   WALL THICKNESS      (R)   THICKNESS OF THE CHAMBER WALLS (cms)
	//"                             (DEFAULTS TO 0.273)
	//"
	//"   CAVITY RADIUS       (R)   OUTER RADIUS OF THE CAVITY (cms)
	//"
	//"   CAVITY LENGHT       (R)   LENGTH OF THE CAVITY (cms) (DEFAULTS TO 0.2)
	//"
	//"   ELECTRODE RADIUS    (R)   RADIUS OF THE ELECTRODE (DEFAULTS TO 0.)
	//"
	//"   WALL MATERIAL       (C)   WALL MATERIAL  (clear enough?)
	//"
	//"CARD CAV1-B     (ONLY IF   ELECTRODE RADIUS > 0.0)
	//"
	//"   ELECTRODE MATERIAL  (C)   ELECTRODE MATERIAL
	//"
	//"         ---------------------------------------------------
	//"
	//"CARD CAV1-C     (ONLY IF   METHOD OF INPUT= groups or = individual)
	//"
	//"   NUMBER OF CAVITY REGIONS (I)   NUMBER OF GEOMETRICAL ZONES
	//"                                  COMPRISING THE CAVITY.
	//"
	//"
	//"   REGION NUMBERS OF THE CAVITY  (M)
	//"                                  THE ARRAY OF THESE ZONES
	///"                                  (THE CAVITY REGION NUMBERS)
	//;
	//"
	//"*******************************************************************************
	//"                             SOURCE INPUT
	//"                            **************
	//"*******************************************************************************
	//" SOURCE DELIMETERS:    :start source inputs:
	//"                       :stop source inputs:
	//"
	//"FOR ALL SOURCES
	//"                                      Charge of the incident beam
	//"  INCIDENT PARTICLE= electron   (-1)  electrons
	//"                     photon     (0)   photons
	//"                     positron   (1)   positrons
	//"
	//"  (for SOURCE 21,22)     all    (2)  include all of the particles
	//"                                     in the phase space file
	//"                                     [IQIN]
	//"                    charged     (3)  include e+ and e-
	//"
	//"  SOURCE NUMBER                 (I)   number of the source
	//"                                      [ISOURC]
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  0 <<<<<<<<
	//"
	//"     PARALLEL BEAM INCIDENT FROM THE FRONT (+VE Z-AXIS) ""toc:
	//"
	//"
	//"  SOURCE OPTIONS            (M4)  RBEAM, UINC, VINC, WINC
	//"
	//"               RBEAM          radius of parallel beam in cm
	//"                              (defaults to max radius of geometry)
	//"               UINC           incident x-axis direction cosine
	//"               VINC           incident y-axis direction cosine
	//"               WINC           incident z-axis direction cosine
	//"                              NOTE: (UINC,VINC,WINC)
	//"                              get automatically normalized
	//"                              defaults to (0.0,0.0,1.0)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  1 <<<<<<<<
	//"
	//"     POINT SOURCE ON AXIS INCIDENT FROM THE FRONT    ""toc:
	//"
	//"  SOURCE OPTIONS                (M4)  DISTZ, RBEAM, 0, 0
	//"
	//"               DISTZ          distance of the point source from the
	//"                              front of the target in cm (DEFAULT 100.)
	//"               RBEAM          radius of the beam at the front of the
	//"                              target in cm (defaults to MAX radius)
	//"
	//"------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  2 <<<<<<<<
	//"
	//"        BROAD PARALLEL BEAM INCIDENT FROM FRONT (+VE Z-AXIS) ""toc:
	//"               WITH UNIT AREA BEAM AND LARGE SCORING AREA
	//"
	//"  SOURCE OPTIONS          (M4)  0, 0, 0, 0
	//"
	//"------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  3 <<<<<<<<
	//"
	//"     UNIFORM ISOTROPICALLY RADIATING DISK OF FINITE SIZE   ""toc:
	//"            (MUST BE ALLOWED FOR IN THE GEOMETRICAL DEFINITIONS)
	//"
	//"  SOURCE OPTIONS                (M4)  RMINBM, RBEAM, ZSMIN, ZSMAX
	//"
	//"               RMINBM,RBEAM           inner and outer radii of source region
	//"                                      must be inside geometry
	//"               ZSMIN,ZSMAX            min and max z values for source
	//"
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  4 <<<<<<<<
	//"
	//"     FOR CENTRAL AXIS FLUENCE VS BEAM RADIUS      ""toc:
	//"
	//"  SOURCE OPTIONS            (M4)  RCAXIS, 0, 0, 0
	//"
	//"               RCAXIS       radius of central axis scoring zone (cm)
	//"
	//"        NOTE: this source option treats the cylindrical radii input
	//"              above as beam radii. the largest radius must be infinite
	//"              and the phantom must be homogeneous (at least in each layer)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  10 <<<<<<<<
	//"
	//"     PARALLEL BEAM INCIDENT FROM THE SIDE (+VE Y-AXIS)    ""toc:
	//"
	//"  SOURCE OPTIONS           (M4)  XBEAM, ZBEAM, 0, 0
	//"
	//"               XBEAM             half-width of the rectangular beam in cm
	//"                                 (defaults to max radius)
	//"               ZBEAM             half-height of the rectangular beam in cm
	//"                                 (defaults to max)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  11 <<<<<<<<
	//"
	//"     POINT SOURCE INCIDENT FROM THE SIDE     ""toc:
	//"
	//"
	//"  SOURCE OPTIONS                (M4)  DISTRH, XBEAM, ZBEAM, 0
	//"
	//"               DISTRH                 distance of the source from the middle
	//"                                      of the target in cm (defaults to 100.)
	//"               XBEAM                  half-width of the beam at the center of
	//"                                      the target in cm (defaults to max radius)
	//"               ZBEAM                  half-height of the beam at the center of
	//"                                      the target in cm (defaults to max)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  12 <<<<<<<<
	//"
	//"   POINT SOURCE OFF AXIS         ""toc:
	//"
	//"  SOURCE OPTIONS                (M4)  DISTRH, DISTZ, 0, 0
	//"
	//"               DISTRH                 distance of the point source off the
	//"                                      Z-axis.
	//"               DISTZ                  perpendicular distance of the
	//"                                      point source away from the front face.
	//"                                      a negative value is permitted.
	//"
	//"                                      DISTZ > 0
	//"                                      point located in front of front face
	//"
	//"                                      0 > DISTZ > -(ZPLANE(NPLANE)-ZPLANE(1))
	//"                                      point located between front and rear face
	//"
	//"                                      DISTZ < -(ZPLANE(NPLANE)-ZPLANE(1))
	//"                                      point located rear of rear plane
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  13 <<<<<<<<
	//"
	//"        PARALLEL BEAM FROM ANY ANGLE     ""toc:
	//"
	//"  SOURCE OPTIONS                (M4)  UINC, VINC, WINC, 0
	//"
	//"               UINC                   incident x-axis direction cosine
	//"               VINC                   incident y-axis direction cosine
	//"               WINC                   incident z-axis direction cosine
	//"
	//"                 NOTE: (UINC,VINC,WINC) get automatically normalized
	//"                       default is (0.0,0.0,1.0)
	//"
	//"
	//"------------------------------------------------------------------------------
	//"                     >>>>>>>> SOURCE  14 <<<<<<<<
	//"
	//"   POINT SOURCE ON AXIS INCIDENT FROM THE FRONT WITH ALL   ""toc:
	//"    EVENTS INSIDE RMINBM NOT FOLLOWED (A FUDGE FOR COLLIMATOR STUDIES)
	//"
	//"  SOURCE OPTIONS                (M4)  DISTZ, RBEAM, RMINBM, IGNORED
	//"
	//"               DISTZ                  distance of the point source from the
	//"                                      front of the target in cm
	//"                                      (defaults to 100.)
	//"               RBEAM                  radius of the beam at the front of the
	//"                                      target in cm (defaults to max radius)
	//"               RMINBM                 below this radius, all histories are
	//"                                      terminated by the source routines by
	//"                                      giving them zero weight.
	//"                                      The howfar routines must check for this.
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"
	//"
	//"                     >>>>>>>> SOURCE  15 <<<<<<<<
	//"
	//"  POINT SOURCE OFF AXIS. The same as source 12 but uses an alternative
	//"  implementation for sampling points on the surface of the RZ-geomtry. The
	//"  motivation for implementing this source was to check that source 12 is OK
	//"  and to check the effect of varying weights from the source on the
	//"  statistical uncertainty (contrary to source 12, source 16 produces
	//"  essentially  constant weights if the geometry-to-source distance is large
	//"  compared to the geometry dimension, a typical situation for ion chamber
	//"  simulations)
	//"
	//"  SOURCE OPTIONS                (M4)  DIST, ANGLE, IGNORED, IGNORED
	//"
	//"             DIST                     distance of the centre of the geometry
	//"                                      to the source in cm.
	//"             ANGLE                    angle of rotation around the x-axis.
	//"                                      (because of the cylindrical symetry,
	//"                                      rotations around the x-axis and y-axis
	//"                                      are indistinguishable). 0 degrees
	//"                                      corresponds to a source above the front
	//"                                      face (i.e. the same as source 1), 90
	//"                                      degrees to a source from the side
	//"                                      (i.e. the same as source 11).
	//"                                      The source MUST be outside the geometry,
	//"                                      otherwise the initialization routine
	//"                                      will abort execution.
	//"
	//"           Note that if you are not actually rotating about the center of the
	//"           geometry, you must calculate the angle and distance as if you
	//"           were.
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"
	//"
	//"                     >>>>>>>> SOURCE  16 <<<<<<<<
	//"
	//"           EXTENDED (CIRCULAR OR RECTANGULAR)  SOURCE OFF AXIS.
	//"
	//"  SOURCE OPTIONS                (M4)  DIST, ANGLE, TMP1, TMP2
	//"
	//"            DIST                      distance of geometry centre to source
	//"                                      centre in cm.
	//"
	//"            ANGLE                     angle of rotation around the x-axis
	//"                                      (see comments/explanations to source 15)
	//"
	//"            TMP1, if TMP2 <= 0        radius of the source (i.e., the emitting
	//"        or  TMP2, if TMP1 <= 0        position is picked uniformly within the
	//"                                      circle).
	//"
	//"            TMP1 and TMP2, if both    half-sizes of the radiating rectangle
	//"            >= 0                      in x- and y-directions before rotation,
	//"                                      i.e., initially x and y are picked
	//"                                      within the rectangle and z is set to
	//"                                      -DIST + geometry centre. Then a rotation
	//"                                      around the x-axis is performed.
	//"       In all cases the source plane is perpendicular to the line joining
	//"       it to the center of the geometry.   Note that this introduces a
	//"       slight error if the center of your geometry is not the true point
	//"       of rotation.
	//"
	//"       Note: if TEMP1 <= 0 and TEMP2 <= 0, source 16 becomes a point-source
	//"             off-axis, i.e. the same as source 12 and 15.
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  20 <<<<<<<<
	//"
	//"   RADIAL DISTRIBUTION INPUT      ""toc:
	//"
	//"  MODEIN= Local                 (0)   if radial distribution is to be input
	//"                                      locally through the .egs4inp file
	//"        = External              (1)   if the distribution is to be input
	//"                                      via an external file
	//"
	//"                        -----------------------------
	//" ONLY IF MODEIN= Local
	//"
	//"  NRDIST                        (I)   # radial bins in distribution histogram
	//"  RDISTF                        (M)   top of radial bin.
	//"                                      should be values for 1 to NRDIST.
	//"  RPDF                          (M)   Probability of initial particle being
	//"                                      in this bin.
	//"                                      Probability doesn't need to be normalized
	//"                                      but it should be in units cm**-2
	//"                                      Should be values for 1 to NRDIST.
	//"  RDIST IOUTSP= None            (0)   No distribution data in output summary
	//"              = Include         (1)   include distribution data output summary
	//"
	//"                        -----------------------------
	//"  ONLY IF MODEIN= External
	//"
	//"  RDIST FILENAME                (C)   filename(with ext) contains
	//"                                      distribution information
	//"
	//"  RDIST IOUTSP= None            (0)   No distribution data in output summary
	//"              = Include         (1)   include distribution data output summary
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  21 <<<<<<<<
	//"
	//"    FULL BEAM PHASE-SPACE BEAM DATA, INCIDENT ON FRONT FACE    ""toc:
	//"
	//"  SOURCE OPTIONS                (M4)  IMODE, 0, 0, 0
	//"
	//"               IMODE              0=> 7 variables/record: X,Y,U,V,E,WT,LATCH
	//"                                  2=> 8 variables/record: the above + ZLAST
	//"
	//"  FILSPC                        (C)   filename (with ext) contains
	//"                                      phase space information
	//"                                      (maximum of 80 characters)
	//"                                      (assigned to unit 42)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  22 <<<<<<<<
	//"
	//"    FULL BEAM PHASE-SPACE BEAM DATA FROM ANY ANGLE, INSIDE OR OUTSIDE   ""toc:
	//"
	//"    PARTICLES ARE READ IN FROM A BEAM PHASE SPACE and placed on a plane
	//"    described by the SOURCE OPTIONS inputs (see below). Then it is checked
	//"    whether they are already inside the geometry. If yes, the region index
	//"    is determined and the shower intiated. If not, it is checked whether
	//"    the particle trajectory will intersect the geometry (assuming that the
	//"    geometry is surrounded by vacuum). If not, the particle is rejected and
	//"    the next one taken from the phase-space file. If yes, the particle
	//"    is placed on the entry point and the shower is initiated.
	//"
	//"  SOURCE OPTIONS                (M4)  IMODE, DIST, ANGLE, ZOFFSET
	//"
	//"               IMODE              0=> 7 variables/record: X,Y,U,V,E,WT,LATCH
	//"                                  2=> 8 variables/record: the above + ZLAST
	//"               DIST               Perpendicular distance of the phase-space
	//"                                  plane to the point of rotation in cm.
	//"               ANGLE              Angle of rotation in degrees. The rotation
	//"                                  is performed around an axis that is parallel
	//"                                  to the x-axis and passes through the point
	//"                                  (x,y,z)=(0,0,ZOFFSET).
	//"               ZOFFSET            Point of rotation. If |ZOFFSET| > 1e4,
	//"                                  the centre of the geometry is taken as
	//"                                  the point of rotation (but note that
	//"                                  the maximum value allowed by the input
	//"                                  routine is 1e6, so that |ZOFFSET| must
	//"                                  be between 1e4 and 1e6 to use the centre
	//"                                  of the geometry automatically).
	//"
	//"    Examples:
	//"       - to place a phase-space on the upper z-face of the geometry,
	//"         use DIST=0, ANGLE=0, ZOFFSET=zplane(1)
	//"         This is the same as source 21
	//"       - to place a phase space on the lower z-face of the geometry,
	//"         use DIST=0, ANGLE=180, ZOFFSET=zplane(n)
	//"       - to have a phase file incident from, say, 60 degrees with
	//"         a distance to the centre of the geometry of 30 cm, use
	//"         DIST=30, ANGLE=60, ZOFFSET=9999.
	//"   etc.
	//"
	//"  FILSPC                        (C)   filename (with ext) contains
	//"                                      phase space information
	//"                                      (maximum of 80 characters)
	//"                                      (assigned to unit 42)
	//"
	//"*******************************************************************************
	//"
	//"            NOT REQUIRED IF ISOURC=21 (FULL PHASE SPACE OF READ)
	//"
	//"CARD SC2    (Input from ENSRC.MORTRAN)
	//"
	//" ENSRC DELIMETERS:  :start source inputs:
	//"                    :stop source inputs:
	//"
	//"
	//"  INCIDENT ENERGY
	//"        = monoenergetic  (0)  IF MONOENERGETIC BEAM
	//"        = spectrum       (1)  IF ENERGY SPECTRUM TO BE USED
	//"
	//"           ---------------------------------------
	//"
	//"       ONLY IF INCIDENT ENERGY= Monoenergetic:
	//"
	//"                   INCIDENT KINETIC ENERGY(MEV)   (I)
	//"                                   KINETIC ENERGY OF THE INCIDENT BEAM IN MeV
	//"                                   (DEFAULTS TO 1.25)
	//"
	//"           ---------------------------------------
	//"
	//"       ONLY IF INCIDENT ENERGY= Spectrum:
	//"
	//"                   SPEC FILENAME   (C)  FILENAME (WITH EXT)
	//"                                   CONTAINS SPECTRUM INFORMATION
	//"
	//"                                   FILE FORMAT:
	//"                                   TITLE      SPECTRUM TITLE  (80 char)
	//"                                   NENSRC, ENMIN, MODE
	//"                                   NENSRC     # ENERGY BINS IN SPEC. HISTOGRAM
	//"                                   ENMIN      LOWER ENERGY OF FIRST BIN
	//"                                   MODE       =0, assumes cts/bin
	//"                                              =1  assumes cts/MeV
	//"                                   ENSRCD(I),SRCPDF(I)  I=1,NENSRC
	//"                                   TOP OF ENERGY BIN AND PROBABILITY OF
	//"                                   INITIAL PARTICLE BEING IN THIS BIN.
	//"                                   PROBABILITY DOES NOT NEED TO BE NORMALIZED
	//"
	//"                   SPEC IOUTSP
	//"                        = none     (0)  NO SPECTRUM DATA IN OUTPUT SUMMARY
	//"                        = include  (1)  INCLUDE SPECTRUM DATA IN OUTPUT SUMMARY
	//;
	////"*******************************************************************************
	//"
	//"                         MC TRANSPORT PARAMETER
	//////"                         **********************
	//"
	//"  All input associated with selection of various transport parameter
	//"  is not crucial for the execution as there are default values set.
	//"  Therefore, if some of the input options in this section are
	////"  missing/misspelled, this will be ignored and defualt parameter assumed
	//"  As the transport parameter input routine uses get_inputs, a lot
	//"  of error/warning messages may be produced on UNIT 15, though.
	//"  If you don't have the intention of changing default settings,
	//"  simply ignore the error messages.
	//"
	////"  The delimeters are
	//"
	//"               :start mc transport parameter:
	//"               :stop mc transport parameter:
	//"
	////"  You can change this by including the statement
	//"
	//"  REPLACE {$THE_DELIMETER} WITH {'WHATEVER-YOU-LIKE-THE-DELIMETER-TO-BE'}
	////"
	//"  in your input file.
	////"
	//"  Currently, the following options are available (case does not matter):
	//"
	//"       Global ECUT=     Set a global (in all regions) electron transport
	//"                        cut off energy (in MeV). If this imput is missing,
	//"                        AE(medium) will be used.
	//"                        [ ECUT ]
	//"       Global PCUT=     Set a global (in all regions) photon transport
	//"                        cut off energy (in MeV). If this imput is missing,
	//"                        AP(medium) will be used.
	//"                        [ PCUT ]
	//"       Global SMAX=     Set a global (in all regions) maximum step-size
	//"                        restriction for electron transport (in cm).
	//"                        If missing, no geometrical step-size restrictions will
	//"                        be employed. Note that if you use the default
	//"                        EGSnrc electron-step algorithm, no SMAX-restriction
	//"                        is necessary. Option is useful for transport in low
	//"                        density materials (air) when PRESTA behaviour is
	//"                        turned on (see below)
	//"                        [ SMAXIR ]
	//"       ESTEPE=          Set the maximum fractional energy loss per step.
	//"                        Note that this is a global option only, no
	//"                        region-by-region setting is possible. If missing,
	//"                        the defualt is 0.25 (25%)
	//"                        [ ESTEPE ]
	//"       XImax=           Maximum first elastic scattering moment per step.
	//"                        Default is 0.5, NEVER use value greater than 1 as
	//"                        this is beyond the range of MS data available.
	//"                        [ XIMAX ]
	//"       Boundary crossing algorithm=
	//"                        There are two selections possible: EXACT, means
	//"                        the algorithm will cross boundaries in a single
	//"                        scattering (SS) mode, the distance from a boundary
	//"                        at which the transition to SS mode is made is
	//"                        determined by 'Skin depth for BCA' (see below).
	//"                        The second option is PRESTA-I, if selected boundaries
	//"                        will be crossed a la PRESTA, i.e. with lateral
	//"                        correlations turned off and MS forced at boundaries.
	//"                        Default is EXACT.
	//"                        [ bca_algorithm, exact_bca ]
	//"       Skin depth for BCA=
	//"                        Determines the distance from a boundary (in elastic
	//"                        MFP) at which the algorithm will go into single
	//"                        scattering mode (if EXACT boundary crossing) or
	//"                        swith off lateral correlations (if PRESTA-I boundary
	//"                        crossing). Default value is 3 for EXACT or
	//"                        exp(BLCMIN)/BLCMIN for PRESTA-I (see the PRESTA paper
	//"                        for a definition of BLCMIN). Note that if you choose
	//"                        EXACT boundary crossing and set Skin depth for BCA
	//"                        to a very large number (e.g. 1e10), the entire
	//"                        calculation will be in SS mode. If you choose
	//"                        PRESTA-I boundary crossing and make Skin depth for BCA
	//"                        large, you will get default EGS4 behavious (no PRESTA)
	//"                        [ skindepth_for_bca ]
	//"       Electron-step algorithm=
	//"                        PRESTA-II (the default), the name is
	//"                        used for historical reasons
	//"                        or PRESTA-I
	//"                        Determines the algorithm used to take into account
	//"                        lateral and longitudinal correlations in a
	//"                        condensed history step.
	//"                        [ transport_algorithm ]
	//"       Spin effects=    Off, On, default is On
	//"                        Turns off/on spin effects for electron elastic
	//"                        scattering. Spin On is ABSOLUTELY necessary for
	//"                        good backscattering calculations. Will make a
	//"                        even in `well conditioned' situations (e.g. depth
	//"                        dose curves for RTP energy range electrons).
	//"                        [ spin_effects ]
	//"       Brems angular sampling= Simple, KM, default is KM
	//"                        If Simple, use only the leading term of the Koch-Motz
	//"                        distribution to determine the emission angle of
	//"                        bremsstrahlung photons. If On, complete
	//"                        modified Koch-Motz 2BS is used (modifications
	//"                        concern proper handling of kinematics at low energies,
	//"                        makes 2BS almost the same as 2BN at low energies).
	//"                        [ IBRDST ]
	//"       Brems cross sections= BH, NIST, default is BH
	//"                        If BH is selected, the Bethe-Heitler bremsstrahlung
	//"                        cross sections (Coulomb corrected above 50 MeV)
	//"                        will be used. If NIST is selected, the NIST brems
	//"                        cross section data base (which is the basis for
	//"                        the ICRU radiative stopping powers) will be employed.
	//"                        Differences are negligible for E > ,say, 10 MeV,
	//"                        but signifficant in the keV energy range.
	//"       Bound Compton scattering=  On or Off
	//"                        If Off, Compton scattering will be treated with
	//"                        Klein-Nishina, with On Compton scattering is
	//"                        treated in the Impuls approximation. Default is On.
	//"                        Make sure to turn on for low energy applications,
	//"                        not necessary above, say, 1 MeV.
	//"                        [ IBCMP ]
	//"       Pair angular sampling= Off, Simple or KM
	//"                        If off, pairs are set in motion at an angle m/E
	//"                        relative to the photon direction (m is electron rest
	//"                        energy, E the photon energy). Simple turns on
	//"                        the leading term of the angular distribution
	//"                        (this is sufficient for most applications),
	//"                        KM (comes from Koch and Motz) turns on using 2BS
	//"                        from the article by Koch and Motz.
	//"                        Default is Simple, make sure you always use Simple or
	//"                        KM
	//"                        [ IPRDST ]
	//"       Photoelectron angular sampling= Off or On
	//"                        If Off, photo-electrons get the direction of the
	//"                        `mother' photon, with On, Sauter's furmula is
	//"                        used (which is, striktly speaking, valid only for
	//"                        K-shell photo-absorption).
	//"                        If the user has a better approach, replace the macro
	//"                            $SELECT-PHOTOELECTRON-DIRECTION;
	//"                        The only application that
	//"                        I encountered until now where this option made a
	//"                        small difference was a big ion chamber (cavity size
	//"                        comparable with electron range) with high-Z walls
	//"                        in a low energy photon beam.
	//"                        Default is On
	//"                        [ IPHTER ]
	//"       Rayleigh scattering= Off, On
	//"                        If On, turnd on coherent (Rayleigh) scattering.
	//"                        Default is Off. Should be turned on for low energy
	//"                        applications but needs PEGS4 data set with data.
	//"                        [ IRAYLR ]
	//"       Atomic relaxations= Off, On
	//"                        Default is On. The effect of using On is twofold:
	//"                        - In photo-electric absorption events, the element
	//"                          (if material is mixture) and the shell the photon
	//"                          is interacting with are sampled from the appropriate
	//"                          cross seections
	//"                        - Shell vacancies created in photo-absorption events
	//"                          are relaxed via emission of fluorescent X-Rays,
	//"                          Auger and Koster-Cronig electrons.
	//"                         Make sure to turn this option on for low energy
	//"                         applications.
	//"                         [ IEDGFL ]
	//"
	//"       Atomic relaxations, Rayleigh scattering,
	//"       Photoelectron angular sampling and Bound Compton scattering
	//"                         can also be turned On/Off on a region-by-region
	//"                         basis. To do so, put e.g.
	//"
	//"       Atomic relaxations= On in Regions   or
	//"       Atomic relaxations= Off in regions
	//"
	//"                         in your input file. Then use
	//"
	//"       Bound Compton start region=
	//"       Bound Compton stop region=
	//"                or
	//"       Rayleigh start region=
	//"       Rayleigh stop region=
	//"                or
	//"       Relaxations start region=
	//"       Relaxations stop region=
	//"                or
	//"       PE sampling start region=
	//"       PE sampling stop region=
	//"
	//"                         each followed by a lost of of one or more
	//"                         start and stop regions separated by commas.
	//"                         Example:
	//"        Atomic relaxations= On in Regions
	//"        Relaxations start region=  1, 40
	//"        Relaxations stop region=  10, 99
	//"                         will first turn off relaxations everywhere and
	//"                         then turn off in regions 1-10 and 40-99.
	//"                         Note that input is checked against min. and max.
	//"                         region number and ignored if
	//"                         start region < 1 or stop_region > $MXREG or
	//"                         start region > stop region.
	//"
	//"                         ECUT, PCUT and SMAX can also be set on a
	//"                         region-by-region basis. To do so, iclude
	//"                         in your input file
	//"
	//"         Set XXXX=              f_value1, f_value2, ...
	//"         Set XXXX start region= i_value1, i_value2, ...
	//"         Set XXXX stop region=  j_value1, j_value2, ...
	//"
	//"                         where XXXX is ECUT, PCUT or SMAX ,
	//"                         f_value1, f_value2,... are the desired values for XXXX
	//"                         and i_value_i and j_value_i are the start and
	//"                         stop regions.
	//"
	//"*******************************************************************************
	//"
	//"                    VARIANCE REDUCTION
	//"                    ******************
	//"
	//"  Delimeter:       :start variance reduction:
	//"                   :stop variance reduction:
	//"
	//"  ELECTRON RANGE REJECTION
	//"         = off        (0)  No electron range rejection
	//"         = on         (1)  Do electron range rejection.
	//"                           There are 2 components to range rejection.
	//"                           One uses the EGSnrc range rejection below ESAVEIN
	//"                           and terminates any charged particle which cannot get
	//"                           out of its local region.
	//"                           The second component terminates any charged particle
	//"                           which cannot reach the cylinder which encloses the
	//"                           cavity region and any other region of the same
	//"                           material as the cavity.  This cylinder is determined
	//"                           automatically.
	//"                           The parameter ESAVEIN also plays a role (see below)
	//"                           [IREJCT]
	//"
	//"  ESAVEIN             (R)  If ELECTRON RANGE REJECTION is on, discard an
	//"                           electron  when E< ESAVEIN and RANGE < CDIST
	//"                           where CDIST is closest distance to region of
	//"                           interest specified below. This ignores brem
	//"                           losses below ESAVEIN.
	//"                           This parameter must be input even if not used.
	//"                           Note - ESAVEIN is total energy (with 511 keV)
	//"
	//"  RUSSIAN ROULETTE DEPTH      (R)
	//"                           FOR RUSSIAN ROULETTE -
	//"                           AS ANY PHOTON CROSSES THE Z='RUSSIAN ROULETTE DEPTH'
	//"                           PLANE,  RUSSIAN ROULETTE IS PLAYED.
	//"
	//"  RUSSIAN ROULETTE FRACTION   (R)
	//"                           EACH TIME RUSSIAN ROULETTE IS PLAYED, RRF IS THE
	//"                           PROBABILITY OF SURVIVAL.
	//"                           WEIGHT INCREASES BY 1/RRF,  IF IT SURVIVES
	//"
	//"                    ****** IF BOTH ZERO, NO RUSSIAN ROULETTE IS PLAYED ******
	//"
	//"  EXPONENTIAL TRANSFORM C     (R)
	//"                           PARAMETER FOR PATHLENGTH BIASING <0 FOR SHORTENING
	//"                           IF 0.0, NO BIASING DONE
	;
	//"
	//"  PHOTON FORCING
	//"        = Off         (0)    NORMAL PHOTON TRANSPORT (NO FORCING)
	//"        = On          (1)    FORCE PHOTON INTERACTIONS EXPLICITLY
	//"                             MUST SET START AND STOP FORCING IN THIS CASE
	//"
	//"  START FORCING       (I)    NUMBER OF PHOTON INTERACTION/HISTORY AT WHICH
	//"                             TO START FORCING PHOTON INTERACTIONS
	//"
	//"  STOP FORCING AFTER  (I)    NUMBER OF PHOTON INTERACTION/HISTORY AFTER WHICH
	//"                             TO STOP FORCING PHOTON INTERACTIONS
	//"
	//"                             STOP FORCING AFTER > OR = START FORCING
	//"
	//" PHOTON SPLITTING=    (I)    Number of times to split a photon
	//"                             If missing or < 2  => normal transport
	//"                             If >= 2 (allowed only for ifull=0,1), the macro
	//"                             $SELECT-PHOTON-MFP essentially replaces the
	//"                             entire PHOTON routine.
	//"                             This option increses the efficiency of ifull=0
	//"                             caluclations by up to a factor of 3 compared
	//"                             to simple photon forcing (forcing is ignored
	//"                             if n_split > 1). A rule of thumb for good eff. is
	//"                             n_split >= No/(1-exp(-Lambda))
	//"                             where Lambda is approx. number of photon MFP
	//"                             in the geometry of interest and No >= 5.
	//"                             Note that if you use the above, there will
	//"                             be on average approx. No primary interactions per
	//"                             incident photon => reduce the number of
	//"                             histories by this number.
	//"
	//"        The algorithm works as follows:
	//"         * dpmfp_i = -log(1 - (eta+i)/n_split)
	//"           where dpmfp_i is MFP to the next interaction for the i-th sub-photon
	//"           eta is a random number (the same for all n_split sub-photons)
	//"         * Once at the interaction site, the i'th sub-photon produces
	//"           electrons and/or scattered photons. Scattered photons are
	//"           killed with probability 1/n_split, so that, if they survive,
	//"           they have the weight of the original photon. Electrons have the
	//"           weight of 1/n_split of original weight.
	//"         * In any radiative events (brems, annih, annih at rest), photons
	//"           are killed with probability 1/n_split => they have again the
	//"           weight of the photon that initiated the history, if they survive
	//"         * If ifano = 1 (regeneration), all scattered photons are killed
	//"           and split-photons re-generated with probability 1/n_split.
	//"         * If the user wants to calculate Awall in addition to the dose
	//"           (ifull=1), at all interaction sites the primary photon is
	//"           regenerated with probability 1/n_split and marked as such
	//"           (latch=2). It's descending electrons then only score to the
	//"           dose with attenuation and scatter removed. The ratio of the
	//"           real dose to the dose with attenuation and scatter removed
	//"           is per definition Awall. Due to the strong correlation between
	//"           the two doses the uncertainty on Awall is much smaller than
	//"           the uncertainty on the each individual dose.
	//"
	//"   BE AWARE: When this option is turned on, the normal scoring arrays
	//"             in cavrznrc are ignored and scoring is done on a history-
	//"             by-history basis (instead of batch scoring) for the whole
	//"             cavity only => detailed output for each cavity region
	//"             is not done even if requested.

	//"   BE AWARE: If you select too a big splitting number, a stack overflow
	//"             may result. In such cases either use a smaller n_split
	//"             or increase $MXSTACK (currently 2000)
	//"
	//"        Implemented by I.Kawrakow, February 2000.
	//"
	//"-----------------------------------------------------------------------------
	//"
	//" CS ENHANCEMENT FACTOR= (R) Cross section enhancement factor.
	//"                            If this input is missing or set to <= 1, it
	//"                            has no effect on the simulation. But if
	//"                            the enhancement factor is set to > 1, the
	//"                            effect is dramatic: all other user input
	//"                            concerning photon forcing, splitting, exp.
	//"                            transform, etc., is ignored. In addition,
	//"                            the calculation result corresponds  ALWAYS
	//"                            to 'Aatt and Ascat', no matter what the
	//"                            user requested (but only Awall is calculated,
	//"                            not the individual Ascat and Aatt).
	//"                            The algorithm employed is implemented via
	//"                            $RAYLEIGH-CORRECTION and appropriate calls to
	//"                            AUSGAB and works as follows:
	//"       - the photon cross section is multiplied with the enhancement
	//"         factor (say, C).
	//"       - when a photon arrives at an interaction site, it is split into
	//"         an interacting part (fraction 1/C) and a non-interacting part
	//"         (fraction  1-1/C).
	//"       - Electrons set in motion by the interacting fraction are kept
	//"       - Scattered photons and the non-interacting portion of the
	//"         original photon are killed with probabilities 1/C and 1-1/C,
	//"         respectively.
	//"       - Because the goal is to calculate Awall in addition to the real
	//"         cavity dose, killed non-interacting fractions of primary
	//"         photons are left on
	//"         the stack but marked as such so that, when they are transported,
	//"         all their descending alectrons contribute only to the dose
	//"         with attenuation and scatter removed.
	//"       - Scattered photons descending from a 'killed' non-interacting
	//"         fractions are always removed from the stack
	//"       - Because no forcing is done, it is a good idea to set the
	//"         enhancement factor to a big enough number so that
	//"         photons interact several times in the chamber. For instance,
	//"         the attenuation coefficient of a Co60 beam in graphite is ~0.1,
	//"         if the chamber wall thickness is, say, 1 cm, the average
	//"         number of MFP is 0.1. To have on average 5 interactions per
	//"         incident photon one needs C=50 (0.1*50 = 5)
	//"
	//"   BE AWARE: If you select to a big enhancement factor, a stack overflow
	//"             may result. In such cases either use a smaller enhancement
	//"             factor or increase $MXSTACK (currently 2000)
	//"
	//"        Implemented by I.Kawrakow, April 2001.
	//"
	//"******************************************************************************

	public static int $NSWTCH=8;//       "# OF NRC SWITCHES FOR CONTROLLING SCATTERING "
	public static int $NBATCH=10;//        "OUTPUT BATCHES                             "
	public static int JCASE=0;; //"no. of histories per batch"
	public static int $NCASEMIN=100;//   "min. no. of histories                        "
	public static int $MXDATA=1040;//  "MAXIMUM DATA POINTS FOR ANALYSIS (i.e.($MXREG-1))"
	public static int $MAXIT=4;//        "MAX # OF PARAMETERS TO BE SCORED             "
	//"                                (1) PRIMARY EDEP TO GAS                       "
	//"                                (2) SECONDARY EDEP TO GAS                     "
	//"                                (3) UNATTENUATED PRIMARY EDEP TO GAS          "
	//"                                     - PRIMARY EDEP TO GAS                    "
	//"                                (4) (R**2-R0**2)(UNATTENUATED PRIMARY         "
	//"                                     EDEP TO GAS)                             "
	//"                                (5) UNSOURCED UNATTENUATED PRIMARY EDEP TO    "
	//"                                     GAS WITH GAS MATERIAL REPLACED BY WALL   "
	//"                                     MATERIAL                                 "
	//"                                (6) UNSOURCED UNATTENUATED PRIMARY EDEP TO    "
	//"                                     WALL MATERIAL WITH GAS MATERIAL REPLACED "
	//"                                     BY WALL MATERIAL                         "
	public static int $MAXCMPTS=$MAXIT;//  "FOR THE GRID OUTPUTS"
	public static int $MAX_SC_PLANES=1;//"required to use phase space macros"

	public static int NCASE=0;
	public static int NCASEO=0;
	public static int NCASET=0;
	public static double[][] AMASS;//($MAXZREG,$MAXRADII),
	public static double TMCPUO=0.0;
	public static double TIMMAX=0.0;
	public static double STATLM=0.0;//EIN,
	public static int[] MEDSAV;//($MXREG),
	public static double[] RHOSAV;//($MXREG),
	public static int IDAT=0;
	public static int IDOPES=0;
	public static int IRESTART=0;
	public static int NNREADO=0;
	public static int datcount=0;
	//"AMASS(IZ,IX)  MASS OF ZONE WITH COORDINATES (IZ,IX)
	//"TMCPUO        CPU TIME USED IN PREVIOUS SESSIONS
	//"TIMMAX        MAXIMUM ALLOWED CPU HOURS FOR A GIVEN CALCULATION
	//"STATLM        TARGET STATISTICS IN CAVITY USED FOR AN EARLY EXIT
	//"EIN           KINETIC ENERGY OF THE EXTERNAL BEAM
	//"ISUMCV(NREG)  THE ARRAY OF ZONES COMPRISING THE CAVITY REGION
	//"MEDSAV(NREG)  SAVES MEDIUM NUMBERS FOR CORRELATION SCORING
	//"RHOSAV(NREG)  SAVES DENSITIES  FOR CORRELATION SCORING
	//"IDAT          = 0 STORE DATA ARRAYS FOR RE-USE
	//"              = 1 DON'T STORE THEM
	//"IDOPES        = 1 INCLUDES PHOTOELECTRON ANGLE SELECTION
	//"              = 0 DOES NOT INCLUDE PHOTOELECTRON ANGLE SELECTION
	//"NCASE         NUMBER OF HISTORIES REMAINING TO BE DONE
	//"NCASEO        NUMBER OF HISTORIES DONE IN PREVIOUS SESSIONS
	//"NCASET        NUMBER OF HISTORIES ALREADY DONE
	//"NNREADO       TOTAL NO. OF PARTICLES READ FROM PHSP SOURCE IN PREVIOUS
	//"              RUNS
	//"IRESTART        = 0 => INITIAL RUN
	//"              = 1 => RESTARTED RUN
	//"              = 3 => DATA ANALYSIS ONLY
	//"              = 4 => READ STARTING RANDOM NUMBERS FROM A FILE
	//"              = 5 => analyse previous parallel runs

	public static double RRZ=0.0;//
	public static double RRCUT=0.0;
	public static boolean RUSROU=false;
	//"RRZ      COORDINATE OF PLANE AT WHICH RUSSIAN ROULETTE IS PLAYED
	//"RRCUT    SURVIVAL PROBABILITY AFTER CROSSING THE PLANE
	//"RUSROU   = .FALSE. => RUSSIAN ROULETTE WILL NOT BE PLAYED
	//"         = .TRUE.  => RUSSIAN ROULETTE WILL BE PLAYED
	public static double cav_dose=0.0;//     => total dose in cavity
	public static double cav_dose0=0.0;//    => primary dose in cavity
	public static double cav_dose1=0.0;//    => primary dose corrected for attenuation in cavity
	public static double cav_dose2=0.0;//    => secondary dose in cavity
	public static double cav2_dose=0.0;//    => total dose squared in cavity
	public static double cav2_dose0=0.0;//   => primary dose squared in cavity
	public static double cav2_dose1=0.0;//   => unattenuated primary dose squared in cavity
	public static double cav2_dose2=0.0;//   => secondary dose squared in cavity
	//"             => eventually, they hold the uncertainty in their respective doses
	public static double cav_dosec=0.0;//    => correlation total - primary(unatttenuated) in cavity
	public static double cav_dosec01=0.0;//  => correlation primary - primary (unattenuated) in cavity
	public static double cav_dosec02=0.0;//  => correlation secondary - primary in cavity
	public static double SCSTP=0.0;//        => total no. of charged particle steps
	public static double SCSTP2=0.0;//       => total no. of charged particle steps squared--eventually holds
	//"                uncertainty in SCSTP
	public static double SCCSTP=0.0;//       => no. of charged particle steps in cavity
	public static double SCCSTP2=0.0;//      => no. of charged particle steps in cavity squared--eventually
	//"                holds uncertainty in SCCSTP
	public static double[][][] SCDOSE;//(IZ,IX,IT)      => dose in cavity voxel IZ,IX:
	//"                        IT=1 -- total
	//"                        IT=2 -- primary
	//"                        IT=3 -- primary corrected for attenuation
	//"                        IT=4 -- secondary
	public static double[][][] SCDOSE2;//(IZ,IX,IT)     => dose in cavity voxel IZ,IX squared.  IT same as above.
	//"                        Eventually holds uncertainties in respective doses.
	public static double[][][] SCDOSE_COV;//(IZ,IX,IT)  => correlation in cavity voxel IZ,IX for:
	//"                        IT=1 -- total and primary (unattenuated)
	//"                        IT=2 -- primary and primary (unattenuated)
	//"                        IT=3 -- secondary and primary

	public static double PIISTP=0.0;//                => no. of PRESTA-II steps from previous runs
	//"last_case             => last primary history to score dose in cavity
	public static int SCSTP_LAST=0;//            => last primary history to score charged particle step
	public static int SCCSTP_LAST=0;//           => last primary history to score charged particle step
	//"                         in cavity
	public static int[][] SCDOSE_LAST;//(IZ,IX)    => last primary history to score dose in cavity
	//"                         voxel IZ,IX.
	public static double tmp_dose=0.0;//              => temporay arrays for scoring the different
	public static double tmp_dose0=0.0;//                dose components in the cavity
	public static double tmp_dose1=0.0;//
	public static double tmp_dose2=0.0;//
	public static double corr_02=0.0;//               => correlation primary - secondary cavity doses
	public static double SCSTP_TMP=0.0;//             => temp. variable for scoring total no. of charged
	//"                         particle steps
	public static double SCCSTP_TMP=0.0;//            => temp. variable for scoring total no. of charged
	//"                         particle steps in cavity
	public static double[][][] SCDOSE_TMP;//(IZ,IX,IT)  => temp. variable for scoring dose in cavity voxel
	//"                        IZ,IX.  IT same as for SCDOSE.
	//"cs_enhance            => cross-section enhancement factor

	public static int MXNP=0;//    MAXIMUM LEVEL TO WHICH THE STACK OF DAUGHTER PARTICLES FROM AN
	//"        INCIDENT PARTICLE RISES (STACK MAY INCLUDE INCIDENT PARTICLE)
	//"IFULL   = 0 JUST CALCULATE TOTAL DOSE AND THAT DUE TO STOPPERS
	//"            AND DISCARDS (THE DEFAULT)
	//"        = 1 ABOVE PLUS Aatt, Ascat
	public static int ISTORE=0;//  = 0 DO NOT STORE THE INITIAL RANDOM NUMBERS (THE DEFAULT)
	//"        = 1 STORE THE INITIAL RANDOM NUMBER FOR THE LAST HISTORY
	//"        = 2 STORE THE INITIAL RANDOM NUMBER FOR ALL HISTORIES
	//"            THAT DEPOSIT ENERGY IN THE CAVITY
	//"        = 3 STORE ALL THE INITIAL RANDOM NUMBERS
	public static int IWATCH=0;//  = 0 FOR NORMAL OUTPUT (THE DEFAULT)
	//"        = 1 OUTPUT ON EVERY DISCRETE INTERACTION
	//"        = 2 OUTPUT ON EVERY ELECTRON/PHOTON STEP AS WELL
	//"        = 3 PRINTS OUT ONLY WHEN ENERGY IS DEPOSITED
	//"        = 4 PRINTS OUT FILE FOR GRAPHICS
	public static int IOOPTN=0;//  = 0 SHORT OUTPUT (THE DEFAULT) -JUST CAVITY SUMMARY
	//"            AND THE MATERIAL GRID
	//"        = 1 ABOVE PLUS OUTPUT GRID
	public static int IOUTSP=0;//  = 0 NO SPECTRUM DATA IN OUTPUT SUMMARY
	//"        = 1 INCLUDE SPECTRUM DATA IN OUTPUT SUMMARY
	public static int IFANO=0;//   = 0 NO PHOTON REGENERATION
	//"        = 1 PHOTONS REGENERATED AFTER THEY HAVE INTERACTED
	//"        = 2 NO PHOTON REGENERATION, ELECTRONS FROM CAVITY WALL ARE ELIMINATED

	public static final String NULLs="";
//-------------INPUTS-------------------------------------------
	public static String TITLEs="";
	//@ IWATCH
	public static final int IWATCH_OFF=0;
	public static final int IWATCH_INTERACTIONS=1;
	public static final int IWATCH_STEPS=2;
	public static final int IWATCH_DEPOSITED=3;
	public static final int IWATCH_GRAPH=4;
	//@ STORE INITIAL RANDOM NUMBERS
	public static final int ISTORE_NO=0;
	public static final int ISTORE_LAST=1;
	public static final int ISTORE_ALL_DEPOSITED=2;
	public static final int ISTORE_ALL=3;
	//@ IRESTART
	public static final int IRESTART_FIRST=0;
	public static final int IRESTART_RESTART=1;
	public static final int IRESTART_ANALYZE=2;
	public static final int IRESTART_START_RNS=3;
	public static final int IRESTART_PARALLEL=4;
	//@ OUTPUT OPTIONS
	public static final int IOOPTN_SHORT=0;
	public static final int IOOPTN_DETAIL=1;
	//@ STORE DATA ARRAYS
	public static final int IDAT_YES=0;
	public static final int IDAT_NO=1;
	//@ NUMBER OF HISTORIES
	public static final int NCASE_MIN=1;
	public static final int NCASE_MAX=461168600;//4.611686e18;
	public static final int NCASE_DEFAULT=20000;
	//@ MAX CPU HOURS ALLOWED
	public static final double TIMMAX_MIN=0.0;
	public static final double TIMMAX_MAX=1000.0;
	public static final double TIMMAX_DEFAULT=999.0;
	//@ IFULL
	public static final int IFULL_DOSE_AND_STOPPERS=0;
	public static final int IFULL_AATT_AND_ASCAT=1;
	//@ STATISTICAL ACCURACY SOUGHT
	public static final double STATLM_MIN=0.0;
	public static final double STATLM_MAX=100.0;
	public static final double STATLM_DEFAULT=0.0;
	//@ PHOTON REGENERATION
	public static final int IFANO_NO=0;
	public static final int IFANO_YES=1;
	public static final int IFANO_NO_ELECTRONS_FROM_WALL=2;
	//@ INITIAL RANDOM NO. SEEDS
	public static final int RANLUX_LEVEL_MIN=0;
	public static final int RANLUX_LEVEL_MAX=4;
	public static final int RANLUX_LEVEL_DEFAULT=EGS4.$DEFAULT_LL;//see egs4
	public static final int RANLUX_SEED_MIN=1;
	public static final int RANLUX_SEED_MAX=1073741824;//default seed is set in egs4.setdefaults!!
	public static final int RANLUX_SEED_DEFAULT=999999;
	public static final int RANMAR_SEED_MIN=1;
	public static final int RANMAR_SEED_MAX=30081;
	public static final int RANMAR_SEED_DEFAULT=9373;
	public static int jrng1=0;public static int jrng2=0;
	//@TRANSPORT PARAMS:
	public static double ecut=0.0;
	public static double pcut=0.0;
	public static double smax=0.0;
	public static boolean setEcutRegion=false;
    public static boolean setPcutRegion=false;
    public static boolean setSmaxRegion=false;
	public static final double ECUT_MIN=0.0;//" ECUT "
	public static final double ECUT_MAX=1.E15;
	public static final double PCUT_MIN=0.0;//" PCUT "
	public static final double PCUT_MAX=1.E15;
	public static final double SMAXIR_MIN=0.0;//" SMAX "
	public static final double SMAXIR_MAX=1.E15;
	public static int nEcut=1;//number of data
	public static double[] Ecut;//=new double[EGS4.$MXREG];
	public static int[] startEcutRegion;
	public static int[] stopEcutRegion;
	public static int nPcut=1;//number of data
	public static double[] Pcut;//=new double[EGS4.$MXREG];
	public static int[] startPcutRegion;
	public static int[] stopPcutRegion;
	public static int nSmax=1;//number of data
	public static double[] Smax;//=new double[EGS4.$MXREG];
	public static int[] startSmaxRegion;
	public static int[] stopSmaxRegion;

	public static boolean setIncohRegion=false;
	public static boolean setCohRegion=false;
	public static boolean setRelaxRegion=false;
	public static boolean setPeRegion=false;
	public static int nIncoh=1;//number of data
	public static int[] Incoh;//=new int[EGS4.$MXREG];
	public static int[] startIncohRegion;
	public static int[] stopIncohRegion;
	public static int nCoh=1;//number of data
	public static int[] Coh;//=new int[EGS4.$MXREG];
	public static int[] startCohRegion;
	public static int[] stopCohRegion;
	public static int nRelax=1;//number of data
	public static int[] Relax;//=new int[EGS4.$MXREG];
	public static int[] startRelaxRegion;
	public static int[] stopRelaxRegion;
	public static int nPe=1;//number of data
	public static int[] Pe;//=new int[EGS4.$MXREG];
	public static int[] startPeRegion;
	public static int[] stopPeRegion;

		//" Incoherent (Compton) scattering "
	public static int incoh=0;
	public static final int incoh_OFF=0;
	public static final int incoh_ON=1;
	public static final int incoh_ON_IN_REGIONS=2;
	public static final int incoh_OFF_IN_REGIONS=3;
		//" Radiative corrections for Compton scattering "
	public static final int radc_OFF=0;
	public static final int radc_ON=1;
		//" Coherent (Rayleigh) scattering "
	public static int coh=0;
	public static final int coh_OFF=0;
	public static final int coh_ON=1;
	public static final int coh_ON_IN_REGIONS=2;
	public static final int coh_OFF_IN_REGIONS=3;
		//" Atomic Relaxations "
	public static int relax=0;
	public static final int relax_OFF=0;
	public static final int relax_ON=1;
	public static final int relax_ON_IN_REGIONS=2;
	public static final int relax_OFF_IN_REGIONS=3;
		//" Photoelectron angular sampling "
	public static int pe=0;
	public static final int pe_ang_OFF=0;
	public static final int pe_ang_ON=1;
	public static final int pe_ang_ON_IN_REGIONS=2;
	public static final int pe_ang_OFF_IN_REGIONS=3;
		//" Bremsstrahlung angular sampling "
	public static final int brems_ang_SIMPLE=0;
	public static final int brems_ang_KM=1;
		//" Bremsstrahlung cross sections "
	public static final int brems_cross_BH=0;
	public static final int brems_cross_NIST=1;
		//" Pair angular sampling "
	public static final int pair_ang_OFF=0;
	public static final int pair_ang_SIMPLE=1;
	public static final int pair_ang_KM=2;
	public static final int pair_ang_UNIFORM=3;
	public static final int pair_ang_BLEND=4;
		//" Pair cross sections "
	public static final int pair_cross_BH=0;
	public static final int pair_cross_NRC=1;
		//" Triplet production "
	public static final int triplet_OFF=0;
	public static final int triplet_ON=1;
		//" Spin effects  		"
	public static int ispin=0;
	public static final int spin_OFF=0;
	public static final int spin_ON=1;
		//" Electron impact ionization "
	public static final int eii_OFF=0;
	public static final int eii_ON=1;

	public static final double ESTEPE_MIN=1.0E-5;//" ESTEPE "
	public static final double ESTEPE_MAX=1.0;
	public static final double XIMAX_MIN=0.0;//" XIMAX "
	public static final double XIMAX_MAX=1.0;
		//" BCA "
	public static final int BCA_EXACT=0;
	public static final int BCA_PRESTA_I=1;

	public static final double Skindepth_MIN=-1.0;//" Skindepth "
	public static final double Skindepth_MAX=1.0e15;
		//" Electron-step algorithm "
	public static final int estep_alg_PRESTA_I=1;
	public static final int estep_alg_PRESTA_II=0;

	//Variance reduction
	public static final int irejct_OFF=0;
	public static final int irejct_ON=1;
	public static double ESAVEIN=0.0;
	public static final double ESAVEIN_MIN=0.0;
	public static final double ESAVEIN_MAX=1.0e30;
	public static final double ESAVEIN_DEFAULT=0.0;
	public static final double cs_enhance_MIN=0.0;
	public static final double cs_enhance_MAX=1.0e6;
	public static final double cs_enhance_DEFAULT=0.5;
	public static final double RRDEPTH_MIN=-1.0e30;
	public static final double RRDEPTH_MAX=1.0e30;
	public static final double RRDEPTH_DEFAULT=0.0;
	public static final double RRFRACTION_MIN=-1.0e30;
	public static final double RRFRACTION_MAX=1.0e30;
	public static final double RRFRACTION_DEFAULT=0.0;
	public static final double EXPC_MIN=-1.0e30;
	public static final double EXPC_MAX=1.0e30;
	public static final double EXPC_DEFAULT=0.0;
	public static int IFARCE=0;
	public static final int IFARCE_ON=1;
	public static final int IFARCE_OFF=0;
	public static final int NFMIN_MIN=0;
	public static final int NFMIN_DEFAULT=1;
	public static int NFMIN_MAX=0;
	public static final int NFMAX_MIN=0;
	public static final int NFMAX_DEFAULT=1;
	public static int NFMAX_MAX=0;
	public static int phsplitt=1;
	public static final int phsplitt_MIN=1;
	public static final int phsplitt_DEFAULT=1;
	public static int phsplitt_MAX=0;

	//###########################################
	public static int IHSTRY=0;//=> counter for total no. of histories
	public static double EI=0.0;
	public static double EKMAX=0.0;
	public static double DEPTH=0.0;
	public static double VOLUME=0.0;
	public static double RLOW2=0.0;
	public static int ISUMX=0;

	public static double FMASSC=0.0;
	public static double FMASS=0.0;
	public static double AINFLU_CURRENT=0.0;	//-------
	public static int last_case=0;
	public static int n_photons = 0;
	public static int n_electrons = 0;
	public static double sumE_photons= 0.0;
	public static double sumE2_photons= 0.0;
	public static double sumE_electrons= 0.0;
	public static double sumE2_electrons = 0.0;
	public static int IBTCH=0;
	public static double dcav_current=0.0;
	public static double dcavun=0.0;
	public static int IDECAV=0;
	public static int ipass=0;//not usded
	public static double SCORE_NORM_NUM=0.0;
	public static double SCORE_TEMP=0.0;

	public static int IRL=0;
	public static int MEDNUM=0;

	//========================
	public static boolean createOutputFile=false;
	public static boolean putInFile=false;//internal var defining when and what to print
	private String filename="";
	FileWriter sigfos;
//-------------end INPUTS-------------------------------------------
	//REPLACE {$INVDIM} WITH {1000}    "DIMENSION CONTROLS GRID SIZE FOR INVERSE"
	//CDFINV($INVDIM,2)

	public static boolean is_finished = false;

	/**
	 * Constructor.
	 */
	public cavrz()
	{
		//createOutputFile=false;
		createOutputFile=true;
		putInFile=true;
		Calendar cal=Calendar.getInstance();
		String fs=cal.get(Calendar.YEAR)+"_"+cal.get(Calendar.MONTH)+"_"+cal.get(Calendar.DAY_OF_MONTH)+
		"_"+cal.get(Calendar.HOUR)+"_"+cal.get(Calendar.MINUTE)+"_"+cal.get(Calendar.SECOND)+".txt";
		filename=fs;//will be 2005_10_25_14_30_56.txt

		init();
	}

	/**
	 * Perform basic initialization and RUN the Monte Carlo engine.
	 */
	private void init()
	{
		if(createOutputFile)
		{
			try
			{
				sigfos= new FileWriter(filename);
			}
			catch(Exception ex){}
		}

		EGS4.startSimulationTime=System.currentTimeMillis();

		EGS4.setMXMED(5);//"MAX # OF MEDIA                               "
		EGS4.setMXREG(50);//"#REGIONS, $MAXRADII*($MAXZPLANE-1)+1(VAC)    "
		EGS4.setMXSTACK(4000);//"NEED HUGE STACK FOR CORRELATIONS+splitting   "
		EGS4.setMXRANGE(500);  //"for range arrays used in range_rejection()"
		//--variable init
		EGS4Geom.$MAXZREG=20;//"MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4Geom.$MAXRADII=8;//"MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4SrcEns.$MXRDIST=1000;
		EGS4SrcEns.$NENSRC=300;
		EGS4Geom.$NVALUE=100;
		AMASS=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];//($MAXZREG,$MAXRADII),
		MEDSAV=new int[EGS4.$MXREG];;//($MXREG),
		RHOSAV=new double[EGS4.$MXREG];;//($MXREG),
		Ecut=new double[EGS4.$MXREG];
		startEcutRegion=new int[EGS4.$MXREG];
		stopEcutRegion=new int[EGS4.$MXREG];
		Pcut=new double[EGS4.$MXREG];
		startPcutRegion=new int[EGS4.$MXREG];
		stopPcutRegion=new int[EGS4.$MXREG];
		Smax=new double[EGS4.$MXREG];
		startSmaxRegion=new int[EGS4.$MXREG];
		stopSmaxRegion=new int[EGS4.$MXREG];
		Incoh=new int[EGS4.$MXREG];
		startIncohRegion=new int[EGS4.$MXREG];
		stopIncohRegion=new int[EGS4.$MXREG];
		Coh=new int[EGS4.$MXREG];
		startCohRegion=new int[EGS4.$MXREG];
		stopCohRegion=new int[EGS4.$MXREG];
		Relax=new int[EGS4.$MXREG];
		startRelaxRegion=new int[EGS4.$MXREG];
		stopRelaxRegion=new int[EGS4.$MXREG];
		Pe=new int[EGS4.$MXREG];
		startPeRegion=new int[EGS4.$MXREG];
		stopPeRegion=new int[EGS4.$MXREG];
		//---------LOCAL VAR---------------
		SCDOSE=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE2=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
       	SCDOSE_COV=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][3];
       	SCDOSE_LAST=new int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];
       	SCDOSE_TMP=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
       	//($MAXZREG,$MAXRADII,$MAXIT)
       	EGS4Grid.CDSTBL=new String[EGS4.$MXREG];
       	EGS4Grid.CTRTBL=new String[EGS4.$MXREG];
       	EGS4Grid.CABSRB=new String[EGS4.$MXREG];
       	EGS4Grid.CAVTRACK=new String[EGS4.$MXREG];
		//--interface->LINK
		EGS4.eq=this;//pass the printing mode
		EGS4Core.eq=this;//pass the printing mode
		EGS4Macro.eq=this;//pass the printing mode
		EGS4SrcEns.eq=this;//pass the printing mode
		EGS4Geom.eq=this;//pass the printing mode
		EGS4Grid.eq=this;//pass the printing mode
		//--Macro and defaults param:
		EGS4.egs_set_defaults();//first default then:
		EGS4.RandomUse=1;//use EGS ranlux or ranmar generators
		EGS4.ranluxB=false;//use ranmar!!
		EGS4.ispmfp=EGS4.iCavity;//select photon mfp
		EGS4.iurd=EGS4.iCavity;//user range discard
		EGS4.iraycorr=EGS4.iCavity;//rayleigh correction
		EGS4.isemfp=EGS4.iCavity;//$SELECT-ELECTRON-MFP
		EGS4.iGeom=EGS4.iCavity;//1=RZ GEOM
		EGS4Macro.ismfpfpb=EGS4.iCavity;//select mfp parallel beam
		EGS4Macro.irange_rej=EGS4.iCavity;//range rejection

		EGS4Grid.$MAXCMPTS=$MAXCMPTS;
			//-----------------
		EGS4.USER_CONTROLS_TSTEP_RECURSION=EGS4.iCavity;//no effect here
		EGS4.hatchindex=EGS4.iCavity;//no effect here
			//-----------------
		EGS4.iprint=2;//summary
		//inputs
		is_finished = true;

		//"******************************************************************************
		//"
		//"                       *** SECTION 1 ***
		//"
		//"------------------------------------------------------------------------------
		//"
		//"READ INPUTS AND CALCULATE ONE-TIME ONLY CONSTANTS
		//"
		//"------------------------------------------------------------------------------

		Calendar cal=Calendar.getInstance();
		Date d=cal.getTime();
		EGS4.seqStr="=================================================================================";
      	EGS4.seqStr=" ********************CAVITY: DEMO APPLICATION************************************";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY.";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" IT IS INTENDED FOR USE IN CALCULATING QUANTITIES OF INTEREST FOR THICK-WALLED ION CHAMBERS ";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" EXPOSED TO PHOTON BEAMS. IT ALSO MAY BE USED SIMPLY TO SCORE DOSE IN A CYLINDRICAL GEOMETRY.";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" 					"+d.toString();
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" ********************************************************************************";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
		inputs();if(EGS4.STOPPROGRAM){return;}

		//"OPEN A FILE FOR STORING OR READING RANDOM NUMBERS"
		if( ISTORE > 0 )//NOT ALLOWED
		{// "We want to store the rng state in a file"
		 //   rng_unit = egs_open_file(2,0,1,'.egsrns');
		}
		/*
		ELSE IF( irestart = 4 )//NOT ALLOWED
		[
		    rng_unit = egs_open_datfile(2,0,1,'.egsrns');
		]*/
        /*
		IF(IRESTART.EQ.5)//NOT ALLOWED
		[

		    call egs_combine_runs(combine_results,'.egsdat');

		    NBATCH=0;      "DON'T WANT IT TO RUN ANY HISTORIES"
		    NCASET=NCASEO; "To prevent a wrong normalization if some of the "
		                   "parallel runs not available, IK, Jan 21 1999"
		] "end of IRESTART = 5, DISTRIBUTED POST-PROCESSING"
		ELSE [
		*/
		    if(NCASE/$NBATCH==0){NCASE=$NBATCH;}
		    JCASE=NCASE/$NBATCH; NCASE=JCASE*$NBATCH;//"NUMBER OF HISTORIES PER BATCH
		//]
		MXNP=0; //"reset the maximum stack indicator"
		IHSTRY=NCASEO; //"reset the number of histories counter"

		//"set up the broad parallel beam defaults"
		if (EGS4SrcEns.ISOURC == 2)
		{
			EGS4Geom.NR=1;EGS4Geom.RCYL[1]=1000.;
			EGS4Geom.nreg=EGS4Geom.NZ+1;
			EGS4Geom.CYRAD2[0]=EGS4Geom.RCYL[1]*EGS4Geom.RCYL[1];
		}
		//"set up ausgab calls"
		for(int J=1;J<=5;J++)
		{
			EGS4.iausfl[J-1]=1;//IAUSFL(J)=1;
		}
		for(int J=6;J<=25;J++)
		{
			EGS4.iausfl[J-1]=0;
		} //"NORMAL EXECUTION"

		if(EGS4Macro.IFULL==1)
		{
		    //"these flags are the minimum set needed to identify primary and secondary"
		    EGS4.iausfl[7] =1; //"After BREMSSTRAHLUNG"
		    EGS4.iausfl[13]=1; //"After ANNIHILATION IN FLIGHT"
		    EGS4.iausfl[14]=1; //"After ANNIHILATION AT REST"
		    EGS4.iausfl[18]=1; //"After COMPTON"
		    EGS4.iausfl[20]=1; //"After Photo"
		    EGS4.iausfl[24]=1; //"After Rayleigh"
		}
		if(IFANO == 1)
		{
		    //"AUSGAB will be responsible for making sure that the beam is not"
		    //"attenuated and getting rid of the scattered photons."
		    EGS4.iausfl[15] = 1; //"Before pair"
		    EGS4.iausfl[17] = 1; //"Before Compton"
		    EGS4.iausfl[18] = 1; //"After Compton"
		    EGS4.iausfl[19] = 1; //"Before photoelectric"
		    EGS4.iausfl[20] = 1; //"After photoelectric"
		    EGS4.iausfl[23] = 1; //"Before Rayleigh"
		    EGS4.iausfl[24] = 1; //"After Rayleigh"

		    //"AUSGAB will be responsible for throwing away any photons resulting"
		    //"from a primary electron. ie. True equilibtrium requires that all"
		    //"energy deposition be local to the primary interaction site."
		    EGS4.iausfl[7]  = 1; //"After bremsstrahlung"
		    EGS4.iausfl[13] = 1; //"A positron has annihilated in-flight"
		    EGS4.iausfl[14] = 1; //"A positron has annihilated at rest"
		}

		if( EGS4Macro.n_split > 1 )
		{
		    EGS4.iausfl[7]  = 1; //"After bremsstrahlung"
		    EGS4.iausfl[13] = 1; //"A positron has annihilated in-flight"
		    EGS4.iausfl[14] = 1; //"A positron has annihilated at rest"

		    //"With n_split > 1, we don't need the following calls even if ifano = 1"
		    EGS4.iausfl[15] = 0;
		    EGS4.iausfl[17] = 0;
		    EGS4.iausfl[19] = 0;
		    EGS4.iausfl[23] = 0;
		    EGS4.iausfl[24] = 0;
		}

		if(IFANO == 2)
		{

		    //"AUSGAB will be responsible for discarding energy due to electrons set"
		    //"in motion in the wall"
		    EGS4.iausfl[16] = 1; //"After pair"
		    EGS4.iausfl[18] = 1; //"After Compton"
		    EGS4.iausfl[20] = 1; //"After photoelectric"
		}

		if( EGS4Macro.use_enhance )
		{
		    EGS4.iausfl[15] = 1; //"Before pair"
		    EGS4.iausfl[17] = 1; //"Before Compton"
		    EGS4.iausfl[18] = 1; //"After Compton"
		    EGS4.iausfl[19] = 1; //"Before photoelectric"
		    EGS4.iausfl[20] = 1; //"After photoelectric"
		    EGS4.iausfl[23] = 1; //"Before Rayleigh"
		    EGS4.iausfl[24] = 1; //"After Rayleigh"
		    EGS4.iausfl[7]  = 1; //"After bremsstrahlung"
		    EGS4.iausfl[13] = 1; //"A positron has annihilated in-flight"
		    EGS4.iausfl[14] = 1; //"A positron has annihilated at rest"
		}

		for(int j=1;j<=28;j++)
		{
  			if( EGS4.iausfl[j-1]!=0 )
  			{//write(6,'(i3,$)') j;
  			    int jj=j-1;
		      	EGS4.seqStr=" AUSGAB index call: "+jj;
			  	if(EGS4.iprint>2)
		  			printSequence(EGS4.seqStr);
			}
		}

		//"HATCH CALL PREPARATION AND EXECUTION"
		EGS4.DUNIT=1; //"SET LENGTH UNITS TO CMS"
		HATCH();if(EGS4.STOPPROGRAM){return;}

		if( EGS4Macro.irejct == 1 )
		{
			EGS4Macro.initialize_range_rejection();
		}

		if(EGS4SrcEns.MONOEN==0 && EGS4SrcEns.ISOURC!=21
		&& EGS4SrcEns.ISOURC!=22)
		{//"MONOENERGETIC INPUT BEAM"
		    if(EGS4SrcEns.iqin==0)
		    {
				EI=EGS4SrcEns.ein;
			}
			else
			{
				EI=EGS4SrcEns.ein+EGS4.RM;
			}
		    EKMAX=EGS4SrcEns.ein; //"MAXIMUM KINETIC ENERGY"
		}
		else if(EGS4SrcEns.MONOEN==1)
		{// "ENERGY SPECTRUM"
		    EGS4SrcEns.ENSRC1();// "NORMALIZE THE ENERGY DISTRIBUTION"
		    EKMAX=EGS4SrcEns.ENSRCD[EGS4SrcEns.NENSRC];// "MAXIMUM KINETIC ENERGY IN THE SPECTRUM"
		}
		else if(EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22)
		{// "phase space input"
		    //EKMAX=EKSRCM;//NOT ALLOWED
		}
		else { EKMAX=0; }// " <------------ fixme"

		//"CHECK THAT THE DATA FILE HAD DATA OVER THE ENERGY RANGE REQUIRED"
		for(int I=1;I<=EGS4.NMED;I++)
		{
		    if((EKMAX>EGS4.UP[I-1])||(EKMAX>EGS4.UE[I-1]-EGS4.RM))
		    {
				EGS4.STOPPROGRAM=true;
		      	EGS4.seqStr=" FOR MEDIUM: "+I+"  INCIDENT ENERGY="+
		      	EGS4.format(EKMAX,10,true)+" MeV";
	  			printSequence(EGS4.seqStr);
		      	EGS4.seqStr=" IS GREATER THAN COVERED BY DATA FILE WHERE UP,UE="+
		      	EGS4.format(EGS4.UP[I-1],10,true)+" "+
		      	EGS4.format(EGS4.UE[I-1],10,true)+" MeV";
	  			printSequence(EGS4.seqStr);
				EGS4.seqStr=" EXECUTION WILL BE TERMINATED!";
				printSequence(EGS4.seqStr);

				//return;
			}
		 }// "END OF LOOP OVER MEDIA"
		 if(EGS4.STOPPROGRAM){return;}

		//"CALCULATE THE MASS OF EACH ZONE (AREAL MASS FOR ISOURC=2 OR 4)
		for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)
		{
		    DEPTH=EGS4Geom.ZPLANE[IZ]-EGS4Geom.ZPLANE[IZ-1];
		    for(int IX=1;IX<=EGS4Geom.NR;IX++)
		    {
		        //$GET-IRL(IZ,IX);
		        IRL=EGS4Geom.GET_IRL(IZ,IX);
		        MEDNUM=EGS4.MED[IRL-1];
		        if(MEDNUM!=0)
		        {
		            if((EGS4SrcEns.ISOURC==2)||
		            (EGS4SrcEns.ISOURC==4)){VOLUME=DEPTH;}
		            else
		            {
		                if(IX==1){RLOW2=0.0;}
		                else{RLOW2=EGS4Geom.CYRAD2[IX-2];}
		                VOLUME=Math.PI*DEPTH*(EGS4Geom.CYRAD2[IX-1]-RLOW2);
					}
		            AMASS[IZ-1][IX-1]=EGS4.RHOR[IRL-1]*VOLUME;
				}
		        else{AMASS[IZ-1][IX-1]=0.0;}
		     }//"END OF IX LOOP"
		 }//"END OF IZ LOOP"

		 //"CALCULATE ONE-TIME-ONLY CONSTANTS FOR SOURCE"
		 EGS4SrcEns.SRCOTO();//(WEIGHT);
		 if((EGS4SrcEns.IFPB==0)&&(EGS4Macro.IFORCE!=0)&&
		 (EGS4SrcEns.iqin==0)&&(EGS4SrcEns.MONOEN==0))
		 {
		     EGS4Macro.SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM();
		 }

		 //"INITIALIZE DATA ARRAYS FOR FLUORESCENT X-RAYS IF NEEDED"
		 ISUMX=0;
		 for(int JJ=1;JJ<=EGS4Geom.nreg;JJ++)
		 {
			 ISUMX=ISUMX+EGS4.iedgfl[JJ-1]; //"NON-ZERO IF X-RAYS ANYWHERE"
		 }
		 if(ISUMX!=0)
		 {
			 EGS4.EDGSET(2,EGS4Geom.nreg);
		 }
		 //"NOTE THE ABOVE WILL PRODUCE LOTS OF EXTRA OUTPUT AND SHOULD BE"
		 //"CLEANED UP"
		ISUMRY(); //"PRINT THE SUMMARY OF INPUTS"

		//"******************************************************************************
		//"
		//"                       *** SECTION 2 ***
		//"
		//"------------------------------------------------------------------------------
		//"
		//"LOOP THROUGH THE NUMBER OF HISTORIES. CALCULATE CONSTANTS THAT MAY CHANGE FOR
		//"EACH HISTORY AND DO THE SIMULATION
		//"
		//"------------------------------------------------------------------------------

		//"WRITE THE HEADER"
		//WRITE(IOUT,100) TITLE; call egs_fdate(iout); write(iout,*);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="                   EXECUTION INFORMATION";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		//WRITE(IOUT,200);WRITE(6,200); "PRINT HEADER FOR EXECUTION MESSAGES"

		//"PRINT EXECUTION MODE"
		//if(IRESTART ==  0){WRITE(6,201);WRITE(IOUT,201);}
		//ELSEIF(IRESTART.EQ.1)[
		//    WRITE(6,202) NCASE,NCASEO;
		//    write(6,'(21x,a,$)') 'New RNG state: ';
		//    $SHOW-RNG-STATE(6); write(6,*);
		//    write(iout,'(21x,a,$)') 'New RNG state: ';
		//    $SHOW-RNG-STATE(iout); write(iout,*);
		//]
		//ELSEIF(IRESTART.EQ.3)[WRITE(6,204);WRITE(IOUT,204);GO TO :END-SIM:;]
		//ELSEIF(IRESTART.EQ.4)[WRITE(6,205);WRITE(IOUT,205);]
		//ELSEIF(IRESTART.EQ.5)[WRITE(6,206);WRITE(IOUT,206);GO TO :END-SIM:;]

		//"Initialize IWATCH routine"
		if(IWATCH != 0) EGS4.WATCH(-99,IWATCH);

		//"SET CLOCK AT THE BEGINNING OF SIMULATIONS"
		//$INITIALIZE_ELAPSED_CPU_TIME;
		//$SET_ELAPSED_CPUTIME(CPUT1);
		//$INITIALIZE_ELAPSED_TOTAL_TIME;
		//ETIMETOT=0;
		//TIMEB=0;
		//NETADJ=0;

		//"dcav_old = 0.0;"

		//"Calculate the cavity mass for in-flight dose display"
		FMASSC=0.0; //"TOTAL CAVITY MASS"
		for(int IX=1;IX<=EGS4Geom.NR;IX++)
		{
		    for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)
		    {
		        //$GET-IRL(IZ,IX);
		        IRL=EGS4Geom.GET_IRL(IZ, IX);
		        if(EGS4Geom.ntrack[IRL-1]==1)
		        {
		            FMASS=AMASS[IZ-1][IX-1];
		            FMASSC=FMASSC+FMASS; //"SUM THE CAVITY MASS USED LATER"
		        }
		    }
		}
		//"end of calculation of cavity mass"

		//"Initialize variables for the cs_enhance scoring "
		tmp_dose = 0; tmp_dose1 = 0;
		last_case = 0;EGS4SrcEns.NHSTRY = 0;

		//"Open file for data storage, if requested "
		//"The file is opened in the temporary working directory"
		//IF( idat = 0 ) data_unit = egs_open_file(4,0,1,'.egsdat');

		n_photons=0; n_electrons = 0;
		sumE_photons=0.0; sumE2_photons=0.0; sumE_electrons=0.0; sumE2_electrons=0.0;

		//"IK: New parallel processing implementation. Only used if there is a
		//"    working C compiler.
	//	#ifdef HAVE_C_COMPILER;
	//	;
	//	/last_dose,last2_dose/ = 0;  n_tot = ncaseo;
	//	first_time = .true.; is_finished = .false.;

	//	:start_parallel_loop:;

	//	IF( n_parallel > 0 ) [  "Job is part of a parallel run "

	//	    last_dose = cav_dose - last_dose + tmp_dose;
	//	    last2_dose = cav2_dose - last2_dose + tmp_dose*tmp_dose;
	//	    tmpf = fmassc*ainflu/ncaset/1.602e-10;
	//	    part_dose = last_dose/tmpf; part2_dose = last2_dose/(tmpf*tmpf);
	//	    last_dose = cav_dose + tmp_dose;
	//	    last2_dose = cav2_dose + tmp_dose*tmp_dose;
	//	    call egs_pjob_control(ncase,n_run,n_left,n_tot,part_dose,part2_dose,
	//	                          current_result, current_uncertainty);
	//	    IF( n_run = 0 ) [
	//	        write(6,'(//a,a//)') '****** No histories left in job control file',
	//	                      '       => end simulation';
	//	        goto :END-SIM:;
	//	    ]
	//	    IF( statlm > 0 & current_uncertainty < statlm ) [
	//	        write(6,'(//a,a//)') '****** Desired uncertainty reached',
	//	                      '       => end simulation';
	//	        goto :END-SIM:;
	//	    ]
	//	    jcase = n_run/$NBATCH;
	//	    IF( jcase < 1 ) [ jcase = 1; n_run = jcase*$NBATCH; ]
	//	    IF( first_time ) [
	//	        first_time = .false.; n_last = n_run;
	//	        write(6,'(//a,i12,a//)') '****** Running ',n_run,' histories';
	//	    ]
	//	    ELSE [
	//	        write(6,'(//a,i12,a)') '***** Finished ',n_last,' histories';
	//	        write(6,'(/a/,20x,1pe11.4,a,0pf5.2,a/,a,i12,a//)')
	//	  '      current result including previous runs and other parallel jobs: ',
	//	         current_result, ' +/- ',current_uncertainty,' %',
	//	  '      will run another ',n_run,' histories';
	//	    ]
	//	]
	//	#endif;


		//"Output batches. Statistical analysis is done after each batch. Execution
		//"stops if the desired statistical accuracy is obtained or there is not enough
		//"time to do another batch.
		boolean ENDSIM=false;
		for(int IBATCH=1;IBATCH<=$NBATCH;IBATCH++)
		{
			ENDSIM=false;
			long startTime=System.currentTimeMillis();
		    IBTCH=IBATCH;
		    if(IBATCH == 1)
		    {
				EGS4.seqStr=" BATCH"+EGS4.format("",5)+"ELAPSED"+EGS4.format("",3)+
				"time of day"+EGS4.format("",2)+"cavity stats(%)"+EGS4.format("",2)+
				"cav.dose(Gy.cm^2)";
				if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

//				String timePerBatch=EGS4.timeElapsed(startTime);
//				EGS4.seqStr=EGS4.format("",2)+"1"+EGS4.format("",5)+timePerBatch+
//				EGS4.format("",3)+timeday;
//				//"time of day"+EGS4.format("",2)+"cavity stats(%)"+EGS4.format("",2)+
//				//"cav.dose(Gy.cm^2)";
//				if(EGS4.iprint>1)
//					printSequence(EGS4.seqStr);
//
//		        OUTPUT;
//		        (/' BATCH',2X,'ELAPSED',2X,'CPUtime',2X,'RATIO',2X,
//		        'time of day',2X,'cavity stats(%)',2X,'cav.dose(Gy.cm^2)'//
//		        2X,'1',5X,'0.0',6X,'0.0',6X,'0.00',3X,' ',$); call egs_time(6);
//
//		        //"IK: it is annoing that for batch runs we don't see the progress"
//		        //"    info in the log file until the job is finished. This is because"
//		        //"    Fortran uses buffered I/O. The following flushes unit 6 so that"
//		        //"    we can see the progress of the calculation. "
//		        //$FLUSH_UNIT(6);
		    }
		    else
		    {//" not first batch"
		        //$SET_ELAPSED_TOTAL_TIME(TIMEB);
		        // / *
		        //   we don't need this because the current EGSnrc implementation
		        //   takes care of run times > 24 h.

		        //IF(TIMEB < 0.0)["batch past midnight relative to TZERO"
		        //     NETADJ=NETADJ+1;
		        //     TIMEB=TIMEB + 86400.;
		        //     IF(NETADJ<2)[
		        //       OUTPUT; (' Wall clock has gone past 24:00 hrs.'/
		        //       ' Elapsed time adjusted assuming batches took < 1 day',
		        //       ' to complete.');
		        //     ]
		        //]
		        // * /
		        //ETIMETOT = ETIMETOT+TIMEB;
		        //$INITIALIZE_ELAPSED_TOTAL_TIME;"reset to get elapsed time/batch"
		        //$SET_ELAPSED_CPUTIME(CPUT2);
		        //TIMCPU = (CPUT2-CPUT1)*$CONVERSION_TO_SECONDS+$TIME_RESOLUTION;

		        //"*****************************************************************"
		        //OUTPUT IBATCH,ETIMETOT,TIMCPU,ETIMETOT/TIMCPU;
		        //          (1X,I2,F8.1,1X,F8.1,2X,F8.2,3X,' ',$);  call egs_time(6);
		        //"*****************************************************************"

		        //"IK: it is annoing that for batch runs we don't see the progress"
		        //"    info in the log file until the job is finished. This is because"
		        //"    Fortran uses buffered I/O. The following flushes unit 6 so that"
		        //"    we can see the progress of the calculation. "
		        //$FLUSH_UNIT(6);

		        //"Check there is time left for another batch"
		        //BATCHT = TIMCPU/dble(IBATCH-1);//"time per batch so far"
		        //IF(TIMCPU+1.1*BATCHT.GT.TIMMAX*3600.)[
		        //    "not enough time for another batch"
		        //    "print message and exit simulation loop"
		        //    WRITE(IOUT,210) TIMMAX,IBATCH-1,IHSTRY-NCASEO,IHSTRY;
		        //    WRITE(6,210) TIMMAX,IBATCH-1,IHSTRY-NCASEO,IHSTRY;
		        //    "adjust the incident fluence"
		        //    "IK: do it at the end for all possible exits from the
		        //    "shower loop
		        //    "AINFLU = AINFLU*dble(IHSTRY)/dble(NCASET);
		        //    GO TO :END-SIM:;
		        //}
		    }//" end of before batch ne 1 block"

		    for( int ICASE = 1;ICASE<=JCASE;ICASE++)
		    {//"now fill each IS bin"

		            if( EGS4SrcEns.ISOURC != 23 ) IHSTRY = IHSTRY+1; //"increment history counter"
		              //" For source 23 ihstry is set in srchst "

		            //:CORRELATION-RESTART:; //"restart here for a correlated history"

		            EGS4Macro.NFTIME = 0; //"reset the photon forced interaction counter"

		            //"retrieve starting random numbers if reading from a file"
		            //IF(IRESTART.EQ.4) [
		            //    $RETRIEVE RNG STATE FROM UNIT rng_unit; "was $RESET-RNG(2);"
		            //]

		            //"store the initial random number, another pass might be needed"
		            //IRNG = 1;
		            //$STORE-RNG(IRNG);     "was $STORE-RNG(0);"

		            //"store initial random #s if requested"
		            //IF(ISTORE =  1)[
		            //    $STORE RNG STATE ON UNIT rng_unit;
		            //]
		            //ELSEIF(ISTORE =  2)[
		            //    "temporarily store the initial random number seed"
		            //    IRNG = 3;
		            //    $STORE-RNG(IRNG);  "was $STORE-RNG(-2);"
		            //    "clear the flag that signals energy deposit in the cavity"
		            //    IDECAV = 0;
		            //]
		            //ELSEIF(ISTORE.EQ.3)[
		            //    "STORE THE INITIAL RANDOM NUMBER SEED"
		            //    $PUT RNG STATE ON UNIT rng_unit;
		            //]

		            //"calculate the source dependent values which change for each
		            //"history these include :
		            //    "entry point into target,
		            //    "initial direction cosines,
		            //    "statistical weight,
		            //    "entry flag(nrcflg)
		            //CALL SRCHST(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
		            EGS4SrcEns.SRCHST();

		            //"calculate the initial energy if a distribution is to be used"
		            if(EGS4SrcEns.MONOEN != 0)
		            {//"if equal to 0, it is monoenergetic"
		                //EGS4EnsSrc.ENSRCH(EIN);   "returns K.E. from distribution"
		                EGS4SrcEns.ein=EGS4SrcEns.ENSRCH();
		                if(EGS4SrcEns.iqin==0)
		                {EI = EGS4SrcEns.ein;}else{EI = EGS4SrcEns.ein+EGS4.RM;}// "total energy"
		                //" there was a check that the data file had data over the energy
		                //"range required, the location of it will eventually be in
		                //"ESRCIN.MOR
		            }
		            //ELSEIF(ISOURC = 21 | ISOURC = 22 | ISOURC = 23 )[ EI = EIN; ]

		            //"Set photon weights if gamma interactions are to be forced in the
		            //"target in the frontal parallel beam case if monoenergetic
		            if((EGS4SrcEns.MONOEN==0)&&
		            (EGS4SrcEns.iqin==0)&&(EGS4Macro.IFORCE==1)&&(EGS4SrcEns.IFPB==0)&&
		            (EGS4SrcEns.ISOURC!=21)&&(EGS4SrcEns.ISOURC!=22))
		            {
		                  int IX=(EGS4SrcEns.irin-2)/EGS4Geom.NZ+1;
		                  EGS4Macro.GWAIT=EGS4Macro.GWATE[IX-1];
		                  EGS4SrcEns.WEIGHT=EGS4Macro.GWAIT;
		            }

		            //"FOR AN INPUT ENERGY SPECTRUM, DETAILED FORCING MACRO IS USED"

		            EGS4.LATCHI=0;

		            if((IWATCH != 0) && (IWATCH != 4))
		            {
						EGS4.seqStr=EGS4.format(1,5)+EGS4.format(EGS4SrcEns.ein,9,true)+
						EGS4.format(EGS4SrcEns.iqin,4)+EGS4.format(EGS4SrcEns.irin,4)+
						EGS4.format(EGS4SrcEns.xin,8,true)+EGS4.format(EGS4SrcEns.yin,8,true)+
						EGS4.format(EGS4SrcEns.zin,8,true)+EGS4.format(EGS4SrcEns.uin,8,true)+
						EGS4.format(EGS4SrcEns.vin,8,true)+EGS4.format(EGS4SrcEns.win,8,true)+
						EGS4.format(EGS4.LATCHI,10)+EGS4.format(EGS4SrcEns.WEIGHT,10,false);
						if(EGS4.iprint>1)
							printSequence(EGS4.seqStr);
		                //OUTPUT 1,EIN,IQIN,IRIN,XIN,YIN,ZIN,UIN,VIN,WIN,LATCHI,WEIGHT;
		                //(/' INITIAL SHOWER VALUES',T36,':',
		                //I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
		            }
		            //"ALL INITIAL SHOWER VARIABLES ARE SET, CALL THE SHOWER ROUTINE"
		            if( EGS4SrcEns.iqin == 0 )
		            {
		                n_photons = n_photons + 1;
		                sumE_photons = sumE_photons + EI;
		                sumE2_photons = sumE2_photons + EI*EI;
		            }
		            else
		            {
		                n_electrons = n_electrons + 1;
		                sumE_electrons = sumE_electrons + EI-EGS4.RM;
		                sumE2_electrons = sumE2_electrons + (EI-EGS4.RM)*(EI-EGS4.RM);
		            }

		            //CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
		            SHOWER();if(EGS4.STOPPROGRAM){return;}

		            //IF(ISTORE.EQ.2)[
		            //    "STORE INITIAL RN SEED IF ENERGY SCORED IN THE CAVITY"
		            //     IF(IDECAV.EQ.1)[ "ENERGY WAS DEPOSITED IN THE CAVITY"
		            //             IRNG = 1;
		            //             $STORE-RNG(IRNG);              "was $STORE-RNG(0);"
		            //             IRNG = 3;
		            //             $RESET-RNG(IRNG);              "was $RESET-RNG(-2);"
		            //             $PUT RNG STATE ON UNIT rng_unit;"was $STORE-RNG(2);"
		            //             IRNG = 1;
		            //             $RESET-RNG(IRNG);              "was $RESET-RNG(0);"
		            //     ]
		            //]

		            //"SIGNAL THE END OF A HISTORY IF WATCH MODE IS SET"
		            //"JAN CHANGE: ZERO STACK ELEMENTS TO PREVENT CROSS TALKING"
		            //"IK: Huh ? $CLEAN-STACK;"
		            //"END OF JAN CHANGE"
		            if(IWATCH>0) EGS4.WATCH(-1,IWATCH);

		    }// "END OF THE ICASE LOOP"

		    //"Succesful completion of a batch. delete the raw data from the last batch"
		    //"and record the new batch only if requested"
		    //if(IDAT == 0)
		    //{

		        //" discard data from previous batches (if any) "
		        //rewind(data_unit);

		    //    TSCSTP=SCSTP+SCSTP_TMP;
		    //    TSCSTP2=SCSTP2+SCSTP_TMP*SCSTP_TMP;
		    //    TSCCSTP=SCCSTP+SCCSTP_TMP;
		    //    TSCCSTP2=SCCSTP2+SCCSTP_TMP*SCCSTP_TMP;

		        //WRITE(data_unit,*)TSCSTP,TSCSTP2,TSCCSTP,TSCCSTP2;

		//"******************"
		//"history by history"
		//"  EMH  March, 2002"
		//"******************"
		    //    tcav_dose  = cav_dose  + tmp_dose;
		    //    tcav2_dose = cav2_dose + tmp_dose*tmp_dose;

		    //    tcav_dose0  = cav_dose0   + tmp_dose0;
		    //    tcav2_dose0 = cav2_dose0 + tmp_dose0*tmp_dose0;

		    //    tcav_dose1  = cav_dose1  + tmp_dose1;
		    //    tcav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1;

		    //    tcav_dose2  = cav_dose2  + tmp_dose2;
		    //    tcav2_dose2 = cav2_dose2 + tmp_dose2*tmp_dose2;

		    //    tcav_dosec   = cav_dosec   + tmp_dose*tmp_dose1;
		    //    tcav_dosec01 = cav_dosec01 + tmp_dose0*tmp_dose1;
		    //    tcav_dosec02 = cav_dosec02 + tmp_dose0*tmp_dose2;

		    //    write(data_unit,*) tcav_dose, tcav_dose0, tcav_dose1, tcav_dose2;
		    //    write(data_unit,*) tcav2_dose,tcav2_dose0,tcav2_dose1,tcav2_dose2;
		    //    write(data_unit,*) tcav_dosec,tcav_dosec01,tcav_dosec02;

		//"*****************"

		    //    IF(NSUMCV>1)["store data for individual cavity regions"

		    //       DO IZ=1,NZ[
		    //           DO IX=1,NR[
		    //              $GET-IRL(IZ,IX);
		    //              IF(NTRACK(IRL).EQ.1)[
		    //               DO IT=1,4[
		    //                TSCDOSE(IZ,IX,IT)=SCDOSE(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT);
		    //                TSCDOSE2(IZ,IX,IT)=SCDOSE2(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT)*
		    //                                                      SCDOSE_TMP(IZ,IX,IT);
		    //               ]
		    //               TSCDOSE_COV(IZ,IX,1)=SCDOSE_COV(IZ,IX,1)+
		    //                                   SCDOSE_TMP(IZ,IX,1)*SCDOSE_TMP(IZ,IX,3);
		    //               TSCDOSE_COV(IZ,IX,2)=SCDOSE_COV(IZ,IX,2)+
		    //                                   SCDOSE_TMP(IZ,IX,2)*SCDOSE_TMP(IZ,IX,3);
		    //               TSCDOSE_COV(IZ,IX,3)=SCDOSE_COV(IZ,IX,3)+
		    //                                   SCDOSE_TMP(IZ,IX,2)*SCDOSE_TMP(IZ,IX,4);
		    //               WRITE(data_unit,*)
		    //                 (TSCDOSE(IZ,IX,IT),TSCDOSE2(IZ,IX,IT),IT=1,4);
		    //               WRITE(data_unit,*)(TSCDOSE_COV(IZ,IX,IT),IT=1,3);
		    //              ]
		    //           ]
		    //       ]
		    //    ]

		    //]"end of conditional data storage"

		    //$SET_ELAPSED_CPUTIME(CPUT2);
		    //TIMCPU=$CONVERSION_TO_SECONDS*(CPUT2-CPUT1)+TMCPUO;
		    //IF(IDAT.EQ.0)[
		    //    $PUT RNG STATE ON UNIT data_unit;
		    //    WRITE(data_unit,*) IHSTRY,TIMCPU,NNREAD,PIISTP+count_pII_steps;
		    //]
		//"******************"
		//"history by history"
		//"    EMH March 2002"
		//"******************"
		    if(EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22)
		    {
		     // AINFLU_CURRENT=
		     //       EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/EGS4SrcEns.dble(NCASE_PHSP)*NINCSRC;
		    }
		    else if( EGS4SrcEns.ISOURC==23 ) { AINFLU_CURRENT=IHSTRY; }
		    else
		    {
				AINFLU_CURRENT=EGS4SrcEns.AINFLU*EGS4SrcEns.dble(IHSTRY)/EGS4SrcEns.dble(EGS4SrcEns.NCASET);
			}
		    dcav_current = cav_dose*1.602E-10/(FMASSC*AINFLU_CURRENT);
		    //double mevtojoule=1.60218E-13;//1MeV=1,60218 · 10-13 J
		    //and MASS is in g not KG!!!=>E-10!!
		    // // *
		    //IF( ISOURC=23 ) [
		    //    dcavun = (cav2_dose*nhstry-cav_dose*cav_dose)/(nhstry-1);
		    //]
		    //ELSE [
		    //    dcavun = (cav2_dose*ihstry-cav_dose*cav_dose)/(ihstry-1);
		    //]
		    // * /

		    dcavun = (cav2_dose*IHSTRY-cav_dose*cav_dose)/(IHSTRY-1);
		    if( dcavun > 0 ) { dcavun = Math.sqrt(dcavun); };
		    dcavun = dcavun*1.602E-10/(FMASSC*AINFLU_CURRENT);
		    dcavun = 100*dcavun/dcav_current;

//#####################################################################################
				String timePerBatch=EGS4.timeElapsedShort(startTime);

				Calendar call=Calendar.getInstance();
				String timeday=call.get(Calendar.HOUR)+":"+
				call.get(Calendar.MINUTE)+":"+call.get(Calendar.SECOND);

				EGS4.seqStr=EGS4.format("",2)+EGS4.format(IBATCH,3)+EGS4.format("",5)+timePerBatch+
				EGS4.format("",3)+timeday+EGS4.format("",6)+EGS4.format(dcavun,7,true)+
				EGS4.format("",10)+EGS4.format(dcav_current,11,false);
				if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

		        //OUTPUT;
		        //(/' BATCH',2X,'ELAPSED',2X,'CPUtime',2X,'RATIO',2X,
		        //'time of day',2X,'cavity stats(%)',2X,'cav.dose(Gy.cm^2)'//
		        //2X,'1',5X,'0.0',6X,'0.0',6X,'0.00',3X,' ',$); call egs_time(6);

		        //"IK: it is annoing that for batch runs we don't see the progress"
		        //"    info in the log file until the job is finished. This is because"
		        //"    Fortran uses buffered I/O. The following flushes unit 6 so that"
		        //"    we can see the progress of the calculation. "
		        //$FLUSH_UNIT(6);
//#####################################################################################
//OUTPUT DCAVUN,DCAV_CURRENT;('+',5X,F5.2,10X,1PE11.4,5X,1PE10.3);

		    if(dcavun<=STATLM && STATLM>0.)
		    {
		      //"we have reached the desired statistics, print a message and exit"
				EGS4.seqStr=" DESIRED STATISTICAL ACCURACY OBTAINED.";
				if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
				EGS4.seqStr=" STATS IN CAVITY= "+EGS4.format(dcavun,5,true)+"%"+
				" AFTER "+EGS4.format(IBTCH,2)+" BATCHES";
				if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

		      //WRITE(6,230)dcavun,IBTCH;WRITE(IOUT,230)dcavun,IBTCH;
		      //GO TO :END-SIM:;
		      ENDSIM=true;
		      break;
		    }
		}// "END OF SIMULATIONS"

		//#ifdef HAVE_C_COMPILER;
		//;
		//IF( n_parallel > 0 ) [ goto :start_parallel_loop:; ]
//
		//#endif;

		//"PRINT INSUFFICIENT STATS WARNING"
		if(!ENDSIM)
		{
			EGS4.seqStr=" DESIRED STATISTICAL ACCURACY OF "+EGS4.format(STATLM,5,true)+"%"+
			" NOT REACHED!";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" STATS IN CAVITY= "+EGS4.format(dcavun,5,true)+" %"+
			" AFTER "+EGS4.format(IBTCH,2)+" BATCHES";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

			//WRITE(IOUT,240) STATLM,DCAVUN,IBTCH;WRITE(6,240) STATLM,DCAVUN,IBTCH;
		}

		//:END-SIM:;

		SCSTP=SCSTP+SCSTP_TMP; SCSTP2=SCSTP2+SCSTP_TMP*SCSTP_TMP;
		SCCSTP=SCCSTP+SCCSTP_TMP; SCCSTP2=SCCSTP2+SCCSTP_TMP*SCCSTP_TMP;

		//"******************"
		//"history by history"
		//"    EMH March 2002"
		//"******************"
		//"we have to add the temporary scoring"
		//"variables for the last history      "
		cav_dose  = cav_dose  + tmp_dose; cav2_dose = cav2_dose + tmp_dose*tmp_dose;

		cav_dose0  = cav_dose0   + tmp_dose0;
		cav2_dose0 = cav2_dose0 + tmp_dose0*tmp_dose0;

		cav_dose1  = cav_dose1  + tmp_dose1;
		cav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1;

		cav_dose2  = cav_dose2  + tmp_dose2;
		cav2_dose2 = cav2_dose2 + tmp_dose2*tmp_dose2;

		cav_dosec   = cav_dosec   + tmp_dose*tmp_dose1;
		cav_dosec01 = cav_dosec01 + tmp_dose0*tmp_dose1;
		cav_dosec02 = cav_dosec02 + tmp_dose0*tmp_dose2;

		if(EGS4Geom.NSUMCV>1)
		{//"store data for individual cavity regions"

		    for (int IZ=1;IZ<=EGS4Geom.NZ;IZ++)
		    {
				for(int IX=1;IX<=EGS4Geom.NR;IX++)
				{
		        	//$GET-IRL(IZ,IX);
		        	IRL=EGS4Geom.GET_IRL(IZ,IX);

		        	if(EGS4Geom.ntrack[IRL-1]==1)
		        	{
		            	for(int IT=1;IT<=4;IT++)
		            	{
		            	    SCDOSE[IZ-1][IX-1][IT-1]=SCDOSE[IZ-1][IX-1][IT-1]+
		            	    SCDOSE_TMP[IZ-1][IX-1][IT-1];
		            	    SCDOSE2[IZ-1][IX-1][IT-1]=SCDOSE2[IZ-1][IX-1][IT-1]+
		            	    SCDOSE_TMP[IZ-1][IX-1][IT-1]*SCDOSE_TMP[IZ-1][IX-1][IT-1];
		            	}
		            	SCDOSE_COV[IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]+
		                        SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][2];
		            	SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]+
		                        SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][2];
		            	SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]+
		                        SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][3];
		        	}
		     	}
			}
		}

		EGS4.seqStr=" FINAL RANDOM NUMBER STATE:";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.SHOW_RNG_STATE();
		//write(6,'(/a,$)') '****** FINAL RANDOM NUMBER STATE:';
		//$SHOW-RNG-STATE(6); write(6,'(a)') '  ******';
		//write(iout,'(/a,$)') '********* FINAL RANDOM NUMBER STATE:';
		//$SHOW-RNG-STATE(iout); write(iout,'(a)') '  *********';
		//$SET_ELAPSED_TOTAL_TIME(TIMEB);
		//ETIMETOT=ETIMETOT+TIMEB;
		//$SET_ELAPSED_CPUTIME(CPUT2);
		//TIMCPU=(CPUT2-CPUT1)*$CONVERSION_TO_SECONDS+TMCPUO;
		//IF(IRESTART=3)["just analyzing data--no elapsed time"
		//  WRITE(IOUT,250)TMCPUO,TMCPUO/3600;
		//  WRITE(6,250)TMCPUO,TMCPUO/3600;
		//  WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		//  WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		//]
		//ELSEIF(IRESTART=5)["output time results for parallel runs"
		//  WRITE(IOUT,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT);
		//  WRITE(6,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT);
		//  WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		//  WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		//]
		//ELSE[
		//  IF(TMCPUO.EQ.0)["This is first run"
		//   RATIO=ETIMETOT/TIMCPU;
		//   WRITE(IOUT,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO;
		//   WRITE(6,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO;
		//  ]
		//  ELSE[ "There was previous run, but don't have elapsed time for it"
		//   RATIO = ETIMETOT/(TIMCPU-TMCPUO); "ratio this run only"
		//   WRITE(IOUT,261) ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO;
		//   WRITE(6,261) ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO;
		//  ]
		//  IF(IHSTRY.NE.0 .AND. TIMCPU.NE.0.0) ["this should always happen"
		//   WRITE(IOUT,280) TIMCPU/dble(IHSTRY),3600.*dble(IHSTRY)/TIMCPU;
		//  ]
		//]



		//"******************************************************************************
		//"
		//"                       *** SECTION 3 ***
		//"
		//"------------------------------------------------------------------------------
		//
		//"STATISTICAL AND OTHER DATA HANDLING AND CALL THE OUTPUT SUMMARY ROUTINE"
		//
		//"------------------------------------------------------------------------------

		//:STATS-ANAL:;

		FMASSC=0.0; //"total cavity mass"
		for(int IX=1;IX<=EGS4Geom.NR;IX++)
		{
			for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ
			{
			    //$GET-IRL(IZ,IX);
			    IRL=EGS4Geom.GET_IRL(IZ,IX);
			    if(EGS4Geom.ntrack[IRL-1]==1)
			    {
			        FMASS=AMASS[IZ-1][IX-1];
			        FMASSC=FMASSC+FMASS;
			    }
			}
		}

		//IF(ISOURC = 21|ISOURC = 22)[
		//   "normalize dose to number of incident particles from primary source
		//   AINFLU = dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*NINCSRC;
		//   "we estimate the total number of particles from the primary source "
		//   "(original non-phase space source) by taking the ratio of the total "
		//   "number of particles read from the phase source in this simulation "
		//   "to the total number of particles in the phase space source and multiply "
		//   "this by the number of particles from the primary source that were used "
		//   "to obtain this phase space  source."
		//   OUTPUT AINFLU, NNREAD, NNREAD-IHSTRY, NCASE_PHSP, NINCSRC;
		//   (/' Corresponding number of particles in original BEAM simulation =', F13.0/
		//     ' Based on reading', I13, ' particles from the phase space file,'/
		//     ' rejecting',I13,' particles (due to wrong charge, missing the'/
		//     ' geometry or being multiple passers),'/
		//     ' and  file having', I13, ' particles in it.'/
		//     ' The phase space file was generated by', F13.0, ' initial particles'/);
		//   WRITE(IOUT,:tgh1:) AINFLU, NNREAD, NNREAD-IHSTRY, NCASE_PHSP, NINCSRC;
		//   :tgh1: FORMAT(/' Corresponding number of particles in original BEAM',
		//     ' simulation =', F13.0/
		//     ' Based on reading', I13, ' particles from the phase space file,'/
		//     ' rejecting',I13,' particles (due to wrong charge, missing the'/
		//     ' geometry or being multiple passers),'/
		//     ' and phase space file having', I13, ' particles in it'/
		//     ' The phase space file was generated by', F13.0, ' initial particles'/);
		//   SCORE_NORM_NUM=AINFLU;
		//]
		//ELSEIF( ISOURC = 23 ) [
		//    // *
		//    //AINFLU = NHSTRY; SCORE_NORM_NUM=NHSTRY;
		//    // * /
		//    AINFLU = IHSTRY; SCORE_NORM_NUM=IHSTRY;
		//    OUTPUT AINFLU;
		//    (/'Source 23: number of particles in BEAM source: ',F13.0/);
		//    IF( n_photons > 0 ) [
		//        sumE_photons = sumE_photons/n_photons;
		//        sumE2_photons = sumE2_photons/n_photons;
		//        sumE2_photons = sumE2_photons - sumE_photons*sumE_photons;
		//        IF( sumE2_photons > 0 ) sumE2_photons = sqrt(sumE2_photons/n_photons);
		//        OUTPUT n_photons,sumE_photons,sumE2_photons;
		//        ('   Number of photons:      ',i10/
		//         '   Average energy:         ',f10.5,' +/- ',f10.5,/);
		//    ] ELSE [ OUTPUT; ('   Number of photons:        0'); ]
		//    IF( n_electrons > 0 ) [
		//        sumE_electrons = sumE_electrons/n_electrons;
		//        sumE2_electrons = sumE2_electrons/n_electrons;
		//        sumE2_electrons = sumE2_electrons - sumE_electrons*sumE_electrons;
		//        IF( sumE2_electrons > 0 ) [
		//            sumE2_electrons = sqrt(sumE2_electrons/n_electrons);
		//        ]
		//        OUTPUT n_electrons,sumE_electrons,sumE2_electrons;
		//        ('   Number of electrons:    ',i10/
		//         '   Average energy:         ',f10.5,' +/- ',f10.5,/);
		//    ] ELSE [ OUTPUT; ('   Number of electrons:      0'); ]

		//]
		//ELSE[
		   //"IK: adjust incident fluence
		   EGS4SrcEns.AINFLU =
		   	EGS4SrcEns.AINFLU*EGS4SrcEns.dble(IHSTRY)/EGS4SrcEns.dble(EGS4SrcEns.NCASET);
		   SCORE_NORM_NUM=EGS4SrcEns.dble(IHSTRY);
		//]

		//"******************"
		//"history by history"
		//"    EMH March 2002"
		//"******************"

//REPLACE {$ANALYZE(#,#:#)} WITH {;
//"Macro to analyze uncertainty:"
//"{P1}{P2}=scoring array (eg SCDOSE(IDZ,IDX,ITDOSE))"
//"{P3}=quantity to normalize by (eg incident no. of particles)"
//"Calculates the uncertainty on {P1}{P2}/{P3}.  The "
//"uncertainty is stored in {P1}2{P2} and is expressed as a percentage of"
//"{P1}{P2}/{P3} (max 99.9%).  Note that you must define the REAL*8 variable"
//"SCORE_TEMP in any subroutine where this macro is used.  This macro"
//"is only used in the analysis of no. of steps."

//SCORE_TEMP={P1}{P2}/{P3};
//{P1}2{P2}={P1}2{P2}/{P3};
//{P1}2{P2}=({P1}2{P2}-SCORE_TEMP*SCORE_TEMP)/({P3}-1);
//IF({P1}2{P2}>=0.) {P1}2{P2}= SQRT({P1}2{P2});
//IF(SCORE_TEMP~=0.)[
//    {P1}2{P2}= MIN({P1}2{P2}/SCORE_TEMP*100.,99.9D00);
//]
//ELSE[
//    {P1}2{P2}=99.9D00;
//]

		//$ANALYZE(SCOMEG, :SCORE_NORM_NUM);
		SCORE_TEMP=EGS4SrcEns.SCOMEG/SCORE_NORM_NUM;
		EGS4SrcEns.SCOMEG2=EGS4SrcEns.SCOMEG2/SCORE_NORM_NUM;
		EGS4SrcEns.SCOMEG2=
			(EGS4SrcEns.SCOMEG2-SCORE_TEMP*SCORE_TEMP)/(SCORE_NORM_NUM-1);
		if(EGS4SrcEns.SCOMEG2>=0.) EGS4SrcEns.SCOMEG2= Math.sqrt(EGS4SrcEns.SCOMEG2);
		if(SCORE_TEMP!=0.)
		{
		    EGS4SrcEns.SCOMEG2= Math.min(EGS4SrcEns.SCOMEG2/SCORE_TEMP*100.,99.9);
		}
		else
		{
		    EGS4SrcEns.SCOMEG2=99.9;
		}


		EGS4SrcEns.SCOMEG = EGS4SrcEns.SCOMEG/EGS4SrcEns.dble(IHSTRY);//"Corrected, IK May 4 1999"

		//OUTPUT SCOMEG,SCOMEG2;(/' OMEG =',1PE12.3,'(',0PF5.1,'%)'/);
		EGS4.seqStr=" OMEG ="+EGS4.format(EGS4SrcEns.SCOMEG,12,false)+"("+
			EGS4.format(EGS4SrcEns.SCOMEG2,5)+"%)";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);


		//$ANALYZE(SCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP=SCSTP/SCORE_NORM_NUM;
		SCSTP2=SCSTP2/SCORE_NORM_NUM;
		SCSTP2=
			(SCSTP2-SCORE_TEMP*SCORE_TEMP)/(SCORE_NORM_NUM-1);
		if(SCSTP2>=0.) SCSTP2= Math.sqrt(SCSTP2);
		if(SCORE_TEMP!=0.)
		{
		    SCSTP2= Math.min(SCSTP2/SCORE_TEMP*100.,99.9);
		}
		else
		{
		    SCSTP2=99.9;
		}

		//$ANALYZE(SCCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP=SCCSTP/SCORE_NORM_NUM;
		SCCSTP2=SCCSTP2/SCORE_NORM_NUM;
		SCCSTP2=
			(SCCSTP2-SCORE_TEMP*SCORE_TEMP)/(SCORE_NORM_NUM-1);
		if(SCCSTP2>=0.) SCCSTP2= Math.sqrt(SCCSTP2);
		if(SCORE_TEMP!=0.)
		{
		    SCCSTP2= Math.min(SCCSTP2/SCORE_TEMP*100.,99.9);
		}
		else
		{
		    SCCSTP2=99.9;
		}

		//"first, estimate uncertainties for entire cavity"
		cav2_dose = (cav2_dose*SCORE_NORM_NUM - cav_dose*cav_dose)/
		                (SCORE_NORM_NUM-1);
		if( cav2_dose > 0 ) cav2_dose = Math.sqrt(cav2_dose);

		if( EGS4Macro.IFULL == 1 )
		{

		        cav2_dose0 = (cav2_dose0*SCORE_NORM_NUM - cav_dose0*cav_dose0)/
		                     (SCORE_NORM_NUM-1);
		        if( cav2_dose0 > 0 ) cav2_dose0 = Math.sqrt(cav2_dose0);
		        cav2_dose1 = (cav2_dose1*SCORE_NORM_NUM - cav_dose1*cav_dose1)/
		                     (SCORE_NORM_NUM-1);
		        if( cav2_dose1 > 0 ) cav2_dose1 = Math.sqrt(cav2_dose1);
		        cav2_dose2 = (cav2_dose2*SCORE_NORM_NUM - cav_dose2*cav_dose2)/
		                     (SCORE_NORM_NUM-1);
		        if( cav2_dose2 > 0 ) cav2_dose2 = Math.sqrt(cav2_dose2);

		        corr_02=(cav_dosec02*SCORE_NORM_NUM-cav_dose0*cav_dose2)/
		                (SCORE_NORM_NUM-1);
		        corr_02=cav2_dose0*cav2_dose0+cav2_dose2*cav2_dose2+ 2*corr_02;
		        if (corr_02 > 0) corr_02 = Math.sqrt(corr_02);

		        cav_dosec   = (cav_dosec*SCORE_NORM_NUM   - cav_dose*cav_dose1)/
		                      (SCORE_NORM_NUM-1);
		        cav_dosec01 = (cav_dosec01*SCORE_NORM_NUM - cav_dose0*cav_dose1)/
		                      (SCORE_NORM_NUM-1);
		        cav_dosec02 = (cav_dosec02*SCORE_NORM_NUM - cav_dose0*cav_dose2)/
		                      (SCORE_NORM_NUM-1);
		}

		cav_dose  = cav_dose*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav_dose0 = cav_dose0*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav_dose1 = cav_dose1*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav_dose2 = cav_dose2*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);

		cav2_dose  = cav2_dose*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav2_dose0 = cav2_dose0*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav2_dose1 = cav2_dose1*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		cav2_dose2 = cav2_dose2*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);

		corr_02    = corr_02*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);

		cav_dosec   =
		cav_dosec*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC));
		cav_dosec01 =
		cav_dosec01*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC));
		cav_dosec02 =
		cav_dosec02*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC));

		if(EGS4Geom.NSUMCV>1)
		{//"multiple cavity regions, analyze quantities in each region"

		 //   "FOR ISOURC=4 WE NEED THE DATA FOR CIRCLES, NOT RINGS, SO ADD IT UP"
		 //   "THIS SHOULD ONLY BE USED IF THE CAVITY HAS AN INFINITE DIAMETER"
		    if((EGS4SrcEns.ISOURC==4)&&(EGS4Geom.NR>1))
		    {
		       for(int IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[
		       {
		            for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ[
		            {
		               //$GET-IRL(IZ,IX);
		               IRL=EGS4Geom.GET_IRL(IZ,IX);
		               if(EGS4Geom.ntrack[IRL-1]==1)
		               {
		                   for(int IT=1;IT<=4;IT++)
		                   {
		                    //    "IK: this must be a bug. for ix=1 ix-1=0 and the"
		                    //    "    array is not defined"
		                      SCDOSE[IZ-1][IX-1][IT-1]=SCDOSE[IZ-1][IX-1][IT-1]+
		                      	SCDOSE[IZ-1][IX-2][IT-1];
		                      SCDOSE2[IZ-1][IX-1][IT-1]=SCDOSE2[IZ-1][IX-1][IT-1]+
		                      	SCDOSE2[IZ-1][IX-2][IT-1];
		                      if(IT<4)
		                      {
		                        SCDOSE_COV[IZ-1][IX-1][IT-1]=SCDOSE_COV[IZ-1][IX-1][IT-1]+
		                                       SCDOSE_COV[IZ-1][IX-2][IT-1];
		                      }
		                  }
		              }
		           }
		      }
		    }

		    for(int IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[
		    {
		        for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ[
		        {
		            //$GET-IRL(IZ,IX);
		            IRL=EGS4Geom.GET_IRL(IZ,IX);
		            if(EGS4Geom.ntrack[IRL-1]==1)//NTRACK(IRL).EQ.1)[
		            {
		                SCDOSE2[IZ-1][IX-1][0]=(SCDOSE2[IZ-1][IX-1][0]*SCORE_NORM_NUM -
		                     SCDOSE[IZ-1][IX-1][0]*SCDOSE[IZ-1][IX-1][0])/(SCORE_NORM_NUM-1);
		                if( SCDOSE2[IZ-1][IX-1][0] > 0 )
		                	SCDOSE2[IZ-1][IX-1][0]=Math.sqrt(SCDOSE2[IZ-1][IX-1][0]);
		                if(EGS4Macro.IFULL==1)
		                {
		                  SCDOSE2[IZ-1][IX-1][1]=(SCDOSE2[IZ-1][IX-1][1]*SCORE_NORM_NUM -
		                     SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][1])/(SCORE_NORM_NUM-1);
		                  if( SCDOSE2[IZ-1][IX-1][1] > 0 )
		                  	SCDOSE2[IZ-1][IX-1][1]=Math.sqrt(SCDOSE2[IZ-1][IX-1][1]);
		                  SCDOSE2[IZ-1][IX-1][2]=(SCDOSE2[IZ-1][IX-1][2]*SCORE_NORM_NUM -
		                    SCDOSE[IZ-1][IX-1][2]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1);
		                  if( SCDOSE2[IZ-1][IX-1][2] > 0 )
		                  	SCDOSE2[IZ-1][IX-1][2]=Math.sqrt(SCDOSE2[IZ-1][IX-1][2]);
		                  SCDOSE2[IZ-1][IX-1][3]=(SCDOSE2[IZ-1][IX-1][3]*SCORE_NORM_NUM -
		                    SCDOSE[IZ-1][IX-1][3]*SCDOSE[IZ-1][IX-1][3])/(SCORE_NORM_NUM-1);
		                  if( SCDOSE2[IZ-1][IX-1][3] > 0 )
		                  	SCDOSE2[IZ-1][IX-1][3]=Math.sqrt(SCDOSE2[IZ-1][IX-1][3]);
		                  //"now calculate the covariances"

		                  SCDOSE_COV[IZ-1][IX-1][0]=
		                  (SCDOSE_COV[IZ-1][IX-1][0]*SCORE_NORM_NUM -
		                  SCDOSE[IZ-1][IX-1][0]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1);

		                  SCDOSE_COV[IZ-1][IX-1][1]=
		                  (SCDOSE_COV[IZ-1][IX-1][1]*SCORE_NORM_NUM -
		                  SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1);

		                  SCDOSE_COV[IZ-1][IX-1][2]=
		                  (SCDOSE_COV[IZ-1][IX-1][2]*SCORE_NORM_NUM -
                          SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][3])/(SCORE_NORM_NUM-1);
		                }

		                //"now normalize quantities and convert to dose"
		                FMASS=AMASS[IZ-1][IX-1];

		                SCDOSE[IZ-1][IX-1][0]=
		                SCDOSE[IZ-1][IX-1][0]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE[IZ-1][IX-1][1]=
		                SCDOSE[IZ-1][IX-1][1]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE[IZ-1][IX-1][2]=
		                SCDOSE[IZ-1][IX-1][2]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE[IZ-1][IX-1][3]=
		                SCDOSE[IZ-1][IX-1][3]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);

		                SCDOSE2[IZ-1][IX-1][0]=
		                SCDOSE2[IZ-1][IX-1][0]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE2[IZ-1][IX-1][1]=
		                SCDOSE2[IZ-1][IX-1][1]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE2[IZ-1][IX-1][2]=
		                SCDOSE2[IZ-1][IX-1][2]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		                SCDOSE2[IZ-1][IX-1][3]=
		                SCDOSE2[IZ-1][IX-1][3]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);

		                SCDOSE_COV[IZ-1][IX-1][0]=
		                SCDOSE_COV[IZ-1][IX-1][0]*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS));
		                SCDOSE_COV[IZ-1][IX-1][1]=
		                SCDOSE_COV[IZ-1][IX-1][1]*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS));
		                SCDOSE_COV[IZ-1][IX-1][2]=
		                SCDOSE_COV[IZ-1][IX-1][2]*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS));
		            }
		        }
		    }
		}

		OSUMRY(); //"PRINT THE OUTPUT SUMRY"

		//:END-OF-RUN:;

		//;"******************************************************************************
		//"
		//"                       *** SECTION 4 ***
		//"
		//"------------------------------------------------------------------------------
		//"
		//"THE CONCLUSION"
		//"
		//"------------------------------------------------------------------------------

		//:END:;
		//OUTPUT; (//'END OF RUN',9X,' ',$); call egs_fdate(6);
		//OUTPUT; (//);
		//write(iout,'(/a,$)') 'END OF RUN          '; call egs_fdate(iout);
		//write(iout,'(////)');

		//call egs_finish;

		//#ifdef HAVE_C_COMPILER;
		//;
		//IF( n_parallel > 0 & ~is_finished ) [
		//    call egs_pjob_finish(n_job);
		///    IF( n_job = 0 ) [
		//        is_finished = .true.;
		//        call egs_combine_runs(combine_results,'.egsdat');
		//        NCASET=NCASEO;  IHSTRY=NCASET;
		//        CALL SRCOTO(WEIGHT);
		//        goto :STATS-ANAL:;
		//    ]
		//]
		//#endif;

		EGS4SrcEns.SRCEND();//do nothing

		//$CALL_EXIT(0);

		//"FORTRAN FORMAT STATEMENTS. FORMAT STATEMENT N## IS FIRST USED IN SECTION N."
		//%I0
		//100  FORMAT(80A1//'Calculation using CAVRZnrc(EGSnrc) '$VERSION' ',
		//             /' ON '$MACHINE' ',T55,' ',$);
		//200  FORMAT(//,79('*')/
		//            // ,T20,'EXECUTION INFORMATION AND WARNING MESSAGES'/
		//            // ,79('*')/
		//            //'USING CAVRZnrc(EGSnrc) '$VERSION' ON '$MACHINE' ');
		//201  FORMAT(/'********* NEW INPUT FILE *********'/);
		//202  FORMAT(/'********* RESTARTED INPUT FILE ********* '/
		//             ' ',10X,I12,' NEW + ',I12,' OLD HISTORIES');
		//204  FORMAT(/'********* DATA ANALYSIS ONLY *********'/);
		//205  FORMAT(/'********* RANDOM NUMBERS READ FROM FILE *********'/);
		//206  FORMAT(/' ********* ANALYZING RESULTS FROM PARALLEL RUNS *******'/);
		//210  FORMAT(/'********* NOT ENOUGH TIME TO FINISH WITHIN',
		//       ' LIMIT OF',F8.2,' HOURS',I5,' BATCHES USED********'/
		//       ' ',I12,' HISTORIES RUN, ',I12,' HISTORIES ANALYZED'//);
		//230  FORMAT(/'DESIRED STATISTICAL ACCURACY OBTAINED.'/
		//            ' STATS IN CAVITY= ',F5.2,'%',
		//            ' AFTER ',I2,' BATCHES');
		//240  FORMAT(/'*********DESIRED STATISTICAL ACCURACY OF ',F5.2,'%',
		//            ' NOT REACHED*********'/
		//            ' STATS IN CAVITY= ',F5.2,' % AFTER ',I2,' BATCHES');
		//250  FORMAT(/' FOR OLD RUN:'/
		//             ' ----------- '/
		//             ' Total cputime =',F8.1,'s (=',F5.2,' hr)');
		//255  FORMAT(/' FOR PARALLEL RUNS:'/
		//             ' ----------------- '/
		//             ' On ',I5,' machines '/
		//             ' Total cputime =',F8.1,'s (=',F8.2,' hr), cputime/machine =',
		//              F8.1,'s');
		//260  FORMAT(/'Finished simulations:'/'     time elapsed,cputime',
		//             ',ratio= ',2F8.1,'(=',F5.2,'hr)',F8.2);
		//261  FORMAT(/' Finished: time elapsed this run', F10.1/
		//             '           CPUtime total run    ', F10.1,'(=',F8.2,'hr)'/
		//             '           Ratio ELAPSED/CPU this run:', F8.3);
		//280  FORMAT(/'    CPUtime/history=',F10.5,'  sec.  Histories/hour=',F12.0);
//
		//END; "END OF MAIN ROUTINE-CAVRZnrc"


		EGS4.timeElapsed();
		Calendar call=Calendar.getInstance();
		Date da=call.getTime();
		EGS4.seqStr="End of run:          "+da.toString();
	    if(EGS4.iprint>0)
	    	printSequence(EGS4.seqStr);

//===============================================================
		if(createOutputFile)
		{
			try
			{
				sigfos.close();
			}
			catch (Exception ex){ }

			putInFile=false;
			EGS4.seqStr="Check current directory for results, in text file:"+filename;
			if(EGS4.iprint>0)
    	    	 printSequence(EGS4.seqStr);
		 }
	}

	/**
	 * Where to print runtime information. Interface method.
	 * @param s the String to be printed
	 */
	public void printSequence(String s)
	{
		//write file?
		if(createOutputFile && putInFile)
		{
			try
			{
				sigfos.write(s+" \n");
			}
			catch (Exception ex){ }
		}

		//output data to console
		System.out.println(s);
	}

	//"******************************************************************************
	//"
	//"
	//"                               **********
	//"                               *        *
	//"                               * AUSGAB *
	//"                               *        *
	//"                               **********
	//"
	//"
	//"       An AUSGAB routine to be used with cavrznrc.mortran
	//"
	//"       This routine scores the dose and other ionisation cavity parameters
	//"       in a finite, azimuthally symmetric cylindrical geometry which the
	//"       user defines via plane and radius coordinates. The user must specify
	//"       both the target geometry as well as the planes and radii between
	//"       which the quantities are to be scored. All the geometrical checks for
	//"       crossing 'geometrical' or 'dose' regions are handled by the subroutine
	//"       HOWFAR.
	//"
	//"       FOR IT = 1      the total primary dose
	//"              = 2      the total dose dose less the total primary dose
	//"                       (i.e. the scatter fraction for calculating Ascat)
	//"              = 3      the total primary unattenuated dose less the
	//"                       total primary dose (for calculating Aatt)
	//"              = 4      the total primary unattenuated dose with the source
	//"                       distribution factored out less the total primary
	//"                       unattenuated dose (for calculating Apn)
	//"              = 5      the total primary unattenuated dose to the cavity
	//"                       gas with the source distribution factored out with
	//"                       the electron tranport taking place in the chamber
	//"                       with the cavity filled with wall material (for
	//"                       calculating Afl and the stopping power ratio)
	//"              = 6      the total primary unattenuated dose to the wall
	//"                       material with the source distribution factored out with
	//"                       the electron tranport taking place in the chamber
	//"                       with the cavity filled with wall material (for
	//"                       calculating the stopping power ratio)
	//"
	//"
	//;"******************************************************************************

	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.<br>
	 * This routine scores the dose and other ionisation cavity parameters in a finite, azimuthally symmetric cylindrical geometry which the 
	 * user defines via plane and radius coordinates. The user must specify both the target geometry as well as the planes and radii between 
	 * which the quantities are to be scored. All the geometrical checks for crossing 'geometrical' or 'dose' regions are handled by the subroutine 
	 *  HOWFAR.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG)
	{

		//$IMPLICIT-NONE;

		double FTMP=0.;
		int ip=0;

		//$INTEGER IRL,IX,IZ,IQL,LATCHL,IARG,IDUMMY;
		//$REAL WTL,FDUMMY,xsi;

		//;COMIN/
		//ELECIN,EPCONT,GEOM,MEDIA,PHOTIN,RUSROU,SCORE,SOURCE,STACK,USEFUL,USER,
		// RANDOM,BOUNDS/;
		double xsi=0.0;
		//"STACK OVERFLOW CHECK"
		MXNP=Math.max(MXNP,EGS4.NP);//"keep track of how deep stack is"
		                  //"MXNP is not output but it should be"

		if(IWATCH > 0) {EGS4.WATCH(IARG,IWATCH);} //"signal watch routine if active"]

		//"check if particle is leaving the transport geometry"
		int IRL=EGS4.IR[EGS4.NP-1];  //"local region number"
		if(IRL ==  1) return; //"outside the chamber"

		//"OBTAIN FREQUENTLY USED LOCAL VARIABLES"
		//$GET-IX-IZ(IRL);   //"local plane and radius numbers"
		int IX=EGS4Geom.GET_IX(IRL);
		int IZ=EGS4Geom.GET_IZC(IRL);

		int IQL=EGS4.IQ[EGS4.NP-1];        //"local charge variable"
		double WTL=EGS4.WT[EGS4.NP-1];        //"local weight variable"
		int LATCHL=EGS4.LATCH[EGS4.NP-1];  //"LATCHL=0 for primaries, 1 otherwise"
//#########################################################
		if( EGS4Macro.use_enhance || EGS4Macro.n_split > 1 )
		{
		    if( IARG < 5 )
		    {
		        if( EGS4.EDEP > 0 && WTL > 0 && EGS4Geom.ntrack[IRL-1] == 1 )
		        {
		            //"If we use cross section enhancement or photon splitting,"
		            //"all energy scoring is done here and the rest of ausgab is ignored"
		            //"We use the technique proposed by the PENELOPE group for scoring "
		            //"the energy deposition. This results in a much better estimate   "
		            //"of the uncertainty                                              "

		            FTMP = WTL*EGS4.EDEP;
		            if( EGS4SrcEns.NHSTRY == last_case )
		            {
		                //" Still the same history scoring into the cavity => update    "
		                //" temporary variables                                         "
		                if( LATCHL != 2 ) tmp_dose = tmp_dose + FTMP;
		                if( LATCHL != 3 ) tmp_dose1 = tmp_dose1 + FTMP;
		            }
		            else
		            {
		                //" A new history scoring into the cavity. "
		                last_case = EGS4SrcEns.NHSTRY;
		                cav_dose = cav_dose + tmp_dose;
		                cav2_dose = cav2_dose + tmp_dose*tmp_dose;
		                cav_dose1 = cav_dose1 + tmp_dose1;
		                cav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1;
		                cav_dosec = cav_dosec + tmp_dose*tmp_dose1;
		                if( LATCHL != 2 ) { tmp_dose = FTMP; } else { tmp_dose = 0.; }
		                if( LATCHL != 3 ) { tmp_dose1 = FTMP;} else { tmp_dose1 = 0.; }
		            }
		        }
		        return;
		    }
		}
//#########################################################
		if( EGS4Macro.use_enhance )
		{// "If we use cross section enhancement, all scoring "
		 //                   " is done here and the rest of ausgab is ignored  "

		    if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23)
		    {
		        //"A pair/Compton/photoelectric/Rayleigh event is about to take place"
		        //"As we have increased the photon cross section by a factor of      "
		        //"cs_enhance, we must split the photon into a scattering portion    "
		        //"(1/cs_enhance) and a nor-scattering portion (1-1/cs_enhance)      "
		        //"Start with placing an identical photon on the stack               "
		        EGS4.NP = EGS4.NP + 1;
		        if(EGS4.NP + 1 > EGS4.$MXSTACK)
		        {
				    EGS4.STOPPROGRAM=true;
	      		    EGS4.seqStr=" ***************************************************"+"  \n"+
		  		    " Calculation with CS-enhancement: unable to boost stack."+"  \n"+
		        	" ***************************************************";
	  		    	if(EGS4.iprint>0)
		  			   printSequence(EGS4.seqStr);
					return;
		            //OUTPUT;
		            //( ' Calculation with CS-enhancement: unable to boost stack.'/
		            //    ' Stopping.'/ 1x,80('*')/);
		            //stop;
		        }
		        //$TRANSFER PROPERTIES TO (np) FROM (np - 1);
				EGS4.X[EGS4.NP-1]=EGS4.X[EGS4.NP-2];
				EGS4.Y[EGS4.NP-1]=EGS4.Y[EGS4.NP-2];
				EGS4.Z[EGS4.NP-1]=EGS4.Z[EGS4.NP-2];
				EGS4.IR[EGS4.NP-1]=EGS4.IR[EGS4.NP-2];
				EGS4.WT[EGS4.NP-1]=EGS4.WT[EGS4.NP-2];
				EGS4.DNEAR[EGS4.NP-1]=EGS4.DNEAR[EGS4.NP-2];
				EGS4.LATCH[EGS4.NP-1]=EGS4.LATCH[EGS4.NP-2];

		        EGS4.E[EGS4.NP-1]  =  EGS4.E[EGS4.NP-2];
		        EGS4.U[EGS4.NP-1]  =  EGS4.U[EGS4.NP-2];
		        EGS4.V[EGS4.NP-1]  =  EGS4.V[EGS4.NP-2];
		        EGS4.W[EGS4.NP-1]  =  EGS4.W[EGS4.NP-2];
		        EGS4.IQ[EGS4.NP-1] = EGS4.IQ[EGS4.NP-2];
		        if( EGS4.LATCH[EGS4.NP-2]  != 2 )
		        {
		          //" This is either a primary photon that has not yet been attenuated "
		          //" away or a scattered photon. Let's decide what to do with the     "
		          //" unscattered fraction of that photon (which is at np-1)           "
		            xsi=EGS4.random01();
		            if( EGS4Macro.cs_enhance*xsi < 1. )
		            {//  " The photon doesn't survive. "
		                if( EGS4.LATCH[EGS4.NP-2] == 3 )
		                {// " It is a scattered photon => kill it"
		                    EGS4.WT[EGS4.NP-2] = 0.0;
		                    EGS4.DNEAR[EGS4.NP-2] = -1.;
		                }
		                else
		                {// " This is a primary => mark it as attenuated.         "
		                 //      " From now on, all descendents of this photon will    "
		                 //      " only contribute to the cavity dose with attenuation "
		                 //      " and scatter removed                                 "
		                    EGS4.LATCH[EGS4.NP-2] = 2;
		                }
		            }
		        }
		        //"Adjust the weight of to be scattered photon"
		        EGS4.WT[EGS4.NP-1] = EGS4.WT[EGS4.NP-1]/EGS4Macro.cs_enhance;
		        return;
		    }

		    if( IARG == 18 || IARG == 20 || IARG == 24 ||
		         //" A Compton/photo-absorption/Rayleigh event just occured"
		        IARG == 7 || IARG == 13 || IARG == 14 )
		    {
		         //" A bremas/annihilation event just occured"
		         //" All scattered photons and photons originating in brems/annihilation"
		         //" events contribute to the scattered dose. But because all of them   "
		         //" have a small weight (initial weight/cs_enhance), we will play      "
		         //" Russian Roulette with them, using 1/cs_enhance as a sirvivng       "
		         //" probability. If they survive, their weight will become equal to the"
		         //" intial photon weight. In addition, we have to set their latch to 3 "
		         //" so that they and rheir descendents only contribute to the scattered"
		         //" dose.                                                              "
		        xsi=EGS4.random01(); xsi = xsi*EGS4Macro.cs_enhance;
		        for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)
		        {
		            if( EGS4.IQ[ip-1] == 0 )
		            {
		                if( EGS4.LATCH[ip-1] == 2 )
		                {// "that's a descendent of a photon that"
		                 //                     "has been attenuated away => kill it"
		                    EGS4.WT[ip-1] = 0.0; EGS4.DNEAR[ip-1] = -1.;
		                }
		                else
		                {
		                  if( EGS4.E[ip-1] >= EGS4.PCUT[IRL-1] )
		                  {
		                    if( xsi < 1. )
		                    {
								EGS4.LATCH[ip-1] = 3;
								EGS4.WT[ip-1] = EGS4.WT[ip-1]*EGS4Macro.cs_enhance;
							}
		                    else
		                    {
								EGS4.WT[ip-1] = 0.; EGS4.DNEAR[ip-1] = -1.;
							}
		                  }
		                  else { EGS4.LATCH[ip-1] = 3; }
		                   //" i.e. we don't need the Russian Roulette for photons below"
		                   //" threshold because they will be discarded and their energy"
		                   //" deposited locally anyway                                 "
		                }
		            }
		        }
		        return;
		    }
		    return;
		}
//#########################################################

		if( EGS4Macro.n_split > 1 )
		{

		    if( IARG == 7 || IARG == 13 || IARG == 14 )
		    {
		        if( EGS4Macro.iifano == 1 )
		        {
		            for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)
		            {
		                if( EGS4.IQ[ip-1] == 0 )
		                { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; }
		            }
		            return;
		        }
		        for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [
		        {
		            if( EGS4.IQ[ip-1] == 0 )
		            {
		                xsi=EGS4.random01();
		                if( xsi*EGS4Macro.n_split > 1)
		                { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; }
		                else
		                { EGS4.LATCH[ip-1] = 3;EGS4.WT[ip-1]=EGS4.WT[ip-1]*EGS4Macro.n_split; }
		            }
		        }
		        return;
		    }
		}
//#########################################################


//REPLACE {$SCORE(#,#:#)} WITH {;
//"Scoring macro used in AUSGAB for quantities other than DOSE and KERMA"
//"{P1}{P2}=scoring array (eg SCSTP)"
//"{P3}=quantity to be scored (eg 1)"

//"If the (primary) history number, NHSTRY, is the same as the history"
//"that last scored in this array, {P1}_LAST{P2}, then {P3} is added"
//"to a temporary array, {P1}_TMP{P2}.  Otherwise, we add"
//"{P1}_TMP{P2} to {P1}{P2}, {P1}_TMP{P2}*{P1}_TMP{P2} to {P1}2{P2},"
//"set {P1}_TMP{P2}={P3}, and set {P1}_LAST{P2}=NHSTRY."
//"This scoring method allows us to calculate  uncorrelated value"
//"of {P1}2{P2} which is then used to calculate the uncertainty"
//"in {P1}{P2}.  In cavrznrc, this macro is only used for counting"
//"no. of charged particle steps."

//IF(NHSTRY={P1}_LAST{P2})[
//  {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
//]
//ELSE[
//  {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
//  {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
//  {P1}_TMP{P2}={P3};
//  {P1}_LAST{P2}=NHSTRY;
//]
//;
//}


		if(IARG == 0)
		{//"ABOUT TO TRANSPORT A PARTICLE"
		    if(IQL!=0)
		    {
		        if(LATCHL == 0)
		        {// "COUNT PRIMARY CHARGED PARTICLES ONLY"
		                //$SCORE(SCSTP, :1);"COUNT CHARGED PARTICLE STEPS TAKEN"
					if(EGS4SrcEns.NHSTRY==SCSTP_LAST)//IF(NHSTRY={P1}_LAST{P2})[
					{
					//  {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
					  SCSTP_TMP=SCSTP_TMP + 1.;
					}
					else
					{
					//  {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
					  SCSTP=SCSTP+SCSTP_TMP;
					//  {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
					  SCSTP2=SCSTP2 + SCSTP_TMP*SCSTP_TMP;
					//  {P1}_TMP{P2}={P3};
					  SCSTP_TMP=1.;
					//  {P1}_LAST{P2}=NHSTRY;
					  SCSTP_LAST=EGS4SrcEns.NHSTRY;
					}

		            if(EGS4Geom.ntrack[IRL-1] == 1)
		            {
						//$SCORE(SCCSTP, :1);
						if(EGS4SrcEns.NHSTRY==SCCSTP_LAST)//IF(NHSTRY={P1}_LAST{P2})[
						{
						//  {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
						  SCCSTP_TMP=SCCSTP_TMP + 1.;
						}
						else
						{
						//  {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
						  SCCSTP=SCCSTP+SCCSTP_TMP;
						//  {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
						  SCCSTP2=SCCSTP2 + SCCSTP_TMP*SCCSTP_TMP;
						//  {P1}_TMP{P2}={P3};
						  SCCSTP_TMP=1.;
						//  {P1}_LAST{P2}=NHSTRY;
						  SCCSTP_LAST=EGS4SrcEns.NHSTRY;
						}

				        //"WRITE(*,*)' sccstp ',SCCSTP;"
		             }
		        }//"COUNT STEPS IN CAVITY REGION"
		    }
		    else
		    {//"PHOTON STEP - PLAY RUSSIAN ROULETTE?"
		        if(RUSROU&&(EGS4.W[EGS4.NP-1]>0.0))
		        {//"YES, PLAY IF CROSSES RRZ "
		            if((EGS4.Z[EGS4.NP-1]<=RRZ)&&
		            (EGS4.Z[EGS4.NP-1]+EGS4.USTEP*EGS4.W[EGS4.NP-1]>=RRZ))
		            {//"CROSSES"
		                xsi=EGS4.random01();
		                if(xsi<RRCUT)
		                {//"PARTICLE SURVIVES"
		                    EGS4.WT[EGS4.NP-1]=WTL/RRCUT;
						}
		                else
		                {//"DISCARD PARTICLE ON NEXT CALL TO HOWFAR"
		                    EGS4.WT[EGS4.NP-1]=0.0;
					    }
		             } //"END TEST IF CROSSES RUSSIAN ROULETTE PLANE"
		        } //"END TEST FOR PLAYING RUSSIAN ROULETTE"
		    }//"END TEST FOR PHOTON STEP"
	  }// "END TEST FOR IARG = 0"

		if (IFANO == 1)
		{
		    if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23)
		    {
		        //"A pair/Compton/photoelectric/Rayleigh event is about to take place"
		        EGS4.NP = EGS4.NP + 1; //"Boost the stack"
		        if(EGS4.NP + 1 > EGS4.$MXSTACK)
		        {
				    EGS4.STOPPROGRAM=true;
	      		    EGS4.seqStr=" ***************************************************"+"  \n"+
		  		    " Fano calculation unable to boost stack."+"  \n"+
		        	" ***************************************************";
	  		    	if(EGS4.iprint>0)
		  			   printSequence(EGS4.seqStr);
					return;
		        }
		        //"Create an identical photon"
		        //$TRANSFER PROPERTIES TO (np) FROM (np - 1);
				EGS4.X[EGS4.NP-1]=EGS4.X[EGS4.NP-2];
				EGS4.Y[EGS4.NP-1]=EGS4.Y[EGS4.NP-2];
				EGS4.Z[EGS4.NP-1]=EGS4.Z[EGS4.NP-2];
				EGS4.IR[EGS4.NP-1]=EGS4.IR[EGS4.NP-2];
				EGS4.WT[EGS4.NP-1]=EGS4.WT[EGS4.NP-2];
				EGS4.DNEAR[EGS4.NP-1]=EGS4.DNEAR[EGS4.NP-2];
				EGS4.LATCH[EGS4.NP-1]=EGS4.LATCH[EGS4.NP-2];

		        EGS4.E[EGS4.NP-1]  =  EGS4.E[EGS4.NP-2];
		        EGS4.U[EGS4.NP-1]  =  EGS4.U[EGS4.NP-2];
		        EGS4.V[EGS4.NP-1]  =  EGS4.V[EGS4.NP-2];
		        EGS4.W[EGS4.NP-1]  =  EGS4.W[EGS4.NP-2];
		        EGS4.IQ[EGS4.NP-1] = EGS4.IQ[EGS4.NP-2];

		        return;
		    }

		    //"Throw away any scattered photons from the primary interaction site."
		    //" Now there is a stack pointer NPold which points to the particle "
		    //" befor the last discrete interaction. This change was necessary "
		    //" for the implementation of atomic relaxations. So, check all particles"
		    //" between NPold and NP and discard photons "
		    //"IF ( (iarg = 18 & NP > NPold)" " Compton has occured"
		    if (  IARG == 18//               " Compton has occured"
		        || IARG == 20 //              " After photo-absorption "
		        || IARG == 24 //              " After Rayleigh "
		        || IARG == 7  //              " After brems "
		        || IARG == 13 //              " After annihilation "
		        || IARG == 14)//              " After annihilation at rest "
		    {
		        for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)
		        {
		            if( EGS4.IQ[ip-1] == 0 )
		            { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; }
		        }
		    }
		}
		if (IFANO == 2)
		{
		    if (  IARG == 16   //            " After pair production "
		        || IARG == 18  //             " After Compton "
		        || IARG == 20) //             " After photo absorption "
		    {
		       if ( EGS4Geom.ntrack[EGS4.IR[EGS4.NP-1]-1] == 0)
		       {
		          for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP
		          {
		             if( EGS4.IQ[ip-1] == 0 )//if( iq(ip) ~= 0 )
		             { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; }//{ wt(ip) = 0; e(ip) = 0; }
		          }
		       }
		    }
		}

		//"SCORE THE ENERGY AS REQUIRED FOR THE DIFFERENT MODES"
		//"****************************************************"

		if((EGS4.EDEP!=0.0)&&
		(WTL>0.0)&&(IARG<5)&&(EGS4Geom.ntrack[IRL-1]==1))
		{
		 //"ENERGY HAS BEEN DEPOSITED IN THE CAVITY REGION"
		 //"SCORE PRIMARY AND SECONDARY ENERGY DEPOSITED"
		    FTMP=EGS4.WT[EGS4.NP-1]*EGS4.EDEP;

		//"***************************************************************************"
		//"                                                                           "
		//" Implementation of a history by history scoring scheme for the cavity dose "
		//" taken from Iwan's splitting implementation, EMH, March 2002               "
		//"                                                                           "
		//" For the meaning of the variables see the COMIN/SCORE block                "
		//"                                                                           "
		//"         expmfp = exp(di) - 1 (see macro $SELECT-PHOTON-MFP)               "
		//"                                                                           "
		//"***************************************************************************"
		    if( EGS4SrcEns.NHSTRY == last_case )
		    {
		            //" Still the same history scoring into the cavity => update    "
		            //" temporary variables                                         "
		            tmp_dose = tmp_dose + FTMP;
		            if( LATCHL == 0 )
		            {
		              tmp_dose0 = tmp_dose0 + FTMP;
		              tmp_dose1 = tmp_dose1 + FTMP*(EGS4Macro.EXPMFP+1);
		            }
		            else
		            {
		              tmp_dose2 = tmp_dose2 + FTMP;
		            }
		    }
		    else
		    {
		            //" A new history scoring into the cavity. "
		            last_case = EGS4SrcEns.NHSTRY;

		            cav_dose  = cav_dose + tmp_dose;
		            cav2_dose = cav2_dose + tmp_dose*tmp_dose;

		            cav_dose0  = cav_dose0 + tmp_dose0;
		            cav2_dose0 = cav2_dose0 + tmp_dose0*tmp_dose0;

		            cav_dose1  = cav_dose1 + tmp_dose1;
		            cav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1;

		            cav_dose2  = cav_dose2 + tmp_dose2;
		            cav2_dose2 = cav2_dose2 + tmp_dose2*tmp_dose2;

		            cav_dosec = cav_dosec + tmp_dose*tmp_dose1;
		            cav_dosec01 = cav_dosec01 + tmp_dose0*tmp_dose1;
		            cav_dosec02 = cav_dosec02 + tmp_dose0*tmp_dose2;

		            tmp_dose = FTMP;

		            if( LATCHL == 0 )
		            {
		                tmp_dose0 = FTMP ;
		                tmp_dose1 = FTMP*(EGS4Macro.EXPMFP+1) ;
		                tmp_dose2 = 0.0;
		            }
		            else
		            {
		                tmp_dose0 = 0.0;
		                tmp_dose1 = 0.0;
		                tmp_dose2 = FTMP;
		            }
			}

		  IDECAV=1;

		  if(EGS4Geom.NSUMCV>1)
		  {//"calculate quantities in individual cavity regions"
		   //            "do it the same as for the overall cavity above"

		        if(EGS4SrcEns.NHSTRY==SCDOSE_LAST[IZ-1][IX-1])
		        {
		             SCDOSE_TMP[IZ-1][IX-1][0]=SCDOSE_TMP[IZ-1][IX-1][0]+FTMP;
		             if(LATCHL==0)
		             {//"primary dose"
		                SCDOSE_TMP[IZ-1][IX-1][1]=SCDOSE_TMP[IZ-1][IX-1][1]+FTMP;
		                SCDOSE_TMP[IZ-1][IX-1][2]=SCDOSE_TMP[IZ-1][IX-1][2]+FTMP*(1+EGS4Macro.EXPMFP);
		             }
		             else
		             {//"secondary dose"
		                SCDOSE_TMP[IZ-1][IX-1][3]=SCDOSE_TMP[IZ-1][IX-1][3]+FTMP;
		             }
		        }
		        else
		        {
		             SCDOSE_LAST[IZ-1][IX-1]=EGS4SrcEns.NHSTRY;

		             SCDOSE[IZ-1][IX-1][0]=SCDOSE[IZ-1][IX-1][0]+SCDOSE_TMP[IZ-1][IX-1][0];
		             SCDOSE2[IZ-1][IX-1][0]=SCDOSE2[IZ-1][IX-1][0]+
		             	SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][0];

		             SCDOSE[IZ-1][IX-1][1]=SCDOSE[IZ-1][IX-1][1]+SCDOSE_TMP[IZ-1][IX-1][1];
		             SCDOSE2[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]+
		             	SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][1];

		             SCDOSE[IZ-1][IX-1][2]=SCDOSE[IZ-1][IX-1][2]+SCDOSE_TMP[IZ-1][IX-1][2];
		             SCDOSE2[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX-1][2]+
		             	SCDOSE_TMP[IZ-1][IX-1][2]*SCDOSE_TMP[IZ-1][IX-1][2];

		             SCDOSE[IZ-1][IX-1][3]=SCDOSE[IZ-1][IX-1][3]+SCDOSE_TMP[IZ-1][IX-1][3];
		             SCDOSE2[IZ-1][IX-1][3]=SCDOSE2[IZ-1][IX-1][3]+
		             	SCDOSE_TMP[IZ-1][IX-1][3]*SCDOSE_TMP[IZ-1][IX-1][3];

		             SCDOSE_COV[IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]+
		             	SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][2];
		             SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]+
		             	SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][2];
		             SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]+
		             	SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][3];

		             SCDOSE_TMP[IZ-1][IX-1][0]=FTMP;

		             if(LATCHL==0)
		             {
		                SCDOSE_TMP[IZ-1][IX-1][1]=FTMP;
		                SCDOSE_TMP[IZ-1][IX-1][2]=FTMP*(1+EGS4Macro.EXPMFP);
		                SCDOSE_TMP[IZ-1][IX-1][3]= 0.0;
		             }
		             else
		             {
		                SCDOSE_TMP[IZ-1][IX-1][1]=0.0;
		                SCDOSE_TMP[IZ-1][IX-1][2]=0.0;
		                SCDOSE_TMP[IZ-1][IX-1][3]=FTMP;
		             }
		        }

		  }

		} //"END OF ENERGY DEPOSITED IN THE CAVITY"

		//"SET FLAG FOR SECONDARY INTERACTIONS"
		//"***********************************"

		if((EGS4Macro.IFULL>0)&&(IARG>5)&&(LATCHL == 0))
		{
		//"ONLY IF PRIMARY PARTICLES HAVE INTERACTED DISCRETELY                 "
		//"IF A SECONDARY PARTICLE IS CREATED ON THE SECOND PASS, GIVE IT A ZERO"
		//"WEIGHT SO THAT HOWFAR WILL DISCARD IT.                               "

		    if( IARG == 7 )
		    {//  "brem has occured"

		        if( IQL == 0 ) { EXCHANGE_STACK(EGS4.NP,EGS4.NP-1); }
		        EGS4.LATCH[EGS4.NP-2] = 1; //" Flag the photon as a secondary"
		        if( ipass >= 1 ) { EGS4.WT[EGS4.NP-2] = 0.; }// "To save time in correlation runs"

		    }
		    else if( IARG == 18 )
		    {// "Compton has occured, with binding effects"
		    //"taken into account, 0, 1, or more particles" "may have resulted"

		            for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [
		            {
		                if( EGS4.IQ[ip-1] == 0 )
		                {
		                    EGS4.LATCH[ip-1] = 1;
		                    if( ipass >=1 ) { EGS4.WT[ip-1] = 0.; }// "To save time"
		                }
		            }
		        //" NP = NPold means the interactiopn has been rejected and thus "
		        //" the emerging (unscattered) photon  is still a primary "
		    }
		    else if( IARG == 9 )
		    {// "Moller has occured. For now there is only"
		     //  "one secondary. When impact ionization is implemented"
		     //                     "the following should be changed"
		        if( EGS4.E[EGS4.NP-1] < EGS4.E[EGS4.NP-2] ) { EXCHANGE_STACK(EGS4.NP,EGS4.NP-1); }
		    }
		    else if( IARG == 13 || IARG == 14 )
		    {//"Annihilation, flag the photons"
		        EGS4.LATCH[EGS4.NP-1] = 1;
		        EGS4.LATCH[EGS4.NP-2] = 1;
		        if( ipass >= 1 ) { EGS4.WT[EGS4.NP-1] = 0.; EGS4.WT[EGS4.NP-2] = 0.; }
		    }
		    else if (IARG == 20)
		    {
		        for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [
		        {
		            if( EGS4.IQ[ip-1] == 0 )
		            {
		               EGS4.LATCH[ip-1] = 1;
		               if( ipass >= 1 ) { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; }
		            }
		        }
		    }
		    else if( IARG == 24 )
		    {
		        EGS4.LATCH[EGS4.NP-1] = 1;
		        if( ipass >= 1 ) { EGS4.WT[EGS4.NP-1] = 0.; }
		    }
		}

		return;
	}//"END OF AUSGAB"

	//"The following is the $CALL-HOWNEAR macro for PRESTA-II

	//REPLACE {$CALL-HOWNEAR(#);} WITH {
	//    ;
	//    "write(6,'(2i3,4e15.8)') np,ir(np),dnear(np), "
	//    "   sqrt(x(np)*x(np)+y(np)*y(np)),z(np),tustep; "
	//    IF( dnear(np) < tustep ) [
	//        call hownear({P1},x(np),y(np),z(np),ir(np));
	//        " write(6,*) ' --> new dnear: ',{P1}; "
	//    ]
	//    ELSE [ {P1} = dnear(np); ]
	//}
	//"*********************************************************************"
	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR()
	{
		//"Subroutine arguments
		//$REAL
		//	    tperp, "nearest distance to any boundary (output)
		//	    tustep,
		//	    x,     "x-position of the particle (input)
		//	    y,     "y-position of the particle (input)
		//	    z;     "z-position of the particle (input)

		//$INTEGER
		//	    ir     "region number of the particle

		//"Local variables
		double r=0.0;
		int ix=0;// "current cylindrical radius number
		int iz=0;//  "current planar slab number

		double z=EGS4.Z[EGS4.NP-1];
		double y=EGS4.Y[EGS4.NP-1];
		double x=EGS4.X[EGS4.NP-1];
		int ir=EGS4.IR[EGS4.NP-1];
	//FROM MACRO=>ELSE [ {P1} = dnear(np); ]->not necessary, kind of varred(skip calc.)!@@@@@@@@@@@
		if( EGS4.DNEAR[EGS4.NP-1] >= EGS4.TUSTEP)
		{
			EGS4.tperp=EGS4.DNEAR[EGS4.NP-1];
			return;
		}
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		ix = (ir - 2)/EGS4Geom.NZ + 1;//NZ=planar zones=NPLANE-1!!!NR=nr of Radii
		iz = ir - 1 - EGS4Geom.NZ*(ix - 1);
		r = Math.sqrt(x*x + y*y);
		EGS4.tperp = EGS4.min(z - EGS4Geom.ZPLANE[iz-1],
							  EGS4Geom.ZPLANE[iz] - z,
							  EGS4Geom.RCYL[ix] - r);
		if(ix != 1)
		{
			EGS4.tperp = Math.min(EGS4.tperp,r - EGS4Geom.RCYL[ix - 1]);
		}

			//"IF( tperp < -1e-6 ) [
			//"  OUTPUT IHSTRY, tperp; (/' Error in HOWNEAR: IHSTRY=',I12,
			//"                      '  tperp negative =',1PE10.3);
			//"  OUTPUT ir,x,y,z,r; (' ir x y z r=', I5, 4F15.7);
			//"  OUTPUT iz,ix;(' depth,radial regions iz,ix=',2I10);
			//"  OUTPUT zplane(iz),zplane(iz+1);( ' upper/lower z planes=',2F15.7);
			//"  IF( ix > 1) [
			//"    OUTPUT rcyl(ix-1),rcyl(ix); (' Radial boundaries=',2F15.7);
			//"  ]
			//"  ELSE [
			//"    OUTPUT rcyl(ix); (' Radial boundary of 1st region =', F15.7);
			//"  ]
			//"]

		return;


	}//"end of subroutine HOWNEAR"

	//;"******************************************************************************
	//"
	//"                               **********
	//"                               *        *
	//"                               * HOWFAR *
	//"                               *        *
	//"                               **********
	//"
	//"       A GENERAL PURPOSE CYLINDRICAL GEOMETRY ROUTINE FOR USE WITH THE EGS4
	//"       CODE SYSTEM ADAPTED FOR USE WITH CAVRZnrc.
	//"
	//"       FOR PARTICLE NP ON THE STACK IN REGION IR(NP), THIS ROUTINE
	//"       DETERMINES IF THE PARTICLE CAN GO A DISTANCE USTEP WITHOUT CHANGING
	//"       ZONES. IF USTEP CAUSES A ZONE CROSSING, IT IS REDUCED TO PLACE IT ON
	//"       THE BOUNDRY AND IRNEW IS SET TO THE ZONE NUMBER ON THE FAR SIDE OF
	//"       THE BOUNDARY. IF IR(NP) IS 1 THEN THE PARTICLE HAS ESCAPED THE REGION
	//"       OF INTEREST AND THE HISTORY IS TERMINATED.(IDISC IS SET TO 1.)
	//"
	//"
	//"
	//"
	//"       SOME VARIABLES
	//"       ==============
	//;"
	//"OUTEND =       .TRUE.  =>      PARTICLE MAY TRANSMIT OR BACKSCATTER OUT ENDS
	//"       =       .FALSE. =>      PARTICLE STAYS WITHIN THE END BOUNDARIES
	//"OUTSID =       .TRUE.  =>      PARTICLE MAY TRANSMIT OUT THE SIDES
	//"       =       .FALSE. =>      PARTICLE STAYS WITHIN THE SIDE BOUNDARY
	//"IRL    =       STARTING REGION NUMBER THE PARTICLE IS IN
	//"IZ     =       STARTING PLANAR ZONE NUMBER THE PARTICLE IS IN.
	//"               THE PARTICLE IS BETWEEN ZPLANE(IZ) AND ZPLANE(IZ+1).
	//"IX     =       STARTING CYLINDRICAL ZONE NUMBER THE PARTICLE IS IN.
	//"               THE PARTICLE IS BETWEEN RCYL(IX-1) AND RCYL(IX).
	//"
	//"       COMMON/GEOM/
	//"               ZPLANE(IZ)      Z VALUES OF PLANES
	//"                               1<=IZ<=NZ+1
	//"               RCYL(IRR)       RADII OF CYLINDERS
	//"                               1<=IRR<=NR
	//"               CYRAD2(IRR)     =RCYL(IRR)**2
	//"               NZ              # PLANAR GEOMETRICAL ZONES (NPLANE-1)
	//"                               ZONE(I) IS BETWEEN ZPLANE(I) AND ZPLANE(I+1)
	//"               NR              # CYLINDRICAL GEOMETRICAL ZONES
	//"                               ZONE(I) IS BETWEEN RCYL(I-1) AND RCYL(I)
	//"               NREG            TOTAL # GEOMETRICAL ZONES =NR*NZ +1
	//"                                       +1 FOR VACUUM ENVELOPE
	//%E     "cavrznrc.mortran"
	//"       DEFINITIONS OF REGION NUMBER, PLANAR ZONE, CYLINDRICAL ZONE
	//"       ===========================================================
	//"               Z AXIS RUNS ACROSS PAGE SHOWN AS .......
	//"
	//"
	//"                                       1
	///"
	//"
	//"       --------------------------------------------------------- RCYL(NR)
	//"       |(NR-1) |(NR-1) |(NR-1) |    . . . .    | NR*NZ | NR*NZ |    IX=NR
	//"       | *NZ+2 | *NZ+3 | *NZ+4 |               |       |   +1  |
	//"       --------------------------------------------------------- RCYL(NR-1)
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       --------------------------------------------------------- RCYL(2)
	//"       |  NZ+2 |  NZ+3 |  NZ+4 |    . . . .    |  2NZ  | 2NZ+1 |    IX=2
	//"       --------------------------------------------------------- RCYL(1)
	//"..1....|...2...|...3...|...4...|...............|...NZ..|..NZ+1.|....IX=1..1..
	//;"       ---------------------------------------------------------
	//"       |       |       |       |    . . . .    |       |       |
	//"       ---------------------------------------------------------
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       |   .   |   .   |   .   |               |   .   |   .   |
	//"       ---------------------------------------------------------
	//"       |       |       |       |    . . . .    |       |       |
	//"       |       |       |       |               |       |       |
	//"       ---------------------------------------------------------
	//"         IZ=1    IZ=2    IZ=3                   IZ=NZ-1  IZ=NZ
	//"
	//"                                       1
	//"
	//"
	//"
	//"       VERSION 1       ADAPTED FROM CAVITY HOWFAR          06/84  ERIC FOX
	//"       VERSION 2       THE SUBROUTINE CALLS TO PLANES AND  10/87  AFB
	//"                       CYLINDER HAVE BEEN REPLACED BY MACROS
	//"                       TO SPEED THINGS UP
	//"
	//"
	//"******************************************************************************

	/**
	 * The following is a general specification of HOWFAR: 
	 * Given a particle at (X,Y,Z) in region IR and going in direction 
	 * (U,V,W), this routine answers the question, can the particle go 
	 * a distance USTEP without crossing a boundary? If yes, it merely returns; 
	 * If no, it sets USTEP=distance to boundary in the current 
	 * direction and sets IRNEW to the region number on the far side of the boundary. 
	 * Interface method.
	 */
	public void HOWFAR()
	{
		//$IMPLICIT-NONE;

		//"MACRO USED LOCALLY TO CHANGE REGIONS, ADJUST USTEP, AND EXIT"
		//REPLACE {$SET NEW REGION(#,#);} WITH
		//{
		//IF({P1}.LE.USTEP)[USTEP={P1};IRNEW={P2};]RETURN;
		//}

		//"Debug macro, usually not used, to verify that region is correct on entry"
		//"    It must be used after IX and IZ have been determined"
		//"    REAL Rdebug must be declared"
		//REPLACE {$Check_region;} WITH
		//{
		//Rdebug = SQRT(X(NP)**2 + Y(NP)**2);
		//IF(IX = 1)["we are in inner radial region"
		//   IF(Rdebug > (RCYL(IX)+1.E-6) | (Z(NP) < ZPLANE(IZ) & Z(NP) > ZPLANE(IZ+1)))[
		//      OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX),
		//             ZPLANE(IZ),ZPLANE(IZ+1);
		//      (/' Error in HOWFAR, particle in wrong region'/
		//       '      IRL, IX,IZ =', 3I7/
		//       '      X,Y,Z,R=', 4F15.7/
		//       '      RCYL(1)=', F15.7/
		//       '      ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		//   ]
		//]"end of IX=1 block"
		//ELSE [
		//   IF((Rdebug > (RCYL(IX)+1E-6) & Rdebug < (RCYL(IX-1) -1.E-6)) |
		//                   (Z(NP) < ZPLANE(IZ) & Z(NP) > ZPLANE(IZ+1)))[
		//      OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX-1), RCYL(IX),
		//             ZPLANE(IZ),ZPLANE(IZ+1);
		//      (/' Error in HOWFAR, particle in wrong region'/
		//       '      IRL, IX,IZ =', 3I7/
		//       '      X,Y,Z,R=', 4F15.7/
		//       '      RCYL(IX-1), RCYL(IX)=', 2F15.7/
		//       '      ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		//   ]
		//]
		//;}

		//LOGICAL OUTEND,OUTSID;

		//$INTEGER IRL,IX,IZ,IHITP,IHITC,IZNEW,IXNEW;
		//$REAL WL,TPLANE,U1,V1,A,TCYL,X1,Y1,B,B2,C,COUT,CIN,RAD;
		//"$REAL Rdebug; " "Uncomment if using $Check_region;"

		//"First set idisc and irnew "
		EGS4.IDISC = 0; EGS4.IRNEW = EGS4.IR[EGS4.NP-1];

		//"DISCARD ZERO WEIGHT PARTICLES"
		if(EGS4.WT[EGS4.NP-1] == 0.0){EGS4.IDISC=1;return;}

		//"INITIALLY ASSUME PARTICLE STAYS IN THE TARGET"
		//boolean OUTEND=false;
		//boolean OUTSID=false;

		int IRL=EGS4.IR[EGS4.NP-1];// "LOCAL REGION NUMBER"

		//"DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY"
		if(IRL == 1){EGS4.IDISC=1;return;}

		//$GET-IX-IZ(IRL); //"GET PLANAR AND CYLINDRICAL ZONES NUMBERS"
        int IX=EGS4Geom.GET_IX(IRL);
        int IZ=EGS4Geom.GET_IZC(IRL);

		//"Following commented out usually"
		//"$Check_region;"
		//"write(6,*);"
		//"write(6,*) 'howfar: ',ix,iz,ustep,nz,nr;"

		//$PLANES(IZ,IZ+1,IHITP,TPLANE,ustep);"GET DISTANCE TO PLANE"
		//        "IHITP  =  1 => HITS GREATER Z PLANE"
		//        "       =  0 => MISSES BOTH PLANES"
		//        "       = -1 => HITS LESSER Z PLANE"
		EGS4Geom.ustep=EGS4.USTEP;
		EGS4Geom.PLANES(IZ,IZ+1);
		EGS4.USTEP=EGS4Geom.ustep;
		//$CYLNDR(IX,IHITC,TCYL,ustep);"GET DISTANCE TO CYLINDER"
		//"       IHITC   =  1 => HITS OUTER CYLINDER"
		//"               =  0 => MISSES BOTH CYLINDERS"
		//"               = -1 => HITS INNER CYLINDER"
		EGS4Geom.ustep=EGS4.USTEP;
		EGS4Geom.CYLNDR(IX);
		EGS4.USTEP=EGS4Geom.ustep;

		return;
	}//"END OF SUBROUTINE HOWFAR"

	/**
	 * Gather all media data required for this simulation.
	 */
	private void HATCH()
	{
		Calendar cal=Calendar.getInstance();
		Date d=cal.getTime();
		EGS4.seqStr="Start of run:         "+d.toString();
	    if(EGS4.iprint>0)
	    	printSequence(EGS4.seqStr);

		EGS4.HATCH();
	}

	/**
	 * Start the shower, i.e. the actual simulation for electron-photon transport
	 */
	private void SHOWER()
	{
		////CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
		EGS4Core.SHOWER(EGS4SrcEns.iqin,EGS4SrcEns.ein,
		EGS4SrcEns.xin,EGS4SrcEns.yin,EGS4SrcEns.zin,
		EGS4SrcEns.uin,EGS4SrcEns.vin,EGS4SrcEns.win,
		EGS4SrcEns.irin,EGS4SrcEns.WEIGHT);
	}

	/**
	 * Setup input variables.
	 */
	private void inputs()
	{
		//@ TITLE
		TITLEs="cavrznrc_template: 1.25 MeV on graphite pancake chamber";
		//@ IWATCH
		//#off,interactions,steps,deposited,graph;
        //#debug output with increasing detail, graph outputs .gph file for EGS_Windows
        //#if not "off" use very few histories
		IWATCH=IWATCH_OFF;
		//@ STORE INITIAL RANDOM NUMBERS
		//#no,last,all deposited,all;
		//#last: store initial random numbers for last history in .egsrns
		//#all deposited: store initial random numbers for each history that deposits energy
		//     in the cavity
		//#all: store initial random numbers for each history
		ISTORE=ISTORE_NO;
		//@ IRESTART
		//#first,restart,make,analyze,for graphics,parallel;
		//#first: first run
		//#restart: restart of old run (requires .egsdat file)
		//#make: just create an input file and exit
		//#analyze: read in data from .egsdat file and do statistical analysis and output results
		//#         to .egslst
		//#for graphics: read starting random numbers from .egsrns--eg for output to graphics package
		//#parallel: read .egsdat files from parallel jobs (named inputfile_w*), do statistical
        //#          analysis and output to .egslst
		IRESTART=IRESTART_FIRST;
		//@ OUTPUT OPTIONS
		//#short,cavity details;
		//#short: output cavity summary + dose grid
		//#cavity details: above plus details for every zone in cavity
		IOOPTN=IOOPTN_SHORT;
		//@ STORE DATA ARRAYS
		//#yes,no;
		//#yes: output .egsdat file for restarts, parallel post-processing, etc
		IDAT=IDAT_YES;
		//@ NUMBER OF HISTORIES
		//#splits into $STAT statistical batches
		//#must be >=$STAT**2 if IWATCH= Off
		//#can have less than this if IWATCH set to another option
		NCASE=10000;
		//@ MAX CPU HOURS ALLOWED
		//#Will shut down cleanly prior to exceeding this limit, as long as one #batch has completed.
		TIMMAX=90.000;
		//@ IFULL
		//#dose and stoppers,Aatt and Ascat,Ap,Afl and <s>g/w;
		//#dose and stoppers: output total dose plus that due to stoppers and discards
		//#Aatt and Ascat: above plus Aatt, Ascat
		//#Ap: above plus Ap
		//#Afl and <s>g/w: above plus Afl and stopping power ratio gas/water
		EGS4Macro.IFULL=IFULL_AATT_AND_ASCAT;
		//@ STATISTICAL ACCURACY SOUGHT
		//#If 0, goes until number of histories or CPU limit exceeded.
		//#If not zero goes until this uncertainty (in %) is achieved in the peak dose region
		STATLM=0.0;
		//@ PHOTON REGENERATION
		//#no,yes,no electrons from wall;
		//#no: normal calculation
		//#yes: regenerate parent photon after interaction  (used for FANO calculations)
		//#no electrons from wall: photons not regenerated,
		//#							secondary electrons from cavity wall are eliminated
		IFANO=IFANO_NO;
		//@ INITIAL RANDOM NO. SEEDS
		//#With ranmar: these must be between 1 and 30081 (default to 9373)
		//#With ranlux: 1st is luxury level 0->4 allowed but should not be 0
        //#             2nd is seed 1 -> 1073741824
		jrng1=1;jrng2=3;
		//right here, initialize random generator!!!
		init_random_generator();
//###########################################################################################
		//@ METHOD OF INPUT
		//#groups,individual,cavity information:
		//#group: input groups of slabs of equal thickness
		//#individual: input Z of bottom of every slab
		//#cavity information: generate simple geometry from cavity info input in section below.
		//#   If you use this, there are no more inputs in this section
		EGS4Geom.iterseindex=EGS4Geom.iINDIVIDUAL;
		if(EGS4Geom.iterseindex==EGS4Geom.iCAVITY_INFORMATION)//2
		{
			EGS4Geom.WALLTH=0.5;//WALL THICKNESS
			EGS4Geom.CAVRAD=1.0;//CAVITY OUTER RADIUS
			EGS4Geom.CAVLNG=2.0;//CAVITY LENGTH
			EGS4Geom.ELERAD=0.01;//ELECTRODE RADIUS
			EGS4Geom.SLENGHT="H2O_fortran";//WALL MATERIAL
			EGS4Geom.airs="AIR521ICRU_fortran";//ELERAD=0.0
			EGS4Geom.electrods="AL521ICRU_fortran";//ELERAD!=0.0;
		}
		else////ITERSE->0 sau 1
		{
			EGS4Geom.Z_OF_FRONT_FACE=0.0;//#Beginning of first slab
			if(EGS4Geom.iterseindex==EGS4Geom.iGROUPS)//0
			{
				//NSLAB
				EGS4Geom.nNSLAB=2;//2 value
			//#Define a group of 10 slabs with thickness 1 cm
			//#followed by 10 slabs with thickness 2 cm
				EGS4Geom.NSLAB[0]=10;EGS4Geom.NSLAB[1]=10;
				//SLAB THICKNESS
				EGS4Geom.DELTAZ[0]=1.0;EGS4Geom.DELTAZ[1]=2.0;
			}

			if(EGS4Geom.iterseindex==EGS4Geom.iINDIVIDUAL)//1
			{
				//NSLAB
				EGS4Geom.nNSLAB=3;
				//DEPTH BOUNDARIES
				EGS4Geom.ZPLANE[1]=0.3;
				EGS4Geom.ZPLANE[2]=0.5;
				EGS4Geom.ZPLANE[3]=0.8;
			}

			EGS4Geom.nCyl=2;//"number of radial cylinders input"
			//#Radii of cylinders
			EGS4Geom.RCYL[1]=1.0;EGS4Geom.RCYL[2]=1.3;

			//MEDIA=the media in the problem. These must match exactly, including case, one
			//of the media names in the pegs4 data set being used in the problem.
			//#Next we specify which media are in which geometric regions
			//#note that by default all regions contain
			//#medium 1 and which medium to input as 1 should be selected with this in mind.
			EGS4Geom.nMEDIA=2;
			EGS4.MEDIA[0]="170C521ICRU_fortran";
			EGS4.MEDIA[1]="AIR521ICRU_fortran";

			//DESCRIPTION BY:#planes,regions;
			//#planes: use slab and cylinder no.'s to define what medium goes where
			//#regions: use region numbers to define this (region numbers start at 2 and
			//  number from top to bottom of geometry and innermost radius to outermost radius)
			EGS4Geom.DESCRIBE=EGS4Geom.DESCRIBE_REGIONS;//##INPUT DATA

			EGS4Geom.nMEDNUM=1;
			EGS4Geom.MEDNUM[0]=2;//MEDNUM->AIR

			if(EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_REGIONS_DENSITY ||
			   EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_PLANES_DENSITY)
			{
				EGS4Geom.nRHOR=1;
				EGS4Geom.RHOR[0]=0.;
			}
			if(EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_REGIONS ||
			   EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_REGIONS_DENSITY)
			{
				EGS4Geom.nNREGLO=1;
				EGS4Geom.nNREGHI=1;
				EGS4Geom.NREGLO[0]=3;//START REGION
				EGS4Geom.NREGHI[0]=3;//STOP REGION
			}
			if(EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_PLANES ||
			   EGS4Geom.DESCRIBE==EGS4Geom.DESCRIBE_PLANES_DENSITY)
			{
				EGS4Geom.nNZLO=2;
				EGS4Geom.nNZHI=2;
				EGS4Geom.nNRLO=2;
				EGS4Geom.nNRHI=2;
			//#This puts MEDIUM 1 everywhere and then #inserts a small column of MEDIUM 2
			//on the central #axis with radius 1cm and going from Z=10cm #to Z=12cm
				EGS4Geom.NZLO[0]=1;//START ZSLAB
				EGS4Geom.NZLO[1]=10;//START ZSLAB
				EGS4Geom.NZHI[0]=20;//STOP ZSLAB
				EGS4Geom.NZHI[1]=12;//STOP ZSLAB
				EGS4Geom.NRLO[0]=1;//START RING
				EGS4Geom.NRLO[1]=1;//START RING
				EGS4Geom.NRHI[0]=2;//STOP RING
				EGS4Geom.NRHI[1]=1;//STOP RING
			}
		}//ITERSE->0 sau 1
		//cavity inputs:
		EGS4Geom.NSUMCV=1;//NUMBER OF CAVITY REGIONS
		//#this defines the small cylinder of
		//#air in the geometry to be the cavity
		EGS4Geom.ISUMCV[0]=3;//REGION NUMBERS OF THE CAVITY= 3

		//CALL GEOMRZ!!
      	EGS4.seqStr=" *** Reading geometrical info ... ***";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);

		EGS4Geom.GEOMRZ();if(EGS4.STOPPROGRAM){return;}

	  	//INCIDENT PARTICLE= photon       #electron,photon,positron,all,charged;
        //                          #all & charged: only for phase space sources
//###########################################################################################
		EGS4SrcEns.ipart=EGS4SrcEns.ipart_photon;
		//  SOURCE NUMBER= 1                #0,1,2,3,4,10,11,12,13,14,15,16,20,21,22,23
		EGS4SrcEns.ISOURC=EGS4SrcEns.point_source_on_axis_incident_from_the_front;//1;
		//SOURCE OPTIONS=  100., 1.3, 0, 0
	//#for source 1: SSD of beam, radius of #    beam on front surface
	    if (EGS4SrcEns.ISOURC!=EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)//20)
	    {
        	EGS4SrcEns.source_option[0]=100.;
        	EGS4SrcEns.source_option[1]=1.3;
        	EGS4SrcEns.source_option[2]=0.;
        	EGS4SrcEns.source_option[3]=0.;

        	if (EGS4SrcEns.ISOURC==EGS4SrcEns.full_phase_space_from_any_angle)//22)
        	{
	        	EGS4SrcEns.source_option[4]=0.;
	        	EGS4SrcEns.source_option[5]=0.;
	        	EGS4SrcEns.source_option[6]=0.;
	        	EGS4SrcEns.source_option[7]=0.;
	        	EGS4SrcEns.source_option[8]=0.;
			}
			else if(EGS4SrcEns.ISOURC==EGS4SrcEns.BEAM_treatment_head_simulation_from_any_angle)//23
			{
				EGS4SrcEns.source_option[4]=0.;
			}
		}
        if (EGS4SrcEns.ISOURC==EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)//20)
        {
        	//source 20 Ex:
			EGS4SrcEns.MODEIN=EGS4SrcEns.MODEIN_LOCAL;
			EGS4SrcEns.NRDIST=1;
			EGS4SrcEns.RDISTF[0]=1.0;
			EGS4SrcEns.RPDF[0]=1.0;
			EGS4SrcEns.RDIST_IOUTSP=EGS4SrcEns.RDIST_IOUTSP_NONE;
		}
		//CALL SRCRZ->"Get source data"
      	EGS4.seqStr=" *** Reading radiation source info ... ***";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
		EGS4SrcEns.SRCRZ();if(EGS4.STOPPROGRAM){return;}
		//source initialization
		EGS4SrcEns.SRCINI();if(EGS4.STOPPROGRAM){return;}

		if (EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22||EGS4SrcEns.ISOURC==23)
		{
			EGS4SrcEns.MONOEN=0;
		}
		//             "no need to input monoen for source 21,22"
		else
		{
			//INCIDENT ENERGY= monoenergetic  #monoenergetic, spectrum;
			EGS4SrcEns.monoindex=EGS4SrcEns.iMONOENERGETIC;
			//INCIDENT KINETIC ENERGY(MEV)= 1.25 #only use for "monoenergetic"
			EGS4SrcEns.ikemev=1.25;

			EGS4SrcEns.ENSRC();if(EGS4.STOPPROGRAM){return;}
		}// "Get data re-source energies"
//###########################################################################################
		//get_transport_parameter();
        setEcutRegion=false;
        setPcutRegion=false;
        setSmaxRegion=false;
		ecut=0.521;
		pcut=0.001;
		smax=1.e10;

		if(!setEcutRegion)
			EGS4.ECUT[0]=ecut;//#Electron cutoff for transport
		if(!setPcutRegion)
			EGS4.PCUT[0]=pcut;//#Photon cutoff for transport
		if(!setSmaxRegion)
			EGS4.SMAXIR[0]=smax;//#Maximum step size in cm (not needed
      					     //#unless old PRESTA algorithm used)
		if(setEcutRegion)
		{
			nEcut=1;//number of data
			Ecut[0]=1.;
			startEcutRegion[0]=1;
			stopEcutRegion[0]=1;
		}
		if(setPcutRegion)
		{
			nPcut=1;//number of data
			Pcut[0]=1.;
			startPcutRegion[0]=1;
			stopPcutRegion[0]=1;
		}
		if(setSmaxRegion)
		{
			nSmax=1;//number of data
			Smax[0]=1.;
			startSmaxRegion[0]=1;
			stopSmaxRegion[0]=1;
		}
		//Bound Compton scattering=  On or Off->IBCMP
		//#######//On means 1 and Off means 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//                                        #Off: Klein-Nishina used for compton
		//                                        #     scattering
		//                                        #On: Impulse approximation used for
		//                                        #    scattering
		//					#It has not been established that  the
		//					#the scoring routines work with this
		//					#option ON.
		incoh=incoh_OFF;
		if(incoh==incoh_OFF || incoh==incoh_ON)
			setIncohRegion=false;
		else
			setIncohRegion=true;
		if(setIncohRegion)
		{
			nIncoh=1;//number of data
			Incoh[0]=incoh_OFF;//0 or 1
			startIncohRegion[0]=1;
			stopIncohRegion[0]=1;
		}
		//Rayleigh scattering=            On->IRAYLR
		//                                        #Off: no coherent scattering
		//                                        #On: simulates coherent scattering
		coh=coh_OFF;
		if(coh==coh_OFF || coh==coh_ON)
			setCohRegion=false;
		else
			setCohRegion=true;
		if(setCohRegion)
		{
			nCoh=1;//number of data
			Coh[0]=coh_OFF;//0 or 1
			startCohRegion[0]=1;
			stopCohRegion[0]=1;
		}
		//Atomic relaxations=             On->IEDGFL
		//					#Note, it has not been verified that
		//					#AUSGAB works with this option ON
		//                                        #On: use correct cross section
		//                                        #  for p.e. events and shell vacancies
		//                                        #  for Compton & p.e. events are relaxed
		//                                        #  via emission of fluorescent X-Rays,
		//                                        #  Auger and Koster-Cronig electrons
		relax=relax_OFF;
		if(relax==relax_OFF || relax==relax_ON)
			setRelaxRegion=false;
		else
			setRelaxRegion=true;
		if(setRelaxRegion)
		{
			nRelax=1;//number of data
			Relax[0]=relax_OFF;//0 or 1
			startRelaxRegion[0]=1;
			stopRelaxRegion[0]=1;
		}
		//Photoelectron angular sampling= On->IPHTER
		//                                        #Off: Photoelectrons get direction of
		//                                        #     photon that creates them
		//                                        #On: Sauter's formula is used
		pe=pe_ang_OFF;
		if(pe==pe_ang_OFF || pe==pe_ang_ON)
			setPeRegion=false;
		else
			setPeRegion=true;
		if(setPeRegion)
		{
			nPe=1;//number of data
			Pe[0]=pe_ang_OFF;//0 or 1
			startPeRegion[0]=1;
			stopPeRegion[0]=1;
		}
		//Brems angular sampling=         On ->IBRDST
		//                                        #Simple: leading term of Koch-Motz
		//                                        #        dist'n used to determine angle
		//                                        #        of bremsstrahlung photons
		//                                        #KM:  Koch-Motz distribution used to
		//                                        #     determine angle -> 2BS(modified)
		EGS4.ibrdst = brems_ang_KM;
		//Pair angular sampling=          On->IPRDST
		//                                     #Simple: use leading term of K-M
		//                                     #        dist'n
		//                                     #KM: use complete Koch and Motz dist'n
		//                                     #Off: angle of pairs is m/E--like old EGS4
		EGS4.iprdst = pair_ang_SIMPLE;
		//ESTEPE=                         0.25->ESTEPE
		//		#Max fractional continuous energy loss
		//		#per step. Use 0.25 unless using
        //      #PRESTA-I
		EGS4.estepe =0.25;
		//XIMAX=                          0.5->XIMAX
		//	#Max first elastic scattering moment
        //	#per step.  Using default.
		EGS4.ximax = 0.5;
		//Boundary crossing algorithm=    exact->bca_algorithm, exact_bca
		//                                 #exact: cross boundaries in single scattering
		//                                 #       mode (distance at which to go into
		//                                 #       single scattering mode determined by
		//                                 #       "Skin depth for BCA"
		//                                 #PRESTA-I: cross boundaries with lateral
		//                                 #          correlations off and force multiple
		//                                 #          scattering mode
		EGS4.bca_algorithm = BCA_EXACT;//exact means =EGS4.$BCA_ALGORITHM_DEFAULT=0;
		//Skin depth for BCA=             3->skindepth_for_bca
		//										#Distance from a boundary (in elastic
		//                                      #MFP) at which the algorithm will go
		//                                      #into single scattering mode (using
        //		                                #default here)
		EGS4.skindepth_for_bca=3.0;
		//Electron-step algorithm=        default->transport_algorithm
		//PRESTA-II (the default),$PRESTA_II = 0;
		//                                        #Determines the algorithm used to take
		//                                        #into account lateral and longitudinal
		//                                        #correlations in a condensed history
		//                                        #step
		EGS4.transport_algorithm = estep_alg_PRESTA_II;//0;
		//Spin effects=                   On->spin_effects
		//                                        #Turns off/on spin effects for electron
		//                                        #elastic scattering. Spin On is
		//                                        #ABSOLUTELY necessary for good
		//                                        #backscattering calculations. Will
		//                                        #make a difference even in `well
		//                                        #conditioned' situations (e.g.  depth
		//                                        #dose curves).
		ispin = spin_ON;
		//Brems cross sections=           BH
		//		                                #BH: Bethe-Heitler cross-sections used
		//                                      #NIST: NIST cross-sections used
		EGS4.ibr_nist=brems_cross_BH;
		//" Pair cross sections "= BH
		EGS4.pair_nrc = pair_cross_BH;
		//Electron impact ionization
		EGS4.eii_flag=eii_OFF;
		//Triplet
		EGS4.itriplet=triplet_OFF;
		//Radiative compton correction
		EGS4.radc_flag=radc_OFF;
//###########################################################################################
		//VARIANCE REDUCTION:
		//ELECTRON RANGE REJECTION
		EGS4Macro.irejct=irejct_ON;
		//                        #On: if charged particle energy is below ESAVEIN
        //                        #    and it cannot get out of current region
        //                        #    with energy > ECUT, the particle is
        //                        #    terminated
		//		#also terminates all electrons which cannot
		//		#reach the cavity under conservative assumptions.
		ESAVEIN=2.0;//#total energy below which range rejection is considered
		EGS4Macro.cs_enhance=1.0;//     #Photon cross section scaling factors
		// RUSSIAN ROULETTE DEPTH=  0.0000 #play Russian Roulette with photons once they
		//                                 #cross this Z plane
		RRZ=0.000;
		// RUSSIAN ROULETTE FRACTION=  0.0000  #probability of photon survival--if this
		//                                     #and #RUSSIAN ROULETTE DEPTH both 0, then
		//                                     #photon Russian Roulette is not played
		RRCUT=0.000;
		//                                #exponential pathlength biasing can be
		//                                #used. See Rogers&Bielajew 1990 review for
		//                                #discussion.  C<0 => pathlength shortening
		//                                #              >0 => pathlength stretching
		//                                #               along z axis both cases
		//                     #CAVRZnrc allows for having the photon cross
		//                     #section scaled to enhance interactions.
		//                     #    If this input is missing or set to <= 1, it
		//                     #    has no effect on the simulation. But if
		//                     #    the enhancement factor is set to > 1, the
		//                     #    effect is dramatic: all other user input
		//                     #    concerning photon forcing, splitting, exp.
		//                     #    transform, etc., is ignored. In addition,
		//                     #    the calculation result corresponds  ALWAYS
		//                     #    to 'Aatt and Ascat', no matter what the
		//                     #    user requested (but only Awall is calculated,
		//                     #    not the individual Ascat and Aatt).
		//                     #    The algorithm employed is implemented via
		//                     #    $RAYLEIGH-CORRECTION and appropriate calls to
		//                     #    AUSGAB. For more detail see the manual and the
		//                     #    header of cavrznrc.mortran
		EGS4Macro.CEXPTR=0.000;
		// PHOTON FORCING= On             #Off (default),On;
		//                                #On: force photons to interact according to
		//                                #    START FORCING and STOP FORCING AFTER inputs
		IFARCE=IFARCE_ON;
		EGS4Macro.NFMIN=1;//#Start forcing at this interaction number
		EGS4Macro.NFMAX=1;//#Number of photon interactions after which
                          //#to stop forcing photon interactions
		// PHOTON SPLITTING= 1            #no. of times to split a photon
		//                                #if < 2-->normal transport
		//                                #overrides PHOTON FORCING if >= 2
		//                                #can only be >= 2 if IFULL= dose and stoppers
		//                                # or if  IFULL= Aatt and Ascat
		phsplitt=1;



		test_inputs();if(EGS4.STOPPROGRAM){return;}
//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################
		//"                        SCORING ARRAY INITIALISATION
		//"                        ****************************
		NCASEO=0;NCASET=0;TMCPUO=0; NNREADO=0; //"SET PREVIOUS RUN COUNTERS"
		if(IRESTART==0||IRESTART==5)
		{// "FRESH START, SET EVERYTHING TO ZERO"
		    EGS4SrcEns.NNREAD=0;
		    SCSTP=0.; SCSTP2=0.; SCSTP_TMP=0.; SCSTP_LAST=0;
		    SCCSTP=0.; SCCSTP2=0.; SCCSTP_TMP=0.; SCCSTP_LAST=0;
		    cav_dose  = 0; cav2_dose  = 0;
		    cav_dose0 = 0; cav2_dose0 = 0;
		    cav_dose1 = 0; cav2_dose1 = 0;
		    cav_dose2 = 0; cav2_dose2 = 0;

		    cav_dosec   = 0;
		    cav_dosec01 = 0;
		    cav_dosec02 = 0;

		    tmp_dose  = 0.0;
		    tmp_dose0 = 0.0;
		    tmp_dose1 = 0.0;
		    tmp_dose2 = 0.0;

		    PIISTP=0.;
		    for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)
		    {
		       for(int IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[
		       {
		           for(int IT=1;IT<=4;IT++)
		           {
		             SCDOSE[IZ-1][IX-1][IT-1]=0.;
		             SCDOSE2[IZ-1][IX-1][IT-1]=0.;
		             SCDOSE_TMP[IZ-1][IX-1][IT-1]=0.;
		             if(IT<4)SCDOSE_COV[IZ-1][IX-1][IT-1]=0.;
					}
		           SCDOSE_LAST[IZ-1][IX-1]=0;
				}
			}
		}
		else if(IRESTART!=4)//NOT ALLOWED HERE
		{
		    //"Restart or stats analysis only, read old data from unit 4"
		    //"open unit 4 as an old file"
		    //OUTPUT;(' About to read the previous .egsdat file');
		    /*
		    my_unit = egs_open_datfile(4,0,1,'.egsdat');
		    READ(my_unit,*) SCSTP,SCSTP2,SCCSTP,SCCSTP2;
		    READ(my_unit,*) cav_dose, cav_dose0, cav_dose1, cav_dose2;
		    READ(my_unit,*) cav2_dose,cav2_dose0,cav2_dose1,cav2_dose2;
		    READ(my_unit,*) cav_dosec,cav_dosec01,cav_dosec02;

		    IF(NSUMCV>1)["get data for individual cavity regions"
		       DO IZ=1,NZ[
		           DO IX=1,NR[
		              $GET-IRL(IZ,IX);
		              IF(NTRACK(IRL).EQ.1)[
		                   READ(my_unit,*)(SCDOSE(IZ,IX,IT),SCDOSE2(IZ,IX,IT),IT=1,4);
		                   READ(my_unit,*)(SCDOSE_COV(IZ,IX,IT),IT=1,3);
		              ]
		           ]
		        ]
		    ]
		    $RETRIEVE RNG STATE FROM UNIT my_unit;
		    READ(my_unit,*,END=:EOFA:)NCASEO,TMCPUO,NNREADO,PIISTP;
		    CLOSE(my_unit);
		    NNREAD=NNREADO;
		    */
		}

		if(IRESTART == 3){NCASE=0;}

		NCASET=NCASE+NCASEO;
		EGS4SrcEns.NCASET=NCASET;

      	EGS4.seqStr=" ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);

	}//inputs

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs()
	{
		if(IWATCH<0 || IWATCH>4){IWATCH=0;}//default
		EGS4Macro.IWATCH=IWATCH;EGS4Geom.IWATCH=IWATCH;
		//if(ISTORE<0 || ISTORE>3){ISTORE=0;}//default
if(ISTORE<0 || ISTORE>0){ISTORE=0;}//default
		//if(IRESTART<0 || IRESTART>4){IRESTART=0;}//default
if(IRESTART<0 || IRESTART>0){IRESTART=0;}//default
		if(IOOPTN<0 || IOOPTN>1){IOOPTN=0;}//default
		if(IDAT<0 || IDAT>1){IDAT=1;}//default
if(IDAT<1 || IDAT>1){IDAT=1;}//default->NO FILE STORAGE
		if(IRESTART == 4)//no parrallelllllllllllll!!
		{
		    IDAT=1;   //"do not store output in this case to avoid biasing"
		    ISTORE=0; //"do not store the starting random numbers either"
		}
		if (NCASE<NCASE_MIN || NCASE>NCASE_MAX){NCASE=NCASE_DEFAULT;}
		if (TIMMAX<TIMMAX_MIN || TIMMAX>TIMMAX_MAX){TIMMAX=TIMMAX_DEFAULT;}
		if(EGS4Macro.IFULL<0 || EGS4Macro.IFULL>1){EGS4Macro.IFULL=0;}//default
		if (STATLM<STATLM_MIN || STATLM>STATLM_MAX){STATLM=STATLM_DEFAULT;}
		if(IFANO<0 || IFANO>2){IFANO=0;}//default
		if(IWATCH==0 && NCASE < $NCASEMIN){NCASE=$NCASEMIN;}
		if( IFANO == 1 )
		{//  "With ifano option turned on, it is a waste of time to"
		//                   "have RAYLEIGH turned on (scattered photon will be killed,"
		//                   "original photon re-created) => turn Rayleigh off"
		    for(int j=1;j<=EGS4.$MXREG;j++) { EGS4.IRAYLR[j-1] = 0; }
	      	EGS4.seqStr=" ******** ifano set => turning off Rayleigh! **** ";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);

		    //OUTPUT; (//' ******** ifano set => turning off Rayleigh! **** '//);
		}
		if(EGS4Geom.iterseindex<0 || EGS4Geom.iterseindex>2)
		{
			EGS4.STOPPROGRAM=true;
	      	EGS4.seqStr=" ******** ERROR in geometrical inputs! **** ";
  			printSequence(EGS4.seqStr);
	  		return;
		}
		if(EGS4Geom.iterseindex==EGS4Geom.iCAVITY_INFORMATION)//2
		{
			if((EGS4Geom.SLENGHT).compareTo(NULLs)==0 || (EGS4Geom.airs).compareTo(NULLs)==0)
			{
				EGS4.STOPPROGRAM=true;
		      	EGS4.seqStr=" ******** ERROR in geometrical inputs! **** ";
	  			printSequence(EGS4.seqStr);
		  		return;
			}
			if(EGS4Geom.ELERAD!=0.0 && (EGS4Geom.electrods).compareTo(NULLs)==0)
			{
				EGS4.STOPPROGRAM=true;
		      	EGS4.seqStr=" ******** ERROR in geometrical inputs! **** ";
	  			printSequence(EGS4.seqStr);
		  		return;
			}
		}

		if(!setEcutRegion)
		{
			if(EGS4.ECUT[0]<ECUT_MIN || EGS4.ECUT[0]>ECUT_MAX)
				EGS4.ECUT[0]=EGS4.$GLOBAL_ECUT;
			for(int I=2;I<=EGS4.$MXREG;I++)
			{ EGS4.ECUT[I-1] = EGS4.ECUT[0];}

		}
		if(!setPcutRegion)
		{
			if(EGS4.PCUT[0]<PCUT_MIN || EGS4.PCUT[0]>PCUT_MAX)
				EGS4.PCUT[0]=EGS4.$GLOBAL_PCUT;
			//"Now set ecut and pcut to the values input by the user"
			for(int I=2;I<=EGS4.$MXREG;I++)
			{EGS4.PCUT[I-1] = EGS4.PCUT[0]; }
		}
		if(!setSmaxRegion)
		{
			if(EGS4.SMAXIR[0]<SMAXIR_MIN || EGS4.SMAXIR[0]>SMAXIR_MAX)
				EGS4.SMAXIR[0]=EGS4.$MAX_SMAX;
			for(int I=2;I<=EGS4.$MXREG;I++)
			{ EGS4.SMAXIR[I-1] = EGS4.SMAXIR[0];}
		}

		if(setEcutRegion)
		{
			for (int i=1; i<=nEcut; i++)
			{
				if(Ecut[i-1]<ECUT_MIN || Ecut[i-1]>ECUT_MAX)
					Ecut[i-1]=EGS4.$GLOBAL_ECUT;
				if(startEcutRegion[i-1]<1 || startEcutRegion[i-1]>EGS4.$MXREG)
					startEcutRegion[i-1]=1;
				if(stopEcutRegion[i-1]<1 || stopEcutRegion[i-1]>EGS4.$MXREG)
					stopEcutRegion[i-1]=1;

				for(int j=startEcutRegion[i-1];j<=stopEcutRegion[i-1];j++)
					EGS4.ECUT[j-1]=Ecut[i-1];
			}
		}
		if(setPcutRegion)
		{
			for (int i=1; i<=nPcut; i++)
			{
				if(Pcut[i-1]<PCUT_MIN || Pcut[i-1]>PCUT_MAX)
					Pcut[i-1]=EGS4.$GLOBAL_PCUT;
				if(startPcutRegion[i-1]<1 || startPcutRegion[i-1]>EGS4.$MXREG)
					startPcutRegion[i-1]=1;
				if(stopPcutRegion[i-1]<1 || stopPcutRegion[i-1]>EGS4.$MXREG)
					stopPcutRegion[i-1]=1;

				for(int j=startPcutRegion[i-1];j<=stopPcutRegion[i-1];j++)
					EGS4.PCUT[j-1]=Pcut[i-1];
			}
		}
		if(setSmaxRegion)
		{
			for (int i=1; i<=nSmax; i++)
			{
				if(Smax[i-1]<SMAXIR_MIN || Smax[i-1]>SMAXIR_MAX)
					Smax[i-1]=EGS4.$MAX_SMAX;
				if(startSmaxRegion[i-1]<1 || startSmaxRegion[i-1]>EGS4.$MXREG)
					startSmaxRegion[i-1]=1;
				if(stopSmaxRegion[i-1]<1 || stopSmaxRegion[i-1]>EGS4.$MXREG)
					stopSmaxRegion[i-1]=1;

				for(int j=startSmaxRegion[i-1];j<=stopSmaxRegion[i-1];j++)
					EGS4.SMAXIR[j-1]=Smax[i-1];
			}
		}


		if(!setIncohRegion)
		{
			EGS4.ibcmp[0] = incoh;
			if(EGS4.ibcmp[0]<0 || EGS4.ibcmp[0]>3)
			{	EGS4.ibcmp[0]=EGS4.$IBCMP_DEFAULT;incoh=EGS4.$IBCMP_DEFAULT;}
		    //"Now set ibcmp for all regions to the value input by the user"
			for(int I=2;I<=EGS4.$MXREG;I++) { EGS4.ibcmp[I-1] = EGS4.ibcmp[0]; }
		}
		if(setIncohRegion)
		{
			for (int i=1; i<=nIncoh; i++)
			{
				if(Incoh[i-1]<0 || Incoh[i-1]>1)//only 2 values
					Incoh[i-1]=EGS4.$IBCMP_DEFAULT;
				if(startIncohRegion[i-1]<1 || startIncohRegion[i-1]>EGS4.$MXREG)
					startIncohRegion[i-1]=1;
				if(stopIncohRegion[i-1]<1 || stopIncohRegion[i-1]>EGS4.$MXREG)
					stopIncohRegion[i-1]=1;

				for(int j=startIncohRegion[i-1];j<=stopIncohRegion[i-1];j++)
					EGS4.ibcmp[j-1]=Incoh[i-1];
			}
		}
		if(!setCohRegion)
		{
			EGS4.IRAYLR[0] = coh;
			if(EGS4.IRAYLR[0]<0 || EGS4.IRAYLR[0]>3)
			{	EGS4.IRAYLR[0]=EGS4.$IRAYLR_DEFAULT;coh=EGS4.$IRAYLR_DEFAULT;}
			//"Now set iraylr for all regions to the value input by the user"
			for(int I=2;I<=EGS4.$MXREG;I++) { EGS4.IRAYLR[I-1] = EGS4.IRAYLR[0]; }
		}
		if(setCohRegion)
		{
			for (int i=1; i<=nCoh; i++)
			{
				if(Coh[i-1]<0 || Coh[i-1]>1)//only 2 values
					Coh[i-1]=EGS4.$IRAYLR_DEFAULT;
				if(startCohRegion[i-1]<1 || startCohRegion[i-1]>EGS4.$MXREG)
					startCohRegion[i-1]=1;
				if(stopCohRegion[i-1]<1 || stopCohRegion[i-1]>EGS4.$MXREG)
					stopCohRegion[i-1]=1;

				for(int j=startCohRegion[i-1];j<=stopCohRegion[i-1];j++)
					EGS4.IRAYLR[j-1]=Coh[i-1];
			}
		}
		if(!setRelaxRegion)
		{
			EGS4.iedgfl[0] = relax;
			if(EGS4.iedgfl[0]<0 || EGS4.iedgfl[0]>3)
			{	EGS4.iedgfl[0]=EGS4.$IEDGFL_DEFAULT;relax=EGS4.$IEDGFL_DEFAULT;}
			//"Now set iedgfl for all regions to the value input by the user"
			for(int I=2;I<=EGS4.$MXREG;I++) { EGS4.iedgfl[I-1] = EGS4.iedgfl[0]; }
		}
		if(setRelaxRegion)
		{
			for (int i=1; i<=nRelax; i++)
			{
				if(Relax[i-1]<0 || Relax[i-1]>1)//only 2 values
					Relax[i-1]=EGS4.$IEDGFL_DEFAULT;
				if(startRelaxRegion[i-1]<1 || startRelaxRegion[i-1]>EGS4.$MXREG)
					startRelaxRegion[i-1]=1;
				if(stopRelaxRegion[i-1]<1 || stopRelaxRegion[i-1]>EGS4.$MXREG)
					stopRelaxRegion[i-1]=1;

				for(int j=startRelaxRegion[i-1];j<=stopRelaxRegion[i-1];j++)
					EGS4.iedgfl[j-1]=Relax[i-1];
			}
		}
		if(!setPeRegion)
		{
			EGS4.iphter[0] = pe;
			if(EGS4.iphter[0]<0 || EGS4.iphter[0]>3)
			{	EGS4.iphter[0]=EGS4.$IPHTER_DEFAULT;pe=EGS4.$IPHTER_DEFAULT;}
			//"Now set iphter for all regions to the value input by the user"
			for(int I=2;I<=EGS4.$MXREG;I++) { EGS4.iphter[I-1] = EGS4.iphter[0]; }
		}
		if(setPeRegion)
		{
			for (int i=1; i<=nPe; i++)
			{
				if(Pe[i-1]<0 || Pe[i-1]>1)//only 2 values
					Pe[i-1]=EGS4.$IPHTER_DEFAULT;
				if(startPeRegion[i-1]<1 || startPeRegion[i-1]>EGS4.$MXREG)
					startPeRegion[i-1]=1;
				if(stopPeRegion[i-1]<1 || stopPeRegion[i-1]>EGS4.$MXREG)
					stopPeRegion[i-1]=1;

				for(int j=startPeRegion[i-1];j<=stopPeRegion[i-1];j++)
					EGS4.iphter[j-1]=Pe[i-1];
			}
		}
		if( (EGS4.ibrdst < 0) || (EGS4.ibrdst > 1 ))
		{
		    EGS4.ibrdst = EGS4.$IBRDST_DEFAULT;
		}
		if( (EGS4.iprdst < 0) || (EGS4.iprdst > 4 ))
		{
		    EGS4.iprdst = EGS4.$IPRDST_DEFAULT;
		}
		if( (EGS4.estepe < ESTEPE_MIN) || (EGS4.estepe > ESTEPE_MAX ))
		{
		    EGS4.estepe = EGS4.$MAX_ELOSS;  //"$MAX-ELOSS is defined in egsnrc.macros at 0.25"
		}
		if(( EGS4.ximax < XIMAX_MIN) || (EGS4.ximax > XIMAX_MAX ))
		{
		    EGS4.ximax = EGS4.$EXACT_BCA_XIMAX; //"$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
		}
		if( (EGS4.bca_algorithm < 0) || (EGS4.bca_algorithm > 1 ))
		{
		    EGS4.bca_algorithm = EGS4.$BCA_ALGORITHM_DEFAULT;
		}
		if((EGS4.skindepth_for_bca < Skindepth_MIN) || (EGS4.skindepth_for_bca > Skindepth_MAX ))
			EGS4.skindepth_for_bca=EGS4.$SKIN_DEPTH_FOR_BCA;
		if( EGS4.bca_algorithm == BCA_EXACT )
		{
		    if( EGS4.skindepth_for_bca <= 0.0 ) EGS4.skindepth_for_bca = EGS4.$SKIN_DEPTH_FOR_BCA;
		}
		if( (EGS4.transport_algorithm < 0) || (EGS4.transport_algorithm > 1 ))
		{
		    EGS4.transport_algorithm = EGS4.$TRANSPORT_ALGORITHM_DEFAULT;
		}
		if (ispin==spin_ON)
			EGS4.spin_effects = true;//so=>
		else
			EGS4.spin_effects = false;//so=>
		if( (EGS4.ibr_nist < 0) || (EGS4.ibr_nist > 1 ))
		{
		    EGS4.ibr_nist = EGS4.$IBR_NIST_DEFAULT;
		}
		if( (EGS4.pair_nrc < 0) || (EGS4.pair_nrc > 1 ))
		{
		    EGS4.pair_nrc = EGS4.$PAIR_NRC_DEFAULT;
		}
		if( (EGS4.eii_flag < 0) || (EGS4.eii_flag > 1 ))
		{
		    EGS4.eii_flag = eii_OFF;//default
		}
		if( (EGS4.itriplet < 0) || (EGS4.itriplet > 1 ))
		{
		    EGS4.itriplet = EGS4.$TRIPLET_DEFAULT;
		}
		if( (EGS4.radc_flag < 0) || (EGS4.radc_flag > 1 ))
		{
		    EGS4.radc_flag = radc_OFF;//default
		}

		if(EGS4Macro.irejct<0 || EGS4Macro.irejct>1){EGS4Macro.irejct=irejct_OFF;}//default
		if(ESAVEIN<ESAVEIN_MIN || ESAVEIN>ESAVEIN_MAX){ESAVEIN=ESAVEIN_DEFAULT;}
		if(EGS4Macro.cs_enhance<cs_enhance_MIN || EGS4Macro.cs_enhance>cs_enhance_MAX)
		{EGS4Macro.cs_enhance=cs_enhance_DEFAULT;}
		if( EGS4Macro.cs_enhance > 1. ) { EGS4Macro.use_enhance = true; }
		else { EGS4Macro.use_enhance = false; }

      	EGS4.seqStr=" *** Reading variance reduction inputs ... ***";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);

      	EGS4.seqStr=" Range rejection is On(1) or Off(0):"+EGS4Macro.irejct;
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);

		if(EGS4Macro.irejct > 0)
		{
	      	EGS4.seqStr=" ESAVEIN cutoff value(total) for range rejection:"+
	      	EGS4.format(ESAVEIN,10,true)+" MeV";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);

		   if (ESAVEIN == 0.0 )
		   {
		      	EGS4.seqStr=" WARNING: Have asked for range rejection but left ESAVEIN=0.0!!";
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);
		      	EGS4.seqStr=" There will be no range rejection!";
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);
		   }
		   for (int i=1;i<=EGS4Geom.nreg;i++)
		   {
			   EGS4.i_do_rr[i-1] = 1;
			   EGS4.e_max_rr[i-1] = ESAVEIN;
		   }
		   //"note  e_max_r is total energy"
		   //"above two arrays needed for each region for EGSnrc RANGE-DISCARD macro"
		}
		if(RRZ<RRDEPTH_MIN || RRZ>RRDEPTH_MAX){RRZ=RRDEPTH_DEFAULT;}
		if(RRCUT<RRFRACTION_MIN || RRCUT>RRFRACTION_MAX){RRCUT=RRFRACTION_DEFAULT;}
		if(EGS4Macro.CEXPTR<EXPC_MIN || EGS4Macro.CEXPTR>EXPC_MAX)
		{EGS4Macro.CEXPTR=EXPC_DEFAULT;}
		RUSROU=false;if(RRZ+RRCUT!=0.0)RUSROU=true;
		if(RUSROU)
		{
	      	EGS4.seqStr=" RUSSIAN ROULETTE WILL BE PLAYED";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
	      	EGS4.seqStr=" RUSSIAN ROULETTE PLANE:"+EGS4.format(RRZ,14,false);
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
	      	EGS4.seqStr=" SURVIVAL PROBABILITY:"+EGS4.format(RRCUT,14,false);
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		else
		{
	      	EGS4.seqStr=" RUSSIAN ROULETTE WILL NOT BE PLAYED";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		if (EGS4Macro.CEXPTR == 0.)
		{
	      	EGS4.seqStr=" NO PATHLENGTH BIASING TO BE DONE";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		else
		{
	      	EGS4.seqStr=" CEXPTR PARAMATER:"+EGS4.format(EGS4Macro.CEXPTR,14,false);
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		if(IFARCE<0 || IFARCE>1){IFARCE=IFARCE_OFF;}//default
		NFMIN_MAX=EGS4Geom.nreg;
		NFMAX_MAX=EGS4Geom.nreg+1;
		if(EGS4Macro.NFMIN<NFMIN_MIN || EGS4Macro.NFMIN>NFMIN_MAX)
		{EGS4Macro.NFMIN=NFMIN_DEFAULT;}
		if(EGS4Macro.NFMAX<NFMAX_MIN || EGS4Macro.NFMAX>NFMAX_MAX)
		{EGS4Macro.NFMAX=NFMAX_DEFAULT;}

		if (IFARCE==0)
		{
		      EGS4Macro.IFORCE=0;EGS4Macro.NFMIN=0;EGS4Macro.NFMAX=0;

	      	EGS4.seqStr=" NO INTERACTION FORCING IS IN EFFECT";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		else if (IFARCE==1)
		{
		    EGS4Macro.IFORCE=1;
		    if (EGS4Macro.NFMAX<EGS4Macro.NFMIN) EGS4Macro.NFMAX=EGS4Macro.NFMIN;

	      	EGS4.seqStr=" FORCED PHOTON INTERACTIONS IN EFFECT FROM "+EGS4Macro.NFMIN+
	      	" TO "+EGS4Macro.NFMAX+" INTERACTIONS ";
		  	if(EGS4.iprint>1)
	  			printSequence(EGS4.seqStr);
		}
		EGS4Macro.IQINC=EGS4SrcEns.iqin; //"NEEDED TO TURN OFF FASTSTEP FOR INCIDENT ELECTRONS"
		            //"WHEN FORCING INTERACTIONS"
		phsplitt_MAX=EGS4.$MXSTACK-2;
		if(phsplitt<phsplitt_MIN || phsplitt>phsplitt_MAX){phsplitt=phsplitt_DEFAULT;}
		EGS4Macro.n_split = phsplitt;
		if( EGS4Macro.n_split > 1 )
		{
		    if( EGS4Macro.IFULL > 1 )
		    {
		      	EGS4.seqStr=" IGNORING INPUT: Photon splitting only for ifull = 0,1! ";
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);

		            EGS4Macro.n_split = 1;
			}
		    else
		    {
		      	EGS4.seqStr=" Calculation with photon splitting, n_split = "+EGS4Macro.n_split;
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);

		            EGS4Macro.iifano = IFANO;
			}
		}

		if( EGS4Macro.use_enhance )
		{
		      	EGS4.seqStr="  Calculation with CS enhancement  ";
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);
		      	EGS4.seqStr="  photon forcing, exp. transform, etc. input will be ignored, IFULL will be set to 1! (i.e. Ascat and Aatt) !";
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);
		      	EGS4.seqStr="  Using cs_enhance = "+EGS4Macro.cs_enhance;
			  	if(EGS4.iprint>1)
		  			printSequence(EGS4.seqStr);

		    EGS4Macro.IFULL = 1; EGS4Macro.IFORCE=0;
		    EGS4Macro.NFMIN=0;EGS4Macro.NFMAX=0; EGS4Macro.n_split = 1;
		}

		//"if we are calculating Aatt and Ascatt, we must have n_split>1,"
		//"cs_enhance on, or photon forcing on"
		if(EGS4Macro.IFULL==1 && EGS4Macro.cs_enhance<=1 &&
			EGS4Macro.IFORCE==0 && EGS4Macro.n_split <=1)
		{
		  EGS4Macro.IFORCE=1;
		  EGS4Macro.NFMIN=1;
		  EGS4Macro.NFMAX=1;

      	EGS4.seqStr=" If you are calculating Aatt, Ascat (IFULL=1), you must be using";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" photon forcing, photon splitting, or cross-section enhancement.";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" Currently, none of these are being used.  Will continue run with";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" photon forcing on and one interaction forced (NFMIN=NFMAX=1).";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
		}

		if(IOOPTN==1 && (EGS4Macro.n_split>1 || EGS4Macro.use_enhance))
		{
		  IOOPTN=0;
      	EGS4.seqStr=" You cannot have cross-section enhancement or photon splitting on ";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
      	EGS4.seqStr=" and still output detailed results for each cavity zone.  IOOPTN reset to 0.";
	  	if(EGS4.iprint>1)
  			printSequence(EGS4.seqStr);
		}

	}

	/**
	 * Initialize random number generator.
	 */
	private void init_random_generator()
	{
		if (EGS4.ranluxB)
		{
			if(jrng1<RANLUX_LEVEL_MIN || jrng1>RANLUX_LEVEL_MAX){jrng1=RANLUX_LEVEL_DEFAULT;}
			if(jrng2<RANLUX_SEED_MIN || jrng2>RANLUX_SEED_MAX){jrng2=RANLUX_SEED_DEFAULT;}
		}
		else
		{
			if(jrng1<RANMAR_SEED_MIN || jrng1>RANMAR_SEED_MAX){jrng1=RANMAR_SEED_DEFAULT;}
			if(jrng2<RANMAR_SEED_MIN || jrng2>RANMAR_SEED_MAX){jrng2=RANMAR_SEED_DEFAULT;}
		}
//IF( i_parallel > 0 ) jrng2 = jrng2 - 1 + i_parallel;//NO PARALLEL JOB USED!!!
		//$INITIALIZE RNG USING jrng1 AND jrng2;
		if(EGS4.ranluxB)
		{
    		EGS4.init_ranlux(jrng1,jrng2);
    		EGS4.ranlux(EGS4.rng_array); EGS4.rng_seed = 1;
		}
		else
		{
  			EGS4.ixx = jrng1; EGS4.jxx = jrng2;
  			EGS4.init_ranmar();
		}
	}

	/**
	 * Print input summary.
	 */
	private void ISUMRY()
	{
		String s="";String s1="";int ll=0;int ioff=0;
		EGS4.seqStr="************************INPUT SUMMARY:*****************************************";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="                   Electron/Photon transport parameter";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Photon transport cutoff(MeV)";ioff=32;
		if(setPcutRegion)
		{
			s=s+EGS4.format("",ioff)+"Set in regions";
		}
		else
		{
    		if( EGS4.PCUT[0] > 1.e-4 )
    		{ioff=30; s=s+EGS4.format("",ioff)+EGS4.format(EGS4.PCUT[0],8,true); }
    		else {s=s+EGS4.format("",ioff)+"AP(medium)";}
		}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Pair angular sampling";ioff=39;
		if(EGS4.iprdst==pair_ang_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(EGS4.iprdst==pair_ang_SIMPLE){s=s+EGS4.format("",ioff)+"SIMPLE";}
		else if(EGS4.iprdst==pair_ang_KM){s=s+EGS4.format("",ioff)+"KM";}
		else if(EGS4.iprdst==pair_ang_UNIFORM){s=s+EGS4.format("",ioff)+"UNIFORM";}
		else if(EGS4.iprdst==pair_ang_BLEND){s=s+EGS4.format("",ioff)+"BLEND";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Pair cross sections";ioff=41;
		if(EGS4.pair_nrc == pair_cross_BH){s=s+EGS4.format("",ioff)+"BH";}
		else if(EGS4.pair_nrc == pair_cross_NRC){s=s+EGS4.format("",ioff)+"NRC";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Triplet production";ioff=42;
		if(EGS4.itriplet==triplet_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(EGS4.itriplet==triplet_ON){s=s+EGS4.format("",ioff)+"ON";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Bound Compton scattering";ioff=36;
		if(incoh==incoh_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(incoh==incoh_ON){s=s+EGS4.format("",ioff)+"ON";}
		else if(incoh==incoh_ON_IN_REGIONS){s=s+EGS4.format("",ioff)+"ON IN REGIONS";}
		else if(incoh==incoh_OFF_IN_REGIONS){s=s+EGS4.format("",ioff)+"OFF IN REGIONS";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Radiative Compton corrections";ioff=31;
		if(EGS4.radc_flag == radc_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(EGS4.radc_flag == radc_ON){s=s+EGS4.format("",ioff)+"ON";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Rayleigh scattering";ioff=41;
		if(coh==coh_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(coh==coh_ON){s=s+EGS4.format("",ioff)+"ON";}
		else if(coh==coh_ON_IN_REGIONS){s=s+EGS4.format("",ioff)+"ON IN REGIONS";}
		else if(coh==coh_OFF_IN_REGIONS){s=s+EGS4.format("",ioff)+"OFF IN REGIONS";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Atomic relaxations";ioff=42;
		if(relax==relax_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(relax==relax_ON){s=s+EGS4.format("",ioff)+"ON";}
		else if(relax==relax_ON_IN_REGIONS){s=s+EGS4.format("",ioff)+"ON IN REGIONS";}
		else if(relax==relax_OFF_IN_REGIONS){s=s+EGS4.format("",ioff)+"OFF IN REGIONS";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Photoelectron angular sampling";ioff=30;
		if(pe==pe_ang_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(pe==pe_ang_ON){s=s+EGS4.format("",ioff)+"ON";}
		else if(pe==pe_ang_ON_IN_REGIONS){s=s+EGS4.format("",ioff)+"ON IN REGIONS";}
		else if(pe==pe_ang_OFF_IN_REGIONS){s=s+EGS4.format("",ioff)+"OFF IN REGIONS";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Electron transport cutoff(MeV)";ioff=30;
		if(setEcutRegion)
		{
			s=s+EGS4.format("",ioff)+"Set in regions";
		}
		else
		{
    		if( EGS4.ECUT[0] > 1.e-4 )
    		{ioff=28; s=s+EGS4.format("",ioff)+EGS4.format(EGS4.ECUT[0],8,true); }
    		else {s=s+EGS4.format("",ioff)+"AE(medium)";}
		}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Bremsstrahlung angular sampling";ioff=29;
		if(EGS4.ibrdst==brems_ang_SIMPLE){s=s+EGS4.format("",ioff)+"SIMPLE";}
		else if(EGS4.ibrdst==brems_ang_KM){s=s+EGS4.format("",ioff)+"KM";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Bremsstrahlung cross sections";ioff=31;
		if(EGS4.ibr_nist == brems_cross_BH){s=s+EGS4.format("",ioff)+"BH";}
		else if(EGS4.ibr_nist == brems_cross_NIST){s=s+EGS4.format("",ioff)+"NIST";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Spin effects";ioff=48;
		if(EGS4.spin_effects){s=s+EGS4.format("",ioff)+"ON";}
		else {s=s+EGS4.format("",ioff)+"OFF";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Electron Impact Ionization";ioff=34;
		if(EGS4.eii_flag==eii_OFF){s=s+EGS4.format("",ioff)+"OFF";}
		else if(EGS4.eii_flag==eii_ON){s=s+EGS4.format("",ioff)+"ON";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Maxium electron step in cm (SMAX)";ioff=27;
		if(setSmaxRegion)
		{
			s=s+EGS4.format("",ioff)+"Set in regions";
		}
		else
		{
    		if( EGS4.SMAXIR[0] > 1.e-4 )
    		{ioff=26; s=s+EGS4.format("",ioff)+EGS4.format(EGS4.SMAXIR[0],6,false); }
    		else {s=s+EGS4.format("",ioff)+"Restriction is off";}
		}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Maximum fractional energy loss/step (ESTEPE)";ioff=16;
		s=s+EGS4.format("",ioff)+EGS4.format(EGS4.estepe,6,true);
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Maximum 1st elastic moment/step (XIMAX)";ioff=21;
		s=s+EGS4.format("",ioff)+EGS4.format(EGS4.ximax,6,true);
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Boundary crossing algorithm";ioff=33;
		if(EGS4.bca_algorithm==BCA_EXACT){s=s+EGS4.format("",ioff)+"EXACT";}
		else if(EGS4.bca_algorithm==BCA_PRESTA_I){s=s+EGS4.format("",ioff)+"PRESTA I";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Skin-depth for boundary crossing (MFP)";ioff=22;
		s=s+EGS4.format("",ioff)+EGS4.format(EGS4.skindepth_for_bca,6,true);
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Electron-step algorithm";ioff=37;
		if(EGS4.transport_algorithm==estep_alg_PRESTA_II){s=s+EGS4.format("",ioff)+"PRESTA II";}
		else if(EGS4.transport_algorithm==estep_alg_PRESTA_I){s=s+EGS4.format("",ioff)+"PRESTA I";}
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="                   MONTE CARLO, TRANSPORT, AND SCATTER CONTROLS";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Max. # of histories to RUN";
		ll=s.length();ll=54-ll;s=s+EGS4.format("",ll);
		EGS4.seqStr=s+EGS4.format(NCASE,12);
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Max # of histories to ANALYZE";
		ll=s.length();ll=54-ll;s=s+EGS4.format("",ll);
		EGS4.seqStr=s+EGS4.format(NCASET,12);
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		s=" Incident Charge";//EGS4.format("",20)+"Incident Charge";
		s1="";
		ll=s.length();ll=68-ll;
		if(EGS4SrcEns.iqin == 0 ) s1="photons";
		if(EGS4SrcEns.iqin == -1) s1="electrons";
		if(EGS4SrcEns.iqin == 1 ) s1="positrons";
		if(EGS4SrcEns.iqin == 2 ) s1="all";
		if(EGS4SrcEns.iqin == 3 ) s1="e- & e+";
		s=s+EGS4.format(s1,ll);
		EGS4.seqStr=s;
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		 if(EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21 &&
		 EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23)
		 {
			 s=" Incident kinetic energy:";
		     ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(EGS4SrcEns.ein,9,true)+" MeV";
		     if(EGS4.iprint>1)
		     	printSequence(EGS4.seqStr);

		 }
		 else
		 {
			 EGS4SrcEns.ENSRCO();
		 }
		 boolean gotopast=false;
		 for(int I=2;I<=EGS4Geom.nreg;I++)
		 {
			 gotopast=false;
			 if( (EGS4.ECUT[I-1] != EGS4.ECUT[1]) ||
			     (EGS4.PCUT[I-1] != EGS4.PCUT[1]) )
			 {
		      //"we failed at least one test, so this means there really are"
		      //"varying ECUTs and these will be printed in the grid if we want them"
		      //"print the first 12 ECUT & PCUT just to be sure"
		      int j = Math.min(12,EGS4Geom.nreg);

			  for(int JJ=2;JJ<=j;JJ++)
			  {
			   EGS4.seqStr="First ECUTs (MeV): "+EGS4.format(EGS4.ECUT[JJ-1],12,true)+"  ,"+
			   "First PCUTs (MeV): "+EGS4.format(EGS4.PCUT[JJ-1],12,true);
			   if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			  }
		      //GO TO :past:
		      gotopast=true;
		      break;
		      }
		  }
		 //"if we get here, they were all the same"
		 if(!gotopast)
		 {
			 s=" GLOBAL ELECTRON TRANSPORT CUT-OFF";
		     ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(EGS4.ECUT[1],9,true)+" MeV";
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s=" GLOBAL PHOTON TRANSPORT CUT-OFF";
		     ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
			 EGS4.seqStr=s+EGS4.format(EGS4.PCUT[1],9,true)+" MeV";
			 if(EGS4.iprint>1)
			  	printSequence(EGS4.seqStr);

		 }
		 //:past:
		 if(EGS4Macro.IFORCE != 0)
		 {
			 s=" Min/max photon step forced";
		     ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(EGS4Macro.NFMIN,6)+"/"+EGS4.format(EGS4Macro.NFMAX,6);
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);

		 }
		 else
		 {
			 s=" Photon force interaction switch";
		     ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+"OFF";
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		 }
		 if(EGS4Macro.irejct > 0)
		 {
			 s=" Range rejection on a region by region basis";
		     //ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s=" Also globally to region between z="+EGS4.format(EGS4Macro.z_cavity_min,10,true)+
			 " &"+EGS4.format(EGS4Macro.z_cavity_max,10,true)+" cm";
		     //ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s="      and inside radius="+EGS4.format(EGS4Macro.r_cavity_max,12,true)+
			 " cm";
		     //ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s=" Range rejection only for electrons < ESAVEIN="+EGS4.format(ESAVEIN,10,true)+
			 " MeV";
		     //ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		 }
		 else
		 {
			 s=" RANGE REJECTION SWITCH";
		     ll=s.length();ll=61-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+"OFF";
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		 }
		 //WRITE(IOUT,260) TIMMAX,STATLM;
		 s=" Maximum cputime allowed";
	     ll=s.length();ll=60-ll;s=s+EGS4.format("",ll);
	     EGS4.seqStr=s+EGS4.format(TIMMAX,6,true)+" hrs";
	     if(EGS4.iprint>1)
		   	printSequence(EGS4.seqStr);
		 s=" Stats in cavity objective";
	     ll=s.length();ll=61-ll;s=s+EGS4.format("",ll);
	     EGS4.seqStr=s+EGS4.format(STATLM,6,true)+" %";
	     if(EGS4.iprint>1)
		   	printSequence(EGS4.seqStr);


		EGS4.seqStr=" Initial RNG state:";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.SHOW_RNG_STATE();

		if(RUSROU)
		{
			//WRITE(IOUT,265)RRZ,RRCUT;
			 s=" RUS ROU FOR PHOTONS CROSSING Z = ";
		     ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(RRZ,10,true)+" cm";
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s=" WITH PROBABILITY OF SURVIVAL:";
		     ll=s.length();ll=59-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(RRCUT,7,true);
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		}

		if(EGS4Macro.CEXPTR!=0)
		{
			 s=" PATHLENGTH EXPONENTIAL TRANSFORMATION";
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
			 s="  VARIABLE FOR FORWARD GOING PHOTNS: ";
		     ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
		     EGS4.seqStr=s+EGS4.format(EGS4Macro.CEXPTR,10,true);
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		}

		if(IFANO == 1)
		{
			 s="  *** REGENERATION REQUESTED (IFANO SET TO 1) ! *** ";
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		}
		else if(IFANO == 2)
		{
			 s="  *** ELECTRONS SET IN MOTION IN WALL WILL BE ELIMINATED (IFANO SET TO 2) ! *** ";
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);
		}
		else
		{
			 s="  *** NO REGENERATION REQUESTED (IFANO SET TO 0) ! *** ";
		     EGS4.seqStr=s;
		     if(EGS4.iprint>1)
			   	printSequence(EGS4.seqStr);

		}

		//"EK0=EIN;"
		//"$PRESTA-INPUT-SUMMARY; OUTPUT THE PRESTA INPUT VARIABLES"
		//"taken out input-summary at upgrade to PRESTA-II"

		//"MATERIAL INPUT SUMMARY"
		//"====================="
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="                   MATERIAL SUMMARY  "+EGS4.NMED+" MATERIALS USED";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="=========================================================================";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr=" # MATERIAL  DENSITY(g/cm**3)"+
		EGS4.format("",6)+"AE(MeV)"+EGS4.format("",4)+"AP(MeV)"+
		EGS4.format("",9)+"UE(MeV)"+EGS4.format("",4)+"UP(MeV)";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr="  - --------  ----------------"+
		EGS4.format("",6)+"-------"+EGS4.format("",4)+"-------"+
		EGS4.format("",9)+"-------"+EGS4.format("",4)+"-------";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);
		for(int I=1;I<=EGS4.NMED;I++)
		{
		    //WRITE(IOUT,310) I,(MEDIA(J,I),J=1,6),RHO(I),AE(I),AP(I),UE(I),UP(I);
		    //310    FORMAT(' ',I1,3X,6A1,4X,1PE10.3,2(7X,0PF9.3,2X,F9.3));
		    String meds="";
		    if(EGS4.MEDIA[I-1].length()>6)
		    {meds=EGS4.MEDIA[I-1].substring(0,6);}
		    else
		    {meds=EGS4.MEDIA[I-1];}
			EGS4.seqStr="  "+EGS4.format(I,3)+EGS4.format("",3)+meds+EGS4.format("",4)+
			EGS4.format(EGS4.RHO[I-1],10,false)+EGS4.format("",6)+EGS4.format(EGS4.AE[I-1],10)+
			EGS4.format("",1)+EGS4.format(EGS4.AP[I-1],10)+EGS4.format("",6)+
			EGS4.format(EGS4.UE[I-1],10)+EGS4.format("",1)+EGS4.format(EGS4.UP[I-1],10);
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
		}

		EGS4Geom.GEOMRZ_ISUMRY();

		EGS4SrcEns.SRCOUT();

		//"       PRINT A GRID OF THE ZONE DEPENDENT VARIABLES"
		//"       ============================================"

		//"Set defaults for non-used variables"
		//CDSTBL(1)='0';CTRTBL(1)='0';CABSRB(1)='0';
		EGS4Grid.CDSTBL[0]=" ";
		EGS4Grid.CTRTBL[0]="0";EGS4Grid.CABSRB[0]="0";
		//"For identifying a cavity region"
		//REPLACE {$IRL} WITH {IZ+1+NZ*(IX-1)}
		for(int IRL=2;IRL<= EGS4Geom.NR*EGS4Geom.NZ+1;IRL++)
		{
		    if (EGS4Geom.ntrack[IRL-1]==1)
		    {
				//EGS4Grid.CAVTRACK[IRL-1]="C";
				EGS4Grid.CDSTBL[IRL-1]="C";
			}
		    else
		    {
				//EGS4Grid.CAVTRACK[IRL-1]=" ";
				EGS4Grid.CDSTBL[IRL-1]=" ";
			}
		}
		//"Make the material grid"
		EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, AMASS, 1, EGS4.ECUT,
		EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);
		//,EGS4Grid.CAVTRACK, EGS4Grid.CTRTBL, EGS4Grid.CABSRB);

	}

	/**
	 * Internally used. Called by AUSGAB.
	 * @param P1 P1
	 * @param P2 P2
	 */
	private void EXCHANGE_STACK(int P1, int P2)
	{
		double FDUMMY=0.0;int IDUMMY=0;
		FDUMMY = EGS4.U[P2-1]; EGS4.U[P2-1] = EGS4.U[P1-1]; EGS4.U[P1-1] = FDUMMY;
		FDUMMY = EGS4.V[P2-1]; EGS4.V[P2-1] = EGS4.V[P1-1]; EGS4.V[P1-1] = FDUMMY;
		FDUMMY = EGS4.W[P2-1]; EGS4.W[P2-1] = EGS4.W[P1-1]; EGS4.W[P1-1] = FDUMMY;
		FDUMMY = EGS4.E[P2-1]; EGS4.E[P2-1] = EGS4.E[P1-1]; EGS4.E[P1-1] = FDUMMY;
		FDUMMY = EGS4.WT[P2-1]; EGS4.WT[P2-1] = EGS4.WT[P1-1]; EGS4.WT[P1-1] = FDUMMY;
		IDUMMY = EGS4.IQ[P2-1]; EGS4.IQ[P2-1] = EGS4.IQ[P1-1]; EGS4.IQ[P1-1] = IDUMMY;
		//"LATCH IS NOW STANDARD"
		IDUMMY = EGS4.LATCH[P2-1]; EGS4.LATCH[P2-1] = EGS4.LATCH[P1-1]; EGS4.LATCH[P1-1] = IDUMMY;
	}

	/**
	 * Print output summary
	 */
	private void OSUMRY()
	{
		String s="";int ll=0;

		//$IMPLICIT-NONE;

		//COMIN/GEOM,IODAT1,IODAT2,PRINTC,SCORE,SOURCE,USER/;
		//COMIN/CH-Steps/;

		int I=0;int IRL=0;int IX=0;int IZ=0;
		double ASCT=0.0;double ASCTUN=0.0;double AATT=0.0;double AATTUN=0.0;
		double AWLL=0.0;double AWLLUN=0.0;double TDAW=0.0;double TDAWUN=0.0;
		double KSCT=0.0;double KATT=0.0;double KWLL=0.0;double KSCTUN=0.0;
		double KATTUN=0.0;double KWLLUN=0.0;

		//"SET UP THE PRINTER"
		//ICHPIN=12; "12 CHARACTERS/INCH"
		//ILPIN=6;   "6 LINES/INCH"
		//IPAGE=0;   "NO PAGE THROW"
		//"CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

		//"WRITE(IOUT,100)TITLE,DATEN,TIMEN; HEADER"

		//IF(ISOURC=21 | ISOURC=22)[
		//    WRITE(IOUT,200) SCSTP,SCSTP2,
		//         SCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		//                       NINCSRC),SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		//]
		//ELSE[

		s="                    # primary charged particle steps";
		ll=s.length();ll=58-ll;s=s+EGS4.format("",ll);

		EGS4.seqStr=s+EGS4.format(SCSTP,10,false)+" +/- "+EGS4.format(SCSTP2,6)+"%";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		s="     # primary charged particle steps/initial history";
		ll=s.length();ll=58-ll;s=s+EGS4.format("",ll);

		EGS4.seqStr=s+EGS4.format(SCSTP/EGS4SrcEns.dble(IHSTRY),10,false)+" +/- "+
		EGS4.format(SCSTP2,6)+"%";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		s="# of presta-II steps/# primary charged particle steps";
		ll=s.length();ll=58-ll;s=s+EGS4.format("",ll);

		EGS4.seqStr=s+EGS4.format((EGS4.count_pII_steps+PIISTP)/SCSTP,10,false)+" +/- "+
		EGS4.format(SCSTP2,6)+"%";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);


		    //WRITE(IOUT,200) SCSTP,SCSTP2,SCSTP/dble(IHSTRY),
		    //              SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		//]
//		200  FORMAT(//' ' ,'                    # primary charged particle steps',T58,
//		             1PE10.3,' +/- ',0PF6.3,'%'/
//		             ' ','     # primary charged particle steps/initial history',T58,
//		             1PE10.3,' +/- ',0PF6.3,'%'/
//		             ' ','# of presta-II steps/# primary charged particle steps',
//		             T58,F10.3,' +/- ',0PF6.3,'%');


		//"PRINT # CHARGED PARTICLE STEPS IN cavity REGION, etc"
		//IF(ISOURC=21 | ISOURC=22)[
		//   WRITE(IOUT,210) SCCSTP,SCCSTP2,
		//          SCCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		//                      NINCSRC),SCCSTP2;
		//]
		//ELSE[
		s="   # primary charged particle steps in cavity region";
		ll=s.length();ll=58-ll;s=s+EGS4.format("",ll);

		EGS4.seqStr=s+EGS4.format(SCCSTP,10,false)+" +/- "+
		EGS4.format(SCCSTP2,6)+"%";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);

		s="    # primary steps in cavity region/initial history";
		ll=s.length();ll=58-ll;s=s+EGS4.format("",ll);

		EGS4.seqStr=s+EGS4.format(SCCSTP/EGS4SrcEns.dble(IHSTRY),10,false)+" +/- "+
		EGS4.format(SCCSTP2,6)+"%";
		if(EGS4.iprint>1)
			printSequence(EGS4.seqStr);


		    //WRITE(IOUT,210) SCCSTP,SCCSTP2,SCCSTP/dble(IHSTRY),SCCSTP2;
		//]
//		210  FORMAT(//' ','   # primary charged particle steps in cavity region'
//		             ,T58,1PE10.3,' +/- ',0PF6.3,'%'/
//		              ' ','    # primary steps in cavity region/initial history'
//		             ,T58,1PE10.3,' +/- ',0PF6.3,'%');

		if( EGS4SrcEns.ISOURC == 15 ) EGS4SrcEns.src15_out();//(iout);

		//"THE CAVITY SUMMARY"
		//"******************"

		if(EGS4Geom.NSUMCV==1)
		{
			EGS4.seqStr="                   SUM OF RESULTS FOR THE CAVITY: 1 REGION";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr="                   ***************************************";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			//WRITE(IOUT,300);
		}
		else
		{
			EGS4.seqStr="                   SUM OF RESULTS FOR THE CAVITY: " +
			EGS4Geom.NSUMCV+" REGION";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr="                   ***************************************";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

			//WRITE(IOUT,301)NSUMCV;
		}

//		300  FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/
//		            ' ',T20,'***************************************');
//		301  FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
//		            ' ',T20,'*****************************************');
		//"******************"
		//"history by history"
		//"    EMH March 2002"
		//"******************"
		cav2_dose = cav2_dose/cav_dose;
		if (EGS4Macro.iifano == 1)
		{
			EGS4.seqStr=" This calculation was performed using regeneration ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" ================================================= ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			s=" D/Awall (grays/incident fluence): ";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(cav_dose,11,false)+" +/- "+
			EGS4.format(100*cav2_dose,6)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		       //write(iout,'(a,t50,1PE11.4,a,0PF6.3,a)')
		       // 'D/Awall (grays/incident fluence): ',cav_dose,' +/-  ',
		       //  100*cav2_dose,'%';
		    return;
		}

		if( EGS4Macro.use_enhance )
		{
			EGS4.seqStr=" This calculation was performed using CS enhancement ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr="  enhancement factor was "+EGS4.format(EGS4Macro.cs_enhance,10,true);
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" ================================================= ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
		}
		else if ( EGS4Macro.n_split > 1 )
		{
			EGS4.seqStr=" This calculation was performed using photon splitting ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr="  splitting number was "+EGS4.format(EGS4Macro.n_split,6);
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" ================================================= ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
		}
		if (IFANO==1)
		{
			EGS4.seqStr=" This calculation was performed using regeneration ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		}
		if (IFANO==2)
		{
			EGS4.seqStr=" This calculation was performed eliminating electrons originating in the cavity wall.";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
		}
		//"THE TOTAL DOSE"
		if(EGS4SrcEns.ISOURC==3||EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22)
		{
			s=" TOTAL DOSE (GRAYS/INCIDENT PARTICLE):";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);

			EGS4.seqStr=s+
			EGS4.format(cav_dose,11,false)+" +/- "+EGS4.format(100*cav2_dose,5)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		   //WRITE(IOUT,311)cav_dose,100*cav2_dose;
		}
		else
		{
			s=" TOTAL DOSE (GRAYS/INCIDENT FLUENCE):";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);

			EGS4.seqStr=s+
			EGS4.format(cav_dose,11,false)+" +/- "+EGS4.format(100*cav2_dose,5)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

//		   WRITE(IOUT,310)cav_dose,100*cav2_dose;
		}
//		310  FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT FLUENCE):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
//		311  FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT PARTICLE):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');

		if(EGS4Macro.IFULL==1)
		{
		    //"CALCULATE Ascat, Aatt, Awall, DOSE/Awall"

		    cav2_dose1 = cav2_dose1/cav_dose1;
		    cav_dosec = cav_dosec/cav_dose/cav_dose1;

		    //"total dose/Awall"
		    TDAW=cav_dose1;
		    TDAWUN=cav2_dose1;

			s="  TOTAL DOSE/Awall:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(TDAW,11,false)+" +/- "+EGS4.format(100*TDAWUN,5)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    //WRITE(IOUT,350)TDAW,TDAWUN*100.;
//		350  FORMAT(' TOTAL DOSE/Awall:',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');

		    cav_dosec = cav2_dose*cav2_dose + cav2_dose1*cav2_dose1 - 2*cav_dosec;
		    if( cav_dosec > 0 ) cav_dosec = Math.sqrt(cav_dosec);
		    //"Awall"
		    AWLL=cav_dose/cav_dose1;
		    AWLLUN=cav_dosec;
		    //WRITE(IOUT,340)AWLL,AWLLUN*100.;
//		340  FORMAT(' Awall:',T50,F8.5,'    +/- ',F5.2,'%');
			s=" Awall:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(AWLL,8,true)+"  +/- "+EGS4.format(100*AWLLUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1 )
		    {
		      corr_02 = corr_02/(cav_dose0+cav_dose2);

		      cav2_dose0 = cav2_dose0/cav_dose0;
		      cav2_dose2 = cav2_dose2/cav_dose2;

		      cav_dosec01 = cav_dosec01/cav_dose0/cav_dose1;
		      cav_dosec02 = cav_dosec02/cav_dose0/cav_dose2;

		      cav_dosec01 = cav2_dose0*cav2_dose0 + cav2_dose1*cav2_dose1 -
		                    2*cav_dosec01;
		      cav_dosec02 = cav2_dose0*cav2_dose0 + cav2_dose2*cav2_dose2 -
		                    2*cav_dosec02;

		      if( cav_dosec01 > 0 ) cav_dosec01 = Math.sqrt(cav_dosec01);
		      if( cav_dosec02 > 0 ) cav_dosec02 = Math.sqrt(cav_dosec02);

		      //"... we had rel.error(x), get now rel.error(1+x)"
		      cav_dosec02 = cav_dosec02/(1 + cav_dose0/cav_dose2);

		      //"Ascatt"
		      ASCT=1+cav_dose2/cav_dose0;
		      ASCTUN=cav_dosec02;
		      //WRITE(IOUT,320)ASCT,ASCTUN*100.;
		      ///		320  FORMAT(' Ascat:',T50,F8.5,'    +/- ',F5.2,'%');
			s=" Ascat:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(ASCT,8,true)+"  +/- "+EGS4.format(100*ASCTUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		      //"Aatt"
		      AATT=cav_dose0/cav_dose1;
		      AATTUN=cav_dosec01;
		      //WRITE(IOUT,330)AATT,AATTUN*100.;
//		330  FORMAT(' Aatt :',T50,F8.5,'    +/- ',F5.2,'%');
			s=" Aatt :";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(AATT,8,true)+"  +/- "+EGS4.format(100*AATTUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);


		      //"Dpr + Dsec"
//		345  FORMAT(' Dprimary + Dsecondary:',T50,1PE11.4,' +/- ',0PF5.2,'%');
		      //WRITE(IOUT,345)cav_dose2+cav_dose0,corr_02*100.;
			s=" Dprimary + Dsecondary:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(cav_dose2+cav_dose0,11,false)+"  +/- "+
			EGS4.format(corr_02*100.,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    }
		} //"END OF IFULL = 1"

		//"OUTPUT Kscat,Katt,Kwall,Kpn,Kfl THE INVERSES OF THE Ai's"
		if(EGS4Macro.IFULL>0)
		{
		    KWLL=1./AWLL;KWLLUN=AWLLUN/(AWLL*AWLL);
			s=" Kwall:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(KWLL,8,true)+"  +/- "+EGS4.format(100*KWLLUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    //WRITE(IOUT,341)KWLL,KWLLUN*100.;
		    //		341  FORMAT(' Kwall:',T50,F8.5,'    +/- ',F5.2,'%');
		    if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1 )
		    {
		      KSCT=1./ASCT;KSCTUN=ASCTUN/(ASCT*ASCT);
			s=" Kscat:";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(KSCT,8,true)+"  +/- "+EGS4.format(100*KSCTUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		      //WRITE(IOUT,321)KSCT,KSCTUN*100.;
		      ///		321  FORMAT(' Kscat:',T50,F8.5,'    +/- ',F5.2,'%');
		      KATT=1./AATT;KATTUN=AATTUN/(AATT*AATT);
		      //WRITE(IOUT,331)KATT,KATTUN*100.;
//		331  FORMAT(' Katt :',T50,F8.5,'    +/- ',F5.2,'%');
			s=" Katt :";
			ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr=s+EGS4.format(KATT,8,true)+"  +/- "+EGS4.format(100*KATTUN,5,true)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    }
		}

		//"THE DETAILED OUTPUT FOR EACH CAVITY ZONE"
		//"****************************************"
		if(IOOPTN==1 && EGS4Geom.NSUMCV>1)
		{
		    //"ONLY IF REQUESTED AND MORE THAN ONE CAVITY ZONE"

		    //WRITE(IOUT,400) NSUMCV;
		    //		400  FORMAT(// ,T15,'DETAILED RESULTS FOR EACH OF THE ',I4,' CAVITY REGIONS'/
///		            ' ',T15,'****************************************************');
			EGS4.seqStr="   DETAILED RESULTS FOR EACH OF THE "+EGS4.format(EGS4Geom.NSUMCV,4)+
			" CAVITY REGIONS";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    //"PRINT THE TABLE HEADER"
		    if(EGS4Macro.IFULL==0)
		    {
		        //WRITE(IOUT,410);
//		410  FORMAT(//'Z# P# C#     Total Dose     '/
///		            ' -- -- -- -------------------');
			EGS4.seqStr=" Z# P# C#     Total Dose     ";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" -- -- -- -------------------";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    }
		    else if(EGS4Macro.IFULL>0)
		    {
		        //"SET UP THE PRINTER"
		        //ICHPIN=16; "16 CHARACTERS/INCH"
		        //ILPIN=6;   "6 LINES/INCH"
		        //IPAGE=0;   "NO PAGE THROW"
		        //"CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

		        //WRITE(IOUT,420);
///		420  FORMAT(//'Z# P# C#',
//		            '     Total Dose   ',
//		            '      Ascat      ',
//		            '       Aatt      ',
//		            '      Awall      ',
//		            '  Total Dose/Awall'/
//		            ' -- -- --',
//		            '     ----------   ',
//		            '      -----      ',
//		            '       ----      ',
//		            '      -----      ',
//		            '  ----------------');
			EGS4.seqStr=" Z# P# C#"+"     Total Dose   "+"      Ascat      "+
			"       Aatt      "+"      Awall      "+"  Total Dose/Awall";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr=" -- -- --"+"     ----------   "+"      -----      "+
			"       ----      "+"      -----      "+"  ----------------";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		    }

		    //"LOOP OVER THE CAVITY ZONES"
		    for( I=1;I<=EGS4Geom.NSUMCV;I++)
		    {
		        IRL=EGS4Geom.ISUMCV[I-1];
		        //$GET-IX-IZ(IRL);
				IX=EGS4Geom.GET_IX(IRL);
				IZ=EGS4Geom.GET_IZC(IRL);

		        //"total dose"
		        if(SCDOSE[IZ-1][IX-1][0]>0.)
		        	SCDOSE2[IZ-1][IX-1][0]=SCDOSE2[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][0];
		        if(EGS4Macro.IFULL==1)
		        {
		            //"CALCULATE Ascat, Aatt, Awall, DOSE/Awall"

		            if(SCDOSE[IZ-1][IX-1][2]>0.)
		            {
		                SCDOSE2[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX-1][2]/SCDOSE[IZ-1][IX-1][2];
		                SCDOSE_COV[IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][0]/
		                                    SCDOSE[IZ-1][IX-1][2];
		            }

		            //"TOTAL DOSE/Awall"
		            TDAW=SCDOSE[IZ-1][IX-1][2];TDAWUN=SCDOSE2[IZ-1][IX-1][2];

		            SCDOSE_COV[IZ-1][IX-1][0] = SCDOSE2[IZ-1][IX-1][0]*SCDOSE2[IZ-1][IX-1][0] +
		                                  SCDOSE2[IZ-1][IX-1][2]*SCDOSE2[IZ-1][IX-1][2] -
		                                  2*SCDOSE_COV[IZ-1][IX-1][0];
		            if(SCDOSE_COV[IZ-1][IX-1][0]>0.)SCDOSE_COV[IZ-1][IX-1][0]=
		                                            Math.sqrt(SCDOSE_COV[IZ-1][IX-1][0]);
		            //"Awall"
		            AWLL=SCDOSE[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][2];
		            AWLLUN=SCDOSE_COV[IZ-1][IX-1][0];

		            if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1 )
		            { //"this condition not strictly necessary because"
		              //"IOOPTN set = 0 if either of these options on"

		               if(SCDOSE[IZ-1][IX-1][1]>0.)
		               {
		                   SCDOSE2[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][1];
		                   SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][1]/
		                                           SCDOSE[IZ-1][IX-1][2];
		                   if(SCDOSE[IZ-1][IX-1][3]>0.)
		                   {
		                     SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]/SCDOSE[IZ-1][IX-1][1]/
		                                           SCDOSE[IZ-1][IX-1][3];
		                   }
		               }
		               if(SCDOSE[IZ-1][IX-1][3]>0.)SCDOSE2[IZ-1][IX-1][3]=SCDOSE2[IZ-1][IX-1][3]/
		                                                      SCDOSE[IZ-1][IX-1][3];

		               SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]*SCDOSE2[IZ-1][IX-1][1]+
		                                   SCDOSE2[IZ-1][IX-1][2]*SCDOSE2[IZ-1][IX-1][2]-
		                                   2*SCDOSE_COV[IZ-1][IX-1][1];
		               SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX-1][1]*SCDOSE2[IZ-1][IX-1][1]+
		                                   SCDOSE2[IZ-1][IX-1][3]*SCDOSE2[IZ-1][IX-1][3]-
		                                   2*SCDOSE_COV[IZ-1][IX-1][2];

		               if(SCDOSE_COV[IZ-1][IX-1][1]>0)SCDOSE_COV[IZ-1][IX-1][1]=
		                                              Math.sqrt(SCDOSE_COV[IZ-1][IX-1][1]);
		               if(SCDOSE_COV[IZ-1][IX-1][2]>0)SCDOSE_COV[IZ-1][IX-1][2]=
		                                              Math.sqrt(SCDOSE_COV[IZ-1][IX-1][2]);

		               //"... we had rel.error(x), get now rel.error(1+x)"
		               if(SCDOSE[IZ-1][IX-1][3]>0.) SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]/
		                                      (1+SCDOSE[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][3]);

		               //"Ascatt"
		               ASCT=1+SCDOSE[IZ-1][IX-1][3]/SCDOSE[IZ-1][IX-1][1];
		               ASCTUN=SCDOSE_COV[IZ-1][IX-1][2];

		               //"Aatt"
		               AATT=SCDOSE[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][2];
		               AATTUN=SCDOSE_COV[IZ-1][IX-1][1];

		            }
		        }
		        if(EGS4Macro.IFULL==0)
		        {
			EGS4.seqStr=" "+EGS4.format(IRL,2)+" "+EGS4.format(IZ,2)+" "+EGS4.format(IX,2)+" "+
			EGS4.format(SCDOSE[IZ-1][IX-1][0],11,false)+" +/-"+EGS4.format(100*SCDOSE2[IZ-1][IX-1][0],5)+"%";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		            //WRITE(IOUT,430)
		            //    IRL,IZ,IX,                     "POSITION"
		            //    SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100.; "TOTAL DOSE"
		        }

//		430  FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,' +/-',0PF5.2,'%');
//		440  FORMAT(' ',I2,2(1X,I2),
//		                1X,1PE11.4,'(',0PF5.2,'%)',
//		                3(F8.5,'(',F7.5,')'),
//		                1X,1PE11.4,'(',0PF5.2,'%)');

		        else if(EGS4Macro.IFULL==1)
		        {
			EGS4.seqStr=" "+EGS4.format(IRL,2)+" "+EGS4.format(IZ,2)+" "+EGS4.format(IX,2)+" "+
			EGS4.format(SCDOSE[IZ-1][IX-1][0],11,false)+"("+EGS4.format(100*SCDOSE2[IZ-1][IX-1][0],5)+"%)"+
			EGS4.format(ASCT,8,true)+"("+EGS4.format(ASCTUN,7,true)+")"+
			EGS4.format(AATT,8,true)+"("+EGS4.format(AATTUN,7,true)+")"+
			EGS4.format(AWLL,8,true)+"("+EGS4.format(AWLLUN,7,true)+")"+
			EGS4.format(TDAW,11,false)+"("+EGS4.format(100*TDAWUN,5)+"%)";
			if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

		            //WRITE(IOUT,440)
		            //    IRL,IZ,IX,                     "POSITION"
		            //    SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100., "TOTAL DOSE"
		            //    ASCT,ASCTUN,                   "Ascat"
		            //    AATT,AATTUN,                   "Aatt"
		            //    AWLL,AWLLUN,                   "Awall"
		            //    TDAW,TDAWUN*100;                   "TOTAL DOSE/Awall"
		        }

		    }// "END OF LOOP OVER CAVITY ZONES"
		}// "END OF DETAILED SUMMARY"

//		"RESET UP THE PRINTER"
//		ICHPIN=12; "12 CHARACTERS/INCH"
//		ILPIN=6;   "6 LINES/INCH"
//		IPAGE=0;   "NO PAGE THROW"
//		"CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

//		RETURN;
//		%I0
//		"FORMATS"
//		200  FORMAT(//' ' ,'                    # primary charged particle steps',T58,
//		             1PE10.3,' +/- ',0PF6.3,'%'/
///		             ' ','     # primary charged particle steps/initial history',T58,
//		             1PE10.3,' +/- ',0PF6.3,'%'/
///		             ' ','# of presta-II steps/# primary charged particle steps',
////		             T58,F10.3,' +/- ',0PF6.3,'%');
///		202   FORMAT(//' ' ,'                 # charged particle steps in run',T58,
///		             1PE10.3,/
///		             ' ','    # charged particle steps in run/initial history',T58,
//		             1PE10.3/
//		             ' ','# of presta-II steps/# primary charged particle steps',
///		             T58,F10.3,' +/- ',0PF6.3,'%');
//		210  FORMAT(//' ','   # primary charged particle steps in cavity region'
//		             ,T58,1PE10.3,' +/- ',0PF6.3,'%'/
//		              ' ','    # primary steps in cavity region/initial history'
//		             ,T58,1PE10.3,' +/- ',0PF6.3,'%');
//		220  FORMAT(// ,T8,'STEP COUNTING RESULTS FOR WALL MATERIAL IN THE CAVITY'/
///		             // ,T8,'# primary charged particle steps',T51,
///		             I12,' +/- ',0PF5.2,'%'/
//		             ' ',T8,'# OF Times mscat switched off',T51,
////		             I12,' +/- ',0PF5.2,'%'/
///		             ' ',T8,'RATIO',T54,F7.3,' +/- ',0PF5.2,'%');
//		230  FORMAT(// ,T8,'# Primary charged particle steps in cavity region'
///		             ,T51,I12,' +/- ',0PF5.2,'%'/
//		             ' ',T8,'# Times mscat switched off in cavity region.'
//		             ,T51,I12,' +/- ',0PF5.2,'%'/
//		             ' ',T8,'Ratio',T56,F7.3,' +/- ',0PF5.2,'%');
///		300  FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/
///		            ' ',T20,'***************************************');
//		301  FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
//		            ' ',T20,'*****************************************');
//		310  FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT FLUENCE):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
//		311  FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT PARTICLE):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
///		320  FORMAT(' Ascat:',T50,F8.5,'    +/- ',F5.2,'%');
///		321  FORMAT(' Kscat:',T50,F8.5,'    +/- ',F5.2,'%');
//		330  FORMAT(' Aatt :',T50,F8.5,'    +/- ',F5.2,'%');
//		331  FORMAT(' Katt :',T50,F8.5,'    +/- ',F5.2,'%');
//		340  FORMAT(' Awall:',T50,F8.5,'    +/- ',F5.2,'%');
//		341  FORMAT(' Kwall:',T50,F8.5,'    +/- ',F5.2,'%');
//		345  FORMAT(' Dprimary + Dsecondary:',T50,1PE11.4,' +/- ',0PF5.2,'%');
//		350  FORMAT(' TOTAL DOSE/Awall:',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
///		360  FORMAT(//'Apn  :',T50,F8.5,'   +/- ',F8.5);
//		361  FORMAT(' Kpn  :',T50,F8.5,'   +/- ',F8.5);
//		370  FORMAT(' TOTAL DOSE/([Apn]*Awall):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
//		380  FORMAT(//'Afl  :',T50,F8.5,'   +/- ',F8.5);
//		381  FORMAT(' Kfl  :',T50,F8.5,'   +/- ',F8.5);
//		390  FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall):',
//		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
//		395  FORMAT(//'<s>g,w  :',T50,F8.5,'   +/- ',F8.5);
///		396  FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall*<s>g,w):',
///		              T50,1PE11.4,' +/-  ',0PF5.2,'%');
//		400  FORMAT(// ,T15,'DETAILED RESULTS FOR EACH OF THE ',I4,' CAVITY REGIONS'/
///		            ' ',T15,'****************************************************');
//		410  FORMAT(//'Z# P# C#     Total Dose     '/
///		            ' -- -- -- -------------------');
///		420  FORMAT(//'Z# P# C#',
//		            '     Total Dose   ',
//		            '      Ascat      ',
//		            '       Aatt      ',
//		            '      Awall      ',
//		            '  Total Dose/Awall'/
//		            ' -- -- --',
//		            '     ----------   ',
//		            '      -----      ',
//		            '       ----      ',
//		            '      -----      ',
//		            '  ----------------');
//		430  FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,' +/-',0PF5.2,'%');
//		440  FORMAT(' ',I2,2(1X,I2),
//		                1X,1PE11.4,'(',0PF5.2,'%)',
//		                3(F8.5,'(',F7.5,')'),
//		                1X,1PE11.4,'(',0PF5.2,'%)');
//		445  FORMAT(' ',I2,2(1X,I2),
//		                1X,1PE11.4,'(',0PF5.2,'%)',
//		                3(F8.5,'(',F7.5,')'),
//		                1X,1PE11.4,'(',0PF5.2,'%)',
//		                F8.5,'(',F7.5,')',
//		                1X,1PE11.4,'(',0PF5.2,'%)');
///		450  FORMAT(' ',I2,2(1X,I2),
//		                1X,1PE11.4,
//		                3(1X,F8.5,1X),
//		                1X,1PE11.4,
//		                1X,F8.5,1X,
//		                1X,1PE11.4,
//		                1X,F8.5,1X,
//		                1X,1PE11.4,
//		                1X,F8.5,1X,
///		                1X,1PE11.4/
////		            ' ','ERRORS = ',
//		                4X,'(',0PF5.2,'%)',
///		                3(1X,'(',F7.5,')'),
//		                4X,'(',0PF5.2,'%)',
///		                1X,'(',F7.5,')',
//		                4X,'(',0PF5.2,'%)',
//		                1X,'(',F7.5,')',
///		                4X,'(',0PF5.2,'%)',
///		                1X,'(',F7.5,')',
//		                4X,'(',0PF5.2,'%)');
//
///		END; "END OF SUBROUTINE OSUMRY"
//

	}
}
