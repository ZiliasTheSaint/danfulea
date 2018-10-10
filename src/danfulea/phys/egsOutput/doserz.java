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
 * This code simulates the passage of an electron or photon beam in a finite, right cylindrical geometry. 
 * It also scores pulse height distributions in an arbitrary volume made up of any number of regions. 
 * The energy deposited within various user defined regions is scored and analyzed statistically following the simulation. 
 * 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 07 NOV. 2005
 */
public class doserz implements EgsQuestion {
	

	// "******************************************************************************
	// "
	// "
	// " ********************
	// " * *
	// " * dosrznrc.mortran *
	// " * *
	// " ********************
	// "
	// "
	// " INTRODUCTION:
	// " This code simulates the passage of an electron or photon beam in a
	// " finite, right cylindrical geometry.
	// "
	// " It also scores pulse height distributions
	// " in an arbitrary volume made up of any number of regions.
	// " There is a write to unit 2 of the calculated spectra/response functions
	// "
	// " The energy deposited within various user defined regions is scored and
	// " analyzed statistically following the simulation.
	// ;
	// " USER CONTROLS:
	// " The user 1) defines the geometry of the target via the input of a
	// " number of planar and cylindrical coordinates which
	// " divide the cylinder into a number of regions, each region
	// " composed of a user specified material.
	// " 2) specifies in which of these regions the dose is to be
	// " scored.
	// " 3) selects the form and degree of detail of the output.
	// " 4) selects either the energy if a monoenergetic beam is to be
	// " used or specifies an energy spectrum consisting of
	// " energy points and corresponding probabibities
	// " 5) selects the source of radiation from amongst parallel and
	// " point sources originating from the side or front,
	// " isotropically radiating sources embedded in regions
	// " in the cylinder and phase space files from BEAM.
	// " 6) selects the number of histories, time limit and statistical
	// " limit. all histories will be run unless time runs out or
	// " the variance calculated in the peak region drops below
	// " the statistical limit.
	// " 7) selects transport controls such as the fractional energy
	// " loss per charged particle step, the maximum step size,
	// " particle energy cutoffs, range rejection parameters.
	// "
	// " SIMULATION:
	// " Each particle emitted from the source is transported along with all its
	// " offspring (with energies greater than the cutoff) using EGS4.
	// "
	// " SCORING:
	// " In each of the defined dose scoring regions, the total dose and the
	// " dose less that due to stopped/discarded particles is scored. In
	// " addition, the dose scored in a region can be broken down into that
	// " due to particles which entered via the front wall, back wall, outer
	// " wall and inner wall. The total number of charged particle steps
	// " taken and the total number taken in the dose scoring region
	// " are counted.
	// "
	// " OUTPUT:
	// " The first section of the output echos the user input.
	// " The second section is is composed of messages and information printed
	// " during the simulation as a means of monitoring the simulation.
	// " The third section details the total doses accumulated with
	// " accompanying statistical uncertainties.
	// " The first and third sections may include grid format summaries if
	// " requested by the user.
	// " Plot files for the xmgr program are also output.
	// ;

	// " DEFINITIONS OF REGION NUMBER, PLANAR ZONE, CYLINDRICAL ZONE
	// " ===========================================================
	// " Z axis = axis of rotation runs across page shown as .......
	// "
	// "
	// " --------------------------------------------------------- RCYL(NR)
	// " |(NR-1) |(NR-1) |(NR-1) | . . . . | NR*NZ | NR*NZ | IX=NR
	// " | *NZ+2 | *NZ+3 | *NZ+4 | | | +1 |
	// " --------------------------------------------------------- RCYL(NR-1)
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " --------------------------------------------------------- RCYL(2)
	// " | NZ+2 | NZ+3 | NZ+4 | . . . . | 2NZ | 2NZ+1 | IX=2
	// " --------------------------------------------------------- RCYL(1)
	// "..1....|...2...|...3...|...4...|...............|...NZ..|..NZ+1.|....IX=1..1..
	// ;" ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " ---------------------------------------------------------
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " | | | | | | |
	// " ---------------------------------------------------------
	// " IZ=1 IZ=2 IZ=3 IZ=NZ-1 IZ=NZ

	// "*******************************************************************************
	// " INPUT/OUTPUT CONTROL INPUT
	// " **************************
	// "*******************************************************************************
	// "
	// " I/O DELIMETERS: :start I/O control:
	// " :stop I/O control:
	// "
	// " IWATCH= off (0) for normal output
	// " = interactions(1) output info on every discrete interaction
	// " = steps (2) output on every electron/photon step as well
	// " = deposited (3) output when energy is deposited as well
	// " = graph (4) outputs file for EGS_Windows graphics
	// " [IWATCH]
	// "
	// " STORE INITIAL RANDOM NUMBERS
	// " = no (0) DO NOT STORE THE INITIAL RANDOM NUMBERS
	// " = last (1) STORE THE INITIAL RANDOM NUMBER FOR THE LAST HISTORY
	// " = all (2) STORE ALL THE INITIAL RANDOM NUMBERS
	// "
	// " IRESTART
	// " = first (0) first run for this input file
	// " = restart (1) restart of a previous run, i.e. add more histories
	// " = analyze (3) just read in the raw data and do the statistical
	// //" analysis( gives no timing - restart 100 histories
	// " to get the same effect and more info)
	// " = start-RNS (4) read starting random numbers from a filE (e.g. FOR
	// " OUTPUT TO A GRAPHICS PACKAGE)
	// " = parallel (5) post-process distributed runs (all files
	// " named <filenamebase>_w#)
	// " [IRESTART ]
	// "
	// " OUTPUT OPTIONS
	// " = short (0) short output -just dose grid(DG)
	// " = dose summary (1) output dose summary only (DS)
	// " = material summary (2) output material summary grid(MG) & DG
	// " = material and dose summary (3) output MG + DS
	// " = long (4) output MG + DS + DG
	// " Note-any time there is a dose summary, there
	// " is also an .egs4dose file, unit 10
	// " [IOOPTN]
	// "
	// " STORE DATA ARRAYS
	// " = yes (0) Store data arrays dor re-use
	// " = no (1) don't store them
	// " [IDAT]
	// "
	// " ELECTRON TRANSPORT
	// " = normal (0) normal electron transport (discrete interactions)
	// " = no interactions (1) no discrete interactions (used for CDSA
	// " calculations but note that special data
	// " sets are also needed to do a full CSDA
	// " calculation. All turning off
	// " interactions does is just that. See use
	// " of IUNRST=2,3,4 PEGS4 data sets for real CSDA)
	// " [ ICSDA]
	// "
	// " For complex geometries you may want to just output
	// " a few regions. These are defined here. Note that
	// " proper transport is done everywhere and this is
	// " NOT related to range rejection.
	// "
	// " DOSE ZBOUND MIN (I) Minimum plane # defining dose region (default=1)
	// " [NZDMIN]
	// " DOSE ZBOUND MAX (I) Maximum plane # defining dose region
	// " [NZDMAX]
	// " DOSE RBOUND MIN (I) Minimum cylinder # defining dose region (default=0)
	// " [NRDMIN]
	// " DOSE RBOUND MAX (I) Maximum cylinder # defining dose region
	// " [NRDMAX]
	// ;
	// "
	// "*******************************************************************************
	// " MONTE CARLO CONTROL INPUT
	// " *************************
	// "*******************************************************************************
	// "
	// " MONTE CARLO DELIMETERS: :start Monte Carlo inputs:
	// " :stop Monte Carlo inputs:
	// "
	// " NUMBER OF HISTORIES (I) # HISTORIES TO RUN
	// " (MIN:100, DEFAULTS TO 20 000,
	// " max=1073741824=2^30 max is a compiler
	// " restriction for 32bit integers
	// " [NCASE]
	// " INITIAL RANDOM NO. SEEDS= INTGER1, INTEGER2
	// " User-code can use RANLUX or RANMAR, depending on selection
	// " Default is RANLUX
	// " RANLUX
	// " INTEGER1 is the luxury level, use 1 to 4, 4 taking longest
	// " default is 1 (set by $DEFAULT-LL in ranlux.macros)
	// " INTEGER2 selects the independent sequence to use, it
	// " can be from 1 to 1073741824 (2**30)
	// " RANMAR
	// " INTEGER1 is a seed between 1 and 31328 (0 =>default 1802)
	// " INTEGER2 is a seed between 1 and 30081 (0 =>default 9937)
	// " Selection of unique INTEGER2 values guarantees independent
	// " sequences.
	// " Note:
	// " After the seeds are first input and used for initialization,
	// " the variables IXXIN and JXXIN are just pointers used by the RNG
	// "
	// " MAX CPU HOURS ALLOWED (I) max CPU time allowed in hours, default=999
	// " [TIMMAX]
	// "
	// " IFULL
	// " = dose and stoppers (0) just calculate total dose and that due
	// " to stoppers and discards.
	// " = entrance regions (1) as well analyze the total dose per
	// " entrance region.
	// " = pulse height distribution (2) score a pulse height distribution in
	// " the volume specified after the
	// " material inputs.
	// " = scatter fraction (3) score the scatter fraction instead of
	// " stoppers. Only for incident photons.
	// " Dose after Compton and for fluorescent
	// " photons if followed.
	// "
	// " STATISTICAL ACCURACY SOUGHT (R) % statistical accuracy of the total
	// " dose in the peak region that is sought
	// " The program executes until this
	// " accuracy is obtained or the CPU time
	// " runs out. If 0, no effect.
	// "
	// " SCORE KERMA
	// " = no (0) do not score kerma
	// " = yes (1) score kerma wherever dose scored and estimate
	// " ratio of dose/kerma using correlated uncertainty
	// " estimate.
	// " This only makes sense for photon beams.
	// " [IKERMA]
	// "
	// ;
	// "*******************************************************************************
	// " CYLINDRICAL GEOMETRY INPUT
	// " **************************
	// "*******************************************************************************
	// "
	// " GEOMRZ DELIMETERS: :start geometrical inputs:
	// " :stop geometrical inputs:
	// "
	// " METHOD OF INPUT
	// " = Groups (0) input groups of slabs of equal thickness
	// " = Individual (1) detailed input of the geometry and media.
	// " [ITERSE]
	// "
	// //"-------------------------------------------------------------------------------
	// "
	// " Information defining depth boundaries along z axis (all dimensions cm)
	// ////"
	// " Only if METHOD OF INPUT= Groups
	// "
	// " Z OF FRONT FACE (R) start of first slab (real)
	// " NSLAB (M) # planar slabs in a group (integers)
	// " SLAB THICKNESS (M) thickness of each slab in the group (reals)
	// "
	// "
	// " Only if METHOD OF INPUT= Individual
	// "
	// " Z OF FRONT FACE (R) start of first slab (real)
	// " DEPTH BOUDARIES (M) geometrical z-plane coordinates (reals)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " Information defining radial boundaries
	// "
	// " RADII (M) radii of cylinders defining the geometry (reals)
	// "
	// "-------------------------------------------------------------------------------
	// " MATERIAL INPUT
	// " **************
	// "
	// " MEDIA (M) material name which must match that in the
	// " pegs4 data set EXACTLY, including case.
	// " 24 characters max per medium, ended by , or ;
	// "
	// " Define which media in which regions, numbering in order given above.
	// "
	// " DESCRIPTION BY= Regions(0) using the geometric region number
	// " = Planes (1) USING THE IX, IZ PLANES
	// " [DESCRIBE]
	// "
	// " If DESCRIPTION BY= Regions
	// "
	// " MEDNUM (M) the material number (integers)
	// " (MEDNUM=0 => vacuum)
	// " START REGION (M) initial geometrical zone(irl) (integers) for
	// " this medium [NREGLO]
	// " STOP REGION (M) final geometrical zone(irl) (integers) for
	// " this medium.[NREGHI]
	// " ( >NREGLO to input more than one zone)
	// " DEFAULTS: MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	// " MEDNUM=1 FOR REGION=2,NREG
	// "
	// " These inputs should be thought of as triplets of
	// " MEDNUM,START and STOP REGIONs which are used
	// " to specify the medium numbers for all regions where
	// " the medium is not the default (medium 1).
	// "
	// " If DESCRIPTION BY= Planes
	// "
	// " MEDNUM (M) the material number (integers)
	// " (MEDNUM=0 => vacuum)
	// " START ZSLAB (M) initial zslab (iz) (integers)
	// " STOP ZSLAB (M) final zslab (iz) (integers)
	// " START RING (M) initial radial ring (ix) (integers)
	// " STOP RING (M) final radial ring (ix) (integers)
	// " DEFAULTS: MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	// " MEDNUM=1 FOR REGION=2,NREG
	// " One must use one type of input or the other so you must decide
	// " which is more convenient for any given case.
	// "
	// "
	// ;
	// //"*******************************************************************************
	// " PULSE HEIGHT DISTRIBUTION INPUTS (ONLY IF IFULL=2)
	// " ********************************
	// "*******************************************************************************
	// "
	// " (ONLY IF IFULL= pulse height distribution)
	// "
	// " PULSE HEIGHT DISTRIBUTION DELIMETERS:
	// " :start pulse height distribution input:
	// " :stop pulse height distribution input:
	// "
	// "
	// " REGION OF SENSITIVE VOLUME (M) Region numbers(IRL) of sensitive volume
	// "
	// "
	// " SLOTE (R) for the pulse height distribution, defines the
	// " output energy bins. for SLOTE > 0.0, use equal size
	// " bins of this width in MeV (this will get increased
	// " by factors of two until the whole spectrum input
	// " gets covered.
	// " SLOTE < 0.0, flags input additional 'TOPS OF ENERGY BINS'
	// "
	// " DELTAE (R) code analyses peak efficiencies using a bin width
	// " of 2*deltae about each peak and two background regions
	// " of width deltae above and below the peak. default
	// " value is 0.005 MeV. (meaningless for electrons and positrons
	// "
	// " IF SLOTE < 0.0 Input tops of individual energy bins for pulse height
	// distn.
	// "
	// " TOPS OF ENERGY BINS (M) lowest energy first, tops of bins.
	// " [BINTOP]
	// "
	// "
	// "*******************************************************************************
	// " SOURCE INPUT(check latest version of srcrz)
	// " **************
	// "*******************************************************************************
	// " SOURCE DELIMETERS: :start source inputs:
	// " :stop source inputs:
	// "
	// "FOR ALL SOURCES
	// " Charge of the incident beam
	// " INCIDENT PARTICLE= electron (-1) electrons
	// " photon (0) photons
	// " positron (1) positrons
	// "
	// " (for SOURCE 21,22,23) all (2) include all of the particles
	// " in the phase space file
	// " [IQIN]
	// " charged (3) include e+ and e-
	// "
	// " SOURCE NUMBER (I) number of the source
	// " [ISOURC]
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 0 <<<<<<<<
	// "
	// "     PARALLEL BEAM INCIDENT FROM THE FRONT (+VE Z-AXIS) //"//"toc:
	// "
	// "
	// " SOURCE OPTIONS (M4) RBEAM, UINC, VINC, WINC
	// "
	// " RBEAM radius of parallel beam in cm
	// " (defaults to max radius of geometry)
	// " UINC incident x-axis direction cosine
	// " VINC incident y-axis direction cosine
	// " WINC incident z-axis direction cosine
	// " NOTE: (UINC,VINC,WINC)
	// " get automatically normalized
	// " defaults to (0.0,0.0,1.0)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 1 <<<<<<<<
	// "
	// "     POINT SOURCE ON AXIS INCIDENT FROM THE FRONT    //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) DISTZ, RBEAM, 0, 0
	// "
	// " DISTZ distance of the point source from the
	// " front of the target in cm (DEFAULT 100.)
	// " RBEAM radius of the beam at the front of the
	// " target in cm (defaults to MAX radius)
	// "
	// "------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 2 <<<<<<<<
	// "
	// "        BROAD PARALLEL BEAM INCIDENT FROM FRONT (+VE Z-AXIS) //"//"toc:
	// " WITH UNIT AREA BEAM AND LARGE SCORING AREA
	// "
	// " SOURCE OPTIONS (M4) 0, 0, 0, 0
	// "
	// "------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 3 <<<<<<<<
	// "
	// "     UNIFORM ISOTROPICALLY RADIATING DISK OF FINITE SIZE   //"//"toc:
	// " (MUST BE ALLOWED FOR IN THE GEOMETRICAL DEFINITIONS)
	// "
	// " SOURCE OPTIONS (M4) RMINBM, RBEAM, ZSMIN, ZSMAX
	// "
	// " RMINBM,RBEAM inner and outer radii of source region
	// " must be inside geometry
	// " ZSMIN,ZSMAX min and max z values for source
	// "
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 4 <<<<<<<<
	// "
	// "     FOR CENTRAL AXIS FLUENCE VS BEAM RADIUS      //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) RCAXIS, 0, 0, 0
	// "
	// " RCAXIS radius of central axis scoring zone (cm)
	// "
	// " NOTE: this source option treats the cylindrical radii input
	// " above as beam radii. the largest radius must be infinite
	// " and the phantom must be homogeneous (at least in each layer)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 10 <<<<<<<<
	// "
	// "     PARALLEL BEAM INCIDENT FROM THE SIDE (+VE Y-AXIS)    //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) XBEAM, ZBEAM, 0, 0
	// "
	// " XBEAM half-width of the rectangular beam in cm
	// " (defaults to max radius)
	// " ZBEAM half-height of the rectangular beam in cm
	// " (defaults to max)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 11 <<<<<<<<
	// "
	// "     POINT SOURCE INCIDENT FROM THE SIDE     //"//"toc:
	// "
	// "
	// " SOURCE OPTIONS (M4) DISTRH, XBEAM, ZBEAM, 0
	// "
	// " DISTRH distance of the source from the middle
	// " of the target in cm (defaults to 100.)
	// " XBEAM half-width of the beam at the center of
	// " the target in cm (defaults to max radius)
	// " ZBEAM half-height of the beam at the center of
	// " the target in cm (defaults to max)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 12 <<<<<<<<
	// "
	// "   POINT SOURCE OFF AXIS         //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) DISTRH, DISTZ, 0, 0
	// "
	// " DISTRH distance of the point source off the
	// " Z-axis.
	// " DISTZ perpendicular distance of the
	// " point source away from the front face.
	// " a negative value is permitted.
	// "
	// " DISTZ > 0
	// " point located in front of front face
	// "
	// " 0 > DISTZ > -(ZPLANE(NPLANE)-ZPLANE(1))
	// " point located between front and rear face
	// "
	// " DISTZ < -(ZPLANE(NPLANE)-ZPLANE(1))
	// " point located rear of rear plane
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 13 <<<<<<<<
	// "
	// "        PARALLEL BEAM FROM ANY ANGLE     //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) UINC, VINC, WINC, 0
	// "
	// " UINC incident x-axis direction cosine
	// " VINC incident y-axis direction cosine
	// " WINC incident z-axis direction cosine
	// "
	// " NOTE: (UINC,VINC,WINC) get automatically normalized
	// " default is (0.0,0.0,1.0)
	// "
	// "
	// "------------------------------------------------------------------------------
	// " >>>>>>>> SOURCE 14 <<<<<<<<
	// "
	// "   POINT SOURCE ON AXIS INCIDENT FROM THE FRONT WITH ALL   //"//"toc:
	// " EVENTS INSIDE RMINBM NOT FOLLOWED (A FUDGE FOR COLLIMATOR STUDIES)
	// "
	// " SOURCE OPTIONS (M4) DISTZ, RBEAM, RMINBM, IGNORED
	// "
	// " DISTZ distance of the point source from the
	// " front of the target in cm
	// " (defaults to 100.)
	// " RBEAM radius of the beam at the front of the
	// " target in cm (defaults to max radius)
	// " RMINBM below this radius, all histories are
	// " terminated by the source routines by
	// " giving them zero weight.
	// " The howfar routines must check for this.
	// "
	// "-------------------------------------------------------------------------------
	// "
	// "
	// "
	// " >>>>>>>> SOURCE 15 <<<<<<<<
	// "
	// " POINT SOURCE OFF AXIS. The same as source 12 but uses an alternative
	// " implementation for sampling points on the surface of the RZ-geomtry.
	// The
	// " motivation for implementing this source was to check that source 12 is
	// OK
	// " and to check the effect of varying weights from the source on the
	// " statistical uncertainty (contrary to source 12, source 16 produces
	// " essentially constant weights if the geometry-to-source distance is
	// large
	// " compared to the geometry dimension, a typical situation for ion chamber
	// " simulations)
	// "
	// " SOURCE OPTIONS (M4) DIST, ANGLE, IGNORED, IGNORED
	// "
	// " DIST distance of the centre of the geometry
	// " to the source in cm.
	// " ANGLE angle of rotation around the x-axis.
	// " (because of the cylindrical symetry,
	// " rotations around the x-axis and y-axis
	// " are indistinguishable). 0 degrees
	// " corresponds to a source above the front
	// " face (i.e. the same as source 1), 90
	// " degrees to a source from the side
	// " (i.e. the same as source 11).
	// " The source MUST be outside the geometry,
	// " otherwise the initialization routine
	// " will abort execution.
	// "
	// " Note that if you are not actually rotating about the center of the
	// " geometry, you must calculate the angle and distance as if you
	// " were.
	// "
	// "-------------------------------------------------------------------------------
	// "
	// "
	// "
	// " >>>>>>>> SOURCE 16 <<<<<<<<
	// "
	// " EXTENDED (CIRCULAR OR RECTANGULAR) SOURCE OFF AXIS.
	// "
	// " SOURCE OPTIONS (M4) DIST, ANGLE, TMP1, TMP2
	// "
	// " DIST distance of geometry centre to source
	// " centre in cm.
	// "
	// " ANGLE angle of rotation around the x-axis
	// " (see comments/explanations to source 15)
	// "
	// " TMP1, if TMP2 <= 0 radius of the source (i.e., the emitting
	// " or TMP2, if TMP1 <= 0 position is picked uniformly within the
	// " circle).
	// "
	// " TMP1 and TMP2, if both half-sizes of the radiating rectangle
	// " >= 0 in x- and y-directions before rotation,
	// " i.e., initially x and y are picked
	// " within the rectangle and z is set to
	// " -DIST + geometry centre. Then a rotation
	// " around the x-axis is performed.
	// " In all cases the source plane is perpendicular to the line joining
	// " it to the center of the geometry. Note that this introduces a
	// " slight error if the center of your geometry is not the true point
	// " of rotation.
	// "
	// " Note: if TEMP1 <= 0 and TEMP2 <= 0, source 16 becomes a point-source
	// " off-axis, i.e. the same as source 12 and 15.
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 20 <<<<<<<<
	// "
	// "   RADIAL DISTRIBUTION INPUT      //"//"toc:
	// "
	// " MODEIN= Local (0) if radial distribution is to be input
	// " locally through the .egs4inp file
	// " = External (1) if the distribution is to be input
	// " via an external file
	// "
	// " -----------------------------
	// " ONLY IF MODEIN= Local
	// "
	// " NRDIST (I) # radial bins in distribution histogram
	// " RDISTF (M) top of radial bin.
	// " should be values for 1 to NRDIST.
	// " RPDF (M) Probability of initial particle being
	// " in this bin.
	// " Probability doesn't need to be normalized
	// " but it should be in units cm**-2
	// " Should be values for 1 to NRDIST.
	// " RDIST IOUTSP= None (0) No distribution data in output summary
	// " = Include (1) include distribution data output summary
	// "
	// " -----------------------------
	// " ONLY IF MODEIN= External
	// "
	// " RDIST FILENAME (C) filename(with ext) contains
	// " distribution information
	// "
	// " RDIST IOUTSP= None (0) No distribution data in output summary
	// " = Include (1) include distribution data output summary
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 21 <<<<<<<<
	// "
	// "    FULL BEAM PHASE-SPACE BEAM DATA, INCIDENT ON FRONT FACE    //"//"toc:
	// "
	// " SOURCE OPTIONS (M4) IMODE, 0, 0, 0
	// "
	// " IMODE 0=> 7 variables/record: X,Y,U,V,E,WT,LATCH
	// " 2=> 8 variables/record: the above + ZLAST
	// "
	// " FILSPC (C) filename (with ext) contains
	// " phase space information
	// " (maximum of 80 characters)
	// " (assigned to unit 42)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " >>>>>>>> SOURCE 22 <<<<<<<<
	// "
	// "    FULL BEAM PHASE-SPACE BEAM DATA FROM ANY ANGLE, INSIDE OR OUTSIDE   //"//"toc:
	// "
	// " PARTICLES ARE READ IN FROM A BEAM PHASE SPACE and placed on a plane
	// " described by the SOURCE OPTIONS inputs (see below). Then it is checked
	// " whether they are already inside the geometry. If yes, the region index
	// " is determined and the shower intiated. If not, it is checked whether
	// " the particle trajectory will intersect the geometry (assuming that the
	// " geometry is surrounded by vacuum). If not, the particle is rejected and
	// " the next one taken from the phase-space file. If yes, the particle
	// " is placed on the entry point and the shower is initiated.
	// "
	// " SOURCE OPTIONS (M4) IMODE, DIST, ANGLE, ZOFFSET
	// "
	// " IMODE 0=> 7 variables/record: X,Y,U,V,E,WT,LATCH
	// " 2=> 8 variables/record: the above + ZLAST
	// " DIST Perpendicular distance of the phase-space
	// " plane to the point of rotation in cm.
	// " ANGLE Angle of rotation in degrees. The rotation
	// " is performed around an axis that is parallel
	// " to the x-axis and passes through the point
	// " (x,y,z)=(0,0,ZOFFSET).
	// " ZOFFSET Point of rotation. If |ZOFFSET| > 1e4,
	// " the centre of the geometry is taken as
	// " the point of rotation (but note that
	// " the maximum value allowed by the input
	// " routine is 1e6, so that |ZOFFSET| must
	// " be between 1e4 and 1e6 to use the centre
	// " of the geometry automatically).
	// "
	// " Examples:
	// " - to place a phase-space on the upper z-face of the geometry,
	// " use DIST=0, ANGLE=0, ZOFFSET=zplane(1)
	// " This is the same as source 21
	// " - to place a phase space on the lower z-face of the geometry,
	// " use DIST=0, ANGLE=180, ZOFFSET=zplane(n)
	// " - to have a phase file incident from, say, 60 degrees with
	// " a distance to the centre of the geometry of 30 cm, use
	// " DIST=30, ANGLE=60, ZOFFSET=9999.
	// " etc.
	// "
	// " FILSPC (C) filename (with ext) contains
	// " phase space information
	// " (maximum of 80 characters)
	// " (assigned to unit 42)
	// "
	// ;
	// "-------------------------------------------------------------------------------
	// "
	// //" >>>>>>>> SOURCE 23 <<<<<<<<
	// "
	// "    BEAM TREATMENT HEAD SIMULATION AS SOURCE INCIDENT FROM AND ANGLE,  //"//"toc:
	// "    INSIDE OR OUTSIDE PHANTOM                                          //"//"toc:
	// "
	// " PARTICLES ARE READ DIRECTLY FROM A BEAM SIMULATION COMPILED AS A
	// " SHARED LIBRARY. Particles are read at the scoring plane in
	// " the BEAM simulation (although no phase space file is scored) and are
	// " tranlated/rotated by the inputs DIST, ANGLE, XOFFSET, YOFFSET, ZOFFSET,
	// " described below. Then it is checked
	// " whether they are already inside the geometry. If yes, the region index
	// " is determined and the shower intiated. If not, it is checked whether
	// " the particle trajectory will intersect the geometry (assuming that the
	// " geometry is surrounded by vacuum). If not, the particle is rejected and
	// //" the next one taken from the BEAM simulation (more histories are run
	// in
	// " the BEAM simulation if required). If yes, the particle
	// " is placed on the entry point and the shower is initiated.
	// "
	// " BEAM CODE (C) The name of the accelerator code being
	// " used as a source including the BEAM_
	// " prefix (ie BEAM_accelname). This code
	// " must have been compiled as a shared
	// " library (see the BEAM manual for more
	// " details) and exist as
	// " libBEAM_accelname.so (for Linux/Unix) or
	// " libBEAM_accelname.dll (for Windows) in
	// " directory $EGS_HOME/bin/config.
	// "
	// " INPUT FILE (C) The name of a working input file
	// " (no .egsinp extension) for
	// " the BEAM code BEAM_accelname. This
	// " input file must specify output of a
	// " phase space file at one scoring plane.
	// " Particles that would have been scored
	// " in the phase space file are extracted
	// " and used as the incident particles in
	// " the DOSXYZ simulation instead. The
	// " input file must exist in the directory
	// " $EGS_HOME/BEAM_accelname.
	// "
	// " PEGS FILE (C) The name of the pegs4 data set (no
	// " .pegs4dat extension) used
	// " by BEAM_accelname with the input file
	// " specified by INPUT FILE. The pegs4
	// " data set must exist in either
	// " $HEN_HOUSE/pegs4/data or in
	// " $EGS_HOME/pegs4/data.
	// "
	// " WEIGHT WINDOW (M2) MIN_WEIGHT_23, MAX_WEIGHT_23
	// "
	// " MIN_WEIGHT_23 Min. weight of particles to use from
	// " the BEAM simulation (defaults to -1E30)
	// " MAX_WEIGHT_23 Max. weight of particles to use from
	// " the BEAM simulation (defaults to 1E30)
	// "
	// " SOURCE OPTIONS (M5) DIST, ANGLE, ZOFFSET, XOFFSET, YOFFSET
	// "
	// " DIST Perpendicular distance of the phase-space
	// " plane to the point of rotation in cm.
	// " ANGLE Angle of rotation in degrees. The rotation
	// " is performed around an axis that is parallel
	// " to the x-axis and passes through the point
	// " (x,y,z)=(0,0,ZOFFSET).
	// " ZOFFSET Point of rotation. If |ZOFFSET| > 1e4,
	// " the centre of the geometry is taken as
	// " the point of rotation (but note that
	// " the maximum value allowed by the input
	// " routine is 1e6, so that |ZOFFSET| must
	// " be between 1e4 and 1e6 to use the centre
	// " of the geometry automatically).
	// " XOFFSET,YOFFSET X and Y offset of scoring plane in BEAM
	// " simulation (cm). Offsets are applied before
	// " rotating the source.
	// "
	// " Examples:
	// " - to have BEAM simulation incident on the upper z-face of the geometry,
	// " use DIST=0, ANGLE=0, ZOFFSET=zplane(1)
	// " This is the same as source 21
	// " - to have BEAM simulation incident on the lower z-face of the geometry,
	// " use DIST=0, ANGLE=180, ZOFFSET=zplane(n)
	// " - to have BEAM simulation incident from, say, 60 degrees with
	// " a distance to the centre of the geometry of 30 cm, use
	// " DIST=30, ANGLE=60, ZOFFSET=9999.
	// "
	// "
	// "*******************************************************************************
	// "
	// " Source Energy Inputs (not required if ISOURC=21,22,23 - phase space or
	// BEAM
	// " simulation )
	// "
	// " Input from ensrc.mortran
	// "
	// " ENSRC DELIMETERS: :start source inputs:
	// " :stop source inputs:
	// "
	// "
	// " INCIDENT ENERGY
	// " = monoenergetic (0) if monoenergetic beam
	// " = spectrum (1) if energy spectrum to be used
	// "
	// " ---------------------------------------
	// "
	// " If INCIDENT ENERGY= Monoenergetic:
	// "
	// " INCIDENT KINETIC ENERGY(MEV) (I)
	// " kinetic energy of the incident beam in MeV
	// " (defaults to 1.25)
	// "
	// " ---------------------------------------
	// "
	// " If INCIDENT ENERGY= Spectrum:
	// "
	// " SPEC FILENAME (C) filename (with ext)
	// " contains spectrum information
	// "
	// " FILE FORMAT:
	// " TITLE spectrum title (80 char)
	// " NENSRC, ENMIN, MODE
	// " NENSRC # energy bins in spec. histogram
	// " ENMIN lower energy of first bin
	// " MODE =0, assumes cts/bin
	// " =1 assumes cts/MeV
	// " ENSRCD(I),SRCPDF(I) I=1,NENSRC
	// " top of energy bin and probability of
	// " initial particle being in this bin.
	// " probability does not need to be normalized
	// "
	// " SPEC IOUTSP
	// " = none (0) no spectrum data in output summary
	// " = include (1) include spectrum data in output summary

	// "*******************************************************************************
	// "
	// " TRANSPORT CONTROL INPUT
	// " ***********************
	// "
	// " BE AWARE!!!!!!!!!!!!!!!!!!!!
	// " IK changed this section to make things more flexible.
	// " All mc transport parameter are handled in a separate
	// " routine, the delimeters are
	// " :start mc transport parameter:
	// " :stop mc transport parameter:
	// " Input associated with variance reduction,
	// " previously in the transport control section,
	// " is to be put between the delimeters
	// " :start variance reduction:
	// " :stop variance reduction:
	// "
	// "*******************************************************************************
	// "
	// " MC TRANSPORT PARAMETER
	// " **********************
	// "
	// " All input associated with selection of various transport parameter
	// " is not crucial for the execution as there are default values set.
	// " Therefore, if some of the input options in this section are
	// " missing/misspelled, this will be ignored and defualt parameter assumed
	// " As the transport parameter input routine uses get_inputs, a lot
	// " of error/warning messages may be produced on UNIT 15, though.
	// " If you don't have the intention of changing default settings,
	// " simply ignore the error messages.
	// "
	// " The delimeters are
	// "
	// " :start mc transport parameter:
	// " :stop mc transport parameter:
	// "
	// " Currently, the following options are available (case does not matter):
	// "
	// " Global ECUT= Set a global (in all regions) electron transport
	// " cut off energy (in MeV). If this imput is missing,
	// " AE(medium) will be used.
	// " [ ECUT ]
	// " Global PCUT= Set a global (in all regions) photon transport
	// " cut off energy (in MeV). If this imput is missing,
	// " AP(medium) will be used.
	// " [ PCUT ]
	// " Global SMAX= Set a global (in all regions) maximum step-size
	// " restriction for electron transport (in cm).
	// " If missing, no geometrical step-size restrictions will
	// " be employed. Note that if you use the default
	// " EGSnrc electron-step algorithm, no SMAX-restriction
	// " is necessary. Option is useful for transport in low
	// " density materials (air) when PRESTA behaviour is
	// " turned on (see below)
	// " [ SMAXIR ]
	// " ESTEPE= Set the maximum fractional energy loss per step.
	// " Note that this is a global option only, no
	// " region-by-region setting is possible. If missing,
	// " the defualt is 0.25 (25%)
	// " [ ESTEPE ]
	// " XImax= Maximum first elastic scattering moment per step.
	// " Default is 0.5, NEVER use value greater than 1 as
	// " this is beyond the range of MS data available.
	// " [ XIMAX ]
	// " Boundary crossing algorithm=
	// " There are two selections possible: EXACT, means
	// " the algorithm will cross boundaries in a single
	// " scattering (SS) mode, the distance from a boundary
	// " at which the transition to SS mode is made is
	// " determined by 'Skin depth for BCA' (see below).
	// " The second option is PRESTA-I, if selected boundaries
	// " will be crossed a la PRESTA, i.e. with lateral
	// " correlations turned off and MS forced at boundaries.
	// " Default is EXACT.
	// " [ bca_algorithm, exact_bca ]
	// " Skin depth for BCA=
	// " Determines the distance from a boundary (in elastic
	// " MFP) at which the algorithm will go into single
	// " scattering mode (if EXACT boundary crossing) or
	// " swith off lateral correlations (if PRESTA-I boundary
	// " crossing). Default value is 3 for EXACT or
	// " exp(BLCMIN)/BLCMIN for PRESTA-I (see the PRESTA paper
	// " for a definition of BLCMIN). Note that if you choose
	// " EXACT boundary crossing and set Skin depth for BCA
	// " to a very large number (e.g. 1e10), the entire
	// " calculation will be in SS mode. If you choose
	// " PRESTA-I boundary crossing and make Skin depth for BCA
	// " large, you will get default EGS4 behavious (no PRESTA)
	// " [ skindepth_for_bca ]
	// " Electron-step algorithm=
	// " PRESTA-II (the default), the name is
	// " used for historical reasons
	// " or PRESTA-I
	// " Determines the algorithm used to take into account
	// " lateral and longitudinal correlations in a
	// " condensed history step.
	// " [ transport_algorithm ]
	// " Spin effects= Off, On, default is Off
	// " Turns off/on spin effects for electron elastic
	// " scattering. Spin On is ABSOLUTELY necessary for
	// " good backscattering calculations. Will make a
	// " even in `well conditioned' situations (e.g. depth
	// " dose curves for RTP energy range electrons).
	// " [ spin_effects ]
	// " Brems angular sampling= Simple, KM, default is KM
	// " If Simple, use only the leading term of the Koch-Motz
	// " distribution to determine the emission angle of
	// " bremsstrahlung photons. If On, complete
	// " modified Koch-Motz 2BS is used (modifications
	// " concern proper handling of kinematics at low energies,
	// " makes 2BS almost the same as 2BN at low energies).
	// " [ IBRDST ]
	// " Brems cross sections= BH, NIST, default is BH
	// " If BH is selected, the Bethe-Heitler bremsstrahlung
	// " cross sections (Coulomb corrected above 50 MeV)
	// " will be used. If NIST is selected, the NIST brems
	// " cross section data base (which is the basis for
	// " the ICRU radiative stopping powers) will be employed.
	// " Differences are negligible for E > ,say, 10 MeV,
	// " but signifficant in the keV energy range.
	// " Bound Compton scattering= On or Off
	// " If Off, Compton scattering will be treated with
	// " Klein-Nishina, with On Compton scattering is
	// " treated in the Impuls approximation. Default is Off???
	// " Make sure to turn on for low energy applications,
	// " not necessary above, say, 1 MeV.
	// " [ IBCMP ]
	// " Pair angular sampling= Off, Simple or KM
	// " If off, pairs are set in motion at an angle m/E
	// " relative to the photon direction (m is electron rest
	// " energy, E the photon energy). Simple turns on
	// " the leading term of the angular distribution
	// " (this is sufficient for most applications),
	// " KM (comes from Koch and Motz) turns on using 2BS
	// " from the article by Koch and Motz.
	// " Default is Simple, make sure you always use Simple or
	// " KM
	// " [ IPRDST ]
	// " Photoelectron angular sampling= Off or On
	// " If Off, photo-electrons get the direction of the
	// " `mother' photon, with On, Sauter's furmula is
	// " used (which is, striktly speaking, valid only for
	// " K-shell photo-absorption).
	// " If the user has a better approach, replace the macro
	// " $SELECT-PHOTOELECTRON-DIRECTION;
	// " The only application that
	// " I encountered until now where this option made a
	// " small difference was a big ion chamber (cavity size
	// " comparable with electron range) with high-Z walls
	// " in a low energy photon beam.
	// " Default is Off
	// " [ IPHTER ]
	// " Rayleigh scattering= Off, On
	// " If On, turnd on coherent (Rayleigh) scattering.
	// " Default is Off. Should be turned on for low energy
	// " applications
	// " [ IRAYLR ]
	// " Atomic relaxations= Off, On
	// " Default is Off. The effect of using On is twofold:
	// " - In photo-electric absorption events, the element
	// " (if material is mixture) and the shell the photon
	// " is interacting with are sampled from the appropriate
	// " cross seections
	// " - Shell vacancies created in photo-absorption events
	// " are relaxed via emission of fluorescent X-Rays,
	// " Auger and Koster-Cronig electrons.
	// " Make sure to turn this option on for low energy
	// " applications.
	// " [ IEDGFL ]
	// "
	// " Atomic relaxations, Rayleigh scattering,
	// " Photoelectron angular sampling and Bound Compton scattering
	// " can also be turned On/Off on a region-by-region
	// " basis. To do so, put e.g.
	// "
	// " Atomic relaxations= On in Regions or
	// " Atomic relaxations= Off in regions
	// "
	// " in your input file. Then use
	// "
	// " Bound Compton start region=
	// " Bound Compton stop region=
	// " or
	// " Rayleigh start region=
	// " Rayleigh stop region=
	// " or
	// " Relaxations start region=
	// " Relaxations stop region=
	// " or
	// " PE sampling start region=
	// " PE sampling stop region=
	// "
	// " each followed by a lost of of one or more
	// " start and stop regions separated by commas.
	// " Example:
	// " Atomic relaxations= On in Regions
	// " Relaxations start region= 1, 40
	// " Relaxations stop region= 10, 99
	// " will first turn off relaxations everywhere and
	// " then turn on in regions 1-10 and 40-99.
	// " Note that input is checked against min. and max.
	// " region number and ignored if
	// " start region < 1 or stop_region > $MXREG or
	// " start region > stop region.
	// "
	// " ECUT, PCUT and SMAX can also be set on a
	// " region-by-region basis. To do so, iclude
	// " in your input file
	// "
	// " Set XXXX= f_value1, f_value2, ...
	// " Set XXXX start region= i_value1, i_value2, ...
	// " Set XXXX stop region= j_value1, j_value2, ...
	// "
	// " where XXXX is ECUT, PCUT or SMAX ,
	// " f_value1, f_value2,... are the desired values for XXXX
	// " and i_value_i and j_value_i are the start and
	// " stop regions.
	// "
	// "*******************************************************************************
	// "
	// " VARIANCE REDUCTION
	// " ******************
	// "
	// " Delimeter: :start variance reduction:
	// " :stop variance reduction:
	// "
	// " BREM SPLITTING
	// " = Off (0) no bremsstrahlung splitting
	// " = On (1) there is bremsstrahlung spliting
	// "
	// " NUMBER OF BREMS PER EVENT
	// " (I) number of brems / event if splitting on
	// "
	// " CHARGED PARTICLE RUSSIAN ROULETTE
	// " = Off (0) Do not play Russian Roulette with charged particles
	// " = On (1) Play Russian Roulette with charged particles with
	// " probability of survival=PROB_RR=1/nbr_split.
	// " [I_PLAY_RR]
	// "
	// " ELECTRON RANGE REJECTION
	// " = off (0) No electron range rejection
	// " = on (1) Do electron range rejection. All charged
	// " particles without enough range to get out
	// " of their current region have their
	// " history terminated. This uses EGSnrc internal
	// " range rejection and takes no time to test.
	// " The parameter ESAVEIN also plays a role (see below)
	// " [IREJCT]
	// "
	// " We could/should? reinstitute the old approach by searching all
	// " regions outside the region of interest, finding the one with the
	// " greatest range and then use that range as a test against the distance
	// " to the region of interest. Some of the coding is left in place.
	// " This has been implemented in CAVRZnrc
	// "
	// " ESAVEIN (R) If ELECTRON RANGE REJECTION is on, discard an
	// " electron when E< ESAVEIN and RANGE < CDIST
	// " where CDIST is closest distance to region of
	// " interest specified below. This ignores brem
	// " losses below ESAVEIN.
	// " This parameter must be input even if not used.
	// " ESAVEIN is a total energy.
	// "
	// "
	// "========================================================================
	// " FOLLOWING IS NOT USED - leave in until decide what to do
	// "========================================================================
	// " RANGE REJECTION PARAMETERS
	// " (M5) CRANGE(1,1),CRANGE(2,1),ERANGE,CRANGE(1,2),CRANGE(2,2)
	// "
	// " Coefficients of a 2-piece logarithmic fit to the electron
	// " CSDA range used to discard electrons that cannot reach the
	// " range rejection region.
	// " The electron residual range is computed in HOWFAR as
	// " RANGE=EXP(CRANGE(1,I)+CRANGE(2,I)*ELKE)/RHO
	// " where i=1 for EKE>ERANGE, ELSE I=2.
	// " default values for crange,erange are for carbon
	// " with ERANGE=200keV FROM B&S '64. The range should always
	// " be overestimated with the closest fit being at the energy
	// " where most of the electrons are expected to be.
	// " NB There is only one medium used and hence it must be the
	// " one with the longest range in the problem.
	// "
	// " If ELECTRON RANGE REJECTION= on
	// "
	// " RANGE REJECTION MINIMUM PLANE (I) min. Z boundary for range rejection
	// "
	// " RANGE REJECTION MAXIMUM PLANE (I) max. Z boundary for range rejection
	// "
	// " RANGE REJECTION MINIMUM RADIUS (I) min radial boundary for range
	// rejection
	// "
	// " RANGE REJECTION MAXIMUM RADIUS (I) max radial boundary for range
	// rejection
	// "
	// " The above boundaries define the ``tracking region'' and if
	// " a charged particles range is such that it cannot get to
	// " this region, the history is terminated if E<ESAVEIN.
	// "========================================================================
	// "========================================================================
	// "
	// "
	// " RUSSIAN ROULETTE DEPTH (R)
	// " for russian roulette -
	// " as any photon crosses the Z='RUSSIAN ROULETTE DEPTH'
	// " plane, russian roulette is played.
	// " [RRZ]
	// "
	// " RUSSIAN ROULETTE FRACTION (R)
	// " Each time russian roulette is played, RRF IS THE
	// " probability of survival.
	// " weight increases by 1/RRF, if it survives
	// " [RRCUT]
	// "
	// " ****** IF BOTH ZERO, NO RUSSIAN ROULETTE IS PLAYED ******
	// "
	// " EXPONENTIAL TRANSFORM C (R)
	// " Parameter for pathlength biasing <0 FOR SHORTENING
	// " If 0.0, no biasing done.
	// " Review chapter discusses in detail.
	// " [CEXPTR]

	// "
	// " PHOTON FORCING
	// " = Off (0) normal photon transport (no forcing)
	// " = On (1) force photon interactions explicitly
	// " must set START and STOP FORCING in this case
	// " [IFORCE]
	// "
	// " START FORCING (I) number of photon interaction/history at which
	// " to start forcing photon interactions
	// " This input is required even if forcing off
	// " [NFMIN]
	// " STOP FORCING AFTER (I) number of photon interaction/history after which
	// " to stop forcing photon interactions
	// " [NFMAX]
	// " STOP FORCING AFTER must be >= START FORCING
	// " This input is required even if forcing off
	// "
	// " CS ENHANCEMENT FACTOR (R) can scale the photon cross section by this
	// " factor in a specified set of regions.
	// " From 1 to 10,000. with a default of 200 if
	// " it is on at all.
	// "
	// " CS ENHANCEMENT START REGION (M)
	// " CS ENHANCEMENT STOP REGION (M)
	// " Photon cross section scalled in these defined
	// " sets of regions.
	// " From 0 to NREG. Defaults to region 1 which means
	// " no enhancement since this is outside the geometry.
	// "
	// "

	// //"*******************************************************************************
	// "
	// " PLOT CONTROL INPUTS
	// " *******************
	// "*******************************************************************************
	// "
	// " PLOT CONTROL DELIMETERS: :start plot control:
	// " :stop plot control:
	// "
	// "
	// " PLOTTING
	// " = Off (0) no plots or plot files to be prepared
	// " = On (1) plotting to be prepared or printed
	// "
	// "
	// " ONLY IF PLOTTING= On
	// "
	// " LINE PRINTER OUTPUT
	// " = Off (0) don't plot in egs4lst file
	// " = On (1) do plot in egs4lst file
	// "
	// " ONLY IF PLOTTING= On
	// "
	// " EXTERNAL PLOTTER OUTPUT
	// " = Off (0) don't prepare plot files for xmgr
	// " = On (1) prepare xmgr input files
	// "
	// " ONLY IF EXTERNAL PLOTTER OUTPUT= On
	// " EXTERNAL PLOT TYPE
	// " = Point (1) point plot in xmgr file
	// " = Histogram (2) histogram plot in xmgr file
	// " = Both (3) both point plot and histogram
	// "
	// "
	// " ONLY IF PLOTTING= On
	// " PLOT RADIAL REGION IX (M) radial regions to plot vs depth
	// " (= 0 for no plots)
	// "
	// " ONLY IF PLOTTING= On
	// " PLOT PLANAR REGION IZ (M) planar slab numbers to plot vs radius
	// " (= 0 for no plots)
	// "
	// "******************************************************************************//"
	// " END OF INPUTS
	// " *************
	// "
	// "*******************************************************************************
	// "*******************************************************************************

	public static int $EBIN = 500;// "MAX NUMBER OF BINS IN PULSE HEIGHT DISTN     "
	// Note if you use large values of EBIN so that "
	// 2*EBIN > MXDATA (below), then use 2*EBIN "
	public static int $NSWTCH = 8;// "# OF NRC SWITCHES FOR CONTROLLING SCATTERING "

	public static int $NBATCH = 10;// "OUTPUT BATCHES                             "
	public static int JCASE = 0;; // "no. of histories per batch"
	public static int $NCASEMIN = 100;// "min. no. of histories                        "
	public static int $MXDATA;// =1040;//
								// "MAXIMUM DATA POINTS FOR ANALYSIS (i.e.($MXREG-1))"
	public static int $MAXIT = 7;// "MAX # OF PARAMETERS TO BE SCORED             "
	// "                                (1) TOTAL DOSE                                "
	// "                                (2) STOPPERS AND DISCARDS DOSE                "
	// "                                (3) TOTAL DOSE FROM FRONT WALL                "
	// "                                (4) TOTAL DOSE FROM SIDE WALL                 "
	// "                                (5) TOTAL DOSE FROM BACK WALL                 "
	// "                                (6) TOTAL DOSE FROM INSIDE WALL               "
	// "                                (7) TOTAL DOSE SCORED IN REGION NSRCRG        "
	// "                                    DUE TO PARTICLES CREATED IN THE           "
	// "                                    REGION (ONLY WHEN ISOURC = 3)             "

	public static int $MAXCMPTS = $MAXIT;// "FOR THE GRID OUTPUTS"
	public static int $MAX_SC_PLANES = 1;// "required to use phase space macros"
	public static int $MAXBRSPLIT = 200;// "MAX BREM SPLITTING NUMBER"
	// REPLACE {$PRESTA-II-INPUTS;} WITH { call prestaII_inputs; }
	public static int NCASE = 0;
	public static int NCASEO = 0;
	public static int NCASET = 0;
	public static double[][] AMASS;// ($MAXZREG,$MAXRADII),
	public static double TMCPUO = 0.0;
	public static double TIMMAX = 0.0;
	public static double STATLM = 0.0;// EIN,
	public static int IDAT = 0;
	public static int IRESTART = 0;
	public static int DATCOUNT = 0;
	// "AMASS(IZ,IX) MASS OF ZONE WITH COORDINATES (IZ,IX)
	// "TMCPUO CPU TIME USED IN PREVIOUS SESSIONS
	// "TIMMAX MAXIMUM ALLOWED CPU HOURS FOR A GIVEN CALCULATION
	// "STATLM TARGET STATISTICS IN CAVITY USED FOR AN EARLY EXIT
	// "EIN KINETIC ENERGY OF THE EXTERNAL BEAM
	// "ISUMCV(NREG) THE ARRAY OF ZONES COMPRISING THE CAVITY REGION
	// "IDAT = 0 STORE DATA ARRAYS FOR RE-USE
	// " = 1 DON'T STORE THEM
	// " = 0 DOES NOT INCLUDE PHOTOELECTRON ANGLE SELECTION
	// "NCASE NUMBER OF HISTORIES REMAINING TO BE DONE
	// "NCASEO NUMBER OF HISTORIES DONE IN PREVIOUS SESSIONS
	// "NCASET NUMBER OF HISTORIES ALREADY DONE
	// "IRESTART = 0 => INITIAL RUN
	// " = 1 => RESTARTED RUN
	// " = 3 => DATA ANALYSIS ONLY
	// " = 4 => READ STARTING RANDOM NUMBERS FROM A FILE
	// " = 5 => analyse previous parallel runs
	public static double RRZ = 0.0;//
	public static double RRCUT = 0.0;
	public static boolean RUSROU = false;
	// "RRZ COORDINATE OF PLANE AT WHICH RUSSIAN ROULETTE IS PLAYED
	// "RRCUT SURVIVAL PROBABILITY AFTER CROSSING THE PLANE
	// "RUSROU = .FALSE. => RUSSIAN ROULETTE WILL NOT BE PLAYED
	// " = .TRUE. => RUSSIAN ROULETTE WILL BE PLAYED

	// ACCUMULATES ENERGY DEPOSITED (ENERGY DEPOSITED^2);EVENTUALLY HOLDS DOSE
	// (UNCERTAINTY IN DOSE)
	public static double[][][] SCDOSE;
	public static double[][][] SCDOSE2;
	// ACCUMULATES ENERGY OF ELECTRONS CREATED (ENERGY OF ELECTRONS
	// " CREATED^2), EVENTUALLY HOLDS KERMA (UNCERTAINTY IN KERMA)
	public static double[][][] SCKERMA;
	public static double[][][] SCKERMA2;
	// ACCUMULATES (ENERGY DEPOSITED)*(ENERGY OF ELECTRONS CREATED)
	// " SUMMED OVER ALL PRIMARY HISTORIES, EVENTUALLY HOLDS UNCERTAINTY
	// " IN DOSE/KERMA
	public static double[][][] SCDOSEtoKERMA2;
	// ACCUMULATES COUNTS (COUNTS^2) IN EACH ENERGY BIN OF PULSE
	// " HEIGHT SENSITIVE REGION, EVENTUALLY HOLDS PULSE HEIGHT DISTN
	// " (UNCERTAINTY IN PULSE HEIGHT DISTN)
	public static double[] SCPDST = new double[$EBIN];
	public static double[] SCPDST2 = new double[$EBIN];
	// ACCUMULATES SUM OF COUNTS (COUNTS^2) IN ALL BINS <= CURRENT BIN
	// " IN PULSE HEIGHT REGION, EVENTUALLY HOLDS CUMULATIVE PULSE
	// " HEIGHT DISTN (UNCERTAINTY IN CUMULATIVE PULSE HEIGHT DISTN)
	public static double[] SCPCUM = new double[$EBIN];
	public static double[] SCPCUM2 = new double[$EBIN];
	// ACCUMULATES ALL COUNTS (COUNTS^2) IN ALL BINS IN PULSE HEIGHT
	// " REGION, EVENTUALLY HOLDS TOTAL OF ALL COUNTS (UNCERTAINTY ON
	// " TOTAL) AND IS USED TO NORMALIZE SCPDST, SCPCUM AND SCDFEP
	public static double SCPTOT = 0.0;
	public static double SCPTOT2 = 0.0;
	// ACCUMULATES TOTAL ENERGY (ENERGY^2) DEPOSITED IN PULSE HEIGHT
	// " REGION, EVENTUALL HOLDS ENERGY DEPOSITED IN PULSE HEIGHT REGION
	// " (UNCERTAINTY ON ENERGY IN PULSE HEIGHT REGION) PER INCIDENT
	// " PARTICLE
	public static double SCPHEN = 0.0;
	public static double SCPHEN2 = 0.0;
	// ACCUMULATES BACKGROUND COUNTS (COUNTS^2) FOR EACH PEAK IN PULSE
	// " HEIGHT DISTRIBUTION, WILL BE SUBTRACTED FROM SCDFEP BELOW
	public static double[] SCDFBK = new double[4];
	public static double[] SCDFBK2 = new double[4];
	// ACCUMULATES COUNTS (COUNTS^2) IN EACH ENERGY PEAK IN PULSE
	// " HEIGHT DISTRIBUTION
	public static double[] SCDFEP = new double[4];
	public static double[] SCDFEP2 = new double[4];
	// ACCUMULATES DIFFERENCE (DIFFERENCE^2) BETWEEN PEAK AND
	// " BACKGROUND COUNTS FOR EACH PEAK IN PULSE HEIGHT DISTRIBUTION,
	// " GETS USED TO DETERMINE COVARIANCE BETWEEN SCPTOT AND
	// " (SCDFEP-SCDFBK) SO WE CAN CALCULATE THE UNCERTAINTY ON THE
	// " FINAL VALUE OF SCDFEP
	public static double[] SCDFDIFF = new double[4];
	public static double[] SCDFDIFF2 = new double[4];
	// ACCUMULATES TOTAL PARTICLE STEPS (PARTICLE STEPS^2) IN PHANTOM,
	// " EVENTUALLY, SCSTP2 HOLDS UNCERTAINTY ON THIS NUMBER
	public static double SCSTP = 0.0;
	public static double SCSTP2 = 0.0;
	// ACCUMULATES TOTAL PARTICLE STEPS (PARTICLE STEPS^2) IN DOSE
	// " REGION OF PHANTOM, EVENTUALLY, SCDSTP2 HOLDS UNCERTAINTY ON
	// " THIS NUMBER
	public static double SCDSTP = 0.0;
	public static double SCDSTP2 = 0.0;
	public static double PIISTP = 0.0;// HOLDS NO. OF PRESTA-II STEPS FROM
										// PREVIOUS RUNS
	public static int SCSTP_LAST = 0;// LAST PRIMARY HISTORY TO SCORE A PARTICLE
										// STEP IN PHANTOM
	// LAST PRIMARY HISTORY TO SCORE A STEP IN THE DOSE REGION OF THE
	// " PHANTOM
	public static int SCDSTP_LAST = 0;
	// LAST PRIMARY HISTORY TO DEPOSIT ENERGY IN EACH DOSE REGION OF
	// " THE PHANTOM
	public static int[][][] SCDOSE_LAST;// ($MAXZREG,$MAXRADII,$MAXIT),
	// LAST PRIMARY HISTORY TO DEPOSIT KERMA IN EACH REGION OF
	// " PHANTOM
	public static int[][][] SCKERMA_LAST;// ($MAXZREG,$MAXRADII,$MAXIT),
	// SECOND LAST PRIMARY HISTORY TO DEPOSIT KERMA IN EACH REGION
	// " OF PHANTOM. USED TO CALCULATE SCDOSE*SCKERMA FOR COVARIANCE
	// " IN UNCERTAINTY ON DOSE/KERMA
	public static int[][][] SCKERMA_LASTOLD;// ($MAXZREG,$MAXRADII,$MAXIT),
	// LAST PRIMARY HISTORY TO SCORE ENERGY IN PULSE HEIGHT REGION
	public static int SCPDST_LAST = 0;
	// COUNTER FOR TOTAL NUMBER OF HISTORIES SUCCESSFULLY SIMULATED
	public static int IHSTRY = 0;
	// ENERGY LIMITS FOR FOUR PEAK AREA AND BACKGROUNDS
	public static double[][] DFEN = new double[4][4];// (4,4),
	// ENERGY DEPOSITED IN SENSITIVE VOLUME FOR CURRENT HISTORY ONLY
	public static double PHENER = 0.0;
	// stores value of WT(1) from last history for scoring phd
	public static double WT1OLD = 0.0;
	// TOPS OF ENERGY BINS FOR PULSE HEIGHT DISTRIBUTION
	public static double[] BINTOP = new double[$EBIN];// ($EBIN),
	// >0, WIDTH OF PULSE HEIGHT DISTRIBUTION ENERGY BINS
	// " ELSE FLAG TO USE BINTOP
	public static double SLOTE = 0.0;
	// WIDTH OF ENERGY BINS USED FOR PEAK EFFICIENCIES
	public static double DELTAE = 0.0;
	// ACCUMULATES ENERGY DEPOSITED IN EACH REGION BY CURRENT PRIMARY
	// " HISTORY
	public static double[][][] SCDOSE_TMP;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ACCUMULATES ENERGY OF ELECTRONS CREATED IN EACH REGION BY
	// " CURRENT PRIMARY HISTORY
	public static double[][][] SCKERMA_TMP;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ENERGY OF ELECTRONS CREATED IN EACH REGION BY LAST PRIMARY
	// " HISTORY. USED TO CALCULATE SCDOSE*SCKERMA FOR COVARIANCE
	// " IN UNCERTAINTY ON DOSE/KERMA
	public static double[][][] SCKERMA_TMPOLD;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ACCUMULATES CHARGED PARTICLE STEPS IN CURRENT PRIMARY HISTORY
	public static double SCSTP_TMP = 0.0;
	// ACCUMULATES CHARGED PARTICLE STEPS IN DOSE REGION IN CURRENT
	// " PRIMARY HISTORY
	public static double SCDSTP_TMP = 0.0;
	// MAXIMUM LEVEL TO WHICH THE STACK OF DAUGHTER PARTICLES FROM AN
	// " INCIDENT PARTICLE RISES (STACK MAY INCLUDE INCIDENT PARTICLE)
	public static int MXNP = 0;
	// "IFULL = 0 JUST CALCULATE TOTAL DOSE AND THAT DUE TO STOPPERS
	// " AND DISCARDS (THE DEFAULT)
	// " = 1 ABOVE ANALYSE WHERE THE DOSE IS COMING FROM
	// " = 2 IFULL = 0 SCORING PLUS PULSE HEIGHT DISTRIBUTIONS
	// " = 3 SCORE THE SCATTER FRACTION INSTEAD OF STOPPERS

	// "ISTORE = 0 DO NOT STORE THE INITIAL RANDOM NUMBERS (THE DEFAULT)
	// " = 1 STORE THE INITIAL RANDOM NUMBER FOR LAST HISTORY
	// " = 2 STORE INITIAL RANDOM NUMBERS FOR ALL HISTORIES
	public static int ISTORE = 0;
	// "IKERMA = 0 DO NOT SCORE KERMA
	// " = 1 SCORE KERMA
	public static int IKERMA = 0;
	// "IWATCH = 0 FOR NORMAL OUTPUT (THE DEFAULT)
	// " = 1 OUTPUT ON EVERY DISCRETE INTERACTION
	// " = 2 OUTPUT ON EVERY ELECTRON/PHOTON STEP AS WELL
	// " = 3 PRINTS OUT ONLY WHEN ENERGY IS DEPOSITED
	// " = 4 PRINTS OUT FILE FOR GRAPHICS
	public static int IWATCH = 0;
	// "IOOPTN = 0 SHORT OUTPUT (THE DEFAULT) -JUST DOSE GRID(DG)
	// " = 1 OUTPUT DOSE SUMMARY ONLY (DS)
	// " = 2 OUTPUT MATERIAL SUMMARY GRID(MG) + DG
	// " = 3 OUTPUT MG + DS
	// " = 4 OUTPUT MG + DS + DG
	public static int IOOPTN = 0;
	// "IOUTSP = 0 NO ENERGY INPUT SPECTRUM DATA IN OUTPUT SUMMARY
	// " = 1 INCLUDE ENERGY INOUT SPECTRUM DATA IN OUTPUT SUMMARY
	public static int IOUTSP = 0;
	// FLAG ARRAY FOR EACH REGION, NON-ZERO ONLY IF PULSE HEIGHT
	// " DISTRIBUTION WANTED IN THIS GEOMETRIC REGION
	public static int[] IPHR;// ($MXREG),
	// NUMBER OF ENERGY BINS FOR PULSE HEIGHT DISTRIBUTION
	public static int MAXBIN = 0;
	public static int NCOMPT = 0;
	// public static int IFANO=0;// = 0 NO PHOTON REGENERATION
	// " = 1 PHOTONS REGENERATED AFTER THEY HAVE INTERACTED
	// " = 2 NO PHOTON REGENERATION, ELECTRONS FROM CAVITY WALL ARE ELIMINATED

	// "CDFINV INVERSE OF THE CUMULATIVE ENERGY PROBABILITY DISTRIBUTION
	// FUNCTION
	public static int $INVDIM = 1000;
	public static double[][] CDFINV = new double[$INVDIM][2];
	public static final String NULLs = "";
	public static final String BLANK = " ";
	public static final String DCHAR = "D";
	public static final String TCHAR = "T";
	public static final String ACHAR = "A";
	// -------------INPUTS-------------------------------------------
	public static String TITLEs = "";
	// @ IWATCH
	public static final int IWATCH_OFF = 0;
	public static final int IWATCH_INTERACTIONS = 1;
	public static final int IWATCH_STEPS = 2;
	public static final int IWATCH_DEPOSITED = 3;
	public static final int IWATCH_GRAPH = 4;
	// @ STORE INITIAL RANDOM NUMBERS
	public static final int ISTORE_NO = 0;
	public static final int ISTORE_LAST = 1;
	public static final int ISTORE_ALL_DEPOSITED = 2;
	public static final int ISTORE_ALL = 3;
	// @ IRESTART
	public static final int IRESTART_FIRST = 0;
	public static final int IRESTART_RESTART = 1;
	public static final int IRESTART_ANALYZE = 2;
	public static final int IRESTART_START_RNS = 3;
	public static final int IRESTART_PARALLEL = 4;
	// @ OUTPUT OPTIONS
	public static final int IOOPTN_SHORT = 0;
	public static final int IOOPTN_DOSE_SUMMARY = 1;
	public static final int IOOPTN_MATERIAL_SUMMARY = 2;
	public static final int IOOPTN_MATERIAL_AND_DOSE_SUMMARY = 3;
	public static final int IOOPTN_LONG = 4;
	// @ ELECTRON TRANSPORT
	public static final int ETRANS_NORMAL = 0;
	public static final int ETRANS_NO_INTERACTION = 1;
	// @ STORE DATA ARRAYS
	public static final int IDAT_YES = 0;
	public static final int IDAT_NO = 1;
	// @ NUMBER OF HISTORIES
	public static final int NCASE_MIN = 1;
	public static final int NCASE_MAX = 461168600;// 4.611686e18;
	public static final int NCASE_DEFAULT = 20000;
	// @ MAX CPU HOURS ALLOWED
	public static final double TIMMAX_MIN = 0.0;
	public static final double TIMMAX_MAX = 1000.0;
	public static final double TIMMAX_DEFAULT = 999.0;
	// @ IFULL
	public static final int IFULL_DOSE_AND_STOPPERS = 0;
	// public static final int IFULL_AATT_AND_ASCAT=1;
	public static final int IFULL_ENTRANCE_REGIONS = 1;
	public static final int IFULL_PULSE_HEIGHT_DISTRIBUTION = 2;
	public static final int IFULL_SCATTER_FRACTION = 3;
	public static final int IFULL_OFMET_Fricke = 4;
	// @ STATISTICAL ACCURACY SOUGHT
	public static final double STATLM_MIN = 0.0;
	public static final double STATLM_MAX = 100.0;
	public static final double STATLM_DEFAULT = 0.0;
	// @ PHOTON REGENERATION
	public static final int IFANO_NO = 0;
	public static final int IFANO_YES = 1;
	public static final int IFANO_NO_ELECTRONS_FROM_WALL = 2;
	// @ SCORE KERMA
	public static final int KERMA_YES = 1;
	public static final int KERMA_NO = 0;
	// @ INITIAL RANDOM NO. SEEDS
	public static final int RANLUX_LEVEL_MIN = 0;
	public static final int RANLUX_LEVEL_MAX = 4;
	public static final int RANLUX_LEVEL_DEFAULT = EGS4.$DEFAULT_LL;// see egs4
	public static final int RANLUX_SEED_MIN = 1;
	public static final int RANLUX_SEED_MAX = 1073741824;// default seed is set
															// in
															// egs4.setdefaults!!
	public static final int RANLUX_SEED_DEFAULT = 999999;
	public static final int RANMAR_SEED_MIN = 1;
	public static final int RANMAR_SEED_MAX = 30081;
	public static final int RANMAR_SEED_DEFAULT = 9373;
	public static int jrng1 = 0;
	public static int jrng2 = 0;
	// @TRANSPORT PARAMS:
	public static double ecut = 0.0;
	public static double pcut = 0.0;
	public static double smax = 0.0;
	public static boolean setEcutRegion = false;
	public static boolean setPcutRegion = false;
	public static boolean setSmaxRegion = false;
	public static final double ECUT_MIN = 0.0;// " ECUT "
	public static final double ECUT_MAX = 1.E15;
	public static final double PCUT_MIN = 0.0;// " PCUT "
	public static final double PCUT_MAX = 1.E15;
	public static final double SMAXIR_MIN = 0.0;// " SMAX "
	public static final double SMAXIR_MAX = 1.E15;
	public static int nEcut = 1;// number of data
	public static double[] Ecut;// =new double[EGS4.$MXREG];
	public static int[] startEcutRegion;
	public static int[] stopEcutRegion;
	public static int nPcut = 1;// number of data
	public static double[] Pcut;// =new double[EGS4.$MXREG];
	public static int[] startPcutRegion;
	public static int[] stopPcutRegion;
	public static int nSmax = 1;// number of data
	public static double[] Smax;// =new double[EGS4.$MXREG];
	public static int[] startSmaxRegion;
	public static int[] stopSmaxRegion;

	public static boolean setIncohRegion = false;
	public static boolean setCohRegion = false;
	public static boolean setRelaxRegion = false;
	public static boolean setPeRegion = false;
	public static int nIncoh = 1;// number of data
	public static int[] Incoh;// =new int[EGS4.$MXREG];
	public static int[] startIncohRegion;
	public static int[] stopIncohRegion;
	public static int nCoh = 1;// number of data
	public static int[] Coh;// =new int[EGS4.$MXREG];
	public static int[] startCohRegion;
	public static int[] stopCohRegion;
	public static int nRelax = 1;// number of data
	public static int[] Relax;// =new int[EGS4.$MXREG];
	public static int[] startRelaxRegion;
	public static int[] stopRelaxRegion;
	public static int nPe = 1;// number of data
	public static int[] Pe;// =new int[EGS4.$MXREG];
	public static int[] startPeRegion;
	public static int[] stopPeRegion;

	// " Incoherent (Compton) scattering "
	public static int incoh = 0;
	public static final int incoh_OFF = 0;
	public static final int incoh_ON = 1;
	public static final int incoh_ON_IN_REGIONS = 2;
	public static final int incoh_OFF_IN_REGIONS = 3;
	// " Radiative corrections for Compton scattering "
	public static final int radc_OFF = 0;
	public static final int radc_ON = 1;
	// " Coherent (Rayleigh) scattering "
	public static int coh = 0;
	public static final int coh_OFF = 0;
	public static final int coh_ON = 1;
	public static final int coh_ON_IN_REGIONS = 2;
	public static final int coh_OFF_IN_REGIONS = 3;
	// " Atomic Relaxations "
	public static int relax = 0;
	public static final int relax_OFF = 0;
	public static final int relax_ON = 1;
	public static final int relax_ON_IN_REGIONS = 2;
	public static final int relax_OFF_IN_REGIONS = 3;
	// " Photoelectron angular sampling "
	public static int pe = 0;
	public static final int pe_ang_OFF = 0;
	public static final int pe_ang_ON = 1;
	public static final int pe_ang_ON_IN_REGIONS = 2;
	public static final int pe_ang_OFF_IN_REGIONS = 3;
	// " Bremsstrahlung angular sampling "
	public static final int brems_ang_SIMPLE = 0;
	public static final int brems_ang_KM = 1;
	// " Bremsstrahlung cross sections "
	public static final int brems_cross_BH = 0;
	public static final int brems_cross_NIST = 1;
	// " Pair angular sampling "
	public static final int pair_ang_OFF = 0;
	public static final int pair_ang_SIMPLE = 1;
	public static final int pair_ang_KM = 2;
	public static final int pair_ang_UNIFORM = 3;
	public static final int pair_ang_BLEND = 4;
	// " Pair cross sections "
	public static final int pair_cross_BH = 0;
	public static final int pair_cross_NRC = 1;
	// " Triplet production "
	public static final int triplet_OFF = 0;
	public static final int triplet_ON = 1;
	// " Spin effects  		"
	public static int ispin = 0;
	public static final int spin_OFF = 0;
	public static final int spin_ON = 1;
	// " Electron impact ionization "
	public static final int eii_OFF = 0;
	public static final int eii_ON = 1;

	public static final double ESTEPE_MIN = 1.0E-5;// " ESTEPE "
	public static final double ESTEPE_MAX = 1.0;
	public static final double XIMAX_MIN = 0.0;// " XIMAX "
	public static final double XIMAX_MAX = 1.0;
	// " BCA "
	public static final int BCA_EXACT = 0;
	public static final int BCA_PRESTA_I = 1;

	public static final double Skindepth_MIN = -1.0;// " Skindepth "
	public static final double Skindepth_MAX = 1.0e15;
	// " Electron-step algorithm "
	public static final int estep_alg_PRESTA_I = 1;
	public static final int estep_alg_PRESTA_II = 0;

	// Variance reduction
	public static final int irejct_OFF = 0;
	public static final int irejct_ON = 1;
	public static double ESAVEIN = 0.0;
	public static final double ESAVEIN_MIN = 0.0;
	public static final double ESAVEIN_MAX = 1.0e30;
	public static final double ESAVEIN_DEFAULT = 0.0;
	public static final double cs_enhance_MIN = 0.0;
	public static final double cs_enhance_MAX = 1.0e6;
	public static final double cs_enhance_DEFAULT = 0.5;
	public static final double RRDEPTH_MIN = -1.0e30;
	public static final double RRDEPTH_MAX = 1.0e30;
	public static final double RRDEPTH_DEFAULT = 0.0;
	public static final double RRFRACTION_MIN = -1.0e30;
	public static final double RRFRACTION_MAX = 1.0e30;
	public static final double RRFRACTION_DEFAULT = 0.0;
	public static final double EXPC_MIN = -1.0e30;
	public static final double EXPC_MAX = 1.0e30;
	public static final double EXPC_DEFAULT = 0.0;
	public static int IFARCE = 0;
	public static final int IFARCE_ON = 1;
	public static final int IFARCE_OFF = 0;
	public static final int NFMIN_MIN = 0;
	public static final int NFMIN_DEFAULT = 1;
	public static int NFMIN_MAX = 0;
	public static final int NFMAX_MIN = 0;
	public static final int NFMAX_DEFAULT = 1;
	public static int NFMAX_MAX = 0;
	public static int phsplitt = 1;
	public static final int phsplitt_MIN = 1;
	public static final int phsplitt_DEFAULT = 1;
	public static int phsplitt_MAX = 0;

	// ###########################################
	public static double EI = 0.0;
	public static double EKMAX = 0.0;
	public static double DEPTH = 0.0;
	public static double VOLUME = 0.0;
	public static double RLOW2 = 0.0;
	public static int ISUMX = 0;

	public static double FMASSC = 0.0;
	public static double FMASS = 0.0;
	public static double AINFLU_CURRENT = 0.0; // -------
	public static int last_case = 0;
	public static int n_photons = 0;
	public static int n_electrons = 0;
	public static double sumE_photons = 0.0;
	public static double sumE2_photons = 0.0;
	public static double sumE_electrons = 0.0;
	public static double sumE2_electrons = 0.0;
	public static int IBTCH = 0;
	public static double dcav_current = 0.0;
	public static double dcavun = 0.0;
	public static int IDECAV = 0;
	public static int ipass = 0;// not usded
	public static double SCORE_NORM_NUM = 0.0;
	public static double SCORE_TEMP = 0.0;

	public static int IRL = 0;
	public static int MEDNUM = 0;
	// ================================dose
	public static int ITMAX = 0;
	// "DECISION = 1 if ustep = 0 and region number changes--this ensures that
	// " LATCH bit gets set to reflect NEWNRC before dose is deposited
	public static int DECISION = 0;
	public static int NDOSE = 0;
	public static int NZDOSE = 0;
	public static int NRDOSE = 0;
	public static int MINZ = 0;
	public static int MAXZ = 0;
	public static int MINR = 0;
	public static int MAXR = 0;
	public static int IZD = 0;
	public static int IXD = 0;
	public static int NZDMIN = 0;
	public static int NZDMAX = 0;
	public static int NRDMIN = 0;
	public static int NRDMAX = 0;
	// "IDSTBL(IRL,1)=DOSE PLANAR SLAB NUMBER IF IRL IS IN DOSE REGION
	// "IDSTBL(IRL,2)=DOSE CYLINDRICAL REGION NUMBER IF IRL IS IN DOSE REGION
	public static int[][] IDSTBL;// ($MXREG,2);
	public static double EK0 = 0.0;
	public static double TDSMAX = 0.0;
	public static int IDSMAX = 0;
	public static double TDOS = 0.0;
	public static double TDOS2 = 0.0;
	public static int IB = 0;
	public static double SCORE_TEMP2 = 0.0;
	public static int IGBUG1 = 0;
	public static int IGBUG2 = 0;

	public static int nREGSOLV = 0;
	public static int[] REGSVOL = new int[1000];// =0;
	public static int nTOPEBIN = 0;
	public static double[] TOPEBIN = new double[1000];
	public static final double SLOTE_DEFAULT = 1.25;
	public static final double SLOTE_MIN = -10000.0;
	public static final double SLOTE_MAX = 10000.0;
	public static final double DELTAE_DEFAULT = 0.005;
	public static final double DELTAE_MIN = 1.e-20;
	public static final double DELTAE_MAX = 100000.0;
	public static final double TOPEBIN_DEFAULT = 1.25;
	public static final double TOPEBIN_MIN = 1.e-20;
	public static final double TOPEBIN_MAX = 100000.0;
	public static int nENHREG = 0;
	public static int[] NENHLO = new int[1000];
	public static int[] NENHHI = new int[1000];
	public static int ics_enhance = 0;
	// ========================
	public static boolean createOutputFile = false;
	public static boolean putInFile = false;// internal var defining when and
											// what to print
	private String filename = "";
	FileWriter sigfos;
	// -------------end INPUTS-------------------------------------------
	// REPLACE {$INVDIM} WITH {1000} "DIMENSION CONTROLS GRID SIZE FOR INVERSE"
	// CDFINV($INVDIM,2)

	public static boolean is_finished = false;

	/**
	 * Constructor.
	 */
	public doserz() {
		// createOutputFile=false;
		createOutputFile = true;
		putInFile = true;
		Calendar cal = Calendar.getInstance();
		String fs = cal.get(Calendar.YEAR) + "_" + cal.get(Calendar.MONTH) + "_"
				+ cal.get(Calendar.DAY_OF_MONTH) + "_" + cal.get(Calendar.HOUR) + "_"
				+ cal.get(Calendar.MINUTE) + "_" + cal.get(Calendar.SECOND) + ".txt";
		filename = fs;// will be 2005_10_25_14_30_56.txt

		init();
	}

	/**
	 * Perform basic initialization and RUN the Monte Carlo engine.
	 */
	private void init() {
		if (createOutputFile) {
			try {
				sigfos = new FileWriter(filename);
			} catch (Exception ex) {
			}
		}

		EGS4Geom.$MAXZREG = 61;// "MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4Geom.$MAXRADII = 60;// "MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4.setMXMED(10);// "MAX # OF MEDIA 		"
		EGS4.setMXREG(EGS4Geom.$MAXRADII * EGS4Geom.$MAXZREG + 1);// "$MAXRADII*$MAXZREG)+1   "
		EGS4.setMXSTACK(4000);// "NEED HUGE STACK FOR CORRELATIONS+splitting   "
		EGS4.setMXRANGE(500); // "for range arrays used in range_rejection()"
		// flush
		EGS4.reset();
		EGS4.startSimulationTime = System.currentTimeMillis();
		// --variable init
		$MXDATA = EGS4Geom.$MAXRADII * EGS4Geom.$MAXZREG;// COMPUTE $MXREG-1
		EGS4SrcEns.$MXRDIST = 1000;
		EGS4SrcEns.$NENSRC = 300;
		EGS4Geom.$NVALUE = 100;
		// ->EGS4Macro.
		// "NEWNRC = 0 IF THE PARTICLE IS ABOUT TO LEAVE THE GEOMETRY
		// "NEWNRC = 10,20,30 OR 40 DEPENDING ON WHETHER THE PARTICLE IN
		// QUESTION IS
		// " ABOUT TO ENTER A NEW REGION VIA THE FRONT WALL, OUTER WALL, BACK
		// " WALL OR INNER WALL RESPECTIVELY
		// "NEWNRC = 50 IF THE PARTICLE IS IN THE SOURCE REGION AND NEVER LEFT

		AMASS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];// ($MAXZREG,$MAXRADII),
		Ecut = new double[EGS4.$MXREG];
		startEcutRegion = new int[EGS4.$MXREG];
		stopEcutRegion = new int[EGS4.$MXREG];
		Pcut = new double[EGS4.$MXREG];
		startPcutRegion = new int[EGS4.$MXREG];
		stopPcutRegion = new int[EGS4.$MXREG];
		Smax = new double[EGS4.$MXREG];
		startSmaxRegion = new int[EGS4.$MXREG];
		stopSmaxRegion = new int[EGS4.$MXREG];
		Incoh = new int[EGS4.$MXREG];
		startIncohRegion = new int[EGS4.$MXREG];
		stopIncohRegion = new int[EGS4.$MXREG];
		Coh = new int[EGS4.$MXREG];
		startCohRegion = new int[EGS4.$MXREG];
		stopCohRegion = new int[EGS4.$MXREG];
		Relax = new int[EGS4.$MXREG];
		startRelaxRegion = new int[EGS4.$MXREG];
		stopRelaxRegion = new int[EGS4.$MXREG];
		Pe = new int[EGS4.$MXREG];
		startPeRegion = new int[EGS4.$MXREG];
		stopPeRegion = new int[EGS4.$MXREG];
		// ---------LOCAL VAR---------------
		SCDOSE = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCKERMA = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCKERMA2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSEtoKERMA2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE_LAST = new int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];// ($MAXZREG,$MAXRADII,$MAXIT),
		SCKERMA_LAST = new int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCKERMA_LASTOLD = new int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE_TMP = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCKERMA_TMP = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCKERMA_TMPOLD = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		// ($MAXZREG,$MAXRADII,$MAXIT)
		IPHR = new int[EGS4.$MXREG];
		IDSTBL = new int[EGS4.$MXREG][2];// ($MXREG,2);
		// @@@@@@@@@@@@@@@@@@@@
		EGS4Grid.CDSTBL = new String[EGS4.$MXREG];
		EGS4Grid.CTRTBL = new String[EGS4.$MXREG];
		EGS4Grid.CABSRB = new String[EGS4.$MXREG];
		// --interface->LINK
		EGS4.eq = this;// pass the printing mode
		EGS4Core.eq = this;// pass the printing mode
		EGS4Macro.eq = this;// pass the printing mode
		EGS4SrcEns.eq = this;// pass the printing mode
		EGS4Geom.eq = this;// pass the printing mode
		EGS4Grid.eq = this;// pass the printing mode
		// --Macro and defaults param:
		EGS4.egs_set_defaults();// first default then:
		EGS4.RandomUse = 1;// use EGS ranlux or ranmar generators
		EGS4.ranluxB = false;// use ranmar!!
		EGS4.ispmfp = EGS4.iDose;// select photon mfp
		EGS4.iurd = 0;// no user range discard
		EGS4.iraycorr = EGS4.iDose;// rayleigh correction
		EGS4.isemfp = EGS4.iDose;// $SELECT-ELECTRON-MFP
		EGS4.iGeom = 1;// EGS4.iCavity;//1=RZ GEOM
		EGS4Macro.ismfpfpb = EGS4.iDose;// select mfp parallel beam
		// EGS4Macro.irange_rej=EGS4.iCavity;//not used range rejection

		EGS4Grid.$MAXCMPTS = $MAXCMPTS;
		// -----------------
		EGS4.USER_CONTROLS_TSTEP_RECURSION = EGS4.iDose;// no effect here
		EGS4.hatchindex = EGS4.iDose;// no effect here
		// -----------------
		EGS4.iprint = 2;// summary
		// inputs
		// is_finished = true;
		// "******************************************************************************
		// "
		// " *** SECTION 1 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "READ INPUTS AND CALCULATE ONE-TIME ONLY CONSTANTS
		// "
		// "------------------------------------------------------------------------------

		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "=================================================================================";
		EGS4.seqStr = " ********************DOSE: DEMO APPLICATION************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " It also scores pulse height distributions";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " in an arbitrary volume made up of any number of regions.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " 					" + d.toString();
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " ********************************************************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		inputs();
		if (EGS4.STOPPROGRAM) {
			return;
		}

		// "OPEN A FILE FOR STORING OR READING RANDOM NUMBERS"
		if (ISTORE > 0)// NOT ALLOWED
		{// "We want to store the rng state in a file"
			// rng_unit = egs_open_file(2,0,1,'.egsrns');
		}
		/*
		 * ELSE IF( irestart = 4 )//NOT ALLOWED [ rng_unit =
		 * egs_open_datfile(2,0,1,'.egsrns'); ]
		 */
		// " This is copied directly from below.  Don't know why it needs "
		// " repeating, but it does.  -- JT "
		if ((EGS4Macro.IFULL == 0) || (EGS4Macro.IFULL == 2)
				|| (EGS4Macro.IFULL == 3)) {
			ITMAX = 2;
		} else {
			ITMAX = $MAXIT;
		}
		/*
		 * IF(IRESTART.EQ.5)//NOT ALLOWED [
		 * 
		 * call egs_combine_runs(combine_results,'.egsdat');
		 * 
		 * NBATCH=0; "DON'T WANT IT TO RUN ANY HISTORIES" NCASET=NCASEO;
		 * "To prevent a wrong normalization if some of the "
		 * "parallel runs not available, IK, Jan 21 1999" ]
		 * "end of IRESTART = 5, DISTRIBUTED POST-PROCESSING" ELSE [
		 */
		if (NCASE / $NBATCH == 0) {
			NCASE = $NBATCH;
		}
		JCASE = NCASE / $NBATCH;
		NCASE = JCASE * $NBATCH;// "NUMBER OF HISTORIES PER BATCH
		// ]
		DECISION = 0; // "this is for howfar and ausgab when IFULL=1"
		MXNP = 0; // "reset the maximum stack indicator"
		IHSTRY = NCASEO; // "reset the number of histories counter"

		EGS4SrcEns.NHSTRY = 0; // "start the no. of primary histories counter at zero"

		// "set up the broad parallel beam defaults"
		if (EGS4SrcEns.ISOURC == 2) {
			EGS4Geom.NR = 1;
			EGS4Geom.RCYL[1] = 1000.;
			EGS4Geom.nreg = EGS4Geom.NZ + 1;
			EGS4Geom.CYRAD2[0] = EGS4Geom.RCYL[1] * EGS4Geom.RCYL[1];
		}

		NDOSE = NZDOSE * NRDOSE; // "NUMBER OF DOSE SCORING REGIONS"
		if (EGS4Macro.IFULL == 1) {
			ITMAX = $MAXIT;
		} else {
			ITMAX = 2;
		}// "# OF DOSE COMPONENTS"

		if (EGS4Macro.irejct == 1) {
			// "GET COORDINATES USED BY HOWFAR FOR RANGE REJECTION"
			EGS4Macro.ZMINR = EGS4Geom.ZPLANE[MINZ - 1];// "LESSER PLANE POSITION"
			EGS4Macro.ZMAXR = EGS4Geom.ZPLANE[MAXZ - 1];// "GREATER PLANE POSITION"
			if (MINR == 0) {
				EGS4Macro.RMINR = 0.0;
			} else {
				EGS4Macro.RMINR = EGS4Geom.RCYL[MINR];
			}// "INNER CYLINDER RADIUS"
			EGS4Macro.RMAXR = EGS4Geom.RCYL[MAXR]; // "OUTER CYLINDER RADIUS"
		}

		// "SET UP TABLES CORRESPONDING TO DOSE AND RANGE REJECTION TRACKING REGIONS"
		// "FOR GEOMETRICAL REGION NUMBER 'IRL' A DOSE SCORING REGION THEN IDSTBL(IRL,1)"
		// "ASSIGNED DOSE COORDINATE IZD AND IDSTBL(IRL,2) ASSIGNED DOSE COORDINATE IXD"
		// "AND CDSTBL(IRL) ASSIGNED 'D'. IF REGION IRL WITHIN TRACKING REGION FOR"
		// "RANGE REJECTION THEN CTRTBL(IRL) ASSIGNED 'T'"
		// "ALSO SET UP TABLE NTRACK WITH ENTRY OF 1 FOR DOSE SCORING ZONE ELSE 0. THIS"
		// "IS REDUNDANT IN VIEW OF CDSTBL AND IS ONLY USED IN MACRO 'COUNT-NOSCAT-IN"
		// "-CAVITY' BUT IS LESS DANGEROUS THAN INTRODUCING A NEW ARRAY TO EGS"
		EGS4Grid.CTRTBL[0] = BLANK;
		EGS4Grid.CDSTBL[0] = BLANK;
		for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
			for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
				// $GET-IRL(IZ,IX); //"DERIVE CORRESPONDING DOSE ZONE NUMBERS"
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				IZD = IZ + 1 - NZDMIN;
				IXD = IX - NRDMIN;
				if ((IZD <= 0) || (IZD > NZDOSE) || (IXD <= 0)
						|| (IXD > NRDOSE)) {
					EGS4Grid.CDSTBL[IRL - 1] = BLANK;
					IDSTBL[IRL - 1][0] = 0;// (IRL,1)=0;
					IDSTBL[IRL - 1][1] = 0;// (IRL,2)=0;
					EGS4Geom.ntrack[IRL - 1] = 0;
				} else {
					EGS4Grid.CDSTBL[IRL - 1] = DCHAR;
					IDSTBL[IRL - 1][0] = IZD;
					IDSTBL[IRL - 1][1] = IXD;
					EGS4Geom.ntrack[IRL - 1] = 1;
				}
				if ((IZ < MINZ) || (IZ >= MAXZ) || (IX <= MINR) || (IX > MAXR)) {
					EGS4Grid.CTRTBL[IRL - 1] = BLANK;
				} else {
					EGS4Grid.CTRTBL[IRL - 1] = TCHAR;
				}
			}
		}

		// "set up ausgab calls"
		for (int J = 1; J <= 5; J++) {
			EGS4.iausfl[J - 1] = 1;// IAUSFL(J)=1;
		}
		for (int J = 6; J <= 25; J++) {
			EGS4.iausfl[J - 1] = 0;
		} // "NORMAL EXECUTION"

		if (EGS4Macro.IFULL == 1 || EGS4Macro.IFULL == 3 || IKERMA == 1) {
			// "need to call ausgab to set flag after photon interactions"
			// "IKERMA=1 means scoring KERMA"
			// /IAUSFL(6),IAUSFL(17),IAUSFL(19),IAUSFL(21)/=1;
			// "17 => call after pair, 19=>call after compton, 21=>call after photoelectric"
			// "for KERMA, rayleigh scatter has no effect"
			EGS4.iausfl[5] = 1; // "AFTER TRANSPORT"
			EGS4.iausfl[16] = 1; // "after pair
			EGS4.iausfl[18] = 1; // "After COMPTON"
			EGS4.iausfl[20] = 1; // "After Photo"
		}
		if (EGS4Macro.IFULL == 4) {
			EGS4.iausfl[7] = 1;// "After bremsstrahlung"
			// "iausfl(14) = 1;" "A positron has annihilated in-flight"
			// "iausfl(15) = 1;" "A positron has annihilated at rest"
		}

		/*
		 * if(IFANO == 1) {
		 * //"AUSGAB will be responsible for making sure that the beam is not"
		 * //"attenuated and getting rid of the scattered photons."
		 * EGS4.iausfl[15] = 1; //"Before pair" EGS4.iausfl[17] = 1;
		 * //"Before Compton" EGS4.iausfl[18] = 1; //"After Compton"
		 * EGS4.iausfl[19] = 1; //"Before photoelectric" EGS4.iausfl[20] = 1;
		 * //"After photoelectric" EGS4.iausfl[23] = 1; //"Before Rayleigh"
		 * EGS4.iausfl[24] = 1; //"After Rayleigh"
		 * 
		 * //"AUSGAB will be responsible for throwing away any photons resulting"
		 * //"from a primary electron. ie. True equilibtrium requires that all"
		 * //"energy deposition be local to the primary interaction site."
		 * EGS4.iausfl[7] = 1; //"After bremsstrahlung" EGS4.iausfl[13] = 1;
		 * //"A positron has annihilated in-flight" EGS4.iausfl[14] = 1;
		 * //"A positron has annihilated at rest" }
		 * 
		 * if( EGS4Macro.n_split > 1 ) { EGS4.iausfl[7] = 1;
		 * //"After bremsstrahlung" EGS4.iausfl[13] = 1;
		 * //"A positron has annihilated in-flight" EGS4.iausfl[14] = 1;
		 * //"A positron has annihilated at rest"
		 * 
		 * //"With n_split > 1, we don't need the following calls even if ifano = 1"
		 * EGS4.iausfl[15] = 0; EGS4.iausfl[17] = 0; EGS4.iausfl[19] = 0;
		 * EGS4.iausfl[23] = 0; EGS4.iausfl[24] = 0; }
		 * 
		 * if(IFANO == 2) {
		 * 
		 * //"AUSGAB will be responsible for discarding energy due to electrons set"
		 * //"in motion in the wall" EGS4.iausfl[16] = 1; //"After pair"
		 * EGS4.iausfl[18] = 1; //"After Compton" EGS4.iausfl[20] = 1;
		 * //"After photoelectric" }
		 * 
		 * if( EGS4Macro.use_enhance ) { EGS4.iausfl[15] = 1; //"Before pair"
		 * EGS4.iausfl[17] = 1; //"Before Compton" EGS4.iausfl[18] = 1;
		 * //"After Compton" EGS4.iausfl[19] = 1; //"Before photoelectric"
		 * EGS4.iausfl[20] = 1; //"After photoelectric" EGS4.iausfl[23] = 1;
		 * //"Before Rayleigh" EGS4.iausfl[24] = 1; //"After Rayleigh"
		 * EGS4.iausfl[7] = 1; //"After bremsstrahlung" EGS4.iausfl[13] = 1;
		 * //"A positron has annihilated in-flight" EGS4.iausfl[14] = 1;
		 * //"A positron has annihilated at rest" }
		 */
		if (EGS4Macro.cs_enhance > 1.0001) {
			// write(6,*) 'flagged all photon intereaction types';
			EGS4.iausfl[15] = 1; // "Before pair"
			EGS4.iausfl[17] = 1; // "Before Compton"
			EGS4.iausfl[18] = 1; // "After Compton"
			EGS4.iausfl[19] = 1; // "Before photoelectric"
			EGS4.iausfl[20] = 1; // "After photoelectric"
			EGS4.iausfl[23] = 1; // "Before Rayleigh"
			EGS4.iausfl[24] = 1; // "After Rayleigh"
		} else {
			for (int j = 1; j <= EGS4.$MXREG; j++) {
				EGS4Macro.iefl[j - 1] = 0;
			}
		}

		for (int j = 1; j <= 28; j++) {
			if (EGS4.iausfl[j - 1] != 0) {// write(6,'(i3,$)') j;
				int jj = j - 1;
				EGS4.seqStr = " AUSGAB index call: " + jj;
				if (EGS4.iprint > 2)
					printSequence(EGS4.seqStr);
			}
		}

		// "HATCH CALL PREPARATION AND EXECUTION"
		EGS4.DUNIT = 1; // "SET LENGTH UNITS TO CMS"
		HATCH();
		if (EGS4.STOPPROGRAM) {
			return;
		}
		/*
		 * if( EGS4Macro.irejct == 1 ) { EGS4Macro.initialize_range_rejection();
		 * }
		 */
		if (EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21
				&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {// "MONOENERGETIC INPUT BEAM"
			if (EGS4SrcEns.iqin == 0) {
				EI = EGS4SrcEns.ein;
			} else {
				EI = EGS4SrcEns.ein + EGS4.RM;
			}
			EKMAX = EGS4SrcEns.ein; // "MAXIMUM KINETIC ENERGY"
		} else if (EGS4SrcEns.MONOEN == 1) {// "ENERGY SPECTRUM"
			EGS4SrcEns.ENSRC1();// "NORMALIZE THE ENERGY DISTRIBUTION"
			EKMAX = EGS4SrcEns.ENSRCD[EGS4SrcEns.NENSRC];// "MAXIMUM KINETIC ENERGY IN THE SPECTRUM"
		} else if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {// "phase space input"
		// EKMAX=EKSRCM;//NOT ALLOWED
		}
		// else { EKMAX=0; }// " <------------ fixme"

		// "CHECK THAT THE DATA FILE HAD DATA OVER THE ENERGY RANGE REQUIRED"
		for (int I = 1; I <= EGS4.NMED; I++) {
			if ((EKMAX > EGS4.UP[I - 1]) || (EKMAX > EGS4.UE[I - 1] - EGS4.RM)) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " FOR MEDIUM: " + I + "  INCIDENT ENERGY="
						+ EGS4.format(EKMAX, 10, true) + " MeV";
				printSequence(EGS4.seqStr);
				EGS4.seqStr = " IS GREATER THAN COVERED BY DATA FILE WHERE UP,UE="
						+ EGS4.format(EGS4.UP[I - 1], 10, true)
						+ " "
						+ EGS4.format(EGS4.UE[I - 1], 10, true) + " MeV";
				printSequence(EGS4.seqStr);
				EGS4.seqStr = " EXECUTION WILL BE TERMINATED!";
				printSequence(EGS4.seqStr);

				// return;
			}
		}// "END OF LOOP OVER MEDIA"
		if (EGS4.STOPPROGRAM) {
			return;
		}

		// "CALCULATE THE MASS OF EACH ZONE (AREAL MASS FOR ISOURC=2 OR 4)
		for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
			DEPTH = EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1];
			for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
				// $GET-IRL(IZ,IX);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				MEDNUM = EGS4.MED[IRL - 1];
				if (MEDNUM != 0) {
					if ((EGS4SrcEns.ISOURC == 2) || (EGS4SrcEns.ISOURC == 4)) {
						VOLUME = DEPTH;
					} else {
						if (IX == 1) {
							RLOW2 = 0.0;
						} else {
							RLOW2 = EGS4Geom.CYRAD2[IX - 2];
						}
						VOLUME = Math.PI * DEPTH
								* (EGS4Geom.CYRAD2[IX - 1] - RLOW2);
					}
					AMASS[IZ - 1][IX - 1] = EGS4.RHOR[IRL - 1] * VOLUME;
				} else {
					AMASS[IZ - 1][IX - 1] = 0.0;
				}
			}// "END OF IX LOOP"
		}// "END OF IZ LOOP"

		// "CALCULATE ONE-TIME-ONLY CONSTANTS FOR SOURCE"
		EGS4SrcEns.SRCOTO();// (WEIGHT);
		if ((EGS4SrcEns.IFPB == 0) && (EGS4Macro.IFORCE != 0)
				&& (EGS4SrcEns.iqin == 0) && (EGS4SrcEns.MONOEN == 0)) {
			EGS4Macro.SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM();
			// -------------
			EGS4Macro.do_fast_step = true;
		}
		// -------------
		else {
			EGS4Macro.do_fast_step = false;
		}

		// "INITIALIZE DATA ARRAYS FOR FLUORESCENT X-RAYS IF NEEDED"
		ISUMX = 0;
		for (int JJ = 1; JJ <= EGS4Geom.nreg; JJ++) {
			ISUMX = ISUMX + EGS4.iedgfl[JJ - 1]; // "NON-ZERO IF X-RAYS ANYWHERE"
		}
		if (ISUMX != 0) {
			EGS4.EDGSET(2, EGS4Geom.nreg);
		}
		// "NOTE THE ABOVE WILL PRODUCE LOTS OF EXTRA OUTPUT AND SHOULD BE"
		// "CLEANED UP"

		EK0 = EGS4SrcEns.ein;
		ISUMRY(); // "PRINT THE SUMMARY OF INPUTS"

		// "******************************************************************************
		// "
		// " *** SECTION 2 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "LOOP THROUGH THE NUMBER OF HISTORIES. CALCULATE CONSTANTS THAT MAY
		// CHANGE FOR
		// "EACH HISTORY AND DO THE SIMULATION
		// "
		// "------------------------------------------------------------------------------

		// "WRITE THE HEADER"
		// WRITE(IOUT,100) TITLE; call egs_fdate(iout); write(iout,*);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   EXECUTION INFORMATION";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,200);WRITE(6,200); "PRINT HEADER FOR EXECUTION MESSAGES"

		// "PRINT EXECUTION MODE"
		// if(IRESTART == 0){WRITE(6,201);WRITE(IOUT,201);}
		// ELSEIF(IRESTART.EQ.1)[
		// WRITE(6,202) NCASE,NCASEO;
		// write(6,'(21x,a,$)') 'New RNG state: ';
		// $SHOW-RNG-STATE(6); write(6,*);
		// write(iout,'(21x,a,$)') 'New RNG state: ';
		// $SHOW-RNG-STATE(iout); write(iout,*);
		// ]
		// ELSEIF(IRESTART.EQ.3)[WRITE(6,204);WRITE(IOUT,204);GO TO :END-SIM:;]
		// ELSEIF(IRESTART.EQ.4)[WRITE(6,205);WRITE(IOUT,205);]
		// ELSEIF(IRESTART.EQ.5)[WRITE(6,206);WRITE(IOUT,206);GO TO :END-SIM:;]

		// "Initialize IWATCH routine"
		if (IWATCH != 0)
			EGS4.WATCH(-99, IWATCH);

		// "SET CLOCK AT THE BEGINNING OF SIMULATIONS"
		// $INITIALIZE_ELAPSED_CPU_TIME;
		// $SET_ELAPSED_CPUTIME(CPUT1);
		// $INITIALIZE_ELAPSED_TOTAL_TIME;
		// ETIMETOT=0;
		// TIMEB=0;
		// NETADJ=0;

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// IF( idat = 0 ) data_unit = egs_open_file(4,0,1,'.egsdat');
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// "Calculate the cavity mass for in-flight dose display"
		/*
		 * FMASSC=0.0; //"TOTAL CAVITY MASS" for(int IX=1;IX<=EGS4Geom.NR;IX++)
		 * { for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++) { //$GET-IRL(IZ,IX);
		 * IRL=EGS4Geom.GET_IRL(IZ, IX); if(EGS4Geom.ntrack[IRL-1]==1) {
		 * FMASS=AMASS[IZ-1][IX-1]; FMASSC=FMASSC+FMASS;
		 * //"SUM THE CAVITY MASS USED LATER" } } }
		 * //"end of calculation of cavity mass"
		 */
		/*
		 * //"Initialize variables for the cs_enhance scoring " tmp_dose = 0;
		 * tmp_dose1 = 0; last_case = 0;EGS4SrcEns.NHSTRY = 0;
		 * 
		 * //"Open file for data storage, if requested "
		 * //"The file is opened in the temporary working directory" //IF( idat
		 * = 0 ) data_unit = egs_open_file(4,0,1,'.egsdat');
		 * 
		 * n_photons=0; n_electrons = 0; sumE_photons=0.0; sumE2_photons=0.0;
		 * sumE_electrons=0.0; sumE2_electrons=0.0;
		 */
		// "IK: New parallel processing implementation. Only used if there is a
		// " working C compiler.
		// #ifdef HAVE_C_COMPILER;
		// ;
		// /last_dose,last2_dose/ = 0; n_tot = ncaseo;
		// first_time = .true.; is_finished = .false.;

		// :start_parallel_loop:;

		// IF( n_parallel > 0 ) [ "Job is part of a parallel run "

		// last_dose = cav_dose - last_dose + tmp_dose;
		// last2_dose = cav2_dose - last2_dose + tmp_dose*tmp_dose;
		// tmpf = fmassc*ainflu/ncaset/1.602e-10;
		// part_dose = last_dose/tmpf; part2_dose = last2_dose/(tmpf*tmpf);
		// last_dose = cav_dose + tmp_dose;
		// last2_dose = cav2_dose + tmp_dose*tmp_dose;
		// call egs_pjob_control(ncase,n_run,n_left,n_tot,part_dose,part2_dose,
		// current_result, current_uncertainty);
		// IF( n_run = 0 ) [
		// write(6,'(//a,a//)') '****** No histories left in job control file',
		// ' => end simulation';
		// goto :END-SIM:;
		// ]
		// IF( statlm > 0 & current_uncertainty < statlm ) [
		// write(6,'(//a,a//)') '****** Desired uncertainty reached',
		// ' => end simulation';
		// goto :END-SIM:;
		// ]
		// jcase = n_run/$NBATCH;
		// IF( jcase < 1 ) [ jcase = 1; n_run = jcase*$NBATCH; ]
		// IF( first_time ) [
		// first_time = .false.; n_last = n_run;
		// write(6,'(//a,i12,a//)') '****** Running ',n_run,' histories';
		// ]
		// ELSE [
		// write(6,'(//a,i12,a)') '***** Finished ',n_last,' histories';
		// write(6,'(/a/,20x,1pe11.4,a,0pf5.2,a/,a,i12,a//)')
		// ' current result including previous runs and other parallel jobs: ',
		// current_result, ' +/- ',current_uncertainty,' %',
		// ' will run another ',n_run,' histories';
		// ]
		// ]
		// #endif;

		// "Output batches. Statistical analysis is done after each batch.
		// Execution
		// "stops if the desired statistical accuracy is obtained or there is
		// not enough
		// "time to do another batch.
		boolean ENDSIM = false;
		for (int IBATCH = 1; IBATCH <= $NBATCH; IBATCH++) {
			ENDSIM = false;
			long startTime = System.currentTimeMillis();
			IBTCH = IBATCH;
			if (IBATCH == 1) {
				EGS4.seqStr = " BATCH" + EGS4.format("", 5) + "ELAPSED"
						+ EGS4.format("", 3) + "time of day"
						+ EGS4.format("", 2) + "peak region stats(%)";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
				// //"IK: it is annoing that for batch runs we don't see the progress"
				// //"    info in the log file until the job is finished. This is because"
				// //"    Fortran uses buffered I/O. The following flushes unit 6 so that"
				// //"    we can see the progress of the calculation. "
				// //$FLUSH_UNIT(6);
			} else {// " not first batch"
			// $SET_ELAPSED_TOTAL_TIME(TIMEB);
			// / *
			// we don't need this because the current EGSnrc implementation
			// takes care of run times > 24 h.

				// IF(TIMEB < 0.0)["batch past midnight relative to TZERO"
				// NETADJ=NETADJ+1;
				// TIMEB=TIMEB + 86400.;
				// IF(NETADJ<2)[
				// OUTPUT; (' Wall clock has gone past 24:00 hrs.'/
				// ' Elapsed time adjusted assuming batches took < 1 day',
				// ' to complete.');
				// ]
				// ]
				// * /
				// ETIMETOT = ETIMETOT+TIMEB;
				// $INITIALIZE_ELAPSED_TOTAL_TIME;"reset to get elapsed time/batch"
				// $SET_ELAPSED_CPUTIME(CPUT2);
				// TIMCPU =
				// (CPUT2-CPUT1)*$CONVERSION_TO_SECONDS+$TIME_RESOLUTION;

				// "*****************************************************************"
				// OUTPUT IBATCH,ETIMETOT,TIMCPU,ETIMETOT/TIMCPU;
				// (1X,I2,F8.1,1X,F8.1,2X,F8.2,3X,' ',$); call egs_time(6);
				// "*****************************************************************"

				// "IK: it is annoing that for batch runs we don't see the progress"
				// "    info in the log file until the job is finished. This is because"
				// "    Fortran uses buffered I/O. The following flushes unit 6 so that"
				// "    we can see the progress of the calculation. "
				// $FLUSH_UNIT(6);

				// "Check there is time left for another batch"
				// BATCHT = TIMCPU/dble(IBATCH-1);//"time per batch so far"
				// IF(TIMCPU+1.1*BATCHT.GT.TIMMAX*3600.)[
				// "not enough time for another batch"
				// "print message and exit simulation loop"
				// WRITE(IOUT,210) TIMMAX,IBATCH-1,IHSTRY-NCASEO,IHSTRY;
				// WRITE(6,210) TIMMAX,IBATCH-1,IHSTRY-NCASEO,IHSTRY;
				// "adjust the incident fluence"
				// "IK: do it at the end for all possible exits from the
				// "shower loop
				// "AINFLU = AINFLU*dble(IHSTRY)/dble(NCASET);
				// GO TO :END-SIM:;
				// }
			}// " end of before batch ne 1 block"

			for (int ICASE = 1; ICASE <= JCASE; ICASE++) {// "now fill each IS bin"

				if (EGS4SrcEns.ISOURC != 23)
					IHSTRY = IHSTRY + 1; // "increment history counter"
					// " For source 23 ihstry is set in srchst "

				// :CORRELATION-RESTART:;
				// //"restart here for a correlated history"

				EGS4Macro.NFTIME = 0; // "reset the photon forced interaction counter"

				// "retrieve starting random numbers if reading from a file"
				// IF(IRESTART.EQ.4) [
				// $RETRIEVE RNG STATE FROM UNIT rng_unit; "was $RESET-RNG(2);"
				// ]

				// "store the initial random number, another pass might be needed"
				// IRNG = 1;
				// $STORE-RNG(IRNG); "was $STORE-RNG(0);"

				// "store initial random #s if requested"
				// IF(ISTORE = 1)[
				// $STORE RNG STATE ON UNIT rng_unit;
				// ]
				// ELSEIF(ISTORE = 2)[
				// "temporarily store the initial random number seed"
				// IRNG = 3;
				// $STORE-RNG(IRNG); "was $STORE-RNG(-2);"
				// "clear the flag that signals energy deposit in the cavity"
				// IDECAV = 0;
				// ]
				// ELSEIF(ISTORE.EQ.3)[
				// "STORE THE INITIAL RANDOM NUMBER SEED"
				// $PUT RNG STATE ON UNIT rng_unit;
				// ]

				// "calculate the source dependent values which change for each
				// "history these include :
				// "entry point into target,
				// "initial direction cosines,
				// "statistical weight,
				// "entry flag(nrcflg)
				// CALL SRCHST(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
				EGS4SrcEns.SRCHST();

				// "calculate the initial energy if a distribution is to be used"
				if (EGS4SrcEns.MONOEN != 0 && EGS4SrcEns.ISOURC != 21
						&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {// "if equal to 0, it is monoenergetic"
				// EGS4EnsSrc.ENSRCH(EIN); "returns K.E. from distribution"
					EGS4SrcEns.ein = EGS4SrcEns.ENSRCH();
					if (EGS4SrcEns.iqin == 0) {
						EI = EGS4SrcEns.ein;
					} else {
						EI = EGS4SrcEns.ein + EGS4.RM;
					}// "total energy"
					// " there was a check that the data file had data over the
					// energy
					// "range required, the location of it will eventually be in
					// "ESRCIN.MOR
				}
				// ELSEIF(ISOURC = 21 | ISOURC = 22 | ISOURC = 23 )[ EI = EIN; ]
				// IF((ISOURC=21 | ISOURC=22 | ISOURC=23) & IFULL=2 &
				// WEIGHT<1.)[
				// OUTPUT WEIGHT; (//' ****WARNING****'/
				// ' A particle of weight ',F10.5,' from the phase space source
				// is being used'/
				// ' Pulse height distribution only makes sense if all particles
				// have'/
				// ' weight=1. Will run simulation with IFULL= dose and
				// stoppers'//);
				// WRITE(1,'(//'' ****WARNING****''/
				// '' A particle of weight '',F10.5,'' from the phase space
				// source is being used''/
				// '' Pulse height distribution only makes sense if all
				// particles have''/
				// '' weight=1. Will run simulation with IFULL= dose and
				// stoppers''//)');
				// IFULL=0;
				// ]

				// "Set photon weights if gamma interactions are to be forced in
				// the
				// "target in the frontal parallel beam case if monoenergetic
				/*
				 * if((EGS4SrcEns.MONOEN==0)&&
				 * (EGS4SrcEns.iqin==0)&&(EGS4Macro.IFORCE
				 * ==1)&&(EGS4SrcEns.IFPB==0)&&
				 * (EGS4SrcEns.ISOURC!=21)&&(EGS4SrcEns.ISOURC!=22)) { int
				 * IX=(EGS4SrcEns.irin-2)/EGS4Geom.NZ+1;
				 * EGS4Macro.GWAIT=EGS4Macro.GWATE[IX-1];
				 * EGS4SrcEns.WEIGHT=EGS4Macro.GWAIT; }
				 */
				if (EGS4Macro.do_fast_step) {
					int IX = (EGS4SrcEns.irin - 2) / EGS4Geom.NZ + 1;
					EGS4Macro.GWAIT = EGS4Macro.GWATE[IX - 1];
					EGS4SrcEns.WEIGHT = EGS4Macro.GWAIT;
				}

				// "FOR AN INPUT ENERGY SPECTRUM, DETAILED FORCING MACRO IS USED"

				EGS4.LATCHI = 0;
				// "SET INITIAL DOSE COMPONENTS"
				if (EGS4Macro.IFULL == 1) {
					EGS4Macro.NEWNRC = EGS4SrcEns.NRCFLG;
					// EGS4.LATCHI=//IBSET(EGS4.LATCHI,EGS4SrcEns.NRCFLG/10);
					EGS4.LATCHI = EGS4.IBSET_LATCHI(EGS4SrcEns.NRCFLG / 10);
				}

				if ((IWATCH != 0) && (IWATCH != 4)) {
					EGS4.seqStr = EGS4.format(1, 5)
							+ EGS4.format(EGS4SrcEns.ein, 9, true)
							+ EGS4.format(EGS4SrcEns.iqin, 4)
							+ EGS4.format(EGS4SrcEns.irin, 4)
							+ EGS4.format(EGS4SrcEns.xin, 8, true)
							+ EGS4.format(EGS4SrcEns.yin, 8, true)
							+ EGS4.format(EGS4SrcEns.zin, 8, true)
							+ EGS4.format(EGS4SrcEns.uin, 8, true)
							+ EGS4.format(EGS4SrcEns.vin, 8, true)
							+ EGS4.format(EGS4SrcEns.win, 8, true)
							+ EGS4.format(EGS4.LATCHI, 10)
							+ EGS4.format(EGS4SrcEns.WEIGHT, 10, false);
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
					// OUTPUT
					// 1,EIN,IQIN,IRIN,XIN,YIN,ZIN,UIN,VIN,WIN,LATCHI,WEIGHT;
					// (/' INITIAL SHOWER VALUES',T36,':',
					// I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
				}
				// "ALL INITIAL SHOWER VARIABLES ARE SET, CALL THE SHOWER ROUTINE"
				/*
				 * if( EGS4SrcEns.iqin == 0 ) { n_photons = n_photons + 1;
				 * sumE_photons = sumE_photons + EI; sumE2_photons =
				 * sumE2_photons + EI*EI; } else { n_electrons = n_electrons +
				 * 1; sumE_electrons = sumE_electrons + EI-EGS4.RM;
				 * sumE2_electrons = sumE2_electrons +
				 * (EI-EGS4.RM)*(EI-EGS4.RM); }
				 */
				// CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
				SHOWER();
				if (EGS4.STOPPROGRAM) {
					return;
				}

				// IF(ISTORE.EQ.2)[
				// "STORE INITIAL RN SEED IF ENERGY SCORED IN THE CAVITY"
				// IF(IDECAV.EQ.1)[ "ENERGY WAS DEPOSITED IN THE CAVITY"
				// IRNG = 1;
				// $STORE-RNG(IRNG); "was $STORE-RNG(0);"
				// IRNG = 3;
				// $RESET-RNG(IRNG); "was $RESET-RNG(-2);"
				// $PUT RNG STATE ON UNIT rng_unit;"was $STORE-RNG(2);"
				// IRNG = 1;
				// $RESET-RNG(IRNG); "was $RESET-RNG(0);"
				// ]
				// ]

				// "SIGNAL THE END OF A HISTORY IF WATCH MODE IS SET"
				// "JAN CHANGE: ZERO STACK ELEMENTS TO PREVENT CROSS TALKING"
				// "IK: Huh ? $CLEAN-STACK;"
				// "END OF JAN CHANGE"
				if (IWATCH > 0)
					EGS4.WATCH(-1, IWATCH);

			}// "END OF THE ICASE LOOP"

			// "Succesful completion of a batch. delete the raw data from the last batch"
			// "and record the new batch only if requested"
			// if(IDAT == 0)
			// {

			// " discard data from previous batches (if any) "
			// rewind(data_unit);

			// TSCSTP=SCSTP+SCSTP_TMP;
			// TSCSTP2=SCSTP2+SCSTP_TMP*SCSTP_TMP;
			// TSCCSTP=SCCSTP+SCCSTP_TMP;
			// TSCCSTP2=SCCSTP2+SCCSTP_TMP*SCCSTP_TMP;

			// WRITE(data_unit,*)TSCSTP,TSCSTP2,TSCCSTP,TSCCSTP2;

			// "******************"
			// "history by history"
			// "  EMH  March, 2002"
			// "******************"
			// tcav_dose = cav_dose + tmp_dose;
			// tcav2_dose = cav2_dose + tmp_dose*tmp_dose;

			// tcav_dose0 = cav_dose0 + tmp_dose0;
			// tcav2_dose0 = cav2_dose0 + tmp_dose0*tmp_dose0;

			// tcav_dose1 = cav_dose1 + tmp_dose1;
			// tcav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1;

			// tcav_dose2 = cav_dose2 + tmp_dose2;
			// tcav2_dose2 = cav2_dose2 + tmp_dose2*tmp_dose2;

			// tcav_dosec = cav_dosec + tmp_dose*tmp_dose1;
			// tcav_dosec01 = cav_dosec01 + tmp_dose0*tmp_dose1;
			// tcav_dosec02 = cav_dosec02 + tmp_dose0*tmp_dose2;

			// write(data_unit,*) tcav_dose, tcav_dose0, tcav_dose1, tcav_dose2;
			// write(data_unit,*)
			// tcav2_dose,tcav2_dose0,tcav2_dose1,tcav2_dose2;
			// write(data_unit,*) tcav_dosec,tcav_dosec01,tcav_dosec02;

			// "*****************"

			// IF(NSUMCV>1)["store data for individual cavity regions"

			// DO IZ=1,NZ[
			// DO IX=1,NR[
			// $GET-IRL(IZ,IX);
			// IF(NTRACK(IRL).EQ.1)[
			// DO IT=1,4[
			// TSCDOSE(IZ,IX,IT)=SCDOSE(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT);
			// TSCDOSE2(IZ,IX,IT)=SCDOSE2(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT)*
			// SCDOSE_TMP(IZ,IX,IT);
			// ]
			// TSCDOSE_COV(IZ,IX,1)=SCDOSE_COV(IZ,IX,1)+
			// SCDOSE_TMP(IZ,IX,1)*SCDOSE_TMP(IZ,IX,3);
			// TSCDOSE_COV(IZ,IX,2)=SCDOSE_COV(IZ,IX,2)+
			// SCDOSE_TMP(IZ,IX,2)*SCDOSE_TMP(IZ,IX,3);
			// TSCDOSE_COV(IZ,IX,3)=SCDOSE_COV(IZ,IX,3)+
			// SCDOSE_TMP(IZ,IX,2)*SCDOSE_TMP(IZ,IX,4);
			// WRITE(data_unit,*)
			// (TSCDOSE(IZ,IX,IT),TSCDOSE2(IZ,IX,IT),IT=1,4);
			// WRITE(data_unit,*)(TSCDOSE_COV(IZ,IX,IT),IT=1,3);
			// ]
			// ]
			// ]
			// ]

			// ]"end of conditional data storage"

			// $SET_ELAPSED_CPUTIME(CPUT2);
			// TIMCPU=$CONVERSION_TO_SECONDS*(CPUT2-CPUT1)+TMCPUO;
			// IF(IDAT.EQ.0)[
			// $PUT RNG STATE ON UNIT data_unit;
			// WRITE(data_unit,*) IHSTRY,TIMCPU,NNREAD,PIISTP+count_pII_steps;
			// ]
			// "******************"
			// "history by history"
			// "    EMH March 2002"
			// "******************"
			/*
			 * IF(IDAT.EQ.0)[
			 * "add unscored portions of _TMP arrays before writing"
			 * SCSTP=SCSTP+SCSTP_TMP; SCSTP2=SCSTP2+SCSTP_TMP*SCSTP_TMP;
			 * SCSTP_TMP=0; SCDSTP=SCDSTP+SCDSTP_TMP;
			 * SCDSTP2=SCDSTP2+SCDSTP_TMP*SCDSTP_TMP; SCDSTP_TMP=0; DO
			 * IZ=1,NZDOSE[ DO IX=1,NRDOSE[ DO IT=1,ITMAX[
			 * SCDOSE(IZ,IX,IT)=SCDOSE(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT);
			 * SCDOSE2(IZ,IX,IT)=SCDOSE2(IZ,IX,IT)+SCDOSE_TMP(IZ,IX,IT)*
			 * SCDOSE_TMP(IZ,IX,IT); IF(IKERMA=1)[
			 * SCKERMA(IZ,IX,IT)=SCKERMA(IZ,IX,IT)+SCKERMA_TMP(IZ,IX,IT);
			 * SCKERMA2(IZ,IX,IT)=SCKERMA2(IZ,IX,IT)+SCKERMA_TMP(IZ,IX,IT)*
			 * SCKERMA_TMP(IZ,IX,IT);
			 * IF(SCDOSE_LAST(IZ,IX,IT)=SCKERMA_LAST(IZ,IX,IT))[
			 * "same history being scored, include these in"
			 * "uncertainty estimate"
			 * SCDOSEtoKERMA2(IZ,IX,IT)=SCDOSEtoKERMA2(IZ,IX,IT)+
			 * SCDOSE_TMP(IZ,IX,IT)*SCKERMA_TMP(IZ,IX,IT); ]
			 * SCKERMA_TMP(IZ,IX,IT)=0.; ] SCDOSE_TMP(IZ,IX,IT)=0.; ] ] ]
			 * "OPEN(UNIT=4,file='fort.4',STATUS='UNKNOWN');"
			 * "IK: the .egsdat file is opened before the shower loop."
			 * "    The only reason for opening and closing the I/O unit is to "
			 * "    overwrite data from previous batches, but we can do this "
			 * "    by rewinding the unit" rewind(data_unit);
			 * WRITE(data_unit,*)SCSTP
			 * ,SCSTP2,SCDSTP,SCDSTP2,PIISTP+count_pII_steps; WRITE(data_unit,*)
			 * (((SCDOSE(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
			 * WRITE(data_unit,*)
			 * (((SCDOSE2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
			 * IF(IKERMA=1)[ WRITE(data_unit,*)
			 * (((SCKERMA(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
			 * WRITE(data_unit,*)
			 * (((SCKERMA2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
			 * WRITE(data_unit,*)
			 * (((SCDOSEtoKERMA2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE
			 * ),IZ=1,NZDOSE); ] IF(IFULL=2)[ IF(PHENER>0.)[
			 * "FIND WHAT BIN WE ARE IN" IF(SLOTE.GT.0.0)[
			 * "EQUAL ENERGY BINS CASE" IB=MIN0(IFIX(PHENER/SLOTE+0.999),$EBIN);
			 * ] ELSE[ IB = MAXBIN; UNTIL((IB.EQ.1).OR.(BINTOP(IB-1).LT.PHENER))
			 * [IB=IB-1;] ]
			 * 
			 * "ACCUMULATE THE PULSE HEIGHT DISTRIBUTION"
			 * SCPDST(IB)=SCPDST(IB)+WT1OLD;
			 * SCPDST2(IB)=SCPDST2(IB)+WT1OLD*WT1OLD; "add to cumulative distn"
			 * DO ICUM=IB,MAXBIN[ SCPCUM(ICUM)=SCPCUM(ICUM)+WT1OLD;
			 * SCPCUM2(ICUM)=SCPCUM2(ICUM)+WT1OLD*WT1OLD; ]
			 * 
			 * IF(IWATCH.EQ.3)[ OUTPUT PHENER,IB,1; (' PULSE HEIGHT ENERGY=',
			 * F10.4,' MeV, IN BIN',I3,' WITH WEIGHT',1PE10.3); ]
			 * 
			 * "NOW SCORE PROBABILITIES FOR COUNTS IN PEAKS" DO IPK=1,4[
			 * "FOR EACH PEAK, F.E., ESCAPES AND 511"
			 * IF((PHENER.GE.DFEN(IPK,2)).AND.(PHENER.LE.DFEN(IPK,3)))[
			 * "IT IS IN THE PEAK" SCDFEP(IPK) = SCDFEP(IPK)+WT1OLD;
			 * SCDFEP2(IPK) = SCDFEP2(IPK)+WT1OLD*WT1OLD;
			 * SCDFDIFF(IPK)=SCDFDIFF(IPK)+WT1OLD;
			 * SCDFDIFF2(IPK)=SCDFDIFF2(IPK)+WT1OLD*WT1OLD; IF(IWATCH.EQ.3)[
			 * OUTPUT IPK;(T50,'IT WAS IN ONE OF THE PEAKS,IPK=',I3/); ] ]
			 * ELSEIF((PHENER.GE.DFEN(IPK,1)).AND.(PHENER.LT.DFEN(IPK,2)))[
			 * "IT IS IN THE BKGD" SCDFBK(IPK)=SCDFBK(IPK)+WT1OLD;
			 * SCDFBK2(IPK)=SCDFBK2(IPK)+WT1OLD*WT1OLD;
			 * SCDFDIFF(IPK)=SCDFDIFF(IPK)-WT1OLD;
			 * SCDFDIFF2(IPK)=SCDFDIFF2(IPK)-WT1OLD*WT1OLD; ] ]"END IPK LOOP"
			 * PHENER=0; ]
			 * WRITE(data_unit,*)(SCPDST(IB),SCPDST2(IB),IB=1,MAXBIN);
			 * WRITE(data_unit,*)(SCPCUM(IB),SCPCUM2(IB),IB=1,MAXBIN);
			 * WRITE(data_unit,*)(SCDFEP(IPK),SCDFEP2(IPK),IPK=1,4);
			 * WRITE(data_unit,*)(SCDFBK(IPK),SCDFBK2(IPK),IPK=1,4);
			 * WRITE(data_unit,*)(SCDFDIFF(IPK),SCDFDIFF2(IPK),IPK=1,4); ]
			 * ]"END OF CONDITIONAL DATA STORAGE"
			 * 
			 * $SET_ELAPSED_CPUTIME(CPUT2);
			 * TIMCPU=$CONVERSION_TO_SECONDS*(CPUT2-CPUT1)+TMCPUO;
			 * IF(IDAT.EQ.0)[ $PUT RNG STATE ON UNIT data_unit;
			 * WRITE(data_unit,*) IHSTRY,TIMCPU,NNREAD; write(data_unit,*)
			 * SCOMEG,SCOMEG2;
			 * "IK: don't need to close, see ebove. CLOSE(UNIT=4);" ]
			 */
			/*
			 * if(EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22) { //
			 * AINFLU_CURRENT= //
			 * EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/
			 * EGS4SrcEns.dble(NCASE_PHSP)*NINCSRC; } else if(
			 * EGS4SrcEns.ISOURC==23 ) { AINFLU_CURRENT=IHSTRY; } else {
			 * AINFLU_CURRENT
			 * =EGS4SrcEns.AINFLU*EGS4SrcEns.dble(IHSTRY)/EGS4SrcEns
			 * .dble(EGS4SrcEns.NCASET); } dcav_current =
			 * cav_dose*1.602E-10/(FMASSC*AINFLU_CURRENT); //double
			 * mevtojoule=1.60218E-13;//1MeV=1,60218 · 10-13 J //and MASS is in
			 * g not KG!!!=>E-10!! // // * //IF( ISOURC=23 ) [ // dcavun =
			 * (cav2_dose*nhstry-cav_dose*cav_dose)/(nhstry-1); //] //ELSE [ //
			 * dcavun = (cav2_dose*ihstry-cav_dose*cav_dose)/(ihstry-1); //] //
			 * * /
			 * 
			 * dcavun = (cav2_dose*IHSTRY-cav_dose*cav_dose)/(IHSTRY-1); if(
			 * dcavun > 0 ) { dcavun = Math.sqrt(dcavun); }; dcavun =
			 * dcavun*1.602E-10/(FMASSC*AINFLU_CURRENT); dcavun =
			 * 100*dcavun/dcav_current;
			 * 
			 * //################################################################
			 * ##################### String
			 * timePerBatch=EGS4.timeElapsedShort(startTime);
			 * 
			 * Calendar call=Calendar.getInstance(); String
			 * timeday=call.get(cal.HOUR)+":"+
			 * call.get(cal.MINUTE)+":"+call.get(cal.SECOND);
			 * 
			 * EGS4.seqStr=EGS4.format("",2)+EGS4.format(IBATCH,3)+EGS4.format(""
			 * ,5)+timePerBatch+
			 * EGS4.format("",3)+timeday+EGS4.format("",6)+EGS4
			 * .format(dcavun,7,true)+
			 * EGS4.format("",10)+EGS4.format(dcav_current,11,false);
			 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
			 * 
			 * //OUTPUT; //(/' BATCH',2X,'ELAPSED',2X,'CPUtime',2X,'RATIO',2X,
			 * //'time of day',2X,'cavity stats(%)',2X,'cav.dose(Gy.cm^2)'//
			 * //2X,'1',5X,'0.0',6X,'0.0',6X,'0.00',3X,' ',$); call egs_time(6);
			 * 
			 * //"IK: it is annoing that for batch runs we don't see the progress"
			 * //
			 * "    info in the log file until the job is finished. This is because"
			 * //
			 * "    Fortran uses buffered I/O. The following flushes unit 6 so that"
			 * //"    we can see the progress of the calculation. "
			 * //$FLUSH_UNIT(6);
			 * //##############################################
			 * ####################################### //OUTPUT
			 * DCAVUN,DCAV_CURRENT;('+',5X,F5.2,10X,1PE11.4,5X,1PE10.3);
			 * 
			 * if(dcavun<=STATLM && STATLM>0.) {
			 * //"we have reached the desired statistics, print a message and exit"
			 * EGS4.seqStr=" DESIRED STATISTICAL ACCURACY OBTAINED.";
			 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
			 * EGS4.seqStr=" STATS IN CAVITY= "+EGS4.format(dcavun,5,true)+"%"+
			 * " AFTER "+EGS4.format(IBTCH,2)+" BATCHES"; if(EGS4.iprint>1)
			 * printSequence(EGS4.seqStr);
			 * 
			 * //WRITE(6,230)dcavun,IBTCH;WRITE(IOUT,230)dcavun,IBTCH; //GO TO
			 * :END-SIM:; ENDSIM=true; break; }
			 */
			// "DO STATISTICAL ANALYSIS ON THE PEAK DOSE REGION TO SEE IF WE QUIT EARLY"
			TDSMAX = 0.0;
			for (IRL = 2; IRL <= EGS4Geom.nreg; IRL++) {
				if (EGS4Grid.CDSTBL[IRL - 1].compareTo(DCHAR) == 0) {
					// IZD=IDSTBL(IRL,1);IXD=IDSTBL(IRL,2);
					IZD = IDSTBL[IRL - 1][0];
					IXD = IDSTBL[IRL - 1][1];
					// FMASS=AMASS(IZD+NZDMIN-1,IXD+NRDMIN);
					FMASS = AMASS[IZD + NZDMIN - 2][IXD + NRDMIN - 1];
					// System.out.println("@@@ "+SCDOSE[IZD-1][IXD-1][0]);
					if ((SCDOSE[IZD - 1][IXD - 1][0] + SCDOSE_TMP[IZD - 1][IXD - 1][0])
							/ FMASS > TDSMAX) {
						// TDSMAX=(SCDOSE(IZD,IXD,1)+SCDOSE_TMP(IZD,IXD,1))/FMASS;
						TDSMAX = (SCDOSE[IZD - 1][IXD - 1][0] + SCDOSE_TMP[IZD - 1][IXD - 1][0])
								/ FMASS;
						IDSMAX = IRL;
					}
				}
			}

			// "NOW DO STATS ON THE PEAK REGION"
			if (TDSMAX > 0.0) {
				// IZD=IDSTBL[IDSMAX][1];IXD=IDSTBL(IDSMAX,2);
				IZD = IDSTBL[IDSMAX - 1][0];
				IXD = IDSTBL[IDSMAX - 1][1];
				// TDOS=SCDOSE(IZD,IXD,1)+SCDOSE_TMP(IZD,IXD,1);
				TDOS = SCDOSE[IZD - 1][IXD - 1][0]
						+ SCDOSE_TMP[IZD - 1][IXD - 1][0];
				// TDOS2=SCDOSE2(IZD,IXD,1)+SCDOSE_TMP(IZD,IXD,1)*SCDOSE_TMP(IZD,IXD,1);
				TDOS2 = SCDOSE2[IZD - 1][IXD - 1][0]
						+ SCDOSE_TMP[IZD - 1][IXD - 1][0]
						* SCDOSE_TMP[IZD - 1][IXD - 1][0];
				// "normalize by incident no. of primary histories--so far"
				if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {
					// SCORE_NORM_NUM=
					// EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*NINCSRC;
				} else {
					SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);
				}
				if (SCORE_NORM_NUM > 1.) {
					TDOS = TDOS / SCORE_NORM_NUM;
					TDOS2 = TDOS2 / SCORE_NORM_NUM;
					TDOS2 = (TDOS2 - TDOS * TDOS) / (SCORE_NORM_NUM - 1);
					if (TDOS2 >= 0.)
						TDOS2 = Math.sqrt(TDOS2);
					TDOS2 = Math.min(TDOS2 / TDOS * 100., 99.9);

					if ((TDOS2 <= STATLM) && (STATLM != 0.0)) {
						// "REACHED OBJECTIVE - PRINT MESSAGE AND JUMP OUT OF SIMULATION LOOP"
						// WRITE(6,230)IDSMAX,TDOS2,IBTCH;
						// 230 FORMAT(/' DESIRED STATISTICAL ACCURACY
						// OBTAINED.'/
						// ' STATS IN PEAK DOSE REGION (REGION ',I3,')= ',F6.3,
						// ' AFTER ',I2,' BATCHES');
						// WRITE(IOUT,230)IDSMAX,TDOS2,IBTCH;
						EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OBTAINED.";
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);
						EGS4.seqStr = " STATS IN PEAK DOSE REGION (REGION "
								+ EGS4.format(IDSMAX, 3)
								+ EGS4.format(TDOS2, 6, true) + "%" + " AFTER "
								+ EGS4.format(IBTCH, 2) + " BATCHES";
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);

						EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU
								* EGS4SrcEns.dble(IHSTRY)
								/ EGS4SrcEns.dble(EGS4SrcEns.NCASET); // "FIX NORM ON EARLY EXIT"
						// GO TO :END-SIM:;
						ENDSIM = true;
						break;
					}
				}
			}
			// OUTPUT TDOS2;(' ',9X,F6.3);
			// System.out.println(" ncaset "+EGS4SrcEns.NCASET+
			// "ainflu "+EGS4SrcEns.AINFLU);
			// EGS4SrcEns.AINFLU=EGS4SrcEns.AINFLU*EGS4SrcEns.dble(IHSTRY)/EGS4SrcEns.dble(EGS4SrcEns.NCASET);
			String timePerBatch = EGS4.timeElapsedShort(startTime);

			Calendar call = Calendar.getInstance();
			String timeday = call.get(Calendar.HOUR) + ":" + call.get(Calendar.MINUTE)
					+ ":" + call.get(Calendar.SECOND);

			EGS4.seqStr = EGS4.format("", 2) + EGS4.format(IBATCH, 3)
					+ EGS4.format("", 5) + timePerBatch + EGS4.format("", 3)
					+ timeday + EGS4.format("", 6)
					+ EGS4.format(TDOS2, 7, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}// "END OF SIMULATIONS"

		// #ifdef HAVE_C_COMPILER;
		// ;
		// IF( n_parallel > 0 ) [ goto :start_parallel_loop:; ]
		//
		// #endif;

		// "PRINT INSUFFICIENT STATS WARNING"
		// WRITE(IOUT,240) STATLM,TDOS2,IBTCH;WRITE(6,240) STATLM,TDOS2,IBTCH;
		// 240 FORMAT(/' *********DESIRED STATISTICAL ACCURACY OF ',F6.3,'%',
		// ' NOT REACHED*********'/
		// ' STATS IN PEAK DOSE REGION= ',F6.3,' % AFTER ',I2,' BATCHES');
		if (!ENDSIM) {
			EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OF "
					+ EGS4.format(STATLM, 6, true) + "%" + " NOT REACHED!";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " STATS IN PEAK DOSE REGION= "
					+ EGS4.format(TDOS2, 6, true) + " %" + " AFTER "
					+ EGS4.format(IBTCH, 2) + " BATCHES";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		// :END-SIM:;
		/*
		 * write(6,'(/a,$)') '********* FINAL RANDOM NUMBER STATE:';
		 * $SHOW-RNG-STATE(6); write(6,'(a)') ' *********'; write(iout,'(/a,$)')
		 * '********* FINAL RANDOM NUMBER STATE:'; $SHOW-RNG-STATE(iout);
		 * write(iout,'(a)') ' *********';
		 * 
		 * $SET_ELAPSED_TOTAL_TIME(TIMEB); ETIMETOT=ETIMETOT+TIMEB;
		 * $SET_ELAPSED_CPUTIME(CPUT2);
		 * TIMCPU=(CPUT2-CPUT1)*$CONVERSION_TO_SECONDS+TMCPUO;
		 * IF(TMCPUO.EQ.0)[RATIO=ETIMETOT/TIMCPU;]ELSE[RATIO=0.;]
		 * "ONLY FOR SINGLE RUN"
		 * IF(IRESTART=3)["just analyzing data--no elapsed time"
		 * WRITE(IOUT,250)TMCPUO,TMCPUO/3600; WRITE(6,250)TMCPUO,TMCPUO/3600;
		 * WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		 * WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; ]
		 * ELSEIF(IRESTART=5)["output time results for parallel runs"
		 * WRITE(IOUT,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT);
		 * WRITE(6,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT);
		 * WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO;
		 * WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; ] ELSE[
		 * IF(RATIO.NE.0.0)["I.E. WE HAVE CORRECT ELAPSED TIME"
		 * WRITE(IOUT,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO;
		 * WRITE(6,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO; ] ELSE[
		 * WRITE(IOUT,270) TIMCPU,TIMCPU/3600.; WRITE(6,270)
		 * TIMCPU,TIMCPU/3600.; ]
		 * IF((IHSTRY.NE.0).AND.(TIMCPU.NE.0.0))["THIS SHOULD ALWAYS HAPPEN"
		 * WRITE(IOUT,280)TIMCPU/dble(IHSTRY),3600.*dble(IHSTRY)/TIMCPU;
		 * WRITE(6,280)TIMCPU/dble(IHSTRY),3600.*dble(IHSTRY)/TIMCPU; ] ]
		 */
		/*
		 * SCSTP=SCSTP+SCSTP_TMP; SCSTP2=SCSTP2+SCSTP_TMP*SCSTP_TMP;
		 * SCCSTP=SCCSTP+SCCSTP_TMP; SCCSTP2=SCCSTP2+SCCSTP_TMP*SCCSTP_TMP;
		 * 
		 * //"******************" //"history by history" //"    EMH March 2002"
		 * //"******************" //"we have to add the temporary scoring"
		 * //"variables for the last history      " cav_dose = cav_dose +
		 * tmp_dose; cav2_dose = cav2_dose + tmp_dose*tmp_dose;
		 * 
		 * cav_dose0 = cav_dose0 + tmp_dose0; cav2_dose0 = cav2_dose0 +
		 * tmp_dose0*tmp_dose0;
		 * 
		 * cav_dose1 = cav_dose1 + tmp_dose1; cav2_dose1 = cav2_dose1 +
		 * tmp_dose1*tmp_dose1;
		 * 
		 * cav_dose2 = cav_dose2 + tmp_dose2; cav2_dose2 = cav2_dose2 +
		 * tmp_dose2*tmp_dose2;
		 * 
		 * cav_dosec = cav_dosec + tmp_dose*tmp_dose1; cav_dosec01 = cav_dosec01
		 * + tmp_dose0*tmp_dose1; cav_dosec02 = cav_dosec02 +
		 * tmp_dose0*tmp_dose2;
		 * 
		 * if(EGS4Geom.NSUMCV>1) {//"store data for individual cavity regions"
		 * 
		 * for (int IZ=1;IZ<=EGS4Geom.NZ;IZ++) { for(int
		 * IX=1;IX<=EGS4Geom.NR;IX++) { //$GET-IRL(IZ,IX);
		 * IRL=EGS4Geom.GET_IRL(IZ,IX);
		 * 
		 * if(EGS4Geom.ntrack[IRL-1]==1) { for(int IT=1;IT<=4;IT++) {
		 * SCDOSE[IZ-1][IX-1][IT-1]=SCDOSE[IZ-1][IX-1][IT-1]+
		 * SCDOSE_TMP[IZ-1][IX-1][IT-1];
		 * SCDOSE2[IZ-1][IX-1][IT-1]=SCDOSE2[IZ-1][IX-1][IT-1]+
		 * SCDOSE_TMP[IZ-1][IX-1][IT-1]*SCDOSE_TMP[IZ-1][IX-1][IT-1]; }
		 * SCDOSE_COV[IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]+
		 * SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][2];
		 * SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]+
		 * SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][2];
		 * SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]+
		 * SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][3]; } } } }
		 * 
		 * EGS4.seqStr=" FINAL RANDOM NUMBER STATE:"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr); EGS4.SHOW_RNG_STATE();
		 * //write(6,'(/a,$)') '****** FINAL RANDOM NUMBER STATE:';
		 * //$SHOW-RNG-STATE(6); write(6,'(a)') ' ******';
		 * //write(iout,'(/a,$)') '********* FINAL RANDOM NUMBER STATE:';
		 * //$SHOW-RNG-STATE(iout); write(iout,'(a)') ' *********';
		 * //$SET_ELAPSED_TOTAL_TIME(TIMEB); //ETIMETOT=ETIMETOT+TIMEB;
		 * //$SET_ELAPSED_CPUTIME(CPUT2);
		 * //TIMCPU=(CPUT2-CPUT1)*$CONVERSION_TO_SECONDS+TMCPUO;
		 * //IF(IRESTART=3)["just analyzing data--no elapsed time" //
		 * WRITE(IOUT,250)TMCPUO,TMCPUO/3600; // WRITE(6,250)TMCPUO,TMCPUO/3600;
		 * // WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; //
		 * WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; //]
		 * //ELSEIF(IRESTART=5)["output time results for parallel runs" //
		 * WRITE(IOUT,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT); //
		 * WRITE(6,255)DATCOUNT,TMCPUO,TMCPUO/3600.,TMCPUO/dble(DATCOUNT); //
		 * WRITE(IOUT,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; //
		 * WRITE(6,280)TMCPUO/dble(IHSTRY),3600.*dble(IHSTRY)/TMCPUO; //]
		 * //ELSE[ // IF(TMCPUO.EQ.0)["This is first run" //
		 * RATIO=ETIMETOT/TIMCPU; //
		 * WRITE(IOUT,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO; //
		 * WRITE(6,260)ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO; // ] // ELSE[
		 * "There was previous run, but don't have elapsed time for it" // RATIO
		 * = ETIMETOT/(TIMCPU-TMCPUO); "ratio this run only" // WRITE(IOUT,261)
		 * ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO; // WRITE(6,261)
		 * ETIMETOT,TIMCPU,TIMCPU/3600.,RATIO; // ] // IF(IHSTRY.NE.0 .AND.
		 * TIMCPU.NE.0.0) ["this should always happen" // WRITE(IOUT,280)
		 * TIMCPU/dble(IHSTRY),3600.*dble(IHSTRY)/TIMCPU; // ] //]
		 */
		// "******************************************************************************
		// "
		// " *** SECTION 3 ***
		// "
		// "------------------------------------------------------------------------------
		//
		// "STATISTICAL AND OTHER DATA HANDLING AND CALL THE OUTPUT SUMMARY ROUTINE"
		//
		// "------------------------------------------------------------------------------

		// :STATS-ANAL:;
		/*
		 * FMASSC=0.0; //"total cavity mass" for(int IX=1;IX<=EGS4Geom.NR;IX++)
		 * { for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ { //$GET-IRL(IZ,IX);
		 * IRL=EGS4Geom.GET_IRL(IZ,IX); if(EGS4Geom.ntrack[IRL-1]==1) {
		 * FMASS=AMASS[IZ-1][IX-1]; FMASSC=FMASSC+FMASS; } } }
		 */
		if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {
			// SCORE_NORM_NUM=EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/
			// EGS4SrcEns.dble(NCASE_PHSP)*NINCSRC;
		} else {
			SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);
		}

		// IF(ISOURC = 21|ISOURC = 22)[
		// "normalize dose to number of incident particles from primary source
		// AINFLU =
		// dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*NINCSRC;
		// "we estimate the total number of particles from the primary source "
		// "(original non-phase space source) by taking the ratio of the total "
		// "number of particles read from the phase source in this simulation "
		// "to the total number of particles in the phase space source and multiply "
		// "this by the number of particles from the primary source that were used "
		// "to obtain this phase space  source."
		// OUTPUT AINFLU, NNREAD, NNREAD-IHSTRY, NCASE_PHSP, NINCSRC;
		// (/' Corresponding number of particles in original BEAM simulation =',
		// F13.0/
		// ' Based on reading', I13, ' particles from the phase space file,'/
		// ' rejecting',I13,' particles (due to wrong charge, missing the'/
		// ' geometry or being multiple passers),'/
		// ' and file having', I13, ' particles in it.'/
		// ' The phase space file was generated by', F13.0, ' initial
		// particles'/);
		// WRITE(IOUT,:tgh1:) AINFLU, NNREAD, NNREAD-IHSTRY, NCASE_PHSP,
		// NINCSRC;
		// :tgh1: FORMAT(/' Corresponding number of particles in original BEAM',
		// ' simulation =', F13.0/
		// ' Based on reading', I13, ' particles from the phase space file,'/
		// ' rejecting',I13,' particles (due to wrong charge, missing the'/
		// ' geometry or being multiple passers),'/
		// ' and phase space file having', I13, ' particles in it'/
		// ' The phase space file was generated by', F13.0, ' initial
		// particles'/);
		// SCORE_NORM_NUM=AINFLU;
		// ]
		// ELSEIF( ISOURC = 23 ) [
		// // *
		// //AINFLU = NHSTRY; SCORE_NORM_NUM=NHSTRY;
		// // * /
		// AINFLU = IHSTRY; SCORE_NORM_NUM=IHSTRY;
		// OUTPUT AINFLU;
		// (/'Source 23: number of particles in BEAM source: ',F13.0/);
		// IF( n_photons > 0 ) [
		// sumE_photons = sumE_photons/n_photons;
		// sumE2_photons = sumE2_photons/n_photons;
		// sumE2_photons = sumE2_photons - sumE_photons*sumE_photons;
		// IF( sumE2_photons > 0 ) sumE2_photons =
		// sqrt(sumE2_photons/n_photons);
		// OUTPUT n_photons,sumE_photons,sumE2_photons;
		// (' Number of photons: ',i10/
		// ' Average energy: ',f10.5,' +/- ',f10.5,/);
		// ] ELSE [ OUTPUT; (' Number of photons: 0'); ]
		// IF( n_electrons > 0 ) [
		// sumE_electrons = sumE_electrons/n_electrons;
		// sumE2_electrons = sumE2_electrons/n_electrons;
		// sumE2_electrons = sumE2_electrons - sumE_electrons*sumE_electrons;
		// IF( sumE2_electrons > 0 ) [
		// sumE2_electrons = sqrt(sumE2_electrons/n_electrons);
		// ]
		// OUTPUT n_electrons,sumE_electrons,sumE2_electrons;
		// (' Number of electrons: ',i10/
		// ' Average energy: ',f10.5,' +/- ',f10.5,/);
		// ] ELSE [ OUTPUT; (' Number of electrons: 0'); ]

		// ]
		// ELSE[
		// "IK: adjust incident fluence
		/*
		 * EGS4SrcEns.AINFLU =
		 * EGS4SrcEns.AINFLU*EGS4SrcEns.dble(IHSTRY)/EGS4SrcEns
		 * .dble(EGS4SrcEns.NCASET); SCORE_NORM_NUM=EGS4SrcEns.dble(IHSTRY);
		 */
		// ]

		// "******************"
		// "history by history"
		// "    EMH March 2002"
		// "******************"

		// REPLACE {$ANALYZE(#,#:#)} WITH {;
		// "Macro to analyze uncertainty:"
		// "{P1}{P2}=scoring array (eg SCDOSE(IDZ,IDX,ITDOSE))"
		// "{P3}=quantity to normalize by (eg incident no. of particles)"
		// "Calculates the uncertainty on {P1}{P2}/{P3}.  The "
		// "uncertainty is stored in {P1}2{P2} and is expressed as a percentage of"
		// "{P1}{P2}/{P3} (max 99.9%).  Note that you must define the REAL*8 variable"
		// "SCORE_TEMP in any subroutine where this macro is used.  This macro"
		// "is only used in the analysis of no. of steps."

		// SCORE_TEMP={P1}{P2}/{P3};
		// {P1}2{P2}={P1}2{P2}/{P3};
		// {P1}2{P2}=({P1}2{P2}-SCORE_TEMP*SCORE_TEMP)/({P3}-1);
		// IF({P1}2{P2}>=0.) {P1}2{P2}= SQRT({P1}2{P2});
		// IF(SCORE_TEMP~=0.)[
		// {P1}2{P2}= MIN({P1}2{P2}/SCORE_TEMP*100.,99.9D00);
		// ]
		// ELSE[
		// {P1}2{P2}=99.9D00;
		// ]

		/*
		 * //$ANALYZE(SCOMEG, :SCORE_NORM_NUM);
		 * SCORE_TEMP=EGS4SrcEns.SCOMEG/SCORE_NORM_NUM;
		 * EGS4SrcEns.SCOMEG2=EGS4SrcEns.SCOMEG2/SCORE_NORM_NUM;
		 * EGS4SrcEns.SCOMEG2=
		 * (EGS4SrcEns.SCOMEG2-SCORE_TEMP*SCORE_TEMP)/(SCORE_NORM_NUM-1);
		 * if(EGS4SrcEns.SCOMEG2>=0.) EGS4SrcEns.SCOMEG2=
		 * Math.sqrt(EGS4SrcEns.SCOMEG2); if(SCORE_TEMP!=0.) {
		 * EGS4SrcEns.SCOMEG2=
		 * Math.min(EGS4SrcEns.SCOMEG2/SCORE_TEMP*100.,99.9); } else {
		 * EGS4SrcEns.SCOMEG2=99.9; }
		 * 
		 * 
		 * EGS4SrcEns.SCOMEG =
		 * EGS4SrcEns.SCOMEG/EGS4SrcEns.dble(IHSTRY);//"Corrected, IK May 4 1999"
		 * 
		 * //OUTPUT SCOMEG,SCOMEG2;(/' OMEG =',1PE12.3,'(',0PF5.1,'%)'/);
		 * EGS4.seqStr=" OMEG ="+EGS4.format(EGS4SrcEns.SCOMEG,12,false)+"("+
		 * EGS4.format(EGS4SrcEns.SCOMEG2,5)+"%)"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 */
		// $ANALYZE(SCSTP, :SCORE_NORM_NUM);
		// $ANALYZE(SCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP = SCSTP / SCORE_NORM_NUM;
		SCSTP2 = SCSTP2 / SCORE_NORM_NUM;
		SCSTP2 = (SCSTP2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (SCSTP2 >= 0.)
			SCSTP2 = Math.sqrt(SCSTP2);
		if (SCORE_TEMP != 0.) {
			SCSTP2 = Math.min(SCSTP2 / SCORE_TEMP * 100., 99.9);
		} else {
			SCSTP2 = 99.9;
		}
		/*
		 * //$ANALYZE(SCCSTP, :SCORE_NORM_NUM);
		 * SCORE_TEMP=SCCSTP/SCORE_NORM_NUM; SCCSTP2=SCCSTP2/SCORE_NORM_NUM;
		 * SCCSTP2= (SCCSTP2-SCORE_TEMP*SCORE_TEMP)/(SCORE_NORM_NUM-1);
		 * if(SCCSTP2>=0.) SCCSTP2= Math.sqrt(SCCSTP2); if(SCORE_TEMP!=0.) {
		 * SCCSTP2= Math.min(SCCSTP2/SCORE_TEMP*100.,99.9); } else {
		 * SCCSTP2=99.9; }
		 */
		// $ANALYZE(SCDSTP, :SCORE_NORM_NUM);
		// $ANALYZE(SCCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP = SCDSTP / SCORE_NORM_NUM;
		SCDSTP2 = SCDSTP2 / SCORE_NORM_NUM;
		SCDSTP2 = (SCDSTP2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (SCDSTP2 >= 0.)
			SCDSTP2 = Math.sqrt(SCDSTP2);
		if (SCORE_TEMP != 0.) {
			SCDSTP2 = Math.min(SCDSTP2 / SCORE_TEMP * 100., 99.9);
		} else {
			SCDSTP2 = 99.9;
		}
		/*
		 * //"first, estimate uncertainties for entire cavity" cav2_dose =
		 * (cav2_dose*SCORE_NORM_NUM - cav_dose*cav_dose)/ (SCORE_NORM_NUM-1);
		 * if( cav2_dose > 0 ) cav2_dose = Math.sqrt(cav2_dose);
		 * 
		 * if( EGS4Macro.IFULL == 1 ) {
		 * 
		 * cav2_dose0 = (cav2_dose0*SCORE_NORM_NUM - cav_dose0*cav_dose0)/
		 * (SCORE_NORM_NUM-1); if( cav2_dose0 > 0 ) cav2_dose0 =
		 * Math.sqrt(cav2_dose0); cav2_dose1 = (cav2_dose1*SCORE_NORM_NUM -
		 * cav_dose1*cav_dose1)/ (SCORE_NORM_NUM-1); if( cav2_dose1 > 0 )
		 * cav2_dose1 = Math.sqrt(cav2_dose1); cav2_dose2 =
		 * (cav2_dose2*SCORE_NORM_NUM - cav_dose2*cav_dose2)/
		 * (SCORE_NORM_NUM-1); if( cav2_dose2 > 0 ) cav2_dose2 =
		 * Math.sqrt(cav2_dose2);
		 * 
		 * corr_02=(cav_dosec02*SCORE_NORM_NUM-cav_dose0*cav_dose2)/
		 * (SCORE_NORM_NUM-1);
		 * corr_02=cav2_dose0*cav2_dose0+cav2_dose2*cav2_dose2+ 2*corr_02; if
		 * (corr_02 > 0) corr_02 = Math.sqrt(corr_02);
		 * 
		 * cav_dosec = (cav_dosec*SCORE_NORM_NUM - cav_dose*cav_dose1)/
		 * (SCORE_NORM_NUM-1); cav_dosec01 = (cav_dosec01*SCORE_NORM_NUM -
		 * cav_dose0*cav_dose1)/ (SCORE_NORM_NUM-1); cav_dosec02 =
		 * (cav_dosec02*SCORE_NORM_NUM - cav_dose0*cav_dose2)/
		 * (SCORE_NORM_NUM-1); }
		 * 
		 * cav_dose = cav_dose*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC); cav_dose0 =
		 * cav_dose0*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC); cav_dose1 =
		 * cav_dose1*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC); cav_dose2 =
		 * cav_dose2*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * 
		 * cav2_dose = cav2_dose*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * cav2_dose0 = cav2_dose0*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * cav2_dose1 = cav2_dose1*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * cav2_dose2 = cav2_dose2*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * 
		 * corr_02 = corr_02*1.602e-10/(EGS4SrcEns.AINFLU*FMASSC);
		 * 
		 * cav_dosec =
		 * cav_dosec*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC))*(1.602e-10
		 * /(EGS4SrcEns.AINFLU*FMASSC)); cav_dosec01 =
		 * cav_dosec01*(1.602e-10/(EGS4SrcEns
		 * .AINFLU*FMASSC))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC)); cav_dosec02
		 * =
		 * cav_dosec02*(1.602e-10/(EGS4SrcEns.AINFLU*FMASSC))*(1.602e-10/(EGS4SrcEns
		 * .AINFLU*FMASSC));
		 * 
		 * if(EGS4Geom.NSUMCV>1)
		 * {//"multiple cavity regions, analyze quantities in each region"
		 * 
		 * //
		 * "FOR ISOURC=4 WE NEED THE DATA FOR CIRCLES, NOT RINGS, SO ADD IT UP"
		 * // "THIS SHOULD ONLY BE USED IF THE CAVITY HAS AN INFINITE DIAMETER"
		 * if((EGS4SrcEns.ISOURC==4)&&(EGS4Geom.NR>1)) { for(int
		 * IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[ { for(int
		 * IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ[ { //$GET-IRL(IZ,IX);
		 * IRL=EGS4Geom.GET_IRL(IZ,IX); if(EGS4Geom.ntrack[IRL-1]==1) { for(int
		 * IT=1;IT<=4;IT++) { //
		 * "IK: this must be a bug. for ix=1 ix-1=0 and the" //
		 * "    array is not defined"
		 * SCDOSE[IZ-1][IX-1][IT-1]=SCDOSE[IZ-1][IX-1][IT-1]+
		 * SCDOSE[IZ-1][IX-2][IT-1];
		 * SCDOSE2[IZ-1][IX-1][IT-1]=SCDOSE2[IZ-1][IX-1][IT-1]+
		 * SCDOSE2[IZ-1][IX-2][IT-1]; if(IT<4) {
		 * SCDOSE_COV[IZ-1][IX-1][IT-1]=SCDOSE_COV[IZ-1][IX-1][IT-1]+
		 * SCDOSE_COV[IZ-1][IX-2][IT-1]; } } } } } }
		 * 
		 * for(int IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[ { for(int
		 * IZ=1;IZ<=EGS4Geom.NZ;IZ++)//DO IZ=1,NZ[ { //$GET-IRL(IZ,IX);
		 * IRL=EGS4Geom.GET_IRL(IZ,IX);
		 * if(EGS4Geom.ntrack[IRL-1]==1)//NTRACK(IRL).EQ.1)[ {
		 * SCDOSE2[IZ-1][IX-1][0]=(SCDOSE2[IZ-1][IX-1][0]*SCORE_NORM_NUM -
		 * SCDOSE[IZ-1][IX-1][0]*SCDOSE[IZ-1][IX-1][0])/(SCORE_NORM_NUM-1); if(
		 * SCDOSE2[IZ-1][IX-1][0] > 0 )
		 * SCDOSE2[IZ-1][IX-1][0]=Math.sqrt(SCDOSE2[IZ-1][IX-1][0]);
		 * if(EGS4Macro.IFULL==1) {
		 * SCDOSE2[IZ-1][IX-1][1]=(SCDOSE2[IZ-1][IX-1][1]*SCORE_NORM_NUM -
		 * SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][1])/(SCORE_NORM_NUM-1); if(
		 * SCDOSE2[IZ-1][IX-1][1] > 0 )
		 * SCDOSE2[IZ-1][IX-1][1]=Math.sqrt(SCDOSE2[IZ-1][IX-1][1]);
		 * SCDOSE2[IZ-1][IX-1][2]=(SCDOSE2[IZ-1][IX-1][2]*SCORE_NORM_NUM -
		 * SCDOSE[IZ-1][IX-1][2]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1); if(
		 * SCDOSE2[IZ-1][IX-1][2] > 0 )
		 * SCDOSE2[IZ-1][IX-1][2]=Math.sqrt(SCDOSE2[IZ-1][IX-1][2]);
		 * SCDOSE2[IZ-1][IX-1][3]=(SCDOSE2[IZ-1][IX-1][3]*SCORE_NORM_NUM -
		 * SCDOSE[IZ-1][IX-1][3]*SCDOSE[IZ-1][IX-1][3])/(SCORE_NORM_NUM-1); if(
		 * SCDOSE2[IZ-1][IX-1][3] > 0 )
		 * SCDOSE2[IZ-1][IX-1][3]=Math.sqrt(SCDOSE2[IZ-1][IX-1][3]);
		 * //"now calculate the covariances"
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][0]= (SCDOSE_COV[IZ-1][IX-1][0]*SCORE_NORM_NUM
		 * - SCDOSE[IZ-1][IX-1][0]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1);
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][1]= (SCDOSE_COV[IZ-1][IX-1][1]*SCORE_NORM_NUM
		 * - SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][2])/(SCORE_NORM_NUM-1);
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][2]= (SCDOSE_COV[IZ-1][IX-1][2]*SCORE_NORM_NUM
		 * - SCDOSE[IZ-1][IX-1][1]*SCDOSE[IZ-1][IX-1][3])/(SCORE_NORM_NUM-1); }
		 * 
		 * //"now normalize quantities and convert to dose"
		 * FMASS=AMASS[IZ-1][IX-1];
		 * 
		 * SCDOSE[IZ-1][IX-1][0]=
		 * SCDOSE[IZ-1][IX-1][0]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE[IZ-1][IX-1][1]=
		 * SCDOSE[IZ-1][IX-1][1]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE[IZ-1][IX-1][2]=
		 * SCDOSE[IZ-1][IX-1][2]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE[IZ-1][IX-1][3]=
		 * SCDOSE[IZ-1][IX-1][3]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * 
		 * SCDOSE2[IZ-1][IX-1][0]=
		 * SCDOSE2[IZ-1][IX-1][0]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE2[IZ-1][IX-1][1]=
		 * SCDOSE2[IZ-1][IX-1][1]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE2[IZ-1][IX-1][2]=
		 * SCDOSE2[IZ-1][IX-1][2]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * SCDOSE2[IZ-1][IX-1][3]=
		 * SCDOSE2[IZ-1][IX-1][3]*1.602e-10/(EGS4SrcEns.AINFLU*FMASS);
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][0]=
		 * SCDOSE_COV[IZ-1][IX-1][0]*(1.602e-10/(EGS4SrcEns
		 * .AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS));
		 * SCDOSE_COV[IZ-1][IX-1][1]=
		 * SCDOSE_COV[IZ-1][IX-1][1]*(1.602e-10/(EGS4SrcEns
		 * .AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS));
		 * SCDOSE_COV[IZ-1][IX-1][2]=
		 * SCDOSE_COV[IZ-1][IX-1][2]*(1.602e-10/(EGS4SrcEns
		 * .AINFLU*FMASS))*(1.602e-10/(EGS4SrcEns.AINFLU*FMASS)); } } } }
		 */
		if (IDAT == 1) {// "add unscored portions to _TMP arrays"
			for (int IT = 1; IT <= ITMAX; IT++) {
				for (int IX = 1; IX <= NRDOSE; IX++) {
					for (int IZ = 1; IZ <= NZDOSE; IZ++) {
						// System.out.println("@@@# "+SCDOSE[IZ-1][IX-1][IT-1]+" iZ "+IZ+" IX "+IX+" IT "+IT);
						SCDOSE[IZ - 1][IX - 1][IT - 1] = SCDOSE[IZ - 1][IX - 1][IT - 1]
								+ SCDOSE_TMP[IZ - 1][IX - 1][IT - 1];
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = SCDOSE2[IZ - 1][IX - 1][IT - 1]
								+ SCDOSE_TMP[IZ - 1][IX - 1][IT - 1]
								* SCDOSE_TMP[IZ - 1][IX - 1][IT - 1];
						if (IKERMA == 1) {
							SCKERMA[IZ - 1][IX - 1][IT - 1] = SCKERMA[IZ - 1][IX - 1][IT - 1]
									+ SCKERMA_TMP[IZ - 1][IX - 1][IT - 1];
							SCKERMA2[IZ - 1][IX - 1][IT - 1] = SCKERMA2[IZ - 1][IX - 1][IT - 1]
									+ SCKERMA_TMP[IZ - 1][IX - 1][IT - 1]
									* SCKERMA_TMP[IZ - 1][IX - 1][IT - 1];
							if (SCDOSE_LAST[IZ - 1][IX - 1][IT - 1] == SCKERMA_LAST[IZ - 1][IX - 1][IT - 1]) {
								SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1]
										+ SCDOSE_TMP[IZ - 1][IX - 1][IT - 1]
										* SCKERMA_TMP[IZ - 1][IX - 1][IT - 1];
							}
						}
					}
				}
			}
			if (EGS4Macro.IFULL == 2) {
				if (PHENER > 0.) {
					// "FIND WHAT BIN WE ARE IN"
					if (SLOTE > 0.0) {
						// "EQUAL ENERGY BINS CASE"
						Double dbll = new Double(PHENER / SLOTE + 1.0);// NOT
																		// ZEROES
						// IB=Math.min(IFIX(PHENER/SLOTE+0.999),$EBIN);
						IB = Math.min(dbll.intValue(), $EBIN);
					} else {
						IB = MAXBIN;
						// UNTIL((IB.EQ.1).OR.(BINTOP(IB-1).LT.PHENER))
						// [IB=IB-1;]
						while (true) {
							if ((IB == 1) || (BINTOP[IB - 2] < PHENER))
								break;
							IB = IB - 1;
						}
					}

					// "ACCUMULATE THE PULSE HEIGHT DISTRIBUTION"
					SCPDST[IB - 1] = SCPDST[IB - 1] + WT1OLD;
					SCPDST2[IB - 1] = SCPDST2[IB - 1] + WT1OLD * WT1OLD;
					// "also add to cumulative distn"
					for (int ICUM = IB; ICUM <= MAXBIN; ICUM++) {
						SCPCUM[ICUM - 1] = SCPCUM[ICUM - 1] + WT1OLD;
						SCPCUM2[ICUM - 1] = SCPCUM2[ICUM - 1] + WT1OLD * WT1OLD;
					}

					if (IWATCH == 3) {
						EGS4.seqStr = " PULSE HEIGHT ENERGY="
								+ EGS4.format(PHENER, 10, true)
								+ " MeV, IN BIN" + EGS4.format(IB, 3)
								+ " WITH WEIGHT" + EGS4.format(1, 10, false);// !!!!!!!!!!!!!!!!!!!!!
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);
					}
					// "NOW SCORE PROBABILITIES FOR COUNTS IN PEAKS"
					for (int IPK = 1; IPK <= 4; IPK++) {
						// "FOR EACH PEAK, F.E., ESCAPES AND 511"
						// if((PHENER>=DFEN[IPK][2]).AND.(PHENER.LE.DFEN(IPK,3)))[
						if ((PHENER >= DFEN[IPK - 1][1])
								&& (PHENER <= DFEN[IPK - 1][2])) {
							// "IT IS IN THE PEAK"
							SCDFEP[IPK - 1] = SCDFEP[IPK - 1] + WT1OLD;
							SCDFEP2[IPK - 1] = SCDFEP2[IPK - 1] + WT1OLD
									* WT1OLD;
							SCDFDIFF[IPK - 1] = SCDFDIFF[IPK - 1] + WT1OLD;
							SCDFDIFF2[IPK - 1] = SCDFDIFF2[IPK - 1] + WT1OLD
									* WT1OLD;
							if (IWATCH == 3) {
								EGS4.seqStr = " 	IT WAS IN ONE OF THE PEAKS,IPK="
										+ EGS4.format(IPK, 3);
								if (EGS4.iprint > 1)
									printSequence(EGS4.seqStr);
								// OUTPUT IPK;(T50,'IT WAS IN ONE OF THE
								// PEAKS,IPK=',I3/);
							}
						}
						// else
						// if((PHENER.GE.DFEN(IPK,1)).AND.(PHENER.LT.DFEN(IPK,2)))[
						else if ((PHENER >= DFEN[IPK - 1][0])
								&& (PHENER < DFEN[IPK - 1][1])) {
							// "IT IS IN THE BKGD"
							SCDFBK[IPK - 1] = SCDFBK[IPK - 1] + WT1OLD;
							SCDFBK2[IPK - 1] = SCDFBK2[IPK - 1] + WT1OLD
									* WT1OLD;
							SCDFDIFF[IPK - 1] = SCDFDIFF[IPK - 1] - WT1OLD;
							SCDFDIFF2[IPK - 1] = SCDFDIFF2[IPK - 1] - WT1OLD
									* WT1OLD;
						}
					}// "END IPK LOOP"
					PHENER = 0.0;
				}
			}
		}

		// "FOR ISOURC=4 WE NEED THE DATA FOR CIRCLES, NOT RINGS, SO ADD IT UP"
		if ((EGS4SrcEns.ISOURC == 4) && (EGS4Geom.NR > 1)) {
			for (int IT = 1; IT <= ITMAX; IT++) {
				for (int IX = 2; IX <= NRDOSE; IX++) {
					for (int IZ = 1; IZ <= NZDOSE; IZ++) {
						SCDOSE[IZ - 1][IX - 1][IT - 1] = SCDOSE[IZ - 1][IX - 1][IT - 1]
								+ SCDOSE[IZ - 1][IX - 1][IT - 1];
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = SCDOSE2[IZ - 1][IX - 1][IT - 1]
								+ SCDOSE2[IZ - 1][IX - 1][IT - 1];
						SCKERMA[IZ - 1][IX - 1][IT - 1] = SCKERMA[IZ - 1][IX - 1][IT - 1]
								+ SCKERMA[IZ - 1][IX - 1][IT - 1];
						SCKERMA2[IZ - 1][IX - 1][IT - 1] = SCKERMA2[IZ - 1][IX - 1][IT - 1]
								+ SCKERMA2[IZ - 1][IX - 1][IT - 1];
					}
				}
			}
		}

		// "AT THIS POINT, SCDOSE CONTAINS THE ENERGY (IN MeV) DEPOSITED"
		// "IN EACH MODE. TO GET THE AVERAGE ENERGY DEPOSITED IN THE       "
		// "PULSE HEIGHT DETECTOR REGION, WE SUM THE IT=1 VALUES FOR ALL         "
		// "REGIONS IN THE DETECTOR. (FOR IFULL=2 ONLY)                          "
		if (EGS4Macro.IFULL == 2) {
			for (int IX = 1; IX <= NRDOSE; IX++) {
				for (int IZ = 1; IZ <= NZDOSE; IZ++) {
					// $GET-IRL(IZ,IX);
					IRL = EGS4Geom.GET_IRL(IZ, IX);
					if (IPHR[IRL - 1] != 0) {
						// "THIS REGION IS IN DETECTOR"
						SCPHEN = SCPHEN + SCDOSE[IZ - 1][IX - 1][0];
						SCPHEN2 = SCPHEN2 + SCDOSE2[IZ - 1][IX - 1][0];
					}// "END TEST FOR INSIDE THE DETECTOR"
				}
			}// "END LOOPS OVER REGIONS"
		}// "END OF IFULL = 2 BLOCK"

		// "STATISTICAL ANALYSES ON THE RAW DATA"

		for (int IT = 1; IT <= ITMAX; IT++) {
			for (int IX = 1; IX <= NRDOSE; IX++) {
				for (int IZ = 1; IZ <= NZDOSE; IZ++) {
					// $ANALYZE(SCDOSE,(IZ,IX,IT):SCORE_NORM_NUM);
					// System.out.println("@@@# "+SCDOSE[IZ-1][IX-1][IT-1]+" scorenum "+SCORE_NORM_NUM);
					SCORE_TEMP = SCDOSE[IZ - 1][IX - 1][IT - 1]
							/ SCORE_NORM_NUM;
					SCDOSE2[IZ - 1][IX - 1][IT - 1] = SCDOSE2[IZ - 1][IX - 1][IT - 1]
							/ SCORE_NORM_NUM;
					SCDOSE2[IZ - 1][IX - 1][IT - 1] = (SCDOSE2[IZ - 1][IX - 1][IT - 1] - SCORE_TEMP
							* SCORE_TEMP)
							/ (SCORE_NORM_NUM - 1);
					if (SCDOSE2[IZ - 1][IX - 1][IT - 1] >= 0.)
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = Math
								.sqrt(SCDOSE2[IZ - 1][IX - 1][IT - 1]);
					if (SCORE_TEMP != 0.) {
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = Math.min(
								SCDOSE2[IZ - 1][IX - 1][IT - 1] / SCORE_TEMP
										* 100., 99.9);
					} else {
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = 99.9;
					}

					if (IKERMA == 1) {
						// $ANALYZE(SCKERMA,(IZ,IX,IT):SCORE_NORM_NUM);
						SCORE_TEMP = SCKERMA[IZ - 1][IX - 1][IT - 1]
								/ SCORE_NORM_NUM;
						SCKERMA2[IZ - 1][IX - 1][IT - 1] = SCKERMA2[IZ - 1][IX - 1][IT - 1]
								/ SCORE_NORM_NUM;
						SCKERMA2[IZ - 1][IX - 1][IT - 1] = (SCKERMA2[IZ - 1][IX - 1][IT - 1] - SCORE_TEMP
								* SCORE_TEMP)
								/ (SCORE_NORM_NUM - 1);
						if (SCKERMA2[IZ - 1][IX - 1][IT - 1] >= 0.)
							SCKERMA2[IZ - 1][IX - 1][IT - 1] = Math
									.sqrt(SCKERMA2[IZ - 1][IX - 1][IT - 1]);
						if (SCORE_TEMP != 0.) {
							SCKERMA2[IZ - 1][IX - 1][IT - 1] = Math.min(
									SCKERMA2[IZ - 1][IX - 1][IT - 1]
											/ SCORE_TEMP * 100., 99.9);
						} else {
							SCKERMA2[IZ - 1][IX - 1][IT - 1] = 99.9;
						}

						// "now analyze the uncertainty on the dose/kerma ratio"
						// "first set SCDOSEtoKERMA2(IZ,IX,IT)=cov(dose,kerma)"
						SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1]
								/ SCORE_NORM_NUM
								- SCDOSE[IZ - 1][IX - 1][IT - 1]
								* SCKERMA[IZ - 1][IX - 1][IT - 1]
								/ (SCORE_NORM_NUM * SCORE_NORM_NUM);
						// "now set SCDOSEtoKERMA2(IZ,IX,IT)=cov(dose,kerma)/"
						// "                                  (dose*kerma)"
						SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1]
								/ (SCDOSE[IZ - 1][IX - 1][IT - 1]
										* SCKERMA[IZ - 1][IX - 1][IT - 1] / (SCORE_NORM_NUM * SCORE_NORM_NUM));
						// "now set SCDOSEtoKERMA2(IZ,IX,IT)=cov(dose,kerma)/"
						// "                                 (dose*kerma)/(N-1)"
						SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1]
								/ (SCORE_NORM_NUM - 1);
						// "now estimate the uncertainty on dose/fluence"
						SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = (SCDOSE2[IZ - 1][IX - 1][IT - 1] / 100.)
								* (SCDOSE2[IZ - 1][IX - 1][IT - 1] / 100.)
								+ (SCKERMA2[IZ - 1][IX - 1][IT - 1] / 100.)
								* (SCKERMA2[IZ - 1][IX - 1][IT - 1] / 100.)
								- 2
								* SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1];
						if (SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] > 0.) {
							SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = 100 * Math
									.sqrt(SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1]);
						}
						if (SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] > 99.9) {
							SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = 99.9;
						}
					}
				}
			}
		}

		// $ANALYZE(SCOMEG, :dble(IHSTRY));
		SCORE_TEMP = EGS4SrcEns.SCOMEG / EGS4SrcEns.dble(IHSTRY);
		EGS4SrcEns.SCOMEG2 = EGS4SrcEns.SCOMEG2 / EGS4SrcEns.dble(IHSTRY);
		EGS4SrcEns.SCOMEG2 = (EGS4SrcEns.SCOMEG2 - SCORE_TEMP * SCORE_TEMP)
				/ (EGS4SrcEns.dble(IHSTRY) - 1);
		if (EGS4SrcEns.SCOMEG2 >= 0.)
			EGS4SrcEns.SCOMEG2 = Math.sqrt(EGS4SrcEns.SCOMEG2);
		if (SCORE_TEMP != 0.) {
			EGS4SrcEns.SCOMEG2 = Math.min(EGS4SrcEns.SCOMEG2 / SCORE_TEMP
					* 100., 99.9);
		} else {
			EGS4SrcEns.SCOMEG2 = 99.9;
		}

		EGS4SrcEns.SCOMEG = EGS4SrcEns.SCOMEG / EGS4SrcEns.dble(IHSTRY);// "Corrected, IK May 4 1999"

		// OUTPUT SCOMEG,SCOMEG2;(/' OMEG =',1PE12.3,'(',0PF5.1,'%)'/);
		EGS4.seqStr = " OMEG =" + EGS4.format(EGS4SrcEns.SCOMEG, 12, false)
				+ "(" + EGS4.format(EGS4SrcEns.SCOMEG2, 5) + "%)";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		// "ANALYSIS OF THE PULSE HEIGHT DISTRIBUTIONS"
		if (EGS4Macro.IFULL == 2) {

			for (IB = 1; IB <= MAXBIN; IB++) {
				SCPTOT = SCPTOT + SCPDST[IB - 1];
				SCPTOT2 = SCPTOT2 + SCPDST2[IB - 1];
			}

			// $ANALYZE(SCPTOT, :SCORE_NORM_NUM);
			SCORE_TEMP = SCPTOT / SCORE_NORM_NUM;
			SCPTOT2 = SCPTOT2 / SCORE_NORM_NUM;
			SCPTOT2 = (SCPTOT2 - SCORE_TEMP * SCORE_TEMP)
					/ (SCORE_NORM_NUM - 1);
			if (SCPTOT2 >= 0.)
				SCPTOT2 = Math.sqrt(SCPTOT2);
			if (SCORE_TEMP != 0.) {
				SCPTOT2 = Math.min(SCPTOT2 / SCORE_TEMP * 100., 99.9);
			} else {
				SCPTOT2 = 99.9;
			}

			// $ANALYZE(SCPHEN, :SCORE_NORM_NUM);
			SCORE_TEMP = SCPHEN / SCORE_NORM_NUM;
			SCPHEN2 = SCPHEN2 / SCORE_NORM_NUM;
			SCPHEN2 = (SCPHEN2 - SCORE_TEMP * SCORE_TEMP)
					/ (SCORE_NORM_NUM - 1);
			if (SCPHEN2 >= 0.)
				SCPHEN2 = Math.sqrt(SCPHEN2);
			if (SCORE_TEMP != 0.) {
				SCPHEN2 = Math.min(SCPHEN2 / SCORE_TEMP * 100., 99.9);
			} else {
				SCPHEN2 = 99.9;
			}

			for (int IPK = 1; IPK <= 4; IPK++) {
				// $ANALYZE(SCDFEP,(IPK):SCORE_NORM_NUM);
				SCORE_TEMP = SCDFEP[IPK - 1] / SCORE_NORM_NUM;
				SCDFEP2[IPK - 1] = SCDFEP2[IPK - 1] / SCORE_NORM_NUM;
				SCDFEP2[IPK - 1] = (SCDFEP2[IPK - 1] - SCORE_TEMP * SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCDFEP2[IPK - 1] >= 0.)
					SCDFEP2[IPK - 1] = Math.sqrt(SCDFEP2[IPK - 1]);
				if (SCORE_TEMP != 0.) {
					SCDFEP2[IPK - 1] = Math.min(SCDFEP2[IPK - 1] / SCORE_TEMP
							* 100., 99.9);
				} else {
					SCDFEP2[IPK - 1] = 99.9;
				}
				// $ANALYZE(SCDFBK,(IPK):SCORE_NORM_NUM);
				SCORE_TEMP = SCDFBK[IPK - 1] / SCORE_NORM_NUM;
				SCDFBK2[IPK - 1] = SCDFBK2[IPK - 1] / SCORE_NORM_NUM;
				SCDFBK2[IPK - 1] = (SCDFBK2[IPK - 1] - SCORE_TEMP * SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCDFBK2[IPK - 1] >= 0.)
					SCDFBK2[IPK - 1] = Math.sqrt(SCDFBK2[IPK - 1]);
				if (SCORE_TEMP != 0.) {
					SCDFBK2[IPK - 1] = Math.min(SCDFBK2[IPK - 1] / SCORE_TEMP
							* 100., 99.9);
				} else {
					SCDFBK2[IPK - 1] = 99.9;
				}

				// "subtract background from peak"
				SCDFEP[IPK - 1] = SCDFEP[IPK - 1] - SCDFBK[IPK - 1];

				// "estimate uncertainty on this subtracted value"
				SCDFEP2[IPK - 1] = (SCDFEP2[IPK - 1] * SCDFEP[IPK - 1]
						/ SCORE_NORM_NUM / 100.)
						* (SCDFEP2[IPK - 1] * SCDFEP[IPK - 1] / SCORE_NORM_NUM / 100.)
						+ (SCDFBK2[IPK - 1] * SCDFBK[IPK - 1] / SCORE_NORM_NUM / 100.)
						* (SCDFBK2[IPK - 1] * SCDFBK[IPK - 1] / SCORE_NORM_NUM / 100.)
						+ 2
						/ (SCORE_NORM_NUM * SCORE_NORM_NUM)
						* (SCDFEP[IPK - 1] * SCDFBK[IPK - 1])
						/ (SCORE_NORM_NUM - 1);
				if (SCDFEP2[IPK - 1] > 0.) {
					SCDFEP2[IPK - 1] = Math.sqrt(SCDFEP2[IPK - 1])
							/ (SCDFEP[IPK - 1] / SCORE_NORM_NUM);
				}
				if (SCDFEP2[IPK - 1] > 0.999) {
					SCDFEP2[IPK - 1] = 0.999;
				}
				// ---------------------------------------------------------
				SCDFEP[IPK - 1] = SCDFEP[IPK - 1] / SCPTOT;
				// ---------------------------------------------------------
				// "now estimate the uncertainty on this quotient"
				SCDFEP2[IPK - 1] = SCDFEP2[IPK - 1]
						* SCDFEP2[IPK - 1]
						+ (SCPTOT2 / 100.)
						* (SCPTOT2 / 100.)
						- 2
						* (SCDFDIFF2[IPK - 1] / SCORE_NORM_NUM - SCDFDIFF[IPK - 1]
								* SCPTOT
								/ ((SCORE_NORM_NUM) * (SCORE_NORM_NUM)))
						/ (SCDFDIFF[IPK - 1] * SCPTOT / ((SCORE_NORM_NUM) * (SCORE_NORM_NUM)))
						/ (SCORE_NORM_NUM - 1);
				if (SCDFEP2[IPK - 1] > 0.)
					SCDFEP2[IPK - 1] = 100 * Math.sqrt(SCDFEP2[IPK - 1]);
				if (SCDFEP2[IPK - 1] > 99.9)
					SCDFEP2[IPK - 1] = 99.9;
			}
			for (int IB = 1; IB <= MAXBIN; IB++) {
				// "save SCPDST2(IB) since it is also equal to SCPDST(IB)*SCPTOT summed"
				// "over all primary histories and will be used later"
				SCORE_TEMP2 = SCPDST2[IB - 1];
				// $ANALYZE(SCPDST,(IB):SCORE_NORM_NUM);
				SCORE_TEMP = SCPDST[IB - 1] / SCORE_NORM_NUM;
				SCPDST2[IB - 1] = SCPDST2[IB - 1] / SCORE_NORM_NUM;
				SCPDST2[IB - 1] = (SCPDST2[IB - 1] - SCORE_TEMP * SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCPDST2[IB - 1] >= 0.)
					SCPDST2[IB - 1] = Math.sqrt(SCPDST2[IB - 1]);
				if (SCORE_TEMP != 0.) {
					SCPDST2[IB - 1] = Math.min(SCPDST2[IB - 1] / SCORE_TEMP
							* 100., 99.9);
				} else {
					SCPDST2[IB - 1] = 99.9;
				}

				// "now estimate the uncertainty on scpdst(IB)/scptot"
				SCPDST2[IB - 1] = (SCPDST2[IB - 1] / 100.)
						* (SCPDST2[IB - 1] / 100.)
						+ (SCPTOT2 / 100.)
						* (SCPTOT2 / 100.)
						- 2
						* (SCORE_TEMP2 / SCORE_NORM_NUM - SCPDST[IB - 1]
								* SCPTOT / (SCORE_NORM_NUM * SCORE_NORM_NUM))
						/ (SCPDST[IB - 1] * SCPTOT / (SCORE_NORM_NUM * SCORE_NORM_NUM))
						/ (SCORE_NORM_NUM - 1);

				if (SCPDST2[IB - 1] > 0.)
					SCPDST2[IB - 1] = 100 * Math.sqrt(SCPDST2[IB - 1]);
				if (SCPDST2[IB - 1] > 99.9)
					SCPDST2[IB - 1] = 99.9;
				// ---------------------------------------------------------
				SCPDST[IB - 1] = SCPDST[IB - 1] / SCPTOT;
				// ---------------------------------------------------------
				// "save SCPCUM2(IB) since it is also equal to SCPCUM(IB)*SCPTOT summed"
				// "over all primary histories and will be used later"
				SCORE_TEMP2 = SCPCUM2[IB - 1];
				// $ANALYZE(SCPCUM,(IB):SCORE_NORM_NUM);
				SCORE_TEMP = SCPCUM[IB - 1] / SCORE_NORM_NUM;
				SCPCUM2[IB - 1] = SCPCUM2[IB - 1] / SCORE_NORM_NUM;
				SCPCUM2[IB - 1] = (SCPCUM2[IB - 1] - SCORE_TEMP * SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCPCUM2[IB - 1] >= 0.)
					SCPCUM2[IB - 1] = Math.sqrt(SCPCUM2[IB - 1]);
				if (SCORE_TEMP != 0.) {
					SCPCUM2[IB - 1] = Math.min(SCPCUM2[IB - 1] / SCORE_TEMP
							* 100., 99.9);
				} else {
					SCPCUM2[IB - 1] = 99.9;
				}

				// "now estimate the uncertainty on this quotient"
				SCPCUM2[IB - 1] = (SCPCUM2[IB - 1] / 100.)
						* (SCPCUM2[IB - 1] / 100.)
						+ (SCPTOT2 / 100.)
						* (SCPTOT2 / 100.)
						- 2
						* (SCORE_TEMP2 / SCORE_NORM_NUM - SCPCUM[IB - 1]
								* SCPTOT / (SCORE_NORM_NUM * SCORE_NORM_NUM))
						/ (SCPCUM[IB - 1] * SCPTOT / (SCORE_NORM_NUM * SCORE_NORM_NUM))
						/ (SCORE_NORM_NUM - 1);

				if (SCPCUM2[IB - 1] > 0.)
					SCPCUM2[IB - 1] = 100 * Math.sqrt(SCPCUM2[IB - 1]);
				if (SCPCUM2[IB - 1] > 99.9)
					SCPCUM2[IB - 1] = 99.9;
				// -------------------------------------------------------
				SCPCUM[IB - 1] = SCPCUM[IB - 1] / SCPTOT;
				// -------------------------------------------------------
			}
			// =============================================================================
			SCPTOT = SCPTOT / SCORE_NORM_NUM; // "NORMALIZE TOTAL TO PER HISTORY"
			SCPHEN = SCPHEN / SCORE_NORM_NUM; // "NORMALIZE TO ENERGY PER HISTORY"
			// =============================================================================

		}// "END IFULL=2 BLOCK"

		// "CONVERT DOSE FROM MeV PER REGION PER BATCH TO GRAY PER UNIT INCIDENT FLUENCE"
		if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {
			// "normalize dose to number of incident particles from primary
			// source
			// AINFLU=dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*NINCSRC;
			// "we estimate the total number of particles from the primary source (original"
			// "non-phase space source) by taking the ratio of the total number of"
			// "particles read from the phase source in this simulation to the total number"
			// "of particles in the phase space source and multiply this by the number of"
			// "particles from the primary source that were used to obtain this phase space"
			// "source."
		}
		// ELSEIF(ISOURC=23)[
		// "have to set AINFLU to no. of primary histories
		// AINFLU=dble(IHSTRY);
		// ]

		// "RECALL 1 MeV = 1.602E-06 erg, 1 rad=100 ergs/g, 1 rad=0.01 Gy"
		// "THE UNIT OF DOSE IS Gy-cm**2"
		for (int IT = 1; IT <= ITMAX; IT++) {
			for (int IX = 1; IX <= NRDOSE; IX++) {
				for (int IZ = 1; IZ <= NZDOSE; IZ++) {
					// System.out.println("before "+SCDOSE[IZ-1][IX-1][IT-1]);
					if (SCDOSE[IZ - 1][IX - 1][IT - 1] != 0.0) {
						// FMASS=AMASS(IZ+NZDMIN-1,IX+NRDMIN);
						FMASS = AMASS[IZ + NZDMIN - 2][IX + NRDMIN - 1];
						if (FMASS == 0.0)
							FMASS = 1.0; // "AVOIDS /0 FOR VACUUM"
						SCDOSE[IZ - 1][IX - 1][IT - 1] = SCDOSE[IZ - 1][IX - 1][IT - 1]
								* 1.602E-10 / (FMASS * EGS4SrcEns.AINFLU);
						if (SCKERMA[IZ - 1][IX - 1][IT - 1] != 0) {
							SCKERMA[IZ - 1][IX - 1][IT - 1] = SCKERMA[IZ - 1][IX - 1][IT - 1]
									* 1.602E-10 / (FMASS * EGS4SrcEns.AINFLU);
						}
					}
					// System.out.println("after "+SCDOSE[IZ-1][IX-1][IT-1]+" fmass "+FMASS+" AInflu "+EGS4SrcEns.AINFLU);
				}
			}
		}

		OSUMRY(); // "PRINT THE OUTPUT SUMRY"

		// :END-OF-RUN:;

		// ;"******************************************************************************
		// "
		// " *** SECTION 4 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "THE CONCLUSION"
		// "
		// "------------------------------------------------------------------------------

		// :END:;
		// OUTPUT; (//'END OF RUN',9X,' ',$); call egs_fdate(6);
		// OUTPUT; (//);
		// write(iout,'(/a,$)') 'END OF RUN '; call egs_fdate(iout);
		// write(iout,'(////)');

		// call egs_finish;

		// #ifdef HAVE_C_COMPILER;
		// ;
		// IF( n_parallel > 0 & ~is_finished ) [
		// call egs_pjob_finish(n_job);
		// / IF( n_job = 0 ) [
		// is_finished = .true.;
		// call egs_combine_runs(combine_results,'.egsdat');
		// NCASET=NCASEO; IHSTRY=NCASET;
		// CALL SRCOTO(WEIGHT);
		// goto :STATS-ANAL:;
		// ]
		// ]
		// #endif;

		EGS4SrcEns.SRCEND();// do nothing

		// $CALL_EXIT(0);

		// "FORTRAN FORMAT STATEMENTS. FORMAT STATEMENT N## IS FIRST USED IN SECTION N."
		// %I0
		// 100 FORMAT(80A1//'Calculation using CAVRZnrc(EGSnrc) '$VERSION' ',
		// /' ON '$MACHINE' ',T55,' ',$);
		// 200 FORMAT(//,79('*')/
		// // ,T20,'EXECUTION INFORMATION AND WARNING MESSAGES'/
		// // ,79('*')/
		// //'USING CAVRZnrc(EGSnrc) '$VERSION' ON '$MACHINE' ');
		// 201 FORMAT(/'********* NEW INPUT FILE *********'/);
		// 202 FORMAT(/'********* RESTARTED INPUT FILE ********* '/
		// ' ',10X,I12,' NEW + ',I12,' OLD HISTORIES');
		// 204 FORMAT(/'********* DATA ANALYSIS ONLY *********'/);
		// 205 FORMAT(/'********* RANDOM NUMBERS READ FROM FILE *********'/);
		// 206 FORMAT(/' ********* ANALYZING RESULTS FROM PARALLEL RUNS
		// *******'/);
		// 210 FORMAT(/'********* NOT ENOUGH TIME TO FINISH WITHIN',
		// ' LIMIT OF',F8.2,' HOURS',I5,' BATCHES USED********'/
		// ' ',I12,' HISTORIES RUN, ',I12,' HISTORIES ANALYZED'//);
		// 230 FORMAT(/'DESIRED STATISTICAL ACCURACY OBTAINED.'/
		// ' STATS IN CAVITY= ',F5.2,'%',
		// ' AFTER ',I2,' BATCHES');
		// 240 FORMAT(/'*********DESIRED STATISTICAL ACCURACY OF ',F5.2,'%',
		// ' NOT REACHED*********'/
		// ' STATS IN CAVITY= ',F5.2,' % AFTER ',I2,' BATCHES');
		// 250 FORMAT(/' FOR OLD RUN:'/
		// ' ----------- '/
		// ' Total cputime =',F8.1,'s (=',F5.2,' hr)');
		// 255 FORMAT(/' FOR PARALLEL RUNS:'/
		// ' ----------------- '/
		// ' On ',I5,' machines '/
		// ' Total cputime =',F8.1,'s (=',F8.2,' hr), cputime/machine =',
		// F8.1,'s');
		// 260 FORMAT(/'Finished simulations:'/' time elapsed,cputime',
		// ',ratio= ',2F8.1,'(=',F5.2,'hr)',F8.2);
		// 261 FORMAT(/' Finished: time elapsed this run', F10.1/
		// ' CPUtime total run ', F10.1,'(=',F8.2,'hr)'/
		// ' Ratio ELAPSED/CPU this run:', F8.3);
		// 280 FORMAT(/' CPUtime/history=',F10.5,' sec. Histories/hour=',F12.0);
		//
		// END; "END OF MAIN ROUTINE-CAVRZnrc"

		EGS4.timeElapsed();
		Calendar call = Calendar.getInstance();
		Date da = call.getTime();
		EGS4.seqStr = "End of run:          " + da.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// ===============================================================
		if (createOutputFile) {
			try {
				sigfos.close();
			} catch (Exception ex) {
			}

			putInFile = false;
			EGS4.seqStr = "Check current directory for results, in text file:"
					+ filename;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
	}

	/**
	 * Where to print runtime information. Interface method.
	 * @param s the String to be printed
	 */
	public void printSequence(String s) {
		// write file?
		if (createOutputFile && putInFile) {
			try {
				sigfos.write(s + " \n");
			} catch (Exception ex) {
			}
		}

		// output data to console
		System.out.println(s);
	}

	// "******************************************************************************
	// "
	// "
	// " **********
	// " * *
	// " * AUSGAB *
	// " * *
	// " **********
	// "
	// "
	// " AN AUSGAB ROUTINE TO BE USED WITH DOSRZnrc.mortran
	// "
	// " Called with IARG = -1 after each history is over in order to
	// " score things for the pulse height distribution when IFULL = 2.
	// "
	// " For IFULL = 2 (passed in COMIN SCORE), the parameter phener keeps
	// " track of the energy deposited in the sensitive volume defined by non-
	// " zero elements in the array IPHR($MXREG) (passed IN SCORE).
	// "
	// " This routine scores the dose in a finite, azimuthally symmetric
	// " cylindrical geometry which the user defines via plane and radial
	// " coordinates. The user must specify both the target geometry as well
	// " as the planes and radii between which the dose is to be scored. All
	// " the geometrical checks for crossing 'geometrical' or 'dose' regions
	// " are handled by the subroutine HOWFAR.
	// ;
	// " FOR IT = 1 total dose is scored
	// " = 2 dose less stopped/discarded particles is scored
	// " For IFULL = 3 the scattered dose for incident photons
	// " = 3 dose due to particles entering the dose region from
	// " the front wall
	// " = 4 dose due to particles entering the dose region from
	// " the side wall
	// " = 5 dose due to particles entering the dose region from
	// " the back wall
	// " = 6 dose due to particles entering the dose region from
	// " the inside wall
	// " = 7 dose due to particles originates from within an
	// " isotropically radiating disk buried in the geometry
	// " which have not yet strayed outside the source region
	// "
	// " Some of the logic
	// " Bit 6 of latch is set for all photons scattered after Compton
	// " Bit 7 of LATCH is set for all photons after photoeffect
	// " Bit 8 of LATCH is set when
	// "
	// "
	// "******************************************************************************

	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG) {

		// $IMPLICIT-NONE;

		double FTMP = 0.;
		int ip = 0;

		// $INTEGER IRL,IX,IZ,IQL,LATCHL,IARG,IDUMMY;
		// $REAL WTL,FDUMMY,xsi;
		// IRL,IZD,IXD,IQL,ip,IX,IZ,IGEOM,IB,IPK,I,II,ICUM
		// ;COMIN/
		// ELECIN,EPCONT,GEOM,MEDIA,PHOTIN,RUSROU,SCORE,SOURCE,STACK,USEFUL,USER,
		// RANDOM,BOUNDS/;
		// IGBUG1=0;IGBUG2=0;
		double xsi = 0.0;
		double R1 = 0.0;
		double aux1 = 0.0;
		int IRL = 0;
		int IZD = 0;
		int IXD = 0;
		int IQL = 0;
		int IGEOM = 0;
		int IX = 0;
		//int IZ = 0;
		int IB = 0;
		// "STACK OVERFLOW CHECK"
		MXNP = Math.max(MXNP, EGS4.NP);// "keep track of how deep stack is"
		// "MXNP is not output but it should be"
		/*
		 * if(IWATCH > 0) {EGS4.WATCH(IARG,IWATCH);}
		 * //"signal watch routine if active"]
		 * 
		 * //"check if particle is leaving the transport geometry" int
		 * IRL=EGS4.IR[EGS4.NP-1]; //"local region number" if(IRL == 1) return;
		 * //"outside the chamber"
		 * 
		 * //"OBTAIN FREQUENTLY USED LOCAL VARIABLES" //$GET-IX-IZ(IRL);
		 * //"local plane and radius numbers" int IX=EGS4Geom.GET_IX(IRL); int
		 * IZ=EGS4Geom.GET_IZC(IRL);
		 * 
		 * int IQL=EGS4.IQ[EGS4.NP-1]; //"local charge variable" double
		 * WTL=EGS4.WT[EGS4.NP-1]; //"local weight variable" int
		 * LATCHL=EGS4.LATCH[EGS4.NP-1]; //"LATCHL=0 for primaries, 1 otherwise"
		 * //######################################################### if(
		 * EGS4Macro.use_enhance || EGS4Macro.n_split > 1 ) { if( IARG < 5 ) {
		 * if( EGS4.EDEP > 0 && WTL > 0 && EGS4Geom.ntrack[IRL-1] == 1 ) {
		 * //"If we use cross section enhancement or photon splitting,"
		 * //"all energy scoring is done here and the rest of ausgab is ignored"
		 * //"We use the technique proposed by the PENELOPE group for scoring "
		 * //"the energy deposition. This results in a much better estimate   "
		 * //"of the uncertainty                                              "
		 * 
		 * FTMP = WTL*EGS4.EDEP; if( EGS4SrcEns.NHSTRY == last_case ) {
		 * //" Still the same history scoring into the cavity => update    "
		 * //" temporary variables                                         " if(
		 * LATCHL != 2 ) tmp_dose = tmp_dose + FTMP; if( LATCHL != 3 ) tmp_dose1
		 * = tmp_dose1 + FTMP; } else {
		 * //" A new history scoring into the cavity. " last_case =
		 * EGS4SrcEns.NHSTRY; cav_dose = cav_dose + tmp_dose; cav2_dose =
		 * cav2_dose + tmp_dose*tmp_dose; cav_dose1 = cav_dose1 + tmp_dose1;
		 * cav2_dose1 = cav2_dose1 + tmp_dose1*tmp_dose1; cav_dosec = cav_dosec
		 * + tmp_dose*tmp_dose1; if( LATCHL != 2 ) { tmp_dose = FTMP; } else {
		 * tmp_dose = 0.; } if( LATCHL != 3 ) { tmp_dose1 = FTMP;} else {
		 * tmp_dose1 = 0.; } } } return; } }
		 * //######################################################### if(
		 * EGS4Macro.use_enhance ) {//
		 * "If we use cross section enhancement, all scoring " //
		 * " is done here and the rest of ausgab is ignored  "
		 * 
		 * if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
		 * //"A pair/Compton/photoelectric/Rayleigh event is about to take place"
		 * /
		 * /"As we have increased the photon cross section by a factor of      "
		 * /
		 * /"cs_enhance, we must split the photon into a scattering portion    "
		 * /
		 * /"(1/cs_enhance) and a nor-scattering portion (1-1/cs_enhance)      "
		 * /
		 * /"Start with placing an identical photon on the stack               "
		 * EGS4.NP = EGS4.NP + 1; if(EGS4.NP + 1 > EGS4.$MXSTACK) {
		 * EGS4.STOPPROGRAM=true;
		 * EGS4.seqStr=" ***************************************************"
		 * +"  \n"+
		 * " Calculation with CS-enhancement: unable to boost stack."+"  \n"+
		 * " ***************************************************";
		 * if(EGS4.iprint>0) printSequence(EGS4.seqStr); return; //OUTPUT; //( '
		 * Calculation with CS-enhancement: unable to boost stack.'/ // '
		 * Stopping.'/ 1x,80('*')/); //stop; } //$TRANSFER PROPERTIES TO (np)
		 * FROM (np - 1); EGS4.X[EGS4.NP-1]=EGS4.X[EGS4.NP-2];
		 * EGS4.Y[EGS4.NP-1]=EGS4.Y[EGS4.NP-2];
		 * EGS4.Z[EGS4.NP-1]=EGS4.Z[EGS4.NP-2];
		 * EGS4.IR[EGS4.NP-1]=EGS4.IR[EGS4.NP-2];
		 * EGS4.WT[EGS4.NP-1]=EGS4.WT[EGS4.NP-2];
		 * EGS4.DNEAR[EGS4.NP-1]=EGS4.DNEAR[EGS4.NP-2];
		 * EGS4.LATCH[EGS4.NP-1]=EGS4.LATCH[EGS4.NP-2];
		 * 
		 * EGS4.E[EGS4.NP-1] = EGS4.E[EGS4.NP-2]; EGS4.U[EGS4.NP-1] =
		 * EGS4.U[EGS4.NP-2]; EGS4.V[EGS4.NP-1] = EGS4.V[EGS4.NP-2];
		 * EGS4.W[EGS4.NP-1] = EGS4.W[EGS4.NP-2]; EGS4.IQ[EGS4.NP-1] =
		 * EGS4.IQ[EGS4.NP-2]; if( EGS4.LATCH[EGS4.NP-2] != 2 ) {
		 * //" This is either a primary photon that has not yet been attenuated "
		 * /
		 * /" away or a scattered photon. Let's decide what to do with the     "
		 * /
		 * /" unscattered fraction of that photon (which is at np-1)           "
		 * xsi=EGS4.random01(); if( EGS4Macro.cs_enhance*xsi < 1. ) {//
		 * " The photon doesn't survive. " if( EGS4.LATCH[EGS4.NP-2] == 3 ) {//
		 * " It is a scattered photon => kill it" EGS4.WT[EGS4.NP-2] = 0.0;
		 * EGS4.DNEAR[EGS4.NP-2] = -1.; } else {//
		 * " This is a primary => mark it as attenuated.         " //
		 * " From now on, all descendents of this photon will    " //
		 * " only contribute to the cavity dose with attenuation " //
		 * " and scatter removed                                 "
		 * EGS4.LATCH[EGS4.NP-2] = 2; } } }
		 * //"Adjust the weight of to be scattered photon" EGS4.WT[EGS4.NP-1] =
		 * EGS4.WT[EGS4.NP-1]/EGS4Macro.cs_enhance; return; }
		 * 
		 * if( IARG == 18 || IARG == 20 || IARG == 24 ||
		 * //" A Compton/photo-absorption/Rayleigh event just occured" IARG == 7
		 * || IARG == 13 || IARG == 14 ) {
		 * //" A bremas/annihilation event just occured"
		 * //" All scattered photons and photons originating in brems/annihilation"
		 * /
		 * /" events contribute to the scattered dose. But because all of them   "
		 * /
		 * /" have a small weight (initial weight/cs_enhance), we will play      "
		 * /
		 * /" Russian Roulette with them, using 1/cs_enhance as a sirvivng       "
		 * /
		 * /" probability. If they survive, their weight will become equal to the"
		 * /
		 * /" intial photon weight. In addition, we have to set their latch to 3 "
		 * /
		 * /" so that they and rheir descendents only contribute to the scattered"
		 * /
		 * /" dose.                                                              "
		 * xsi=EGS4.random01(); xsi = xsi*EGS4Macro.cs_enhance; for(
		 * ip=EGS4.NPold;ip<=EGS4.NP;ip++) { if( EGS4.IQ[ip-1] == 0 ) { if(
		 * EGS4.LATCH[ip-1] == 2 ) {// "that's a descendent of a photon that" //
		 * "has been attenuated away => kill it" EGS4.WT[ip-1] = 0.0;
		 * EGS4.DNEAR[ip-1] = -1.; } else { if( EGS4.E[ip-1] >= EGS4.PCUT[IRL-1]
		 * ) { if( xsi < 1. ) { EGS4.LATCH[ip-1] = 3; EGS4.WT[ip-1] =
		 * EGS4.WT[ip-1]*EGS4Macro.cs_enhance; } else { EGS4.WT[ip-1] = 0.;
		 * EGS4.DNEAR[ip-1] = -1.; } } else { EGS4.LATCH[ip-1] = 3; }
		 * //" i.e. we don't need the Russian Roulette for photons below"
		 * //" threshold because they will be discarded and their energy"
		 * //" deposited locally anyway                                 " } } }
		 * return; } return; }
		 * //#########################################################
		 * 
		 * if( EGS4Macro.n_split > 1 ) {
		 * 
		 * if( IARG == 7 || IARG == 13 || IARG == 14 ) { if( EGS4Macro.iifano ==
		 * 1 ) { for( ip=EGS4.NPold;ip<=EGS4.NP;ip++) { if( EGS4.IQ[ip-1] == 0 )
		 * { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; } } return; } for(
		 * ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [ { if( EGS4.IQ[ip-1]
		 * == 0 ) { xsi=EGS4.random01(); if( xsi*EGS4Macro.n_split > 1) {
		 * EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; } else { EGS4.LATCH[ip-1] =
		 * 3;EGS4.WT[ip-1]=EGS4.WT[ip-1]*EGS4Macro.n_split; } } } return; } }
		 * //#########################################################
		 * 
		 * 
		 * //REPLACE {$SCORE(#,#:#)} WITH {;
		 * //"Scoring macro used in AUSGAB for quantities other than DOSE and KERMA"
		 * //"{P1}{P2}=scoring array (eg SCSTP)"
		 * //"{P3}=quantity to be scored (eg 1)"
		 * 
		 * //"If the (primary) history number, NHSTRY, is the same as the history"
		 * //"that last scored in this array, {P1}_LAST{P2}, then {P3} is added"
		 * //"to a temporary array, {P1}_TMP{P2}.  Otherwise, we add"
		 * //"{P1}_TMP{P2} to {P1}{P2}, {P1}_TMP{P2}*{P1}_TMP{P2} to {P1}2{P2},"
		 * //"set {P1}_TMP{P2}={P3}, and set {P1}_LAST{P2}=NHSTRY."
		 * //"This scoring method allows us to calculate  uncorrelated value"
		 * //"of {P1}2{P2} which is then used to calculate the uncertainty"
		 * //"in {P1}{P2}.  In cavrznrc, this macro is only used for counting"
		 * //"no. of charged particle steps."
		 * 
		 * //IF(NHSTRY={P1}_LAST{P2})[ //{P1}_TMP{P2}={P1}_TMP{P2} + {P3}; //]
		 * //ELSE[ //{P1}{P2}={P1}{P2}+{P1}_TMP{P2}; //{P1}2{P2}={P1}2{P2} +
		 * {P1}_TMP{P2}*{P1}_TMP{P2}; //{P1}_TMP{P2}={P3};
		 * //{P1}_LAST{P2}=NHSTRY; //] //; //}
		 * 
		 * 
		 * if(IARG == 0) {//"ABOUT TO TRANSPORT A PARTICLE" if(IQL!=0) {
		 * if(LATCHL == 0) {// "COUNT PRIMARY CHARGED PARTICLES ONLY"
		 * //$SCORE(SCSTP, :1);"COUNT CHARGED PARTICLE STEPS TAKEN"
		 * if(EGS4SrcEns.NHSTRY==SCSTP_LAST)//IF(NHSTRY={P1}_LAST{P2})[ { //
		 * {P1}_TMP{P2}={P1}_TMP{P2} + {P3}; SCSTP_TMP=SCSTP_TMP + 1.; } else {
		 * // {P1}{P2}={P1}{P2}+{P1}_TMP{P2}; SCSTP=SCSTP+SCSTP_TMP; //
		 * {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2}; SCSTP2=SCSTP2 +
		 * SCSTP_TMP*SCSTP_TMP; // {P1}_TMP{P2}={P3}; SCSTP_TMP=1.; //
		 * {P1}_LAST{P2}=NHSTRY; SCSTP_LAST=EGS4SrcEns.NHSTRY; }
		 * 
		 * if(EGS4Geom.ntrack[IRL-1] == 1) { //$SCORE(SCCSTP, :1);
		 * if(EGS4SrcEns.NHSTRY==SCCSTP_LAST)//IF(NHSTRY={P1}_LAST{P2})[ { //
		 * {P1}_TMP{P2}={P1}_TMP{P2} + {P3}; SCCSTP_TMP=SCCSTP_TMP + 1.; } else
		 * { // {P1}{P2}={P1}{P2}+{P1}_TMP{P2}; SCCSTP=SCCSTP+SCCSTP_TMP; //
		 * {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2}; SCCSTP2=SCCSTP2 +
		 * SCCSTP_TMP*SCCSTP_TMP; // {P1}_TMP{P2}={P3}; SCCSTP_TMP=1.; //
		 * {P1}_LAST{P2}=NHSTRY; SCCSTP_LAST=EGS4SrcEns.NHSTRY; }
		 * 
		 * //"WRITE(*,*)' sccstp ',SCCSTP;" } }//"COUNT STEPS IN CAVITY REGION"
		 * } else {//"PHOTON STEP - PLAY RUSSIAN ROULETTE?"
		 * if(RUSROU&&(EGS4.W[EGS4.NP-1]>0.0)) {//"YES, PLAY IF CROSSES RRZ "
		 * if((EGS4.Z[EGS4.NP-1]<=RRZ)&&
		 * (EGS4.Z[EGS4.NP-1]+EGS4.USTEP*EGS4.W[EGS4.NP-1]>=RRZ)) {//"CROSSES"
		 * xsi=EGS4.random01(); if(xsi<RRCUT) {//"PARTICLE SURVIVES"
		 * EGS4.WT[EGS4.NP-1]=WTL/RRCUT; } else
		 * {//"DISCARD PARTICLE ON NEXT CALL TO HOWFAR" EGS4.WT[EGS4.NP-1]=0.0;
		 * } } //"END TEST IF CROSSES RUSSIAN ROULETTE PLANE" }
		 * //"END TEST FOR PLAYING RUSSIAN ROULETTE"
		 * }//"END TEST FOR PHOTON STEP" }// "END TEST FOR IARG = 0"
		 * 
		 * if (IFANO == 1) { if (IARG == 15 || IARG == 17 || IARG == 19 || IARG
		 * == 23) {
		 * //"A pair/Compton/photoelectric/Rayleigh event is about to take place"
		 * EGS4.NP = EGS4.NP + 1; //"Boost the stack" if(EGS4.NP + 1 >
		 * EGS4.$MXSTACK) { EGS4.STOPPROGRAM=true;
		 * EGS4.seqStr=" ***************************************************"
		 * +"  \n"+ " Fano calculation unable to boost stack."+"  \n"+
		 * " ***************************************************";
		 * if(EGS4.iprint>0) printSequence(EGS4.seqStr); return; }
		 * //"Create an identical photon" //$TRANSFER PROPERTIES TO (np) FROM
		 * (np - 1); EGS4.X[EGS4.NP-1]=EGS4.X[EGS4.NP-2];
		 * EGS4.Y[EGS4.NP-1]=EGS4.Y[EGS4.NP-2];
		 * EGS4.Z[EGS4.NP-1]=EGS4.Z[EGS4.NP-2];
		 * EGS4.IR[EGS4.NP-1]=EGS4.IR[EGS4.NP-2];
		 * EGS4.WT[EGS4.NP-1]=EGS4.WT[EGS4.NP-2];
		 * EGS4.DNEAR[EGS4.NP-1]=EGS4.DNEAR[EGS4.NP-2];
		 * EGS4.LATCH[EGS4.NP-1]=EGS4.LATCH[EGS4.NP-2];
		 * 
		 * EGS4.E[EGS4.NP-1] = EGS4.E[EGS4.NP-2]; EGS4.U[EGS4.NP-1] =
		 * EGS4.U[EGS4.NP-2]; EGS4.V[EGS4.NP-1] = EGS4.V[EGS4.NP-2];
		 * EGS4.W[EGS4.NP-1] = EGS4.W[EGS4.NP-2]; EGS4.IQ[EGS4.NP-1] =
		 * EGS4.IQ[EGS4.NP-2];
		 * 
		 * return; }
		 * 
		 * //"Throw away any scattered photons from the primary interaction site."
		 * //" Now there is a stack pointer NPold which points to the particle "
		 * //" befor the last discrete interaction. This change was necessary "
		 * /
		 * /" for the implementation of atomic relaxations. So, check all particles"
		 * //" between NPold and NP and discard photons "
		 * //"IF ( (iarg = 18 & NP > NPold)" " Compton has occured" if ( IARG ==
		 * 18// " Compton has occured" || IARG == 20 //
		 * " After photo-absorption " || IARG == 24 // " After Rayleigh " ||
		 * IARG == 7 // " After brems " || IARG == 13 // " After annihilation "
		 * || IARG == 14)// " After annihilation at rest " { for(
		 * ip=EGS4.NPold;ip<=EGS4.NP;ip++) { if( EGS4.IQ[ip-1] == 0 ) {
		 * EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; } } } } if (IFANO == 2) { if (
		 * IARG == 16 // " After pair production " || IARG == 18 //
		 * " After Compton " || IARG == 20) // " After photo absorption " { if (
		 * EGS4Geom.ntrack[EGS4.IR[EGS4.NP-1]-1] == 0) { for(
		 * ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP { if( EGS4.IQ[ip-1]
		 * == 0 )//if( iq(ip) ~= 0 ) { EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.;
		 * }//{ wt(ip) = 0; e(ip) = 0; } } } } }
		 * 
		 * //"SCORE THE ENERGY AS REQUIRED FOR THE DIFFERENT MODES"
		 * //"****************************************************"
		 * 
		 * if((EGS4.EDEP!=0.0)&&
		 * (WTL>0.0)&&(IARG<5)&&(EGS4Geom.ntrack[IRL-1]==1)) {
		 * //"ENERGY HAS BEEN DEPOSITED IN THE CAVITY REGION"
		 * //"SCORE PRIMARY AND SECONDARY ENERGY DEPOSITED"
		 * FTMP=EGS4.WT[EGS4.NP-1]*EGS4.EDEP;
		 * 
		 * //
		 * "***************************************************************************"
		 * //
		 * "                                                                           "
		 * //
		 * " Implementation of a history by history scoring scheme for the cavity dose "
		 * //
		 * " taken from Iwan's splitting implementation, EMH, March 2002               "
		 * //
		 * "                                                                           "
		 * //
		 * " For the meaning of the variables see the COMIN/SCORE block                "
		 * //
		 * "                                                                           "
		 * //
		 * "         expmfp = exp(di) - 1 (see macro $SELECT-PHOTON-MFP)               "
		 * //
		 * "                                                                           "
		 * //
		 * "***************************************************************************"
		 * if( EGS4SrcEns.NHSTRY == last_case ) {
		 * //" Still the same history scoring into the cavity => update    "
		 * //" temporary variables                                         "
		 * tmp_dose = tmp_dose + FTMP; if( LATCHL == 0 ) { tmp_dose0 = tmp_dose0
		 * + FTMP; tmp_dose1 = tmp_dose1 + FTMP*(EGS4Macro.EXPMFP+1); } else {
		 * tmp_dose2 = tmp_dose2 + FTMP; } } else {
		 * //" A new history scoring into the cavity. " last_case =
		 * EGS4SrcEns.NHSTRY;
		 * 
		 * cav_dose = cav_dose + tmp_dose; cav2_dose = cav2_dose +
		 * tmp_dose*tmp_dose;
		 * 
		 * cav_dose0 = cav_dose0 + tmp_dose0; cav2_dose0 = cav2_dose0 +
		 * tmp_dose0*tmp_dose0;
		 * 
		 * cav_dose1 = cav_dose1 + tmp_dose1; cav2_dose1 = cav2_dose1 +
		 * tmp_dose1*tmp_dose1;
		 * 
		 * cav_dose2 = cav_dose2 + tmp_dose2; cav2_dose2 = cav2_dose2 +
		 * tmp_dose2*tmp_dose2;
		 * 
		 * cav_dosec = cav_dosec + tmp_dose*tmp_dose1; cav_dosec01 = cav_dosec01
		 * + tmp_dose0*tmp_dose1; cav_dosec02 = cav_dosec02 +
		 * tmp_dose0*tmp_dose2;
		 * 
		 * tmp_dose = FTMP;
		 * 
		 * if( LATCHL == 0 ) { tmp_dose0 = FTMP ; tmp_dose1 =
		 * FTMP*(EGS4Macro.EXPMFP+1) ; tmp_dose2 = 0.0; } else { tmp_dose0 =
		 * 0.0; tmp_dose1 = 0.0; tmp_dose2 = FTMP; } }
		 * 
		 * IDECAV=1;
		 * 
		 * if(EGS4Geom.NSUMCV>1)
		 * {//"calculate quantities in individual cavity regions" //
		 * "do it the same as for the overall cavity above"
		 * 
		 * if(EGS4SrcEns.NHSTRY==SCDOSE_LAST[IZ-1][IX-1]) {
		 * SCDOSE_TMP[IZ-1][IX-1][0]=SCDOSE_TMP[IZ-1][IX-1][0]+FTMP;
		 * if(LATCHL==0) {//"primary dose"
		 * SCDOSE_TMP[IZ-1][IX-1][1]=SCDOSE_TMP[IZ-1][IX-1][1]+FTMP;
		 * SCDOSE_TMP[IZ
		 * -1][IX-1][2]=SCDOSE_TMP[IZ-1][IX-1][2]+FTMP*(1+EGS4Macro.EXPMFP); }
		 * else {//"secondary dose"
		 * SCDOSE_TMP[IZ-1][IX-1][3]=SCDOSE_TMP[IZ-1][IX-1][3]+FTMP; } } else {
		 * SCDOSE_LAST[IZ-1][IX-1]=EGS4SrcEns.NHSTRY;
		 * 
		 * SCDOSE[IZ-1][IX-1][0]=SCDOSE[IZ-1][IX-1][0]+SCDOSE_TMP[IZ-1][IX-1][0];
		 * SCDOSE2[IZ-1][IX-1][0]=SCDOSE2[IZ-1][IX-1][0]+
		 * SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][0];
		 * 
		 * SCDOSE[IZ-1][IX-1][1]=SCDOSE[IZ-1][IX-1][1]+SCDOSE_TMP[IZ-1][IX-1][1];
		 * SCDOSE2[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]+
		 * SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][1];
		 * 
		 * SCDOSE[IZ-1][IX-1][2]=SCDOSE[IZ-1][IX-1][2]+SCDOSE_TMP[IZ-1][IX-1][2];
		 * SCDOSE2[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX-1][2]+
		 * SCDOSE_TMP[IZ-1][IX-1][2]*SCDOSE_TMP[IZ-1][IX-1][2];
		 * 
		 * SCDOSE[IZ-1][IX-1][3]=SCDOSE[IZ-1][IX-1][3]+SCDOSE_TMP[IZ-1][IX-1][3];
		 * SCDOSE2[IZ-1][IX-1][3]=SCDOSE2[IZ-1][IX-1][3]+
		 * SCDOSE_TMP[IZ-1][IX-1][3]*SCDOSE_TMP[IZ-1][IX-1][3];
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]+
		 * SCDOSE_TMP[IZ-1][IX-1][0]*SCDOSE_TMP[IZ-1][IX-1][2];
		 * SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]+
		 * SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][2];
		 * SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]+
		 * SCDOSE_TMP[IZ-1][IX-1][1]*SCDOSE_TMP[IZ-1][IX-1][3];
		 * 
		 * SCDOSE_TMP[IZ-1][IX-1][0]=FTMP;
		 * 
		 * if(LATCHL==0) { SCDOSE_TMP[IZ-1][IX-1][1]=FTMP;
		 * SCDOSE_TMP[IZ-1][IX-1][2]=FTMP*(1+EGS4Macro.EXPMFP);
		 * SCDOSE_TMP[IZ-1][IX-1][3]= 0.0; } else {
		 * SCDOSE_TMP[IZ-1][IX-1][1]=0.0; SCDOSE_TMP[IZ-1][IX-1][2]=0.0;
		 * SCDOSE_TMP[IZ-1][IX-1][3]=FTMP; } }
		 * 
		 * }
		 * 
		 * } //"END OF ENERGY DEPOSITED IN THE CAVITY"
		 * 
		 * //"SET FLAG FOR SECONDARY INTERACTIONS"
		 * //"***********************************"
		 * 
		 * if((EGS4Macro.IFULL>0)&&(IARG>5)&&(LATCHL == 0)) {
		 * //"ONLY IF PRIMARY PARTICLES HAVE INTERACTED DISCRETELY                 "
		 * /
		 * /"IF A SECONDARY PARTICLE IS CREATED ON THE SECOND PASS, GIVE IT A ZERO"
		 * /
		 * /"WEIGHT SO THAT HOWFAR WILL DISCARD IT.                               "
		 * 
		 * if( IARG == 7 ) {// "brem has occured"
		 * 
		 * if( IQL == 0 ) { EXCHANGE_STACK(EGS4.NP,EGS4.NP-1); }
		 * EGS4.LATCH[EGS4.NP-2] = 1; //" Flag the photon as a secondary" if(
		 * ipass >= 1 ) { EGS4.WT[EGS4.NP-2] = 0.; }//
		 * "To save time in correlation runs"
		 * 
		 * } else if( IARG == 18 ) {//
		 * "Compton has occured, with binding effects"
		 * //"taken into account, 0, 1, or more particles" "may have resulted"
		 * 
		 * for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [ { if(
		 * EGS4.IQ[ip-1] == 0 ) { EGS4.LATCH[ip-1] = 1; if( ipass >=1 ) {
		 * EGS4.WT[ip-1] = 0.; }// "To save time" } }
		 * //" NP = NPold means the interactiopn has been rejected and thus "
		 * //" the emerging (unscattered) photon  is still a primary " } else
		 * if( IARG == 9 ) {// "Moller has occured. For now there is only" //
		 * "one secondary. When impact ionization is implemented" //
		 * "the following should be changed" if( EGS4.E[EGS4.NP-1] <
		 * EGS4.E[EGS4.NP-2] ) { EXCHANGE_STACK(EGS4.NP,EGS4.NP-1); } } else if(
		 * IARG == 13 || IARG == 14 ) {//"Annihilation, flag the photons"
		 * EGS4.LATCH[EGS4.NP-1] = 1; EGS4.LATCH[EGS4.NP-2] = 1; if( ipass >= 1
		 * ) { EGS4.WT[EGS4.NP-1] = 0.; EGS4.WT[EGS4.NP-2] = 0.; } } else if
		 * (IARG == 20) { for( ip=EGS4.NPold;ip<=EGS4.NP;ip++)//DO ip=NPold,NP [
		 * { if( EGS4.IQ[ip-1] == 0 ) { EGS4.LATCH[ip-1] = 1; if( ipass >= 1 ) {
		 * EGS4.WT[ip-1] = 0.; EGS4.E[ip-1] = 0.; } } } } else if( IARG == 24 )
		 * { EGS4.LATCH[EGS4.NP-1] = 1; if( ipass >= 1 ) { EGS4.WT[EGS4.NP-1] =
		 * 0.; } } }
		 */
		if (EGS4Macro.ienhance == 1) {// " Option to enhance photon cross section in some region"
		// "write(6,*) 'in ausgab to recreate photon';"
			if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
				// "A pair/Compton/photoelectric/pair event is about to take place"
				EGS4.NP = EGS4.NP + 1; // "Boost the stack"
				if (EGS4.NP > EGS4.$MXSTACK) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Calculation with CS-enhancement: unable to boost stack."
							+ "  \n"
							+ " ***************************************************";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);
					return;
				}
				// "Create an identical photon"
				// $TRANSFER PROPERTIES TO (np) FROM (np - 1);
				EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
				EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
				EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
				EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
				EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
				EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

				EGS4.E[EGS4.NP - 1] = EGS4.E[EGS4.NP - 2];
				EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NP - 2];
				EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NP - 2];
				EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NP - 2];
				EGS4.IQ[EGS4.NP - 1] = EGS4.IQ[EGS4.NP - 2];

				R1 = EGS4.random01();
				aux1 = 1. - 1. / EGS4Macro.cs_enhance_current;
				if (R1 > aux1) {
					EGS4.WT[EGS4.NP - 2] = 0.;
					EGS4.E[EGS4.NP - 2] = 0.;
					EGS4.DNEAR[EGS4.NP - 2] = -1.;
				}
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] * 1.0
						/ EGS4Macro.cs_enhance_current;
				// "write(6,*) ' Creating an unscattered photon!
				// ',ir(np),wt(np),wt(np-1);
				return;
			}// "end of block after photon events"

			// " Play Russian Roulette with scattered photons to avoid transport
			// " of many low weight particles.
			// " Now there is a stack pointer NPold which points to the particle "
			// " before the last discrete interaction. This change was necessary "
			// " for the implementation of atomic relaxations. So, check all particles"
			// " between NPold and NP for RR"
			if ((IARG == 18 // " Compton has occured"
					|| IARG == 20 // " After photo-absorption "
			|| IARG == 24) // " After Rayleigh "
					&& EGS4Macro.ienhance == 1) {
				EGS4Macro.ienhance = 0;
				// "write(6,*) ' iarg = ',iarg,' NP NPold: ',NP,NPold;

				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] == 0) {
						R1 = EGS4.random01();
						if (R1 < 1.0 / EGS4Macro.cs_enhance_current) {
							// "PARTICLE SURVIVES"
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1]
									* EGS4Macro.cs_enhance_current;
							// "write(6,*) ' iarg = ',iarg,' particle has
							// survived ',
							// " wt(ip),ip;
						} else {
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
							EGS4.DNEAR[ip - 1] = -1.;
							// "write(6,*) ' Have killed particle ',ip,
							// " ' after iarg = ',iarg;
						}
					}
				}
				return;
			}

		}

		// "Check if particle is leaving the transport geometry"
		IRL = EGS4.IR[EGS4.NP - 1]; // "local region number"
		if (IRL == 1) {
			if (IWATCH > 0)
				EGS4.WATCH(IARG, IWATCH); // "signal watch routine if active"
			return; // "outside the chamber, howfar will discard"
		}

		// "obtain frequently used local variables"
		IRL = EGS4.IR[EGS4.NP - 1];// IR(NP);
		if (IRL == 1)
			return; // "outside the chamber"
			// IZD=IDSTBL(IRL,1);IXD=IDSTBL(IRL,2); "dose zone coordinates"
		IZD = IDSTBL[IRL - 1][0];
		IXD = IDSTBL[IRL - 1][1];
		IQL = EGS4.IQ[EGS4.NP - 1]; // "local variable"
		// "write(1,*) ' Ausgabe: iarg iq edep = ',iarg,iql,edep
		// IF(NHSTRY={P1}_LAST{P2})[
		// {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
		// ]
		// ELSE[
		// {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
		// {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
		// {P1}_TMP{P2}={P3};
		// {P1}_LAST{P2}=NHSTRY;
		// ]

		if (IARG == 0) {// "about to transport a particle"
			if (IQL != 0) {
				// $SCORE(SCSTP, :1);"count charged particle steps taken"
				if (EGS4SrcEns.NHSTRY == SCSTP_LAST) {
					SCSTP_TMP = SCSTP_TMP + 1;
				} else {
					SCSTP = SCSTP + SCSTP_TMP;
					SCSTP2 = SCSTP2 + SCSTP_TMP * SCSTP_TMP;
					SCSTP_TMP = 1;
					SCSTP_LAST = EGS4SrcEns.NHSTRY;
				}

				if (EGS4Grid.CDSTBL[IRL - 1].compareTo(DCHAR) == 0) {
					// $SCORE(SCDSTP, :1);
					if (EGS4SrcEns.NHSTRY == SCDSTP_LAST) {
						SCDSTP_TMP = SCDSTP_TMP + 1;
					} else {
						SCDSTP = SCDSTP + SCDSTP_TMP;
						SCDSTP2 = SCDSTP2 + SCDSTP_TMP * SCDSTP_TMP;
						SCDSTP_TMP = 1;
						SCDSTP_LAST = EGS4SrcEns.NHSTRY;
					}
				}
				// "count steps in dose region"
			} else {// "photon step - play russian roulette?"
				if (RUSROU && (EGS4.W[EGS4.NP - 1] > 0.0)) {// "yes, play if crosses RRZ "
					if ((EGS4.Z[EGS4.NP - 1] <= RRZ)
							&& (EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
									* EGS4.W[EGS4.NP - 1] >= RRZ)) {// "crosses"
						xsi = EGS4.random01();
						if (xsi < RRCUT) {// "it survives"
							EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] / RRCUT;
						} else {// "discard it on next call to HOWFAR"
							EGS4.WT[EGS4.NP - 1] = 0.0;
						}
					}// "end test if crosses russian roulette plane"
				}// "end test for playing russian roulette"
			}// "end test for photon step"
		}// "end test for IARG = 0"

		if (IWATCH > 0)
			EGS4.WATCH(IARG, IWATCH); // "SIGNAL WATCH ROUTINE IF ACTIVE"

		if (EGS4Macro.IFULL == 1) {
			// "check to see if any electrons created BY PHOTONS"
			// "if so, set bit 5 to 1 (ie dose from inside). If e- created by e-"
			// "(ie moller, bhaba), just pass latch value on.  This is"
			// "our arbitrary definition of where dose comes from."
			if (IARG == 16 || IARG == 18 || IARG == 20) {// " after pair production,compton or p.e event"
				if (EGS4.NP >= EGS4.NPold) {// "needed for p.e. case where russian roulette may have"
				// "eliminated ALL electrons"
					for (int II = EGS4.NPold; II <= EGS4.NP; II++) {
						if (EGS4.IQ[II - 1] != 0) {
							// "first, clear latch bits 1-5"
							for (int I = 1; I <= 5; I++) {
								// EGS4.LATCH[II-1]=IBCLR(LATCH(II),I);
								EGS4.LATCH[II - 1] = EGS4.IBCLR_LATCH(II, I);
							}
							// "now set latch bit 5"
							// LATCH(II)=IBSET(LATCH(II),5);
							EGS4.LATCH[II - 1] = EGS4.IBSET_LATCH(II, 5);
						}
					}
				}
			} else if (EGS4.IQ[EGS4.NP - 1] == 0 && (IARG == 1 || IARG == 2)) {
				// "In the rare case that a photon is terminated "
				// "because E<PCUT or E<AP.  We also want"
				// "this to show as dose originating from"
				// "within the volume"
				// "first, clear latch bits 1-5"
				// DO I=1,5[LATCH(NP)=IBCLR(LATCH(NP),I);]
				// "now set latch bit 5"
				// LATCH(NP)=IBSET(LATCH(NP),5);
				for (int I = 1; I <= 5; I++) {
					EGS4.LATCH[EGS4.NP - 1] = EGS4.IBCLR_LATCH(EGS4.NP, I);
				}
				EGS4.LATCH[EGS4.NP - 1] = EGS4.IBSET_LATCH(EGS4.NP, 5);

			} else if (((IARG == 5) && (EGS4.IRNEW != EGS4.IROLD))
					|| DECISION == 1) {
				// "SET LATCH FOR A PARTICLE ENTERING A NEW REGION DURING A STEP"
				// "used to be below, but found that bit setting was wrong for particles"
				// "ending their track on a region boundary"
				// DO I=1,5[LATCH(NP)=IBCLR(LATCH(NP),I);]
				// LATCH(NP)=IBSET(LATCH(NP),NEWNRC/10);
				for (int I = 1; I <= 5; I++) {
					EGS4.LATCH[EGS4.NP - 1] = EGS4.IBCLR_LATCH(EGS4.NP, I);
				}
				EGS4.LATCH[EGS4.NP - 1] = EGS4.IBSET_LATCH(EGS4.NP,
						EGS4Macro.NEWNRC / 10);

				DECISION = 0;
			}
		} else if (EGS4Macro.IFULL == 3) {// "set latch to score scattered dose separately"
		// "currently set as follows:
		// " scattered dose includes:
		// " any dose from compton scattered photons
		// " any dose from fluorescent photon which is re-absorbed
		// "
		// " This counts as primary any dose after bremsstrahlung from e- Why???
		// " and dose from fluorescent photons when not transported.
		// "

			if (IARG == 18) {// "After Compton. With binding effects and subsequent "
			// "atomic relaxations implemented, there may be 0, 1  "
			// "or more additional particles on the stack. NPold   "
			// "is a current addition to STACK. It is set to NP    "
			// "at the beginning of every scattering routine       "

				if (EGS4.NP > EGS4.NPold || EGS4.i_survived_RR > 0) {// " Flag all photons as scattered "
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
						if (EGS4.IQ[ip - 1] == 0) {
							// latch(ip) = IBSET(latch(ip),6);
							EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 6);
						}
					}
				}
				// " NPold = NP & i_survived_rr=0 after Compton means the interaction"
				// "  was rejected due to bindfing effects                           "
				// " => emerging photon is not scattered => don't flag it            "
			} else if (IARG == 20) {// "A photo-absorption has occured. "
			// "If NPold = NP, no particles with energies above "
			// "the specified thresholds resulted from the "
			// "relaxation cascade. If NPold < NP, there were "
			// "such particles emitted. Check them and flag all"
			// "fluorescent photons as secondaries"
				if (EGS4.NP > EGS4.NPold || EGS4.i_survived_RR > 0) {
					// DO ip=NPold,NP
					// ["The particle at NPold is always the photo-electron"
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
						// IF( iq(ip) = 0 ) [ latch(ip) = IBSET(latch(ip),7); ]
						if (EGS4.IQ[ip - 1] == 0) {
							EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 7);
						}
					}
				}
			}
		}// " end IFULL=3 "

		// "IKERMA=1 option will not work with new photon physics!!!"
		// "I intend to fix it in the future, IK January 1999"
		if (IKERMA == 1) {// "want to score KERMA"
		// "we score kerma for scattered component too if ifull=1"
		// "the kerma is part of scattered kerma if latch of initial photon is not zero"
		/*
		 * if(EGS4SrcEns.NHSTRY={P1}_LAST{P2})[ {P1}_TMP{P2}={P1}_TMP{P2} +
		 * {P3}; ] ELSE[ {P1}{P2}={P1}{P2}+{P1}_TMP{P2}; {P1}2{P2}={P1}2{P2} +
		 * {P1}_TMP{P2}*{P1}_TMP{P2}; IF('{P1}'='SCKERMA')[
		 * SCKERMA_TMPOLD{P2}={P1}_TMP{P2}; SCKERMA_LASTOLD{P2}={P1}_LAST{P2}; ]
		 * ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
		 * {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
		 * SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
		 * SCKERMA_TMPOLD{P2}; ] {P1}_TMP{P2}={P3}; {P1}_LAST{P2}=NHSTRY; ]
		 */
			// if (IARG == 4 &&
			// !EGS4.BTEST_LATCH(LATCH(NP),8))["local energy deposition"
			if (IARG == 4 && !EGS4.BTEST_LATCH(EGS4.NP, 8)) {// "local energy deposition"
			// "include deposited energy as kerma"
			// $SCOREDK(SCKERMA,(IZD,IXD,1):WT(NP)*EDEP);
				if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][0]) {
					SCKERMA_TMP[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0]
							+ EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;
				} else {
					SCKERMA[IZD - 1][IXD - 1][0] = SCKERMA[IZD - 1][IXD - 1][0]
							+ SCKERMA_TMP[IZD - 1][IXD - 1][0];
					SCKERMA2[IZD - 1][IXD - 1][0] = SCKERMA2[IZD - 1][IXD - 1][0]
							+ SCKERMA_TMP[IZD - 1][IXD - 1][0]
							* SCKERMA_TMP[IZD - 1][IXD - 1][0];
					// IF('{P1}'='SCKERMA')[
					SCKERMA_TMPOLD[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0];
					SCKERMA_LASTOLD[IZD - 1][IXD - 1][0] = SCKERMA_LAST[IZD - 1][IXD - 1][0];
					// ]
					// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
					// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
					// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
					// SCKERMA_TMPOLD{P2};
					// ]
					SCKERMA_TMP[IZD - 1][IXD - 1][0] = EGS4.WT[EGS4.NP - 1]
							* EGS4.EDEP;
					SCKERMA_LAST[IZD - 1][IXD - 1][0] = EGS4SrcEns.NHSTRY;
				}

				// if(EGS4Macro.IFULL=3 && (BTEST(LATCH(NP),6) |
				// BTEST(LATCH(NP),7)))[
				if (EGS4Macro.IFULL == 3
						&& (EGS4.BTEST_LATCH(EGS4.NP, 6) || EGS4.BTEST_LATCH(
								EGS4.NP, 7))) {
					// $SCOREDK(SCKERMA,(IZD,IXD,2):WT(NP)*EDEP);
					if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][1]) {
						SCKERMA_TMP[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1]
								+ EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;
					} else {
						SCKERMA[IZD - 1][IXD - 1][1] = SCKERMA[IZD - 1][IXD - 1][1]
								+ SCKERMA_TMP[IZD - 1][IXD - 1][1];
						SCKERMA2[IZD - 1][IXD - 1][1] = SCKERMA2[IZD - 1][IXD - 1][1]
								+ SCKERMA_TMP[IZD - 1][IXD - 1][1]
								* SCKERMA_TMP[IZD - 1][IXD - 1][1];
						// IF('{P1}'='SCKERMA')[
						SCKERMA_TMPOLD[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1];
						SCKERMA_LASTOLD[IZD - 1][IXD - 1][1] = SCKERMA_LAST[IZD - 1][IXD - 1][1];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
						// SCKERMA_TMPOLD{P2};
						// ]
						SCKERMA_TMP[IZD - 1][IXD - 1][1] = EGS4.WT[EGS4.NP - 1]
								* EGS4.EDEP;
						SCKERMA_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
					}

				}
			}
			if (IARG == 16) {// "pair event just occured"
				if (EGS4.NP > EGS4.NPold || EGS4.i_survived_RR > 0) {
					// DO IP=NPold,NP[
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
						// if(IQ[ip-1]!=0 && ~BTEST(LATCH(IP),8))[
						if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
							// $SCOREDK(SCKERMA,(IZD,IXD,1):WT(IP)*(E(IP)-PRM));
							if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][0]) {
								SCKERMA_TMP[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0]
										+ EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
							} else {
								SCKERMA[IZD - 1][IXD - 1][0] = SCKERMA[IZD - 1][IXD - 1][0]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][0];
								SCKERMA2[IZD - 1][IXD - 1][0] = SCKERMA2[IZD - 1][IXD - 1][0]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][0]
										* SCKERMA_TMP[IZD - 1][IXD - 1][0];
								// IF('{P1}'='SCKERMA')[
								SCKERMA_TMPOLD[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0];
								SCKERMA_LASTOLD[IZD - 1][IXD - 1][0] = SCKERMA_LAST[IZD - 1][IXD - 1][0];
								// ]
								// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
								// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
								// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
								// SCKERMA_TMPOLD{P2};
								// ]
								SCKERMA_TMP[IZD - 1][IXD - 1][0] = EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
								SCKERMA_LAST[IZD - 1][IXD - 1][0] = EGS4SrcEns.NHSTRY;
							}
							// if(EGS4Macro.IFULL == 3 & (BTEST(LATCH(IP),6) |
							// BTEST(LATCH(IP),7)))[
							if (EGS4Macro.IFULL == 3
									&& (EGS4.BTEST_LATCH(ip, 6) || EGS4
											.BTEST_LATCH(ip, 7))) {
								// $SCOREDK(SCKERMA,(IZD,IXD,2):WT(IP)*(E(IP)-PRM));
								if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][1]) {
									SCKERMA_TMP[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1]
											+ EGS4.WT[ip - 1]
											* (EGS4.E[ip - 1] - EGS4.PRM);
								} else {
									SCKERMA[IZD - 1][IXD - 1][1] = SCKERMA[IZD - 1][IXD - 1][1]
											+ SCKERMA_TMP[IZD - 1][IXD - 1][1];
									SCKERMA2[IZD - 1][IXD - 1][1] = SCKERMA2[IZD - 1][IXD - 1][1]
											+ SCKERMA_TMP[IZD - 1][IXD - 1][1]
											* SCKERMA_TMP[IZD - 1][IXD - 1][1];
									// IF('{P1}'='SCKERMA')[
									SCKERMA_TMPOLD[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1];
									SCKERMA_LASTOLD[IZD - 1][IXD - 1][1] = SCKERMA_LAST[IZD - 1][IXD - 1][1];
									// ]
									// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
									// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
									// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
									// SCKERMA_TMPOLD{P2};
									// ]
									SCKERMA_TMP[IZD - 1][IXD - 1][1] = EGS4.WT[ip - 1]
											* (EGS4.E[ip - 1] - EGS4.PRM);
									SCKERMA_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
								}
							}
							EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
						}
					}
				}
			}// "end of pair case"
			if (IARG == 18) {// "compton event just occured"
			// "must score kerma for all resultant electrons"
				if (EGS4.NP > EGS4.NPold) {// "compton occurred and we have not cleared the stack with"
				// "russian roulette"
				// DO IP=NPold,NP[
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
						// IF(IQ(IP) ~= 0 & ~BTEST(LATCH(IP),8))
						// ["score kerma for the electron"
						if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
							// $SCOREDK(SCKERMA,(IZD,IXD,1):WT(IP)*(E(IP) -
							// PRM));
							if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][0]) {
								SCKERMA_TMP[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0]
										+ EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
							} else {
								SCKERMA[IZD - 1][IXD - 1][0] = SCKERMA[IZD - 1][IXD - 1][0]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][0];
								SCKERMA2[IZD - 1][IXD - 1][0] = SCKERMA2[IZD - 1][IXD - 1][0]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][0]
										* SCKERMA_TMP[IZD - 1][IXD - 1][0];
								// IF('{P1}'='SCKERMA')[
								SCKERMA_TMPOLD[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0];
								SCKERMA_LASTOLD[IZD - 1][IXD - 1][0] = SCKERMA_LAST[IZD - 1][IXD - 1][0];
								// ]
								// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
								// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
								// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
								// SCKERMA_TMPOLD{P2};
								// ]
								SCKERMA_TMP[IZD - 1][IXD - 1][0] = EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
								SCKERMA_LAST[IZD - 1][IXD - 1][0] = EGS4SrcEns.NHSTRY;
							}
							// if(EGS4Macro.IFULL==3 & (BTEST(LATCH(IP),6) |
							// BTEST(LATCH(IP),7)))[
							if (EGS4Macro.IFULL == 3
									&& (EGS4.BTEST_LATCH(ip, 6) || EGS4
											.BTEST_LATCH(ip, 7))) {
								// $SCOREDK(SCKERMA,(IZD,IXD,2):WT(IP)*(E(IP)-PRM));
								if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][1]) {
									SCKERMA_TMP[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1]
											+ EGS4.WT[ip - 1]
											* (EGS4.E[ip - 1] - EGS4.PRM);
								} else {
									SCKERMA[IZD - 1][IXD - 1][1] = SCKERMA[IZD - 1][IXD - 1][1]
											+ SCKERMA_TMP[IZD - 1][IXD - 1][1];
									SCKERMA2[IZD - 1][IXD - 1][1] = SCKERMA2[IZD - 1][IXD - 1][1]
											+ SCKERMA_TMP[IZD - 1][IXD - 1][1]
											* SCKERMA_TMP[IZD - 1][IXD - 1][1];
									// IF('{P1}'='SCKERMA')[
									SCKERMA_TMPOLD[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1];
									SCKERMA_LASTOLD[IZD - 1][IXD - 1][1] = SCKERMA_LAST[IZD - 1][IXD - 1][1];
									// ]
									// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
									// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
									// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
									// SCKERMA_TMPOLD{P2};
									// ]
									SCKERMA_TMP[IZD - 1][IXD - 1][1] = EGS4.WT[ip - 1]
											* (EGS4.E[ip - 1] - EGS4.PRM);
									SCKERMA_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
								}
							}
							// LATCH(IP)=IBSET(LATCH(IP),8);
							EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
						}
					}
				}
			}// "end of compton case"
			if (IARG == 20) {// "photoelectric event just occured"
			// DO IP=NPold,NP[
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					// IF(IQ(IP)~=0 & ~BTEST(LATCH(IP),8))[
					if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
						// $SCOREDK(SCKERMA,(IZD,IXD,1):WT(IP)*(E(IP) - PRM));
						if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][0]) {
							SCKERMA_TMP[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0]
									+ EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
						} else {
							SCKERMA[IZD - 1][IXD - 1][0] = SCKERMA[IZD - 1][IXD - 1][0]
									+ SCKERMA_TMP[IZD - 1][IXD - 1][0];
							SCKERMA2[IZD - 1][IXD - 1][0] = SCKERMA2[IZD - 1][IXD - 1][0]
									+ SCKERMA_TMP[IZD - 1][IXD - 1][0]
									* SCKERMA_TMP[IZD - 1][IXD - 1][0];
							// IF('{P1}'='SCKERMA')[
							SCKERMA_TMPOLD[IZD - 1][IXD - 1][0] = SCKERMA_TMP[IZD - 1][IXD - 1][0];
							SCKERMA_LASTOLD[IZD - 1][IXD - 1][0] = SCKERMA_LAST[IZD - 1][IXD - 1][0];
							// ]
							// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
							// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
							// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
							// SCKERMA_TMPOLD{P2};
							// ]
							SCKERMA_TMP[IZD - 1][IXD - 1][0] = EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
							SCKERMA_LAST[IZD - 1][IXD - 1][0] = EGS4SrcEns.NHSTRY;
						}
						// IF(IFULL = 3 & (BTEST(LATCH(IP),6) |
						// BTEST(LATCH(IP),7)))[
						if (EGS4Macro.IFULL == 3
								&& (EGS4.BTEST_LATCH(ip, 6) || EGS4
										.BTEST_LATCH(ip, 7))) {
							// $SCOREDK(SCKERMA,(IZD,IXD,2):WT(IP)*(E(IP)-PRM));
							if (EGS4SrcEns.NHSTRY == SCKERMA_LAST[IZD - 1][IXD - 1][1]) {
								SCKERMA_TMP[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1]
										+ EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
							} else {
								SCKERMA[IZD - 1][IXD - 1][1] = SCKERMA[IZD - 1][IXD - 1][1]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][1];
								SCKERMA2[IZD - 1][IXD - 1][1] = SCKERMA2[IZD - 1][IXD - 1][1]
										+ SCKERMA_TMP[IZD - 1][IXD - 1][1]
										* SCKERMA_TMP[IZD - 1][IXD - 1][1];
								// IF('{P1}'='SCKERMA')[
								SCKERMA_TMPOLD[IZD - 1][IXD - 1][1] = SCKERMA_TMP[IZD - 1][IXD - 1][1];
								SCKERMA_LASTOLD[IZD - 1][IXD - 1][1] = SCKERMA_LAST[IZD - 1][IXD - 1][1];
								// ]
								// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
								// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
								// SCDOSEtoKERMA2{P2}=SCDOSEtoKERMA2{P2}+{P1}_TMP{P2}*
								// SCKERMA_TMPOLD{P2};
								// ]
								SCKERMA_TMP[IZD - 1][IXD - 1][1] = EGS4.WT[ip - 1]
										* (EGS4.E[ip - 1] - EGS4.PRM);
								SCKERMA_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
							}
						}
						// LATCH(IP)=IBSET(LATCH(IP),8);
						EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
					}
				}
			}// "end of photoelectric case"
		}// "end of IKERMA = 1, kerma scoring block"

		// "do some basic checks to see if scoring is needed"
		if (IARG >= 5 || EGS4.EDEP == 0)
			return;

		// "score total energy deposited"
		// "============================="

		FTMP = EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;

		if (EGS4Grid.CDSTBL[IRL - 1].compareTo(DCHAR) == 0) {// "IN A DOSE SCORING REGION"
		// "SCORE TOTAL ENERGY DEPOSITED"
		// $SCOREDK(SCDOSE,(IZD,IXD,1):FTMP);
			if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][0]) {
				SCDOSE_TMP[IZD - 1][IXD - 1][0] = SCDOSE_TMP[IZD - 1][IXD - 1][0]
						+ FTMP;
			} else {
				SCDOSE[IZD - 1][IXD - 1][0] = SCDOSE[IZD - 1][IXD - 1][0]
						+ SCDOSE_TMP[IZD - 1][IXD - 1][0];
				SCDOSE2[IZD - 1][IXD - 1][0] = SCDOSE2[IZD - 1][IXD - 1][0]
						+ SCDOSE_TMP[IZD - 1][IXD - 1][0]
						* SCDOSE_TMP[IZD - 1][IXD - 1][0];
				// IF('{P1}'='SCKERMA')[
				// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
				// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
				// ]
				// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
				// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
				if (IKERMA == 1
						&& SCDOSE_LAST[IZD - 1][IXD - 1][0] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][0]) {
					SCDOSEtoKERMA2[IZD - 1][IXD - 1][0] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][0]
							+ SCDOSE_TMP[IZD - 1][IXD - 1][0]
							* SCKERMA_TMPOLD[IZD - 1][IXD - 1][0];
				}
				SCDOSE_TMP[IZD - 1][IXD - 1][0] = FTMP;
				SCDOSE_LAST[IZD - 1][IXD - 1][0] = EGS4SrcEns.NHSTRY;
			}

			// if((EGS4Macro.IFULL == 3) && (BTEST(LATCH(NP),6) |
			// BTEST(LATCH(NP),7)))[
			if ((EGS4Macro.IFULL == 3)
					&& (EGS4.BTEST_LATCH(EGS4.NP, 6) || EGS4.BTEST_LATCH(
							EGS4.NP, 7))) {
				// $SCOREDK(SCDOSE,(IZD,IXD,2):FTMP);
				if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][1]) {
					SCDOSE_TMP[IZD - 1][IXD - 1][1] = SCDOSE_TMP[IZD - 1][IXD - 1][1]
							+ FTMP;
				} else {
					SCDOSE[IZD - 1][IXD - 1][1] = SCDOSE[IZD - 1][IXD - 1][1]
							+ SCDOSE_TMP[IZD - 1][IXD - 1][1];
					SCDOSE2[IZD - 1][IXD - 1][1] = SCDOSE2[IZD - 1][IXD - 1][1]
							+ SCDOSE_TMP[IZD - 1][IXD - 1][1]
							* SCDOSE_TMP[IZD - 1][IXD - 1][1];
					// IF('{P1}'='SCKERMA')[
					// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
					// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
					// ]
					// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
					// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
					if (IKERMA == 1
							&& SCDOSE_LAST[IZD - 1][IXD - 1][1] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][1]) {
						SCDOSEtoKERMA2[IZD - 1][IXD - 1][1] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][1]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][1]
								* SCKERMA_TMPOLD[IZD - 1][IXD - 1][1];
					}
					SCDOSE_TMP[IZD - 1][IXD - 1][1] = FTMP;
					SCDOSE_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
				}

			}
			if (IARG == 0) {
				// "SCORE TOTAL ENERGY DEPOSITED LESS STOPPED/DISCARDED"
				if (EGS4Macro.IFULL != 3) {
					// $SCOREDK(SCDOSE,(IZD,IXD,2):FTMP);
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][1]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][1] = SCDOSE_TMP[IZD - 1][IXD - 1][1]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][1] = SCDOSE[IZD - 1][IXD - 1][1]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][1];
						SCDOSE2[IZD - 1][IXD - 1][1] = SCDOSE2[IZD - 1][IXD - 1][1]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][1]
								* SCDOSE_TMP[IZD - 1][IXD - 1][1];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][1] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][1]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][1] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][1]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][1]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][1];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][1] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][1] = EGS4SrcEns.NHSTRY;
					}

				}
			}
			if ((IWATCH > 1) && (IWATCH != 4)) {
				EGS4.seqStr = " 	Weighted dose deposition  = "
						+ EGS4.format(FTMP, 14) + " MeV. IRL= "
						+ EGS4.format(IRL, 3) + " IARG= "
						+ EGS4.format(IARG, 3);
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				// OUTPUT FTMP,IRL,IARG;
				// (9x,' ***weighted dose deposition = ',1PE14.7,
				// ' MeV. IRL= ',I3, ' IARG= ',I3);
			}
		}

		if (EGS4Macro.IFULL == 1) {
			// "SCORE TOTAL ENERGY INTO BINS ACCORDING TO WHICH WALL THE PARTICLE CAME"
			// "FROM. CORNER SHOTS ATTRIBUTED TO PLANAR WALL : SEE ASSIGNMENTS IN HOWFAR"
			if (EGS4Grid.CDSTBL[IRL - 1].compareTo(DCHAR) == 0) {// "IN A DOSE SCORING REGION"
			// $GET-IX-IZ(IRL);
				IX = EGS4Geom.GET_IX(IRL);
				//IZ = EGS4Geom.GET_IZC(IRL);

				for (int I = 1; I <= 5; I++) {
					// if(BTEST(LATCH(NP),I))
					if (EGS4.BTEST_LATCH(EGS4.NP, I)) {
						IGEOM = I * 10;
						break;// EXIT;
					}
				}
				if (IGEOM == 10) {
					// $SCOREDK(SCDOSE,(IZD,IXD,3):FTMP); "FRONT WALL"
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][2]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][2] = SCDOSE_TMP[IZD - 1][IXD - 1][2]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][2] = SCDOSE[IZD - 1][IXD - 1][2]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][2];
						SCDOSE2[IZD - 1][IXD - 1][2] = SCDOSE2[IZD - 1][IXD - 1][2]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][2]
								* SCDOSE_TMP[IZD - 1][IXD - 1][2];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][2] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][2]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][2] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][2]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][2]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][2];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][2] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][2] = EGS4SrcEns.NHSTRY;
					}

				} else if (IGEOM == 20) {
					// $SCOREDK(SCDOSE,(IZD,IXD,4):FTMP); "OUTSIDE WALL" ]
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][3]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][3] = SCDOSE_TMP[IZD - 1][IXD - 1][3]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][3] = SCDOSE[IZD - 1][IXD - 1][3]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][3];
						SCDOSE2[IZD - 1][IXD - 1][3] = SCDOSE2[IZD - 1][IXD - 1][3]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][3]
								* SCDOSE_TMP[IZD - 1][IXD - 1][3];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][3] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][3]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][3] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][3]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][3]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][3];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][3] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][3] = EGS4SrcEns.NHSTRY;
					}

				} else if (IGEOM == 30) {
					// $SCOREDK(SCDOSE,(IZD,IXD,5):FTMP); "BACK WALL" ]
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][4]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][4] = SCDOSE_TMP[IZD - 1][IXD - 1][4]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][4] = SCDOSE[IZD - 1][IXD - 1][4]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][4];
						SCDOSE2[IZD - 1][IXD - 1][4] = SCDOSE2[IZD - 1][IXD - 1][4]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][4]
								* SCDOSE_TMP[IZD - 1][IXD - 1][4];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][4] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][4]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][4] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][4]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][4]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][4];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][4] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][4] = EGS4SrcEns.NHSTRY;
					}

				} else if (IGEOM == 40) {
					// $SCOREDK(SCDOSE,(IZD,IXD,6):FTMP); "INSIDE WALL"
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][5]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][5] = SCDOSE_TMP[IZD - 1][IXD - 1][5]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][5] = SCDOSE[IZD - 1][IXD - 1][5]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][5];
						SCDOSE2[IZD - 1][IXD - 1][5] = SCDOSE2[IZD - 1][IXD - 1][5]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][5]
								* SCDOSE_TMP[IZD - 1][IXD - 1][5];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][5] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][5]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][5] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][5]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][5]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][5];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][5] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][5] = EGS4SrcEns.NHSTRY;
					}

					if (IX == 1) {// "BUG"
						IGBUG2 = IGBUG2 + 1;
						if (IGBUG2 <= 100) {
							EGS4.seqStr = " 	INSIDE DOSE???. BUG NO."
									+ EGS4.format(IGBUG2, 3) + " IGEOM= "
									+ EGS4.format(IGEOM, 3);
							if (EGS4.iprint > 1)
								printSequence(EGS4.seqStr);
							// "OUTPUT IGBUG2,IGE0M; changed 92/11/09"
							// OUTPUT IGBUG2,IGEOM;
							// (' **** INSIDE DOSE???. BUG NO.',I3,'
							// IGEOM=',I3);
						}
					}
				} else if (IGEOM == 50) {
					// $SCOREDK(SCDOSE,(IZD,IXD,7):FTMP); "INSIDE SOURCE"]
					if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZD - 1][IXD - 1][6]) {
						SCDOSE_TMP[IZD - 1][IXD - 1][6] = SCDOSE_TMP[IZD - 1][IXD - 1][6]
								+ FTMP;
					} else {
						SCDOSE[IZD - 1][IXD - 1][6] = SCDOSE[IZD - 1][IXD - 1][6]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][6];
						SCDOSE2[IZD - 1][IXD - 1][6] = SCDOSE2[IZD - 1][IXD - 1][6]
								+ SCDOSE_TMP[IZD - 1][IXD - 1][6]
								* SCDOSE_TMP[IZD - 1][IXD - 1][6];
						// IF('{P1}'='SCKERMA')[
						// SCKERMA_TMPOLD[IZD-1][IXD-1][0]=SCKERMA_TMP[IZD-1][IXD-1][0];
						// SCKERMA_LASTOLD[IZD-1][IXD-1][0]=SCKERMA_LAST[IZD-1][IXD-1][0];
						// ]
						// ELSEIF('{P1}'='SCDOSE' & IKERMA=1 &
						// {P1}_LAST{P2}=SCKERMA_LASTOLD{P2})[
						if (IKERMA == 1
								&& SCDOSE_LAST[IZD - 1][IXD - 1][6] == SCKERMA_LASTOLD[IZD - 1][IXD - 1][6]) {
							SCDOSEtoKERMA2[IZD - 1][IXD - 1][6] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][6]
									+ SCDOSE_TMP[IZD - 1][IXD - 1][6]
									* SCKERMA_TMPOLD[IZD - 1][IXD - 1][6];
						}
						SCDOSE_TMP[IZD - 1][IXD - 1][6] = FTMP;
						SCDOSE_LAST[IZD - 1][IXD - 1][6] = EGS4SrcEns.NHSTRY;
					}

				} else {// "BUG"
					IGBUG1 = IGBUG1 + 1;
					if (IGBUG1 <= 100) {
						// "OUTPUT IGBUG1,IGE0M; changed 92/11/09"
						// OUTPUT IGBUG1,IGEOM;
						// (' **** LOST REGION. BUG NO.',I3,' IGEOM=',I3);
						EGS4.seqStr = " 	LOST REGION. BUG NO."
								+ EGS4.format(IGBUG1, 3) + " IGEOM= "
								+ EGS4.format(IGEOM, 3);
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);

					}
				}
			}

			// "RECOVER REGION ORIENTATION FLAG IF PARTICLE WILL BE DISCARDED"
			if ((EGS4.NP > 1) && (IARG >= 1) && (IARG <= 3)) {
				for (int I = 1; I <= 5; I++) {
					if (EGS4.BTEST_LATCH(EGS4.NP - 1, I)) {
						EGS4Macro.NEWNRC = I * 10;
						break;// EXIT;
					}
				}
			}

		}// "END OF IFULL=1"

		if (EGS4Macro.IFULL == 2 && IPHR[IRL - 1] != 0) {
			if (EGS4SrcEns.NHSTRY == SCPDST_LAST) {// "same primary history"
			// "keep adding energy in pulse height sensitive region"
				PHENER = PHENER + EGS4.EDEP;
				WT1OLD = EGS4.WT[0];
			} else {
				if (PHENER > 0.) {
					// "we either have a new history depositing or are at the end of a batch/run"
					// "find appropriate bin for pulse height"
					// "distn and add initial particle weight to it"
					// "FIND WHAT BIN WE ARE IN"
					if (SLOTE > 0.0) {
						// "EQUAL ENERGY BINS CASE"
						// IB=Math.min(IFIX(PHENER/SLOTE+0.999),$EBIN);
						Double dbll = new Double(PHENER / SLOTE + 1.0);// NOT
																		// ZEROES
						IB = Math.min(dbll.intValue(), $EBIN);
					} else {
						IB = MAXBIN;
						// UNTIL((IB.EQ.1).OR.(BINTOP(IB-1).LT.PHENER))
						// [IB=IB-1;]
						while (true) {
							if ((IB == 1) || (BINTOP[IB - 2] < PHENER))
								break;
							IB = IB - 1;
						}
					}

					// "ACCUMULATE THE PULSE HEIGHT DISTRIBUTION"
					// "USE WT(1) from history that contributed since we may have forcing"
					// "on for a primary source in the case of a phsp source WT(1) must be 1"
					SCPDST[IB - 1] = SCPDST[IB - 1] + WT1OLD;
					SCPDST2[IB - 1] = SCPDST2[IB - 1] + WT1OLD * WT1OLD;
					// "also add this to the cumulative pulse height distn"
					for (int ICUM = IB; ICUM <= MAXBIN; ICUM++) {
						SCPCUM[ICUM - 1] = SCPCUM[ICUM - 1] + WT1OLD;
						SCPCUM2[ICUM - 1] = SCPCUM2[ICUM - 1] + WT1OLD * WT1OLD;
					}

					if (IWATCH == 3) {
						EGS4.seqStr = " PULSE HEIGHT ENERGY="
								+ EGS4.format(PHENER, 10, true)
								+ " MeV, IN BIN" + EGS4.format(IB, 3)
								+ " WITH WEIGHT" + EGS4.format(1, 10, false);// !!!!!!!!!!!!!!!!!!!!!
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);
						// OUTPUT PHENER,IB,1;
						// (' PULSE HEIGHT ENERGY=',
						// F10.4,' MeV, IN BIN',I3,' WITH WEIGHT',1PE10.3);
					}

					// "NOW SCORE PROBABILITIES FOR COUNTS IN PEAKS"
					for (int IPK = 1; IPK <= 4; IPK++) {
						// "FOR EACH PEAK, F.E., ESCAPES AND 511"
						if ((PHENER >= DFEN[IPK - 1][1])
								&& (PHENER <= DFEN[IPK - 1][2])) {
							// "IT IS IN THE PEAK"
							SCDFEP[IPK - 1] = SCDFEP[IPK - 1] + WT1OLD;
							SCDFEP2[IPK - 1] = SCDFEP2[IPK - 1] + WT1OLD
									* WT1OLD;
							SCDFDIFF[IPK - 1] = SCDFDIFF[IPK - 1] + WT1OLD;
							SCDFDIFF2[IPK - 1] = SCDFDIFF2[IPK - 1] + WT1OLD
									* WT1OLD;
							if (IWATCH == 3) {
								EGS4.seqStr = " 	IT WAS IN ONE OF THE PEAKS,IPK="
										+ EGS4.format(IPK, 3);
								if (EGS4.iprint > 1)
									printSequence(EGS4.seqStr);
								// OUTPUT IPK;(T50,'IT WAS IN ONE OF THE
								// PEAKS,IPK=',I3/);
							}
						} else if ((PHENER >= DFEN[IPK - 1][0])
								&& (PHENER < DFEN[IPK - 1][1])) {
							// "IT IS IN THE BKGD"
							SCDFBK[IPK - 1] = SCDFBK[IPK - 1] + WT1OLD;
							SCDFBK2[IPK - 1] = SCDFBK2[IPK - 1] + WT1OLD
									* WT1OLD;
							SCDFDIFF[IPK - 1] = SCDFDIFF[IPK - 1] - WT1OLD;
							SCDFDIFF2[IPK - 1] = SCDFDIFF2[IPK - 1] - WT1OLD
									* WT1OLD;
						}
					}// "END IPK LOOP"
				}
				SCPDST_LAST = EGS4SrcEns.NHSTRY;
				PHENER = EGS4.EDEP;
				WT1OLD = EGS4.WT[0];
			}
		}

		return;
	}// "END OF AUSGAB"

	// "The following is the $CALL-HOWNEAR macro for PRESTA-II

	// REPLACE {$CALL-HOWNEAR(#);} WITH {
	// ;
	// "write(6,'(2i3,4e15.8)') np,ir(np),dnear(np), "
	// "   sqrt(x(np)*x(np)+y(np)*y(np)),z(np),tustep; "
	// IF( dnear(np) < tustep ) [
	// call hownear({P1},x(np),y(np),z(np),ir(np));
	// " write(6,*) ' --> new dnear: ',{P1}; "
	// ]
	// ELSE [ {P1} = dnear(np); ]
	// }
	// "*********************************************************************"
	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		// "Subroutine arguments
		// $REAL
		// tperp, "nearest distance to any boundary (output)
		// tustep,
		// x, "x-position of the particle (input)
		// y, "y-position of the particle (input)
		// z; "z-position of the particle (input)

		// $INTEGER
		// ir "region number of the particle

		// "Local variables
		double r = 0.0;
		int ix = 0;// "current cylindrical radius number
		int iz = 0;// "current planar slab number

		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		int ir = EGS4.IR[EGS4.NP - 1];
		// FROM MACRO=>ELSE [ {P1} = dnear(np); ]->not necessary, kind of
		// varred(skip calc.)!@@@@@@@@@@@
		// if( EGS4.DNEAR[EGS4.NP-1] >= EGS4.TUSTEP)
		// {
		// EGS4.tperp=EGS4.DNEAR[EGS4.NP-1];
		// return;
		// }
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		ix = (ir - 2) / EGS4Geom.NZ + 1;// NZ=planar zones=NPLANE-1!!!NR=nr of
										// Radii
		iz = ir - 1 - EGS4Geom.NZ * (ix - 1);
		r = Math.sqrt(x * x + y * y);
		EGS4.tperp = EGS4.min(z - EGS4Geom.ZPLANE[iz - 1], EGS4Geom.ZPLANE[iz]
				- z, EGS4Geom.RCYL[ix] - r);
		if (ix != 1) {
			EGS4.tperp = Math.min(EGS4.tperp, r - EGS4Geom.RCYL[ix - 1]);
		}

		// "IF( tperp < -1e-6 ) [
		// " OUTPUT IHSTRY, tperp; (/' Error in HOWNEAR: IHSTRY=',I12,
		// " ' tperp negative =',1PE10.3);
		// " OUTPUT ir,x,y,z,r; (' ir x y z r=', I5, 4F15.7);
		// " OUTPUT iz,ix;(' depth,radial regions iz,ix=',2I10);
		// " OUTPUT zplane(iz),zplane(iz+1);( ' upper/lower z planes=',2F15.7);
		// " IF( ix > 1) [
		// " OUTPUT rcyl(ix-1),rcyl(ix); (' Radial boundaries=',2F15.7);
		// " ]
		// " ELSE [
		// " OUTPUT rcyl(ix); (' Radial boundary of 1st region =', F15.7);
		// " ]
		// "]

		return;

	}// "end of subroutine HOWNEAR"

	// ;"******************************************************************************
	// "
	// " **********
	// " * *
	// " * HOWFAR *
	// " * *
	// " **********
	// "
	// " A GENERAL PURPOSE CYLINDRICAL GEOMETRY ROUTINE FOR USE WITH THE EGS4
	// " CODE SYSTEM ADAPTED FOR USE WITH CAVRZnrc.
	// "
	// " FOR PARTICLE NP ON THE STACK IN REGION IR(NP), THIS ROUTINE
	// " DETERMINES IF THE PARTICLE CAN GO A DISTANCE USTEP WITHOUT CHANGING
	// " ZONES. IF USTEP CAUSES A ZONE CROSSING, IT IS REDUCED TO PLACE IT ON
	// " THE BOUNDRY AND IRNEW IS SET TO THE ZONE NUMBER ON THE FAR SIDE OF
	// " THE BOUNDARY. IF IR(NP) IS 1 THEN THE PARTICLE HAS ESCAPED THE REGION
	// " OF INTEREST AND THE HISTORY IS TERMINATED.(IDISC IS SET TO 1.)
	// "
	// "
	// "
	// "
	// " SOME VARIABLES
	// " ==============
	// ;"
	// "OUTEND = .TRUE. => PARTICLE MAY TRANSMIT OR BACKSCATTER OUT ENDS
	// " = .FALSE. => PARTICLE STAYS WITHIN THE END BOUNDARIES
	// "OUTSID = .TRUE. => PARTICLE MAY TRANSMIT OUT THE SIDES
	// " = .FALSE. => PARTICLE STAYS WITHIN THE SIDE BOUNDARY
	// "IRL = STARTING REGION NUMBER THE PARTICLE IS IN
	// "IZ = STARTING PLANAR ZONE NUMBER THE PARTICLE IS IN.
	// " THE PARTICLE IS BETWEEN ZPLANE(IZ) AND ZPLANE(IZ+1).
	// "IX = STARTING CYLINDRICAL ZONE NUMBER THE PARTICLE IS IN.
	// " THE PARTICLE IS BETWEEN RCYL(IX-1) AND RCYL(IX).
	// "
	// " COMMON/GEOM/
	// " ZPLANE(IZ) Z VALUES OF PLANES
	// " 1<=IZ<=NZ+1
	// " RCYL(IRR) RADII OF CYLINDERS
	// " 1<=IRR<=NR
	// " CYRAD2(IRR) =RCYL(IRR)**2
	// " NZ # PLANAR GEOMETRICAL ZONES (NPLANE-1)
	// " ZONE(I) IS BETWEEN ZPLANE(I) AND ZPLANE(I+1)
	// " NR # CYLINDRICAL GEOMETRICAL ZONES
	// " ZONE(I) IS BETWEEN RCYL(I-1) AND RCYL(I)
	// " NREG TOTAL # GEOMETRICAL ZONES =NR*NZ +1
	// " +1 FOR VACUUM ENVELOPE
	// %E "cavrznrc.mortran"
	// " DEFINITIONS OF REGION NUMBER, PLANAR ZONE, CYLINDRICAL ZONE
	// " ===========================================================
	// " Z AXIS RUNS ACROSS PAGE SHOWN AS .......
	// "
	// "
	// " 1
	// /"
	// "
	// " --------------------------------------------------------- RCYL(NR)
	// " |(NR-1) |(NR-1) |(NR-1) | . . . . | NR*NZ | NR*NZ | IX=NR
	// " | *NZ+2 | *NZ+3 | *NZ+4 | | | +1 |
	// " --------------------------------------------------------- RCYL(NR-1)
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " --------------------------------------------------------- RCYL(2)
	// " | NZ+2 | NZ+3 | NZ+4 | . . . . | 2NZ | 2NZ+1 | IX=2
	// " --------------------------------------------------------- RCYL(1)
	// "..1....|...2...|...3...|...4...|...............|...NZ..|..NZ+1.|....IX=1..1..
	// ;" ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " ---------------------------------------------------------
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " | | | | | | |
	// " ---------------------------------------------------------
	// " IZ=1 IZ=2 IZ=3 IZ=NZ-1 IZ=NZ
	// "
	// " 1
	// "
	// "
	// "
	// " VERSION 1 ADAPTED FROM CAVITY HOWFAR 06/84 ERIC FOX
	// " VERSION 2 THE SUBROUTINE CALLS TO PLANES AND 10/87 AFB
	// " CYLINDER HAVE BEEN REPLACED BY MACROS
	// " TO SPEED THINGS UP
	// "
	// "
	// "******************************************************************************

	/**
	 * The following is a general specification of HOWFAR: 
	 * Given a particle at (X,Y,Z) in region IR and going in direction 
	 * (U,V,W), this routine answers the question, can the particle go 
	 * a distance USTEP without crossing a boundary? If yes, it merely returns; 
	 * If no, it sets USTEP=distance to boundary in the current 
	 * direction and sets IRNEW to the region number on the far side of the boundary. 
	 * Interface method.
	 */
	public void HOWFAR() {
		// $IMPLICIT-NONE;

		// "MACRO USED LOCALLY TO CHANGE REGIONS, ADJUST USTEP, AND EXIT"
		// REPLACE {$SET NEW REGION(#,#);} WITH
		// {
		// IF({P1}.LE.USTEP)[USTEP={P1};IRNEW={P2};]RETURN;
		// }

		// "Debug macro, usually not used, to verify that region is correct on entry"
		// "    It must be used after IX and IZ have been determined"
		// "    REAL Rdebug must be declared"
		// REPLACE {$Check_region;} WITH
		// {
		// Rdebug = SQRT(X(NP)**2 + Y(NP)**2);
		// IF(IX = 1)["we are in inner radial region"
		// IF(Rdebug > (RCYL(IX)+1.E-6) | (Z(NP) < ZPLANE(IZ) & Z(NP) >
		// ZPLANE(IZ+1)))[
		// OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX),
		// ZPLANE(IZ),ZPLANE(IZ+1);
		// (/' Error in HOWFAR, particle in wrong region'/
		// ' IRL, IX,IZ =', 3I7/
		// ' X,Y,Z,R=', 4F15.7/
		// ' RCYL(1)=', F15.7/
		// ' ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		// ]
		// ]"end of IX=1 block"
		// ELSE [
		// IF((Rdebug > (RCYL(IX)+1E-6) & Rdebug < (RCYL(IX-1) -1.E-6)) |
		// (Z(NP) < ZPLANE(IZ) & Z(NP) > ZPLANE(IZ+1)))[
		// OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX-1),
		// RCYL(IX),
		// ZPLANE(IZ),ZPLANE(IZ+1);
		// (/' Error in HOWFAR, particle in wrong region'/
		// ' IRL, IX,IZ =', 3I7/
		// ' X,Y,Z,R=', 4F15.7/
		// ' RCYL(IX-1), RCYL(IX)=', 2F15.7/
		// ' ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		// ]
		// ]
		// ;}

		// LOGICAL OUTEND,OUTSID;

		// $INTEGER IRL,IX,IZ,IHITP,IHITC,IZNEW,IXNEW;
		// $REAL WL,TPLANE,U1,V1,A,TCYL,X1,Y1,B,B2,C,COUT,CIN,RAD;
		// "$REAL Rdebug; " "Uncomment if using $Check_region;"

		// "First set idisc and irnew "
		// EGS4.IDISC = 0; EGS4.IRNEW = EGS4.IR[EGS4.NP-1];

		// "DISCARD ZERO WEIGHT PARTICLES"
		if (EGS4.WT[EGS4.NP - 1] == 0.0) {
			EGS4.IDISC = 1;
			return;
		}

		// "INITIALLY ASSUME PARTICLE STAYS IN THE TARGET"
		boolean OUTEND = false;
		boolean OUTSID = false;
		//int IQL = EGS4.IQ[EGS4.NP - 1];
		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"

		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY OR OF THE REGION IS TOTALLY"
		// "ABSORBING"
		if ((IRL == 1) || (EGS4Grid.CABSRB[IRL - 1].compareTo(ACHAR) == 0)) {
			EGS4.IDISC = 1;
			return;
		}

		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY"
		// if(IRL == 1){EGS4.IDISC=1;return;}

		// $GET-IX-IZ(IRL); //"GET PLANAR AND CYLINDRICAL ZONES NUMBERS"
		int IX = EGS4Geom.GET_IX(IRL);
		int IZ = EGS4Geom.GET_IZC(IRL);

		// "Following commented out usually"
		// "$Check_region;"
		// "write(6,*);"
		// "write(6,*) 'howfar: ',ix,iz,ustep,nz,nr;"
		// "CALCULATE DNEAR"
		// "Note this is same as $CALL-HOWNEAR-FOR-NRCC-CYLINDRICAL-GEOMETRY(DNEAR(NP))"
		// "It is redundant if PRESTA is being used since DNEAR is calculated on"
		// "every step anyway"

		// $PLANES(IZ,IZ+1,IHITP,TPLANE,ustep);"GET DISTANCE TO PLANE"
		// "IHITP  =  1 => HITS GREATER Z PLANE"
		// "       =  0 => MISSES BOTH PLANES"
		// "       = -1 => HITS LESSER Z PLANE"
		// EGS4Geom.ustep=EGS4.USTEP;
		EGS4Geom.PLANES2(IZ, IZ + 1);
		// EGS4.USTEP=EGS4Geom.ustep;
		// $CYLNDR(IX,IHITC,TCYL,ustep);"GET DISTANCE TO CYLINDER"
		// "       IHITC   =  1 => HITS OUTER CYLINDER"
		// "               =  0 => MISSES BOTH CYLINDERS"
		// "               = -1 => HITS INNER CYLINDER"
		// EGS4Geom.ustep=EGS4.USTEP;
		EGS4Geom.CYLNDR2(IX);
		// EGS4.USTEP=EGS4Geom.ustep;

		int IZNEW = IZ + EGS4Geom.ihitp;// IHITP; "GET NEW PLANAR REGION"
		if ((IZNEW < 1) || (IZNEW > EGS4Geom.NZ))
			OUTEND = true; // "FLAG IF LEAVES BY THE ENDS"

		int IXNEW = IX + EGS4Geom.ihitc;// IHITC; "GET NEW CYLINDRICAL REGION"
		if (IXNEW > EGS4Geom.NR)
			OUTSID = true; // "FLAG IF LEAVES BY THE SIDES"

		int NWNRCL = 0;
		//int IMSOFF = 0;
		// "DO MOST PROBABLE CASE FIRST WHERE A PLANE AND A CYLINDER CAN BE HIT"
		if ((EGS4Geom.ihitp != 0) && (EGS4Geom.ihitc != 0)) {
			if (EGS4Geom.tplane < EGS4Geom.tcyl) {// "HITS PLANE FIRST"
				if (OUTEND) {
					NWNRCL = 0;
					// $SET NEW REGION(TPLANE,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tplane <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tplane;// {P1};
						EGS4.IRNEW = 1;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				} else {
					NWNRCL = 20 - 10 * EGS4Geom.ihitp;// IHITP;
					// $SET NEW REGION(TPLANE,IRL+IHITP);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tplane <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tplane;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.ihitp;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				}
			} else if (EGS4Geom.tcyl < EGS4Geom.tplane) {// "HITS CYLINDER FIRST"
				if (OUTSID) {
					NWNRCL = 0;
					// $SET NEW REGION(TCYL,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = 1;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				} else {
					NWNRCL = 30 + 10 * EGS4Geom.ihitc;
					// $SET NEW REGION(TCYL,IRL+NZ*IHITC);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.NZ * EGS4Geom.ihitc;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				}
			} else {// "ODD CASE TCYL=TPLANE:HITS PLANE AND CYLINDER TOGETHER"
				if (OUTEND || OUTSID) {
					NWNRCL = 0;
					// $SET NEW REGION(TCYL,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = 1;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				} else {
					NWNRCL = 20 - 10 * EGS4Geom.ihitp;// IHITP;
					// $SET NEW REGION(TCYL,IRL+IHITP+NZ*IHITC);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.ihitp + EGS4Geom.NZ
								* EGS4Geom.ihitc;// {P2};
						//IMSOFF = 0;
						EGS4Macro.NEWNRC = NWNRCL;
						if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
							DECISION = 1;
					}
					return;
				}
			}
		}

		// "DO ODD CASE-PARTICLE CAN HIT PLANE BUT NOT CYLINDER"
		else if (EGS4Geom.ihitp != 0) {
			if (OUTEND) {
				NWNRCL = 0;
				// $SET NEW REGION(TPLANE,1);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tplane <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tplane;// {P1};
					EGS4.IRNEW = 1;// {P2};
					//IMSOFF = 0;
					EGS4Macro.NEWNRC = NWNRCL;
					if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
						DECISION = 1;
				}
				return;
			} else {
				NWNRCL = 20 - 10 * EGS4Geom.ihitp;
				// $SET NEW REGION(TPLANE,IRL+IHITP);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tplane <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tplane;// {P1};
					EGS4.IRNEW = IRL + EGS4Geom.ihitp;// {P2};
					//IMSOFF = 0;
					EGS4Macro.NEWNRC = NWNRCL;
					if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
						DECISION = 1;
				}
				return;
			}
		}

		// "DO ODD CASE-PARTICLE CAN HIT CYLINDER BUT NOT PLANE"
		else {
			if (OUTSID) {
				NWNRCL = 0;
				// $SET NEW REGION(TCYL,1);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tcyl <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tcyl;// {P1};
					EGS4.IRNEW = 1;// {P2};
					//IMSOFF = 0;
					EGS4Macro.NEWNRC = NWNRCL;
					if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
						DECISION = 1;
				}
				return;

			} else {
				NWNRCL = 30 + 10 * EGS4Geom.ihitc;// IHITC;
				// $SET NEW REGION(TCYL,IRL+NZ*IHITC);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tcyl <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tcyl;// {P1};
					EGS4.IRNEW = IRL + EGS4Geom.NZ * EGS4Geom.ihitc;// {P2};
					//IMSOFF = 0;
					EGS4Macro.NEWNRC = NWNRCL;
					if (EGS4.USTEP == 0 && EGS4.IRNEW != EGS4.IROLD)
						DECISION = 1;
				}
				return;
			}
		}

		// "AT THIS STAGE ALL GEOMETRICAL POSSIBILITIES HAVE BEEN CHECKED AND CONTROL"
		// "HAS ALREADY BEEN TRANSFERRED TO EGS"

		// return;
	}// "END OF SUBROUTINE HOWFAR"

	/**
	 * Gather all media data required for this simulation.
	 */
	private void HATCH() {
		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "Start of run:         " + d.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.HATCH();
	}

	/**
	 * Start the shower, i.e. the actual simulation for electron-photon transport
	 */
	private void SHOWER() {
		// //CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
		EGS4Core.SHOWER(EGS4SrcEns.iqin, EGS4SrcEns.ein, EGS4SrcEns.xin,
				EGS4SrcEns.yin, EGS4SrcEns.zin, EGS4SrcEns.uin, EGS4SrcEns.vin,
				EGS4SrcEns.win, EGS4SrcEns.irin, EGS4SrcEns.WEIGHT);
	}

	/**
	 * Setup input variables.
	 */
	private void inputs() {
		// @ TITLE
		TITLEs = "dosrznrc_template--depth dose in H2O due to Cobalt beam";
		// @ IWATCH
		// #off,interactions,steps,deposited,graph;
		// #debug output with increasing detail, graph outputs .gph file for
		// EGS_Windows
		// #if not "off" use very few histories
		IWATCH = IWATCH_OFF;
		// @ STORE INITIAL RANDOM NUMBERS
		// #no,last,all deposited,all;
		// #last: store initial random numbers for last history in .egsrns
		// #all deposited: store initial random numbers for each history that
		// deposits energy
		// in the cavity
		// #all: store initial random numbers for each history
		ISTORE = ISTORE_NO;
		// @ IRESTART
		// #first,restart,make,analyze,for graphics,parallel;
		// #first: first run
		// #restart: restart of old run (requires .egsdat file)
		// #make: just create an input file and exit
		// #analyze: read in data from .egsdat file and do statistical analysis
		// and output results
		// # to .egslst
		// #for graphics: read starting random numbers from .egsrns--eg for
		// output to graphics package
		// #parallel: read .egsdat files from parallel jobs (named
		// inputfile_w*), do statistical
		// # analysis and output to .egslst
		IRESTART = IRESTART_FIRST;
		// @ OUTPUT OPTIONS
		// #short,cavity details;
		// #short: output cavity summary + dose grid
		// #cavity details: above plus details for every zone in cavity
		IOOPTN = IOOPTN_MATERIAL_SUMMARY;// IOOPTN_SHORT;
		// @ ELECTRON TRANSPORT
		// " = normal (0) normal electron transport (discrete interactions)
		// " = no interactions (1) no discrete interactions (used for CDSA
		// " calculations but note that special data
		// " sets are also needed to do a full CSDA
		// " calculation. All turning off
		// " interactions does is just that. See use
		// " of IUNRST=2,3,4 PEGS4 data sets for real CSDA)
		// " [ ICSDA]
		EGS4Macro.ICSDA = ETRANS_NORMAL;
		// ==============================
		// " DOSE ZBOUND MIN (I) Minimum plane # defining dose region
		// (default=1)
		// " [NZDMIN]
		// " DOSE ZBOUND MAX (I) Maximum plane # defining dose region
		// " [NZDMAX]
		// " DOSE RBOUND MIN (I) Minimum cylinder # defining dose region
		// (default=0)
		// " [NRDMIN]
		// " DOSE RBOUND MAX (I) Maximum cylinder # defining dose region
		// " [NRDMAX]
		NZDMIN = 1;
		NZDMAX = 61;
		NRDMIN = 0;
		NRDMAX = 60;
		// ===========================
		// @ STORE DATA ARRAYS
		// #yes,no;
		// #yes: output .egsdat file for restarts, parallel post-processing, etc
		IDAT = IDAT_NO;
		// @ NUMBER OF HISTORIES
		// #splits into $STAT statistical batches
		// #must be >=$STAT**2 if IWATCH= Off
		// #can have less than this if IWATCH set to another option
		NCASE = 10000;
		// @ MAX CPU HOURS ALLOWED
		// #Will shut down cleanly prior to exceeding this limit, as long as one
		// #batch has completed.
		TIMMAX = 90.000;
		// @ IFULL
		// #dose and stoppers,Aatt and Ascat,Ap,Afl and <s>g/w;
		// #dose and stoppers: output total dose plus that due to stoppers and
		// discards
		// #Aatt and Ascat: above plus Aatt, Ascat
		// #Ap: above plus Ap
		// #Afl and <s>g/w: above plus Afl and stopping power ratio gas/water
		EGS4Macro.IFULL = IFULL_PULSE_HEIGHT_DISTRIBUTION;// IFULL_AATT_AND_ASCAT;
		// @ STATISTICAL ACCURACY SOUGHT
		// #If 0, goes until number of histories or CPU limit exceeded.
		// #If not zero goes until this uncertainty (in %) is achieved in the
		// peak dose region
		STATLM = 0.0;
		// @ PHOTON REGENERATION
		// #no,yes,no electrons from wall;
		// #no: normal calculation
		// #yes: regenerate parent photon after interaction (used for FANO
		// calculations)
		// #no electrons from wall: photons not regenerated,
		// # secondary electrons from cavity wall are eliminated
		// IFANO=IFANO_NO;
		// @SCORE KERMA
		IKERMA = KERMA_YES;
		// @ INITIAL RANDOM NO. SEEDS
		// #With ranmar: these must be between 1 and 30081 (default to 9373)
		// #With ranlux: 1st is luxury level 0->4 allowed but should not be 0
		// # 2nd is seed 1 -> 1073741824
		jrng1 = 1;
		jrng2 = 3;
		// right here, initialize random generator!!!
		init_random_generator();
		// ###########################################################################################
		// @ METHOD OF INPUT
		// #groups,individual,cavity information:
		// #group: input groups of slabs of equal thickness
		// #individual: input Z of bottom of every slab
		// #cavity information: generate simple geometry from cavity info input
		// in section below.
		// # If you use this, there are no more inputs in this section
		EGS4Geom.iterseindex = EGS4Geom.iINDIVIDUAL;
		if (EGS4Geom.iterseindex == EGS4Geom.iCAVITY_INFORMATION)// 2
		{
			EGS4Geom.WALLTH = 0.5;// WALL THICKNESS
			EGS4Geom.CAVRAD = 1.0;// CAVITY OUTER RADIUS
			EGS4Geom.CAVLNG = 2.0;// CAVITY LENGTH
			EGS4Geom.ELERAD = 0.01;// ELECTRODE RADIUS
			EGS4Geom.SLENGHT = "H2O_fortran";// WALL MATERIAL
			EGS4Geom.airs = "AIR521ICRU_fortran";// ELERAD=0.0
			EGS4Geom.electrods = "AL521ICRU_fortran";// ELERAD!=0.0;
		} else// //ITERSE->0 sau 1
		{
			EGS4Geom.Z_OF_FRONT_FACE = 0.0;// #Beginning of first slab
			if (EGS4Geom.iterseindex == EGS4Geom.iGROUPS)// 0
			{
				// NSLAB
				EGS4Geom.nNSLAB = 2;// 2 value
				// #Define a group of 10 slabs with thickness 1 cm
				// #followed by 10 slabs with thickness 2 cm
				EGS4Geom.NSLAB[0] = 10;
				EGS4Geom.NSLAB[1] = 10;
				// SLAB THICKNESS
				EGS4Geom.DELTAZ[0] = 1.0;
				EGS4Geom.DELTAZ[1] = 2.0;
			}

			if (EGS4Geom.iterseindex == EGS4Geom.iINDIVIDUAL)// 1
			{
				// NSLAB
				EGS4Geom.nNSLAB = 3;
				// DEPTH BOUNDARIES
				EGS4Geom.ZPLANE[1] = 0.05;
				EGS4Geom.ZPLANE[2] = 6.35;
				EGS4Geom.ZPLANE[3] = 7.30;
			}

			EGS4Geom.nCyl = 2;// "number of radial cylinders input"
			// #Radii of cylinders
			EGS4Geom.RCYL[1] = 3.15;// 6.3/2
			EGS4Geom.RCYL[2] = 3.65;// 7.3/2

			// MEDIA=the media in the problem. These must match exactly,
			// including case, one
			// of the media names in the pegs4 data set being used in the
			// problem.
			// #Next we specify which media are in which geometric regions
			// #note that by default all regions contain
			// #medium 1 and which medium to input as 1 should be selected with
			// this in mind.
			EGS4Geom.nMEDIA = 2;
			EGS4.MEDIA[0] = "AL521ICRU_fortran";// "170C521ICRU_fortran";
			EGS4.MEDIA[1] = "NAI_Fortran";// "AIR521ICRU_fortran";

			// DESCRIPTION BY:#planes,regions;
			// #planes: use slab and cylinder no.'s to define what medium goes
			// where
			// #regions: use region numbers to define this (region numbers start
			// at 2 and
			// number from top to bottom of geometry and innermost radius to
			// outermost radius)
			EGS4Geom.DESCRIBE = EGS4Geom.DESCRIBE_REGIONS;// ##INPUT DATA

			EGS4Geom.nMEDNUM = 2;
			EGS4Geom.MEDNUM[0] = 1;// Al
			EGS4Geom.MEDNUM[1] = 2;// NAI

			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES_DENSITY) {
				EGS4Geom.nRHOR = 1;
				EGS4Geom.RHOR[0] = 0.;
			}
			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY) {
				EGS4Geom.nNREGLO = 2;
				EGS4Geom.nNREGHI = 2;
				EGS4Geom.NREGLO[0] = 2;// START REGION
				EGS4Geom.NREGHI[0] = 7;// STOP REGION
				EGS4Geom.NREGLO[1] = 3;// START REGION
				EGS4Geom.NREGHI[1] = 3;// STOP REGION
			}
			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES_DENSITY) {
				EGS4Geom.nNZLO = 2;
				EGS4Geom.nNZHI = 2;
				EGS4Geom.nNRLO = 2;
				EGS4Geom.nNRHI = 2;
				// #This puts MEDIUM 1 everywhere and then #inserts a small
				// column of MEDIUM 2
				// on the central #axis with radius 1cm and going from Z=10cm
				// #to Z=12cm
				EGS4Geom.NZLO[0] = 1;// START ZSLAB
				EGS4Geom.NZLO[1] = 10;// START ZSLAB
				EGS4Geom.NZHI[0] = 20;// STOP ZSLAB
				EGS4Geom.NZHI[1] = 12;// STOP ZSLAB
				EGS4Geom.NRLO[0] = 1;// START RING
				EGS4Geom.NRLO[1] = 1;// START RING
				EGS4Geom.NRHI[0] = 2;// STOP RING
				EGS4Geom.NRHI[1] = 1;// STOP RING
			}
		}// ITERSE->0 sau 1
		// cavity inputs:
		EGS4Geom.NSUMCV = 1;// NUMBER OF CAVITY REGIONS
		// #this defines the small cylinder of
		// #air in the geometry to be the cavity
		EGS4Geom.ISUMCV[0] = 3;// REGION NUMBERS OF THE CAVITY= 3

		// CALL GEOMRZ!!
		EGS4.seqStr = " *** Reading geometrical info ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4Geom.GEOMRZ();
		if (EGS4.STOPPROGRAM) {
			return;
		}

		// INCIDENT PARTICLE= photon #electron,photon,positron,all,charged;
		// #all & charged: only for phase space sources
		// #######################################PULSE HEIHHT
		// @REGION OF SENSITIVE VOLUME
		nREGSOLV = 1;
		REGSVOL[0] = 3;
		SLOTE = 0.01;
		DELTAE = 0.005;

		if (SLOTE < 0) {
			nTOPEBIN = 2;
			TOPEBIN[0] = 0.12;
			TOPEBIN[1] = 0.02;
		}
		// ###########################################################################################
		EGS4SrcEns.ipart = EGS4SrcEns.ipart_photon;
		// SOURCE NUMBER= 1 #0,1,2,3,4,10,11,12,13,14,15,16,20,21,22,23
		EGS4SrcEns.ISOURC = EGS4SrcEns.point_source_on_axis_incident_from_the_front;// 1;
		// SOURCE OPTIONS= 100., 1.3, 0, 0
		// #for source 1: SSD of beam, radius of # beam on front surface
		if (EGS4SrcEns.ISOURC != EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)// 20)
		{
			EGS4SrcEns.source_option[0] = 1.5;
			EGS4SrcEns.source_option[1] = 7.3;
			EGS4SrcEns.source_option[2] = 0.;
			EGS4SrcEns.source_option[3] = 1.;

			if (EGS4SrcEns.ISOURC == EGS4SrcEns.full_phase_space_from_any_angle)// 22)
			{
				EGS4SrcEns.source_option[4] = 0.;
				EGS4SrcEns.source_option[5] = 0.;
				EGS4SrcEns.source_option[6] = 0.;
				EGS4SrcEns.source_option[7] = 0.;
				EGS4SrcEns.source_option[8] = 0.;
			} else if (EGS4SrcEns.ISOURC == EGS4SrcEns.BEAM_treatment_head_simulation_from_any_angle)// 23
			{
				EGS4SrcEns.source_option[4] = 0.;
			}
		}
		if (EGS4SrcEns.ISOURC == EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)// 20)
		{
			// source 20 Ex:
			EGS4SrcEns.MODEIN = EGS4SrcEns.MODEIN_LOCAL;
			EGS4SrcEns.NRDIST = 1;
			EGS4SrcEns.RDISTF[0] = 1.0;
			EGS4SrcEns.RPDF[0] = 1.0;
			EGS4SrcEns.RDIST_IOUTSP = EGS4SrcEns.RDIST_IOUTSP_NONE;
		}
		// CALL SRCRZ->"Get source data"
		EGS4.seqStr = " *** Reading radiation source info ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4SrcEns.SRCRZ();
		if (EGS4.STOPPROGRAM) {
			return;
		}
		// source initialization
		EGS4SrcEns.SRCINI();
		if (EGS4.STOPPROGRAM) {
			return;
		}

		if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22
				|| EGS4SrcEns.ISOURC == 23) {
			EGS4SrcEns.MONOEN = 0;
		}
		// "no need to input monoen for source 21,22"
		else {
			// INCIDENT ENERGY= monoenergetic #monoenergetic, spectrum;
			EGS4SrcEns.monoindex = EGS4SrcEns.iMONOENERGETIC;
			// INCIDENT KINETIC ENERGY(MEV)= 1.25 #only use for "monoenergetic"
			EGS4SrcEns.ikemev = 0.662;

			EGS4SrcEns.ENSRC();
			if (EGS4.STOPPROGRAM) {
				return;
			}
		}// "Get data re-source energies"
		// ###########################################################################################
		// get_transport_parameter();
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		ecut = 0.521;
		pcut = 0.001;
		smax = 1.e10;

		if (!setEcutRegion)
			EGS4.ECUT[0] = ecut;// #Electron cutoff for transport
		if (!setPcutRegion)
			EGS4.PCUT[0] = pcut;// #Photon cutoff for transport
		if (!setSmaxRegion)
			EGS4.SMAXIR[0] = smax;// #Maximum step size in cm (not needed
		// #unless old PRESTA algorithm used)
		if (setEcutRegion) {
			nEcut = 1;// number of data
			Ecut[0] = 1.;
			startEcutRegion[0] = 1;
			stopEcutRegion[0] = 1;
		}
		if (setPcutRegion) {
			nPcut = 1;// number of data
			Pcut[0] = 1.;
			startPcutRegion[0] = 1;
			stopPcutRegion[0] = 1;
		}
		if (setSmaxRegion) {
			nSmax = 1;// number of data
			Smax[0] = 1.;
			startSmaxRegion[0] = 1;
			stopSmaxRegion[0] = 1;
		}
		// Bound Compton scattering= On or Off->IBCMP
		// #######//On means 1 and Off means 0
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// #Off: Klein-Nishina used for compton
		// # scattering
		// #On: Impulse approximation used for
		// # scattering
		// #It has not been established that the
		// #the scoring routines work with this
		// #option ON.
		incoh = incoh_ON;
		if (incoh == incoh_OFF || incoh == incoh_ON)
			setIncohRegion = false;
		else
			setIncohRegion = true;
		if (setIncohRegion) {
			nIncoh = 1;// number of data
			Incoh[0] = incoh_OFF;// 0 or 1
			startIncohRegion[0] = 1;
			stopIncohRegion[0] = 1;
		}
		// Rayleigh scattering= On->IRAYLR
		// #Off: no coherent scattering
		// #On: simulates coherent scattering
		coh = coh_OFF;
		if (coh == coh_OFF || coh == coh_ON)
			setCohRegion = false;
		else
			setCohRegion = true;
		if (setCohRegion) {
			nCoh = 1;// number of data
			Coh[0] = coh_OFF;// 0 or 1
			startCohRegion[0] = 1;
			stopCohRegion[0] = 1;
		}
		// Atomic relaxations= On->IEDGFL
		// #Note, it has not been verified that
		// #AUSGAB works with this option ON
		// #On: use correct cross section
		// # for p.e. events and shell vacancies
		// # for Compton & p.e. events are relaxed
		// # via emission of fluorescent X-Rays,
		// # Auger and Koster-Cronig electrons
		relax = relax_ON;
		if (relax == relax_OFF || relax == relax_ON)
			setRelaxRegion = false;
		else
			setRelaxRegion = true;
		if (setRelaxRegion) {
			nRelax = 1;// number of data
			Relax[0] = relax_OFF;// 0 or 1
			startRelaxRegion[0] = 1;
			stopRelaxRegion[0] = 1;
		}
		// Photoelectron angular sampling= On->IPHTER
		// #Off: Photoelectrons get direction of
		// # photon that creates them
		// #On: Sauter's formula is used
		pe = pe_ang_ON;
		if (pe == pe_ang_OFF || pe == pe_ang_ON)
			setPeRegion = false;
		else
			setPeRegion = true;
		if (setPeRegion) {
			nPe = 1;// number of data
			Pe[0] = pe_ang_OFF;// 0 or 1
			startPeRegion[0] = 1;
			stopPeRegion[0] = 1;
		}
		// Brems angular sampling= On ->IBRDST
		// #Simple: leading term of Koch-Motz
		// # dist'n used to determine angle
		// # of bremsstrahlung photons
		// #KM: Koch-Motz distribution used to
		// # determine angle -> 2BS(modified)
		EGS4.ibrdst = brems_ang_KM;
		// Pair angular sampling= On->IPRDST
		// #Simple: use leading term of K-M
		// # dist'n
		// #KM: use complete Koch and Motz dist'n
		// #Off: angle of pairs is m/E--like old EGS4
		EGS4.iprdst = pair_ang_SIMPLE;
		// ESTEPE= 0.25->ESTEPE
		// #Max fractional continuous energy loss
		// #per step. Use 0.25 unless using
		// #PRESTA-I
		EGS4.estepe = 0.25;
		// XIMAX= 0.5->XIMAX
		// #Max first elastic scattering moment
		// #per step. Using default.
		EGS4.ximax = 0.5;
		// Boundary crossing algorithm= exact->bca_algorithm, exact_bca
		// #exact: cross boundaries in single scattering
		// # mode (distance at which to go into
		// # single scattering mode determined by
		// # "Skin depth for BCA"
		// #PRESTA-I: cross boundaries with lateral
		// # correlations off and force multiple
		// # scattering mode
		EGS4.bca_algorithm = BCA_EXACT;// exact means
										// =EGS4.$BCA_ALGORITHM_DEFAULT=0;
		// Skin depth for BCA= 3->skindepth_for_bca
		// #Distance from a boundary (in elastic
		// #MFP) at which the algorithm will go
		// #into single scattering mode (using
		// #default here)
		EGS4.skindepth_for_bca = 3.0;
		// Electron-step algorithm= default->transport_algorithm
		// PRESTA-II (the default),$PRESTA_II = 0;
		// #Determines the algorithm used to take
		// #into account lateral and longitudinal
		// #correlations in a condensed history
		// #step
		EGS4.transport_algorithm = estep_alg_PRESTA_II;// 0;
		// Spin effects= On->spin_effects
		// #Turns off/on spin effects for electron
		// #elastic scattering. Spin On is
		// #ABSOLUTELY necessary for good
		// #backscattering calculations. Will
		// #make a difference even in `well
		// #conditioned' situations (e.g. depth
		// #dose curves).
		ispin = spin_ON;
		// Brems cross sections= BH
		// #BH: Bethe-Heitler cross-sections used
		// #NIST: NIST cross-sections used
		EGS4.ibr_nist = brems_cross_BH;
		// " Pair cross sections "= BH
		EGS4.pair_nrc = pair_cross_BH;
		// Electron impact ionization
		EGS4.eii_flag = eii_OFF;
		// Triplet
		EGS4.itriplet = triplet_OFF;
		// Radiative compton correction
		EGS4.radc_flag = radc_OFF;
		// ###########################################################################################
		// VARIANCE REDUCTION:
		// @BREM SPLITTING
		EGS4.nbr_split = 1;// $MAXBRSPLIT
		EGS4.i_play_RR = 0;
		// ELECTRON RANGE REJECTION
		EGS4Macro.irejct = irejct_OFF;
		// #On: if charged particle energy is below ESAVEIN
		// # and it cannot get out of current region
		// # with energy > ECUT, the particle is
		// # terminated
		// #also terminates all electrons which cannot
		// #reach the cavity under conservative assumptions.
		ESAVEIN = 2.0;// #total energy below which range rejection is considered

		// RUSSIAN ROULETTE DEPTH= 0.0000 #play Russian Roulette with photons
		// once they
		// #cross this Z plane
		RRZ = 0.000;
		// RUSSIAN ROULETTE FRACTION= 0.0000 #probability of photon survival--if
		// this
		// #and #RUSSIAN ROULETTE DEPTH both 0, then
		// #photon Russian Roulette is not played
		RRCUT = 0.000;
		// #exponential pathlength biasing can be
		// #used. See Rogers&Bielajew 1990 review for
		// #discussion. C<0 => pathlength shortening
		// # >0 => pathlength stretching
		// # along z axis both cases
		// #CAVRZnrc allows for having the photon cross
		// #section scaled to enhance interactions.
		// # If this input is missing or set to <= 1, it
		// # has no effect on the simulation. But if
		// # the enhancement factor is set to > 1, the
		// # effect is dramatic: all other user input
		// # concerning photon forcing, splitting, exp.
		// # transform, etc., is ignored. In addition,
		// # the calculation result corresponds ALWAYS
		// # to 'Aatt and Ascat', no matter what the
		// # user requested (but only Awall is calculated,
		// # not the individual Ascat and Aatt).
		// # The algorithm employed is implemented via
		// # $RAYLEIGH-CORRECTION and appropriate calls to
		// # AUSGAB. For more detail see the manual and the
		// # header of cavrznrc.mortran
		EGS4Macro.CEXPTR = 0.000;
		// PHOTON FORCING= On #Off (default),On;
		// #On: force photons to interact according to
		// # START FORCING and STOP FORCING AFTER inputs
		IFARCE = IFARCE_OFF;
		EGS4Macro.NFMIN = 1;// #Start forcing at this interaction number
		EGS4Macro.NFMAX = 1;// #Number of photon interactions after which
		// #to stop forcing photon interactions
		// PHOTON SPLITTING= 1 #no. of times to split a photon
		// #if < 2-->normal transport
		// #overrides PHOTON FORCING if >= 2
		// #can only be >= 2 if IFULL= dose and stoppers
		// # or if IFULL= Aatt and Ascat
		// phsplitt=1;

		EGS4Macro.cs_enhance = 1.0;// #Photon cross section scaling factors
		nENHREG = 2;
		NENHLO[0] = 1;
		NENHHI[0] = 1;
		NENHLO[1] = 1;
		NENHHI[1] = 1;
		// CS ENHANCEMENT FACTOR= 1
		// CS ENHANCEMENT START REGION= 1, 1
		// CS ENHANCEMENT STOP REGION= 1, 1

		// PLOTSN();

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			return;
		}
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		/*
		 * //" SCORING ARRAY INITIALISATION //" ****************************
		 * NCASEO=0;NCASET=0;TMCPUO=0; if(IRESTART==0||IRESTART==5) {//
		 * "FRESH START, SET EVERYTHING TO ZERO" EGS4SrcEns.NNREAD=0; SCSTP=0.;
		 * SCSTP2=0.; SCSTP_TMP=0.; SCSTP_LAST=0; SCCSTP=0.; SCCSTP2=0.;
		 * SCCSTP_TMP=0.; SCCSTP_LAST=0; cav_dose = 0; cav2_dose = 0; cav_dose0
		 * = 0; cav2_dose0 = 0; cav_dose1 = 0; cav2_dose1 = 0; cav_dose2 = 0;
		 * cav2_dose2 = 0;
		 * 
		 * cav_dosec = 0; cav_dosec01 = 0; cav_dosec02 = 0;
		 * 
		 * tmp_dose = 0.0; tmp_dose0 = 0.0; tmp_dose1 = 0.0; tmp_dose2 = 0.0;
		 * 
		 * PIISTP=0.; for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++) { for(int
		 * IX=1;IX<=EGS4Geom.NR;IX++)//DO IX=1,NR[ { for(int IT=1;IT<=4;IT++) {
		 * SCDOSE[IZ-1][IX-1][IT-1]=0.; SCDOSE2[IZ-1][IX-1][IT-1]=0.;
		 * SCDOSE_TMP[IZ-1][IX-1][IT-1]=0.;
		 * if(IT<4)SCDOSE_COV[IZ-1][IX-1][IT-1]=0.; } SCDOSE_LAST[IZ-1][IX-1]=0;
		 * } } } else if(IRESTART!=4)//NOT ALLOWED HERE {
		 * //"Restart or stats analysis only, read old data from unit 4"
		 * //"open unit 4 as an old file" //OUTPUT;(' About to read the previous
		 * .egsdat file'); // //my_unit = egs_open_datfile(4,0,1,'.egsdat');
		 * //READ(my_unit,*) SCSTP,SCSTP2,SCCSTP,SCCSTP2; //READ(my_unit,*)
		 * cav_dose, cav_dose0, cav_dose1, cav_dose2; /READ(my_unit,*)
		 * cav2_dose,cav2_dose0,cav2_dose1,cav2_dose2; //READ(my_unit,*)
		 * cav_dosec,cav_dosec01,cav_dosec02;
		 * 
		 * //IF(NSUMCV>1)["get data for individual cavity regions" // DO
		 * IZ=1,NZ[ // DO IX=1,NR[ // $GET-IRL(IZ,IX); // IF(NTRACK(IRL).EQ.1)[
		 * // READ(my_unit,*)(SCDOSE(IZ,IX,IT),SCDOSE2(IZ,IX,IT),IT=1,4); //
		 * READ(my_unit,*)(SCDOSE_COV(IZ,IX,IT),IT=1,3); // ] // ] // ] //]
		 * //$RETRIEVE RNG STATE FROM UNIT my_unit;
		 * //READ(my_unit,*,END=:EOFA:)NCASEO,TMCPUO,NNREADO,PIISTP;
		 * //CLOSE(my_unit); //NNREAD=NNREADO; // }
		 * 
		 * if(IRESTART == 3){NCASE=0;}
		 * 
		 * NCASET=NCASE+NCASEO; EGS4SrcEns.NCASET=NCASET;
		 * 
		 * EGS4.seqStr=" ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 */

		// "CALCULATE THE NUMBER OF DOSE COMPONENTS"
		if ((EGS4Macro.IFULL == 0) || (EGS4Macro.IFULL == 2)
				|| (EGS4Macro.IFULL == 3)) {
			ITMAX = 2;
		} else {
			ITMAX = $MAXIT;
		}

		NCASEO = 0;
		NCASET = 0;
		TMCPUO = 0;

		if (IRESTART == 0 || IRESTART == 5) {// "0-FRESH START, SET EVERYTHING TO"
		// "ZERO OR 5- INITIALIZE ALL ARRAYS TO 0 FOR PARALLEL POST-PROCESSING"

			// NNREAD=0; "have not read any particles from phase space source"

			SCSTP = 0;
			SCSTP2 = 0;
			SCSTP_TMP = 0;
			SCSTP_LAST = 0;
			SCDSTP = 0;
			SCDSTP2 = 0;
			SCDSTP_TMP = 0;
			SCDSTP_LAST = 0;
			PIISTP = 0;

			for (int IT = 1; IT <= ITMAX; IT++) {
				for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
					for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
						SCDOSE[IZ - 1][IX - 1][IT - 1] = 0.0;
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = 0.0;
						SCDOSE_TMP[IZ - 1][IX - 1][IT - 1] = 0.0;
						SCDOSE_LAST[IZ - 1][IX - 1][IT - 1] = 0;
						if (IKERMA == 1) {
							SCKERMA[IZ - 1][IX - 1][IT - 1] = 0.0;
							SCKERMA2[IZ - 1][IX - 1][IT - 1] = 0.0;
							SCKERMA_TMP[IZ - 1][IX - 1][IT - 1] = 0.0;
							SCKERMA_LAST[IZ - 1][IX - 1][IT - 1] = 0;
							SCDOSEtoKERMA2[IZ - 1][IX - 1][IT - 1] = 0.0;
						}
					}
				}
			}
			if (EGS4Macro.IFULL == 2) {
				for (int IB = 1; IB <= MAXBIN; IB++) {
					SCPDST[IB - 1] = 0.0;
					SCPDST2[IB - 1] = 0.0;
					SCPCUM[IB - 1] = 0.0;
					SCPCUM2[IB - 1] = 0.0;
				}
				for (int IPK = 1; IPK <= 4; IPK++) {
					SCDFEP[IPK - 1] = 0.0;
					SCDFEP2[IPK - 1] = 0.0;
					SCDFBK[IPK - 1] = 0.0;
					SCDFBK2[IPK - 1] = 0.0;
					SCDFDIFF[IPK - 1] = 0.0;
					SCDFDIFF2[IPK - 1] = 0.0;
				}
				SCPTOT = 0.0;
				SCPTOT2 = 0.0;
				SCPDST_LAST = 0;
			}
		}// "END OF IRESTART =0 OR 5"

		// ELSEIF(IRESTART.EQ.4)[ "retrieve random numbers from .egsrns file"
		// OUTPUT;(/' Will READ RANDOM NUMBER PARAMETERS FROM UNIT 2:');
		// ]
		/*
		 * ELSE[ "RESTART OR STATS ANALYSIS ONLY, READ OLD DATA FROM UNIT 4"
		 * "OPEN UNIT 4 AS AN OLD FILE" OUTPUT;(/' ***START READING DATA FILE
		 * from PREVIOUS RUN***'/);
		 * "IK: we don't wont ro rely upon a symbolic link to the .egsdat file "
		 * "    to have been made for us before running dosrznrc. "
		 * "    We therefore open the .egsdat file using a file name. "
		 * "    The egs_open_datfile function opens the file named output_file.egsdat"
		 * "    (or input_file.egsdat) where output_file and input_file are the "
		 * "    command line arguments to the -o and -i options."
		 * "OPEN(UNIT=4,file='fort.4',STATUS='OLD');" data_unit =
		 * egs_open_datfile(4,0,1,'.egsdat');
		 * READ(data_unit,*,END=:EOFA:)SCSTP,SCSTP2,SCDSTP,SCDSTP2,PIISTP;
		 * READ(data_unit,*,END=:EOFA:)
		 * (((SCDOSE(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
		 * READ(data_unit,*,END=:EOFA:)
		 * (((SCDOSE2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
		 * IF(IKERMA=1)[ READ(data_unit,*,END=:EOFA:)
		 * (((SCKERMA(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
		 * READ(data_unit,*,END=:EOFA:)
		 * (((SCKERMA2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE);
		 * READ(data_unit,*,END=:EOFA:)
		 * (((SCDOSEtoKERMA2(IZ,IX,IT),IT=1,ITMAX),IX=1,NRDOSE),IZ=1,NZDOSE); ]
		 * IF(IFULL=2)[
		 * READ(data_unit,*,END=:EOFA:)(SCPDST(IB),SCPDST2(IB),IB=1,MAXBIN);
		 * READ(data_unit,*,END=:EOFA:)(SCPCUM(IB),SCPCUM2(IB),IB=1,MAXBIN);
		 * READ(data_unit,*,END=:EOFA:)(SCDFEP(IPK),SCDFEP2(IPK),IPK=1,4);
		 * READ(data_unit,*,END=:EOFA:)(SCDFBK(IPK),SCDFBK2(IPK),IPK=1,4);
		 * READ(data_unit,*,END=:EOFA:)(SCDFDIFF(IPK),SCDFDIFF2(IPK),IPK=1,4); ]
		 * $RETRIEVE RNG STATE FROM UNIT data_unit;
		 * READ(data_unit,*,END=:EOFA:)NCASEO,TMCPUO,NNREAD;
		 * " Read solid angle information as well, IK May 4 1999"
		 * READ(data_unit,*,END=:OMEGA-NOT-THERE:) SCOMEG,SCOMEG2; goto
		 * :OMEGA-THERE:; :OMEGA-NOT-THERE: OUTPUT;(' Warning: Solid angle
		 * information not in data file'); OUTPUT;(' This may cause errors in
		 * absolute normalizations!'); :OMEGA-THERE: CLOSE(data_unit); ]
		 * 
		 * IF(IRESTART .EQ.3)[NCASE=0;]
		 */
		NCASET = NCASE + NCASEO;
		EGS4SrcEns.NCASET = NCASET;

		// :FINISHED: CONTINUE;
		// "************************"
		// "* Check for any errors *"
		// "************************"
		EGS4.seqStr = " ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

	}// inputs

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs() {
		if (IWATCH < 0 || IWATCH > 4) {
			IWATCH = 0;
		}// default
		EGS4Macro.IWATCH = IWATCH;
		EGS4Geom.IWATCH = IWATCH;
		// if(ISTORE<0 || ISTORE>3){ISTORE=0;}//default
		if (ISTORE < 0 || ISTORE > 0) {
			ISTORE = 0;
		}// default
		// if(IRESTART<0 || IRESTART>4){IRESTART=0;}//default
		if (IRESTART < 0 || IRESTART > 0) {
			IRESTART = 0;
		}// default
		if (IOOPTN < 0 || IOOPTN > 4) {
			IOOPTN = 0;
		}// default
		if (EGS4Macro.ICSDA < 0 || EGS4Macro.ICSDA > 1) {
			EGS4Macro.ICSDA = 0;
		}// default
		// ================
		if (NZDMIN < 1 || NZDMIN > EGS4Geom.$MAXZREG) {
			NZDMIN = 1;
		}
		if (NZDMAX < 2 || NZDMAX > EGS4Geom.$MAXZREG + 1) {
			NZDMAX = EGS4Geom.$MAXZREG + 1;
		}
		if (NRDMIN < 0 || NRDMIN > EGS4Geom.$MAXRADII - 1) {
			NRDMIN = 0;
		}
		if (NRDMAX < 1 || NRDMAX > EGS4Geom.$MAXRADII) {
			NRDMAX = EGS4Geom.$MAXRADII;
		}
		NZDOSE = NZDMAX - NZDMIN;
		NRDOSE = NRDMAX - NRDMIN;
		/*
		 * OUTPUT IWATCH,ISTORE,IRESTART ,IOOPTN,IDAT,ICSDA; ( / ' DO NOT
		 * TRACK(0) OR TRACK(>0) EVERY INTERACTION:',T72,I4/ ' DO NOT STORE(0)
		 * OR STORE(1,2) INITIAL RANDOM #s:',T72,I4/ ' FIRST
		 * RUN(0),RESTARTED(1),ANALYZE(3),start-RNS(4),parallel(5):', T72,I4/ '
		 * OUTPUT OPTION (1 THRU 4):',T72,I4/ ' STORE DATA(0) OR
		 * NOT(1):',T72,I4/ ' CSDA CALCULATION(1) OR NOT(0)':,T72,I4/); OUTPUT
		 * ;(/' ******DOSE SCORING BOUNDARIES******'); OUTPUT
		 * NZDMIN,NZDMAX,NRDMIN,NRDMAX; (' NZDMIN,NZDMAX,NRDMIN,NRDMAX:
		 * ',3(I5,','),I5/);
		 */
		// ==============
		if (IDAT < 0 || IDAT > 1) {
			IDAT = 1;
		}// default
		if (IDAT < 1 || IDAT > 1) {
			IDAT = 1;
		}// default->NO FILE STORAGE

		// IF(IDAT.EQ.1)[INEXT=0;]ELSE[INEXT=1;]

		if (IRESTART == 4)// no parrallelllllllllllll!!
		{
			IDAT = 1; // "do not store output in this case to avoid biasing"
			ISTORE = 0; // "do not store the starting random numbers either"
		}
		if (NCASE < NCASE_MIN || NCASE > NCASE_MAX) {
			NCASE = NCASE_DEFAULT;
		}
		if (TIMMAX < TIMMAX_MIN || TIMMAX > TIMMAX_MAX) {
			TIMMAX = TIMMAX_DEFAULT;
		}
		if (EGS4Macro.IFULL < 0 || EGS4Macro.IFULL > 4) {
			EGS4Macro.IFULL = 0;
		}// default
		if (STATLM < STATLM_MIN || STATLM > STATLM_MAX) {
			STATLM = STATLM_DEFAULT;
		}
		// if(IFANO<0 || IFANO>2){IFANO=0;}//default
		if (IWATCH == 0 && NCASE < $NCASEMIN) {
			NCASE = $NCASEMIN;
		}
		// if( IFANO == 1 )
		// {// "With ifano option turned on, it is a waste of time to"
		// "have RAYLEIGH turned on (scattered photon will be killed,"
		// "original photon re-created) => turn Rayleigh off"
		// for(int j=1;j<=EGS4.$MXREG;j++) { EGS4.IRAYLR[j-1] = 0; }
		// EGS4.seqStr=" ******** ifano set => turning off Rayleigh! **** ";
		// if(EGS4.iprint>1)
		// printSequence(EGS4.seqStr);

		// OUTPUT; (//' ******** ifano set => turning off Rayleigh! **** '//);
		// }
		if (IKERMA < 0 || IKERMA > 1) {
			IKERMA = 0;
		}// default

		if (EGS4Geom.iterseindex < 0 || EGS4Geom.iterseindex > 2) {
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
			printSequence(EGS4.seqStr);
			return;
		}
		if (EGS4Geom.iterseindex == EGS4Geom.iCAVITY_INFORMATION)// 2
		{
			if ((EGS4Geom.SLENGHT).compareTo(NULLs) == 0
					|| (EGS4Geom.airs).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				return;
			}
			if (EGS4Geom.ELERAD != 0.0
					&& (EGS4Geom.electrods).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				return;
			}
		}
		if (NZDMAX > EGS4Geom.NZ + 1) {
			NZDMAX = EGS4Geom.NZ + 1;
			// OUTPUT NZDMAX;(/'===> MAX. SCORING PLANE # RESET TO: ',I6);
			EGS4.seqStr = " ===> MAX. SCORING PLANE # RESET TO: "
					+ EGS4.format(NZDMAX, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		if (NZDMIN >= NZDMAX) {
			NZDMIN = NZDMAX - 1;
			// OUTPUT NZDMIN;(/'===> MIN. SCORING PLANE # RESET TO: ',I6);
			EGS4.seqStr = " ===> MIN. SCORING PLANE # RESET TO: "
					+ EGS4.format(NZDMIN, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		if (NRDMAX > EGS4Geom.NR) {
			NRDMAX = EGS4Geom.NR;
			// OUTPUT NRDMAX;(/'===> MAX. SCORING CYLINDER # RESET TO: ',I6);
			EGS4.seqStr = " ===> MAX. SCORING CYLINDER # RESET TO: "
					+ EGS4.format(NRDMAX, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		if (NRDMIN >= NRDMAX) {
			NRDMIN = NRDMAX - 1;
			// OUTPUT NRDMIN;(/'===> MIN. SCORING CYLINDER # RESET TO: ',I6);
			EGS4.seqStr = " ===> MIN. SCORING CYLINDER # RESET TO: "
					+ EGS4.format(NRDMIN, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		NZDOSE = NZDMAX - NZDMIN;
		NRDOSE = NRDMAX - NRDMIN;
		// "-------------------------------------------------------------"
		// PULSE HEIGHT=============
		if (EGS4Macro.IFULL == IFULL_PULSE_HEIGHT_DISTRIBUTION) {
			// OUTPUT;(/' INPUT FOR PULSE HEIGHT DISTRIBUTION'/);
			// "INITIALIZE FLAGS TO NO PULSE HEIGHT DISTRIBUTION IN EACH REGION"
			for (int J = 1; J <= EGS4Geom.nreg; J++) {
				IPHR[J - 1] = 0;
			}
		}
		if (EGS4Macro.IFULL == IFULL_PULSE_HEIGHT_DISTRIBUTION) {
			for (int I = 1; I <= nREGSOLV; I++)
				if (REGSVOL[I - 1] < 0 || REGSVOL[I - 1] > EGS4Geom.nreg) {
					REGSVOL[I - 1] = EGS4Geom.nreg;
				}// default
			if (SLOTE < SLOTE_MIN || SLOTE > SLOTE_MAX) {
				SLOTE = SLOTE_DEFAULT;
			}// default
			if (DELTAE < DELTAE_MIN || DELTAE > DELTAE_MAX) {
				DELTAE = DELTAE_DEFAULT;
			}// default

			for (int J = 1; J <= nREGSOLV; J++) {
				// REGNUM=VALUE(NUM_REGSVOL,J);
				IPHR[REGSVOL[J - 1] - 1] = 1;
				// OUTPUT REGNUM,MED(REGNUM); (/T10,' REGION',I4,' HAS
				// MEDIUM',I3);
				EGS4.seqStr = " 	 REGION" + EGS4.format(REGSVOL[J - 1], 4)
						+ "  HAS MEDIUM"
						+ EGS4.format(EGS4.MED[REGSVOL[J - 1] - 1], 3);
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

			}

			if (SLOTE > 0.0) {
				// OUTPUT SLOTE,SLOTE*dble($EBIN);
				// (/' EQUAL BINS OF',F10.4,' MeV WILL COVER UP TO',F10.3,'
				// MeV');
				EGS4.seqStr = " EQUAL BINS OF" + EGS4.format(SLOTE, 10, true)
						+ " MeV WILL COVER UP TO"
						+ EGS4.format(SLOTE * EGS4SrcEns.dble($EBIN), 10, true)
						+ " MeV";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				// "NOTE THAT LATER, WHEN WE KNOW THE MAXIMUM ENERGY IN THE"
				// "INPUT SPECTRUM, WE WILL INCREASE SLOTE TO MAKE SURE IT WORKS"
			} else {// "SLOTE <= 0.0 means we input energy bins"

				// IVAL=IVAL+1;
				// NUM_TOPEBIN=IVAL;
				// VALUES_SOUGHT(IVAL)='TOPS OF ENERGY BINS';
				// TYPE(IVAL)=1;
				// VALUE_MIN(IVAL)=1.e-20;
				// VALUE_MAX(IVAL)=100000.0;
				// DEFAULT(IVAL)=1.25;
				for (int I = 1; I <= nTOPEBIN; I++)
					if (TOPEBIN[I - 1] < TOPEBIN_MIN
							|| TOPEBIN[I - 1] > TOPEBIN_MAX) {
						TOPEBIN[I - 1] = TOPEBIN_DEFAULT;
					}// default

				// $GET_INPUT(NUM_TOPEBIN);
				// DO J=1, NVALUE(NUM_TOPEBIN) [BINTOP(J)=VALUE(NUM_TOPEBIN,J);]
				for (int J = 1; J <= nTOPEBIN; J++) {
					BINTOP[J - 1] = TOPEBIN[J - 1];
				}
				MAXBIN = nTOPEBIN;// NVALUE(NUM_TOPEBIN);

				if (nTOPEBIN > $EBIN + 1) // (NVALUE(NUM_TOPEBIN) > $EBIN+1)
				{
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " Tried to use more than max number of energy bins";
					printSequence(EGS4.seqStr);
					EGS4.seqStr = " Either increase $EBIN or reduce number of bins";
					printSequence(EGS4.seqStr);

					return;

					// OUTPUT $EBIN;
					// (/' ****Tried to use more than max number of energy bins
					// ',I5/
					// ' Either increase $EBIN or reduce number of bins');
					// STOP;
				}
				// OUTPUT
				// NVALUE(NUM_TOPEBIN),(BINTOP(J),J=1,NVALUE(NUM_TOPEBIN));
				// (/' READ A TOTAL OF',I4,' ENERGY BINS'/ (5F15.4) );
			}

		}

		if (EGS4Macro.IFULL == 3 && EGS4SrcEns.iqin != 0) {// "only allow option 3 for pure photon beam"
		// "                           could modify for full phase space if needed"
		// OUTPUT;
		// (// 1x,70('@')/' Changed IFULL to 0 from 3 since photon beam not
		// input');
			EGS4Macro.IFULL = 0;
			EGS4.seqStr = "  Changed IFULL to 0 from 3 since photon beam not input";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		double EPHTOP = 0.0;
		if (EGS4Macro.IFULL == IFULL_PULSE_HEIGHT_DISTRIBUTION) {
			// "DO SOME CHECKS ON PULSE HEIGHT DISTRIBUTION BINS"
			if (EGS4SrcEns.iqin == 1) {
				EPHTOP = EGS4SrcEns.ein + 1.022;// "INCLUDE ANNIHILATION FOR POSITRONS IN"
			} else {
				EPHTOP = EGS4SrcEns.ein;
			}
			if (SLOTE > 0.0) {
				// UNTIL(SLOTE*dble($EBIN).GT.1.05*EPHTOP)[
				// SLOTE=SLOTE*2.;
				// OUTPUT SLOTE;
				// (/' ***HAVE DOUBLED SLOTE TO',F12.4,
				// ' MeV TO REACH MAXIMUM INPUT ENERGY');
				// ]
				while (true) {
					if (SLOTE * EGS4SrcEns.dble($EBIN) > 1.05 * EPHTOP)
						break;
					SLOTE = SLOTE * 2.;
					EGS4.seqStr = "  HAVE DOUBLED SLOTE TO"
							+ EGS4.format(SLOTE, 12, true)
							+ " MeV TO REACH MAXIMUM INPUT ENERGY";
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
				}
				Double dbb = new Double(1.05 * EPHTOP / SLOTE + 1.0);
				MAXBIN = dbb.intValue();// IFIX(1.05*EPHTOP/SLOTE+0.999);
				EGS4.seqStr = " MAXBIN SET TO" + EGS4.format(MAXBIN, 3)
						+ " TO COVER SPECTRUM";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				// OUTPUT MAXBIN;(/' MAXBIN SET TO',I3,' TO COVER SPECTRUM'/);
			}// "END OF SLOTE>0 BLOCK"
			else {
				if (BINTOP[MAXBIN - 1] <= EPHTOP) {
					EGS4.seqStr = " CHANGED BINTOP(" + EGS4.format(MAXBIN, 3)
							+ "-1)" + " TO REACH"
							+ EGS4.format(EPHTOP, 10, true) + " MeV";
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);

					// OUTPUT MAXBIN,EPHTOP;
					// (/' ***CHANGED BINTOP(',I3,') TO REACH', F10.3,'
					// MeV***');
					BINTOP[MAXBIN - 1] = EPHTOP;
				}
			}// "END SLOTE<=0 BLOCK"

			// "SET UP ENERGY BINS FOR PEAK EFFICIENCY CALCULATION"
			// "THIS ONLY MAKES SENSE IF BEAM IS MONOENERGETIC"
			// DFEN(1,2)=EPHTOP-DELTAE;DFEN(2,2)=DFEN(1,2)-0.511;
			// DFEN(3,2)=DFEN(1,2)-1.022;DFEN(4,2)=0.511-DELTAE;
			DFEN[0][1] = EPHTOP - DELTAE;
			DFEN[1][1] = DFEN[0][1] - 0.511;
			DFEN[2][1] = DFEN[0][1] - 1.022;
			DFEN[3][1] = 0.511 - DELTAE;
			// "I.E. WE SET LOWER ENERGIES FOR FULL ENERGY, SINGLE ESCAPE, DOUBLE"
			// "ESCAPE AND THE 511 KEV LINE"
			for (int IPK = 1; IPK <= 4; IPK++) {
				DFEN[IPK - 1][0] = DFEN[IPK - 1][1] - DELTAE;
				DFEN[IPK - 1][2] = DFEN[IPK - 1][1] + 2 * DELTAE;
				DFEN[IPK - 1][3] = DFEN[IPK - 1][2] + DELTAE;
			}

		}

		// =========================

		if (!setEcutRegion) {
			if (EGS4.ECUT[0] < ECUT_MIN || EGS4.ECUT[0] > ECUT_MAX)
				EGS4.ECUT[0] = EGS4.$GLOBAL_ECUT;
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.ECUT[I - 1] = EGS4.ECUT[0];
			}

		}
		if (!setPcutRegion) {
			if (EGS4.PCUT[0] < PCUT_MIN || EGS4.PCUT[0] > PCUT_MAX)
				EGS4.PCUT[0] = EGS4.$GLOBAL_PCUT;
			// "Now set ecut and pcut to the values input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.PCUT[I - 1] = EGS4.PCUT[0];
			}
		}
		if (!setSmaxRegion) {
			if (EGS4.SMAXIR[0] < SMAXIR_MIN || EGS4.SMAXIR[0] > SMAXIR_MAX)
				EGS4.SMAXIR[0] = EGS4.$MAX_SMAX;
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.SMAXIR[I - 1] = EGS4.SMAXIR[0];
			}
		}

		if (setEcutRegion) {
			for (int i = 1; i <= nEcut; i++) {
				if (Ecut[i - 1] < ECUT_MIN || Ecut[i - 1] > ECUT_MAX)
					Ecut[i - 1] = EGS4.$GLOBAL_ECUT;
				if (startEcutRegion[i - 1] < 1
						|| startEcutRegion[i - 1] > EGS4.$MXREG)
					startEcutRegion[i - 1] = 1;
				if (stopEcutRegion[i - 1] < 1
						|| stopEcutRegion[i - 1] > EGS4.$MXREG)
					stopEcutRegion[i - 1] = 1;

				for (int j = startEcutRegion[i - 1]; j <= stopEcutRegion[i - 1]; j++)
					EGS4.ECUT[j - 1] = Ecut[i - 1];
			}
		}
		if (setPcutRegion) {
			for (int i = 1; i <= nPcut; i++) {
				if (Pcut[i - 1] < PCUT_MIN || Pcut[i - 1] > PCUT_MAX)
					Pcut[i - 1] = EGS4.$GLOBAL_PCUT;
				if (startPcutRegion[i - 1] < 1
						|| startPcutRegion[i - 1] > EGS4.$MXREG)
					startPcutRegion[i - 1] = 1;
				if (stopPcutRegion[i - 1] < 1
						|| stopPcutRegion[i - 1] > EGS4.$MXREG)
					stopPcutRegion[i - 1] = 1;

				for (int j = startPcutRegion[i - 1]; j <= stopPcutRegion[i - 1]; j++)
					EGS4.PCUT[j - 1] = Pcut[i - 1];
			}
		}
		if (setSmaxRegion) {
			for (int i = 1; i <= nSmax; i++) {
				if (Smax[i - 1] < SMAXIR_MIN || Smax[i - 1] > SMAXIR_MAX)
					Smax[i - 1] = EGS4.$MAX_SMAX;
				if (startSmaxRegion[i - 1] < 1
						|| startSmaxRegion[i - 1] > EGS4.$MXREG)
					startSmaxRegion[i - 1] = 1;
				if (stopSmaxRegion[i - 1] < 1
						|| stopSmaxRegion[i - 1] > EGS4.$MXREG)
					stopSmaxRegion[i - 1] = 1;

				for (int j = startSmaxRegion[i - 1]; j <= stopSmaxRegion[i - 1]; j++)
					EGS4.SMAXIR[j - 1] = Smax[i - 1];
			}
		}

		if (!setIncohRegion) {
			EGS4.ibcmp[0] = incoh;
			if (EGS4.ibcmp[0] < 0 || EGS4.ibcmp[0] > 3) {
				EGS4.ibcmp[0] = EGS4.$IBCMP_DEFAULT;
				incoh = EGS4.$IBCMP_DEFAULT;
			}
			// "Now set ibcmp for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.ibcmp[I - 1] = EGS4.ibcmp[0];
			}
		}
		if (setIncohRegion) {
			for (int i = 1; i <= nIncoh; i++) {
				if (Incoh[i - 1] < 0 || Incoh[i - 1] > 1)// only 2 values
					Incoh[i - 1] = EGS4.$IBCMP_DEFAULT;
				if (startIncohRegion[i - 1] < 1
						|| startIncohRegion[i - 1] > EGS4.$MXREG)
					startIncohRegion[i - 1] = 1;
				if (stopIncohRegion[i - 1] < 1
						|| stopIncohRegion[i - 1] > EGS4.$MXREG)
					stopIncohRegion[i - 1] = 1;

				for (int j = startIncohRegion[i - 1]; j <= stopIncohRegion[i - 1]; j++)
					EGS4.ibcmp[j - 1] = Incoh[i - 1];
			}
		}
		if (!setCohRegion) {
			EGS4.IRAYLR[0] = coh;
			if (EGS4.IRAYLR[0] < 0 || EGS4.IRAYLR[0] > 3) {
				EGS4.IRAYLR[0] = EGS4.$IRAYLR_DEFAULT;
				coh = EGS4.$IRAYLR_DEFAULT;
			}
			// "Now set iraylr for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.IRAYLR[I - 1] = EGS4.IRAYLR[0];
			}
		}
		if (setCohRegion) {
			for (int i = 1; i <= nCoh; i++) {
				if (Coh[i - 1] < 0 || Coh[i - 1] > 1)// only 2 values
					Coh[i - 1] = EGS4.$IRAYLR_DEFAULT;
				if (startCohRegion[i - 1] < 1
						|| startCohRegion[i - 1] > EGS4.$MXREG)
					startCohRegion[i - 1] = 1;
				if (stopCohRegion[i - 1] < 1
						|| stopCohRegion[i - 1] > EGS4.$MXREG)
					stopCohRegion[i - 1] = 1;

				for (int j = startCohRegion[i - 1]; j <= stopCohRegion[i - 1]; j++)
					EGS4.IRAYLR[j - 1] = Coh[i - 1];
			}
		}
		if (!setRelaxRegion) {
			EGS4.iedgfl[0] = relax;
			if (EGS4.iedgfl[0] < 0 || EGS4.iedgfl[0] > 3) {
				EGS4.iedgfl[0] = EGS4.$IEDGFL_DEFAULT;
				relax = EGS4.$IEDGFL_DEFAULT;
			}
			// "Now set iedgfl for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.iedgfl[I - 1] = EGS4.iedgfl[0];
			}
		}
		if (setRelaxRegion) {
			for (int i = 1; i <= nRelax; i++) {
				if (Relax[i - 1] < 0 || Relax[i - 1] > 1)// only 2 values
					Relax[i - 1] = EGS4.$IEDGFL_DEFAULT;
				if (startRelaxRegion[i - 1] < 1
						|| startRelaxRegion[i - 1] > EGS4.$MXREG)
					startRelaxRegion[i - 1] = 1;
				if (stopRelaxRegion[i - 1] < 1
						|| stopRelaxRegion[i - 1] > EGS4.$MXREG)
					stopRelaxRegion[i - 1] = 1;

				for (int j = startRelaxRegion[i - 1]; j <= stopRelaxRegion[i - 1]; j++)
					EGS4.iedgfl[j - 1] = Relax[i - 1];
			}
		}
		if (!setPeRegion) {
			EGS4.iphter[0] = pe;
			if (EGS4.iphter[0] < 0 || EGS4.iphter[0] > 3) {
				EGS4.iphter[0] = EGS4.$IPHTER_DEFAULT;
				pe = EGS4.$IPHTER_DEFAULT;
			}
			// "Now set iphter for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.iphter[I - 1] = EGS4.iphter[0];
			}
		}
		if (setPeRegion) {
			for (int i = 1; i <= nPe; i++) {
				if (Pe[i - 1] < 0 || Pe[i - 1] > 1)// only 2 values
					Pe[i - 1] = EGS4.$IPHTER_DEFAULT;
				if (startPeRegion[i - 1] < 1
						|| startPeRegion[i - 1] > EGS4.$MXREG)
					startPeRegion[i - 1] = 1;
				if (stopPeRegion[i - 1] < 1
						|| stopPeRegion[i - 1] > EGS4.$MXREG)
					stopPeRegion[i - 1] = 1;

				for (int j = startPeRegion[i - 1]; j <= stopPeRegion[i - 1]; j++)
					EGS4.iphter[j - 1] = Pe[i - 1];
			}
		}
		if ((EGS4.ibrdst < 0) || (EGS4.ibrdst > 1)) {
			EGS4.ibrdst = EGS4.$IBRDST_DEFAULT;
		}
		if ((EGS4.iprdst < 0) || (EGS4.iprdst > 4)) {
			EGS4.iprdst = EGS4.$IPRDST_DEFAULT;
		}
		if ((EGS4.estepe < ESTEPE_MIN) || (EGS4.estepe > ESTEPE_MAX)) {
			EGS4.estepe = EGS4.$MAX_ELOSS; // "$MAX-ELOSS is defined in egsnrc.macros at 0.25"
		}
		if ((EGS4.ximax < XIMAX_MIN) || (EGS4.ximax > XIMAX_MAX)) {
			EGS4.ximax = EGS4.$EXACT_BCA_XIMAX; // "$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
		}
		if ((EGS4.bca_algorithm < 0) || (EGS4.bca_algorithm > 1)) {
			EGS4.bca_algorithm = EGS4.$BCA_ALGORITHM_DEFAULT;
		}
		if ((EGS4.skindepth_for_bca < Skindepth_MIN)
				|| (EGS4.skindepth_for_bca > Skindepth_MAX))
			EGS4.skindepth_for_bca = EGS4.$SKIN_DEPTH_FOR_BCA;
		if (EGS4.bca_algorithm == BCA_EXACT) {
			if (EGS4.skindepth_for_bca <= 0.0)
				EGS4.skindepth_for_bca = EGS4.$SKIN_DEPTH_FOR_BCA;
		}
		if ((EGS4.transport_algorithm < 0) || (EGS4.transport_algorithm > 1)) {
			EGS4.transport_algorithm = EGS4.$TRANSPORT_ALGORITHM_DEFAULT;
		}
		if (ispin == spin_ON)
			EGS4.spin_effects = true;// so=>
		else
			EGS4.spin_effects = false;// so=>
		if ((EGS4.ibr_nist < 0) || (EGS4.ibr_nist > 1)) {
			EGS4.ibr_nist = EGS4.$IBR_NIST_DEFAULT;
		}
		if ((EGS4.pair_nrc < 0) || (EGS4.pair_nrc > 1)) {
			EGS4.pair_nrc = EGS4.$PAIR_NRC_DEFAULT;
		}
		if ((EGS4.eii_flag < 0) || (EGS4.eii_flag > 1)) {
			EGS4.eii_flag = eii_OFF;// default
		}
		if ((EGS4.itriplet < 0) || (EGS4.itriplet > 1)) {
			EGS4.itriplet = EGS4.$TRIPLET_DEFAULT;
		}
		if ((EGS4.radc_flag < 0) || (EGS4.radc_flag > 1)) {
			EGS4.radc_flag = radc_OFF;// default
		}
		// =================================
		if ((EGS4.nbr_split > $MAXBRSPLIT))
			EGS4.nbr_split = $MAXBRSPLIT;
		// CHARGED PARTICLE RUSSIAN ROULETTE
		if (EGS4.i_play_RR == 1) {
			EGS4.prob_RR = 1. / EGS4SrcEns.dble(EGS4.nbr_split);
		} else {
			EGS4.prob_RR = 1.;
		}
		if (EGS4Macro.IFULL == 2 && EGS4.nbr_split > 1) {// "cannot have this"
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " bremsstrahlung splitting on.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(1,'(//'' ****WARNING****''/
			// ''You cannot calculate a pulse height distribution with''/
			// '' bremsstrahlung splitting on. Will run simulation with''/
			// '' IFULL= dose and stoppers''//)');
			EGS4Macro.IFULL = 0;
		}

		if (EGS4Macro.irejct < 0 || EGS4Macro.irejct > 1) {
			EGS4Macro.irejct = irejct_OFF;
		}// default
		if (ESAVEIN < ESAVEIN_MIN || ESAVEIN > ESAVEIN_MAX) {
			ESAVEIN = ESAVEIN_DEFAULT;
		}
		if (EGS4Macro.cs_enhance < cs_enhance_MIN
				|| EGS4Macro.cs_enhance > cs_enhance_MAX) {
			EGS4Macro.cs_enhance = cs_enhance_DEFAULT;
		}
		// if( EGS4Macro.cs_enhance > 1. ) { EGS4Macro.use_enhance = true; }
		// else { EGS4Macro.use_enhance = false; }

		EGS4.seqStr = " *** Reading variance reduction inputs ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Range rejection is On(1) or Off(0):" + EGS4Macro.irejct;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (EGS4Macro.irejct > 0) {
			EGS4.seqStr = " ESAVEIN cutoff value(total) for range rejection:"
					+ EGS4.format(ESAVEIN, 10, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			if (ESAVEIN == 0.0) {
				EGS4.seqStr = " WARNING: Have asked for range rejection but left ESAVEIN=0.0!!";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
				EGS4.seqStr = " There will be no range rejection!";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}
			for (int i = 1; i <= EGS4Geom.nreg; i++) {
				EGS4.i_do_rr[i - 1] = 1;
				EGS4.e_max_rr[i - 1] = ESAVEIN;
			}
			// "note  e_max_r is total energy"
			// "above two arrays needed for each region for EGSnrc RANGE-DISCARD macro"
		}
		if (RRZ < RRDEPTH_MIN || RRZ > RRDEPTH_MAX) {
			RRZ = RRDEPTH_DEFAULT;
		}
		if (RRCUT < RRFRACTION_MIN || RRCUT > RRFRACTION_MAX) {
			RRCUT = RRFRACTION_DEFAULT;
		}
		if (EGS4Macro.CEXPTR < EXPC_MIN || EGS4Macro.CEXPTR > EXPC_MAX) {
			EGS4Macro.CEXPTR = EXPC_DEFAULT;
		}
		RUSROU = false;
		if (RRZ + RRCUT != 0.0)
			RUSROU = true;
		if (RUSROU) {
			EGS4.seqStr = " RUSSIAN ROULETTE WILL BE PLAYED";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " RUSSIAN ROULETTE PLANE:"
					+ EGS4.format(RRZ, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " SURVIVAL PROBABILITY:"
					+ EGS4.format(RRCUT, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			EGS4.seqStr = " RUSSIAN ROULETTE WILL NOT BE PLAYED";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (EGS4Macro.CEXPTR == 0.) {
			EGS4.seqStr = " NO PATHLENGTH BIASING TO BE DONE";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			EGS4.seqStr = " CEXPTR PARAMATER:"
					+ EGS4.format(EGS4Macro.CEXPTR, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (EGS4Macro.IFULL == 2 && RUSROU) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Russian Roulette on.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
		}
		if (EGS4Macro.IFULL == 2 && EGS4Macro.CEXPTR != 0) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " pathlength biasing.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
		}

		if (IFARCE < 0 || IFARCE > 1) {
			IFARCE = IFARCE_OFF;
		}// default
		NFMIN_MAX = EGS4Geom.nreg;
		NFMAX_MAX = EGS4Geom.nreg + 1;
		if (EGS4Macro.NFMIN < NFMIN_MIN || EGS4Macro.NFMIN > NFMIN_MAX) {
			EGS4Macro.NFMIN = NFMIN_DEFAULT;
		}
		if (EGS4Macro.NFMAX < NFMAX_MIN || EGS4Macro.NFMAX > NFMAX_MAX) {
			EGS4Macro.NFMAX = NFMAX_DEFAULT;
		}

		if (IFARCE == 0) {
			EGS4Macro.IFORCE = 0;
			EGS4Macro.NFMIN = 0;
			EGS4Macro.NFMAX = 0;

			EGS4.seqStr = " NO INTERACTION FORCING IS IN EFFECT";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else if (IFARCE == 1) {
			EGS4Macro.IFORCE = 1;
			if (EGS4Macro.NFMAX < EGS4Macro.NFMIN)
				EGS4Macro.NFMAX = EGS4Macro.NFMIN;

			EGS4.seqStr = " FORCED PHOTON INTERACTIONS IN EFFECT FROM "
					+ EGS4Macro.NFMIN + " TO " + EGS4Macro.NFMAX
					+ " INTERACTIONS ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		if (EGS4Macro.IFULL == 2
				&& (EGS4Macro.NFMAX > EGS4Macro.NFMIN || EGS4Macro.NFMIN > 1)) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " more than 1 interaction forced or if the forced interaction";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " is > the first interaction.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
		}
		ics_enhance = 0;
		for (int jj = 1; jj <= EGS4Geom.nreg; jj++) {
			EGS4Macro.iefl[jj - 1] = 0;
		}
		for (int ii = 1; ii <= nENHREG; ii++) {
			int ics_start = NENHLO[ii - 1];
			int ics_stop = NENHHI[ii - 1];
			for (int jj = ics_start; jj <= ics_stop; jj++) {
				EGS4Macro.iefl[jj - 1] = 1;
			}
		}

		int COUNT = 0;
		for (int jj = 2; jj <= EGS4Geom.nreg; jj++) {
			COUNT = COUNT + EGS4Macro.iefl[jj - 1];
		}
		// "We don't care about region 1 since outside geometry"
		if (COUNT > 0 && (EGS4Macro.cs_enhance > 1.0001)) {// "there is enhancement somewhere"
			EGS4.seqStr = " Cross section enhancement in regions ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			for (int jj = 1; jj <= EGS4Geom.nreg; jj++) {
				EGS4.seqStr = " 		" + EGS4.format(EGS4Macro.iefl[jj - 1], 4);
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}
			EGS4.seqStr = " Cross section enhancement factor: "
					+ EGS4.format(EGS4Macro.cs_enhance, 6, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// OUTPUT; ('Cross section enhancement in regions ');
			// OUTPUT (iefl(jj),jj=1,NREG); (20 I4);
			// OUTPUT cs_enhance;
			// ( ' Cross section enhancement factor: ',T60,F6.1/);
			ics_enhance = 1;
		} else {
			EGS4.seqStr = " No cross section enhancement";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			ics_enhance = 0;
		}

		if (EGS4Macro.IFULL == 2 && ics_enhance == 1) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " cross section enhancement.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4Macro.IFULL = 0;
		}
		/*
		 * EGS4Macro.IQINC=EGS4SrcEns.iqin;
		 * //"NEEDED TO TURN OFF FASTSTEP FOR INCIDENT ELECTRONS"
		 * //"WHEN FORCING INTERACTIONS" phsplitt_MAX=EGS4.$MXSTACK-2;
		 * if(phsplitt<phsplitt_MIN ||
		 * phsplitt>phsplitt_MAX){phsplitt=phsplitt_DEFAULT;} EGS4Macro.n_split
		 * = phsplitt; if( EGS4Macro.n_split > 1 ) { if( EGS4Macro.IFULL > 1 ) {
		 * EGS4
		 * .seqStr=" IGNORING INPUT: Photon splitting only for ifull = 0,1! ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.n_split = 1; } else {
		 * EGS4.seqStr=" Calculation with photon splitting, n_split = "
		 * +EGS4Macro.n_split; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.iifano = IFANO; } }
		 */
		/*
		 * if( EGS4Macro.use_enhance ) {
		 * EGS4.seqStr="  Calculation with CS enhancement  "; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr); EGS4.seqStr=
		 * "  photon forcing, exp. transform, etc. input will be ignored, IFULL will be set to 1! (i.e. Ascat and Aatt) !"
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr="  Using cs_enhance = "+EGS4Macro.cs_enhance;
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.IFULL = 1; EGS4Macro.IFORCE=0;
		 * EGS4Macro.NFMIN=0;EGS4Macro.NFMAX=0; EGS4Macro.n_split = 1; }
		 * 
		 * //"if we are calculating Aatt and Ascatt, we must have n_split>1,"
		 * //"cs_enhance on, or photon forcing on" if(EGS4Macro.IFULL==1 &&
		 * EGS4Macro.cs_enhance<=1 && EGS4Macro.IFORCE==0 && EGS4Macro.n_split
		 * <=1) { EGS4Macro.IFORCE=1; EGS4Macro.NFMIN=1; EGS4Macro.NFMAX=1;
		 * 
		 * EGS4.seqStr=
		 * " If you are calculating Aatt, Ascat (IFULL=1), you must be using";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); EGS4.seqStr=
		 * " photon forcing, photon splitting, or cross-section enhancement.";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); EGS4.seqStr=
		 * " Currently, none of these are being used.  Will continue run with";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); EGS4.seqStr=
		 * " photon forcing on and one interaction forced (NFMIN=NFMAX=1).";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); }
		 * 
		 * if(IOOPTN==1 && (EGS4Macro.n_split>1 || EGS4Macro.use_enhance)) {
		 * IOOPTN=0; EGS4.seqStr=
		 * " You cannot have cross-section enhancement or photon splitting on ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); EGS4.seqStr=
		 * " and still output detailed results for each cavity zone.  IOOPTN reset to 0."
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr); }
		 */
	}

	/**
	 * Initialize random number generator.
	 */
	private void init_random_generator() {
		if (EGS4.ranluxB) {
			if (jrng1 < RANLUX_LEVEL_MIN || jrng1 > RANLUX_LEVEL_MAX) {
				jrng1 = RANLUX_LEVEL_DEFAULT;
			}
			if (jrng2 < RANLUX_SEED_MIN || jrng2 > RANLUX_SEED_MAX) {
				jrng2 = RANLUX_SEED_DEFAULT;
			}
		} else {
			if (jrng1 < RANMAR_SEED_MIN || jrng1 > RANMAR_SEED_MAX) {
				jrng1 = RANMAR_SEED_DEFAULT;
			}
			if (jrng2 < RANMAR_SEED_MIN || jrng2 > RANMAR_SEED_MAX) {
				jrng2 = RANMAR_SEED_DEFAULT;
			}
		}
		// IF( i_parallel > 0 ) jrng2 = jrng2 - 1 + i_parallel;//NO PARALLEL JOB
		// USED!!!
		// $INITIALIZE RNG USING jrng1 AND jrng2;
		if (EGS4.ranluxB) {
			EGS4.init_ranlux(jrng1, jrng2);
			EGS4.ranlux(EGS4.rng_array);
			EGS4.rng_seed = 1;
		} else {
			EGS4.ixx = jrng1;
			EGS4.jxx = jrng2;
			EGS4.init_ranmar();
		}
	}

	/**
	 * Print input summary.
	 */
	private void ISUMRY() {
		String s = "";
		String s1 = "";
		int ll = 0;
		int ioff = 0;
		EGS4.seqStr = "************************INPUT SUMMARY:*****************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   Electron/Photon transport parameter";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Photon transport cutoff(MeV)";
		ioff = 32;
		if (setPcutRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.PCUT[0] > 1.e-4) {
				ioff = 30;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.PCUT[0], 8, true);
			} else {
				s = s + EGS4.format("", ioff) + "AP(medium)";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Pair angular sampling";
		ioff = 39;
		if (EGS4.iprdst == pair_ang_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.iprdst == pair_ang_SIMPLE) {
			s = s + EGS4.format("", ioff) + "SIMPLE";
		} else if (EGS4.iprdst == pair_ang_KM) {
			s = s + EGS4.format("", ioff) + "KM";
		} else if (EGS4.iprdst == pair_ang_UNIFORM) {
			s = s + EGS4.format("", ioff) + "UNIFORM";
		} else if (EGS4.iprdst == pair_ang_BLEND) {
			s = s + EGS4.format("", ioff) + "BLEND";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Pair cross sections";
		ioff = 41;
		if (EGS4.pair_nrc == pair_cross_BH) {
			s = s + EGS4.format("", ioff) + "BH";
		} else if (EGS4.pair_nrc == pair_cross_NRC) {
			s = s + EGS4.format("", ioff) + "NRC";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Triplet production";
		ioff = 42;
		if (EGS4.itriplet == triplet_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.itriplet == triplet_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bound Compton scattering";
		ioff = 36;
		if (incoh == incoh_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (incoh == incoh_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (incoh == incoh_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (incoh == incoh_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Radiative Compton corrections";
		ioff = 31;
		if (EGS4.radc_flag == radc_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.radc_flag == radc_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Rayleigh scattering";
		ioff = 41;
		if (coh == coh_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (coh == coh_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (coh == coh_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (coh == coh_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Atomic relaxations";
		ioff = 42;
		if (relax == relax_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (relax == relax_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (relax == relax_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (relax == relax_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Photoelectron angular sampling";
		ioff = 30;
		if (pe == pe_ang_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (pe == pe_ang_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (pe == pe_ang_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (pe == pe_ang_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron transport cutoff(MeV)";
		ioff = 30;
		if (setEcutRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.ECUT[0] > 1.e-4) {
				ioff = 28;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.ECUT[0], 8, true);
			} else {
				s = s + EGS4.format("", ioff) + "AE(medium)";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bremsstrahlung angular sampling";
		ioff = 29;
		if (EGS4.ibrdst == brems_ang_SIMPLE) {
			s = s + EGS4.format("", ioff) + "SIMPLE";
		} else if (EGS4.ibrdst == brems_ang_KM) {
			s = s + EGS4.format("", ioff) + "KM";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bremsstrahlung cross sections";
		ioff = 31;
		if (EGS4.ibr_nist == brems_cross_BH) {
			s = s + EGS4.format("", ioff) + "BH";
		} else if (EGS4.ibr_nist == brems_cross_NIST) {
			s = s + EGS4.format("", ioff) + "NIST";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Spin effects";
		ioff = 48;
		if (EGS4.spin_effects) {
			s = s + EGS4.format("", ioff) + "ON";
		} else {
			s = s + EGS4.format("", ioff) + "OFF";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron Impact Ionization";
		ioff = 34;
		if (EGS4.eii_flag == eii_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.eii_flag == eii_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maxium electron step in cm (SMAX)";
		ioff = 27;
		if (setSmaxRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.SMAXIR[0] > 1.e-4) {
				ioff = 26;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.SMAXIR[0], 6, false);
			} else {
				s = s + EGS4.format("", ioff) + "Restriction is off";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maximum fractional energy loss/step (ESTEPE)";
		ioff = 16;
		s = s + EGS4.format("", ioff) + EGS4.format(EGS4.estepe, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maximum 1st elastic moment/step (XIMAX)";
		ioff = 21;
		s = s + EGS4.format("", ioff) + EGS4.format(EGS4.ximax, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Boundary crossing algorithm";
		ioff = 33;
		if (EGS4.bca_algorithm == BCA_EXACT) {
			s = s + EGS4.format("", ioff) + "EXACT";
		} else if (EGS4.bca_algorithm == BCA_PRESTA_I) {
			s = s + EGS4.format("", ioff) + "PRESTA I";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Skin-depth for boundary crossing (MFP)";
		ioff = 22;
		s = s + EGS4.format("", ioff)
				+ EGS4.format(EGS4.skindepth_for_bca, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron-step algorithm";
		ioff = 37;
		if (EGS4.transport_algorithm == estep_alg_PRESTA_II) {
			s = s + EGS4.format("", ioff) + "PRESTA II";
		} else if (EGS4.transport_algorithm == estep_alg_PRESTA_I) {
			s = s + EGS4.format("", ioff) + "PRESTA I";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   MONTE CARLO, TRANSPORT, AND SCATTER CONTROLS";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Max. # of histories to RUN";
		ll = s.length();
		ll = 54 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(NCASE, 12);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Max # of histories to ANALYZE";
		ll = s.length();
		ll = 54 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(NCASET, 12);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Incident Charge";// EGS4.format("",20)+"Incident Charge";
		s1 = "";
		ll = s.length();
		ll = 68 - ll;
		if (EGS4SrcEns.iqin == 0)
			s1 = "photons";
		if (EGS4SrcEns.iqin == -1)
			s1 = "electrons";
		if (EGS4SrcEns.iqin == 1)
			s1 = "positrons";
		if (EGS4SrcEns.iqin == 2)
			s1 = "all";
		if (EGS4SrcEns.iqin == 3)
			s1 = "e- & e+";
		s = s + EGS4.format(s1, ll);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21
				&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {
			s = " Incident kinetic energy:";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4SrcEns.ein, 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			EGS4SrcEns.ENSRCO();
		}
		boolean gotopast = false;
		for (int I = 2; I <= EGS4Geom.nreg; I++) {
			gotopast = false;
			if ((EGS4.ECUT[I - 1] != EGS4.ECUT[1])
					|| (EGS4.PCUT[I - 1] != EGS4.PCUT[1])) {
				// "we failed at least one test, so this means there really are"
				// "varying ECUTs and these will be printed in the grid if we want them"
				// "print the first 12 ECUT & PCUT just to be sure"
				int j = Math.min(12, EGS4Geom.nreg);

				for (int JJ = 2; JJ <= j; JJ++) {
					EGS4.seqStr = "First ECUTs (MeV): "
							+ EGS4.format(EGS4.ECUT[JJ - 1], 12, true) + "  ,"
							+ "First PCUTs (MeV): "
							+ EGS4.format(EGS4.PCUT[JJ - 1], 12, true);
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
				}
				// GO TO :past:
				gotopast = true;
				break;
			}
		}
		// "if we get here, they were all the same"
		if (!gotopast) {
			s = " GLOBAL ELECTRON TRANSPORT CUT-OFF";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.ECUT[1], 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " GLOBAL PHOTON TRANSPORT CUT-OFF";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.PCUT[1], 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		// :past:
		if (EGS4Macro.IFORCE != 0) {
			s = " Min/max photon step forced";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.NFMIN, 6) + "/"
					+ EGS4.format(EGS4Macro.NFMAX, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			s = " Photon force interaction switch";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (ics_enhance == 1) {
			// WRITE(IOUT,235) cs_enhance,(j,iefl(j),j=2,NREG);
			// 235 FORMAT(T20,'Cross section enhancement factor of',T56,F8.1/
			// T20,'In regions with a 1:'/
			// (T10, 10('(',I3,',',I1,')')));
			s = " Cross section enhancement factor of";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.cs_enhance, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			s = " No cross section enhancement used";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;// +"OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,236);236 FORMAT(T20,'No cross section enhancement
			// used');
		}

		if (EGS4Macro.irejct > 0) {/*
									 * s=
									 * " Range rejection on a region by region basis"
									 * ;
									 * //ll=s.length();ll=56-ll;s=s+EGS4.format
									 * ("",ll); EGS4.seqStr=s; if(EGS4.iprint>1)
									 * printSequence(EGS4.seqStr);
									 * s=" Also globally to region between z="
									 * +EGS4
									 * .format(EGS4Macro.z_cavity_min,10,true)+
									 * " &"
									 * +EGS4.format(EGS4Macro.z_cavity_max,10
									 * ,true)+" cm";
									 * //ll=s.length();ll=56-ll;s=s
									 * +EGS4.format("",ll); EGS4.seqStr=s;
									 * if(EGS4.iprint>1)
									 * printSequence(EGS4.seqStr);
									 * s="      and inside radius="
									 * +EGS4.format(EGS4Macro
									 * .r_cavity_max,12,true)+ " cm";
									 * //ll=s.length
									 * ();ll=56-ll;s=s+EGS4.format("",ll);
									 * EGS4.seqStr=s; if(EGS4.iprint>1)
									 * printSequence(EGS4.seqStr); s=
									 * " Range rejection only for electrons < ESAVEIN="
									 * +EGS4.format(ESAVEIN,10,true)+ " MeV";
									 * //ll
									 * =s.length();ll=56-ll;s=s+EGS4.format(""
									 * ,ll); EGS4.seqStr=s; if(EGS4.iprint>1)
									 * printSequence(EGS4.seqStr);
									 */
			// 242 FORMAT(' ',T20,'RANGE REJECTION SWITCH',T60,'ON'/
			// T20,' Range rejection for energy <',T56,F9.3,' (MeV)'/
			// T20,' Ranges determined internally from stopping powers');
			s = " RANGE REJECTION SWITCH";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  Range rejection for energy  <";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(ESAVEIN, 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  Ranges determined internally from stopping powers";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			s = " RANGE REJECTION SWITCH";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		for (int I = 1; I <= EGS4Geom.nreg; I++) {
			if (EGS4.IRAYLR[I - 1] == 1) {
				// WRITE(IOUT,244);
				s = " RAYLEIGH SCATTERING INCLUDED'";
				// ll=s.length();ll=61-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;// +"ON";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				break;// EXIT;
			}
		}

		// WRITE(IOUT,260) TIMMAX,STATLM;
		s = " Maximum cputime allowed";
		ll = s.length();
		ll = 60 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(TIMMAX, 6, true) + " hrs";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Stats in PEAK REGION objective";
		ll = s.length();
		ll = 61 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(STATLM, 6, true) + " %";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Initial RNG state:";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.SHOW_RNG_STATE();

		if (EGS4.nbr_split == 1) {
			s = " Bremsstrahlung splitting";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,312);
		} else {
			// WRITE(IOUT,313) nbr_split;
			s = " Bremsstrahlung splitting";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Initially, each bremsstrahlung photon split into ";
			EGS4.seqStr = s + EGS4.format(EGS4.nbr_split, 3) + " photons";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		if (EGS4.i_play_RR == 0) {
			// WRITE(IOUT,314);
			s = " Charged particle Russian Roulette";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			// WRITE(IOUT,315) EGS4.prob_RR;
			s = " Charged particle Russian Roulette";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " With probability of survival = ";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.prob_RR, 9);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		if (RUSROU) {
			// WRITE(IOUT,265)RRZ,RRCUT;
			s = " RUS ROU FOR PHOTONS CROSSING Z = ";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(RRZ, 10, true) + " cm";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " WITH PROBABILITY OF SURVIVAL:";
			ll = s.length();
			ll = 59 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(RRCUT, 7, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		if (EGS4Macro.CEXPTR != 0) {
			s = " PATHLENGTH EXPONENTIAL TRANSFORMATION";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  VARIABLE FOR FORWARD GOING PHOTNS: ";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.CEXPTR, 10, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		/*
		 * if(IFANO == 1) {
		 * s="  *** REGENERATION REQUESTED (IFANO SET TO 1) ! *** ";
		 * EGS4.seqStr=s; if(EGS4.iprint>1) printSequence(EGS4.seqStr); } else
		 * if(IFANO == 2) { s=
		 * "  *** ELECTRONS SET IN MOTION IN WALL WILL BE ELIMINATED (IFANO SET TO 2) ! *** "
		 * ; EGS4.seqStr=s; if(EGS4.iprint>1) printSequence(EGS4.seqStr); } else
		 * { s="  *** NO REGENERATION REQUESTED (IFANO SET TO 0) ! *** ";
		 * EGS4.seqStr=s; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * }
		 */
		// "EK0=EIN;"
		// "$PRESTA-INPUT-SUMMARY; OUTPUT THE PRESTA INPUT VARIABLES"
		// "taken out input-summary at upgrade to PRESTA-II"
		if (EGS4Macro.ICSDA == 0) {
			// WRITE(IOUT,289);
		} else {
			// WRITE(IOUT,290);
		}
		// "PATCH FOR FLUORESCENT X-RAYS"
		/*
		 * ISUMX=0; DO
		 * JJ=1,NREG[ISUMX=ISUMX+IEDGFL(JJ);"NON-ZERO IF X-RAYS ANYWHERE"]
		 * IF(ISUMX.EQ.0)[WRITE(IOUT,307);]
		 * 
		 * ELSE[
		 * 
		 * WRITE(IOUT,306);
		 * 
		 * start_fluor = .false.; DO jj=1,NREG [ IF( IEDGFL(JJ) > 0 & IEDGFL(JJ)
		 * <= $MXELEMENT ) [ IF( ~start_fluor ) [ start_fluor = .true.; i_start
		 * = jj; ] ELSE [ i_stop = jj; ] ] ELSE [ IF( start_fluor ) [
		 * write(iout,'(24x,i5,a,i5)') i_start,' -- ',i_stop; start_fluor =
		 * .false.; ] ] ] IF( start_fluor ) write(iout,'(24x,i5,a,i5)')
		 * i_start,' -- ',i_stop;
		 * 
		 * ]
		 */
		EK0 = EGS4SrcEns.ein;
		// $PRESTA-INPUT-SUMMARY; "OUTPUT THE PRESTA INPUT VARIABLES"

		// "MATERIAL INPUT SUMMARY"
		// "====================="
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   MATERIAL SUMMARY  " + EGS4.NMED
				+ " MATERIALS USED";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " # MATERIAL  DENSITY(g/cm**3)" + EGS4.format("", 6)
				+ "AE(MeV)" + EGS4.format("", 4) + "AP(MeV)"
				+ EGS4.format("", 9) + "UE(MeV)" + EGS4.format("", 4)
				+ "UP(MeV)";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "  - --------  ----------------" + EGS4.format("", 6)
				+ "-------" + EGS4.format("", 4) + "-------"
				+ EGS4.format("", 9) + "-------" + EGS4.format("", 4)
				+ "-------";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		for (int I = 1; I <= EGS4.NMED; I++) {
			// WRITE(IOUT,310)
			// I,(MEDIA(J,I),J=1,6),RHO(I),AE(I),AP(I),UE(I),UP(I);
			// 310 FORMAT(' ',I1,3X,6A1,4X,1PE10.3,2(7X,0PF9.3,2X,F9.3));
			String meds = "";
			if (EGS4.MEDIA[I - 1].length() > 6) {
				meds = EGS4.MEDIA[I - 1].substring(0, 6);
			} else {
				meds = EGS4.MEDIA[I - 1];
			}
			EGS4.seqStr = "  " + EGS4.format(I,3) + EGS4.format("", 3) + meds
					+ EGS4.format("", 4)
					+ EGS4.format(EGS4.RHO[I - 1], 10, false)
					+ EGS4.format("", 6) + EGS4.format(EGS4.AE[I - 1], 10)
					+ EGS4.format("", 1) + EGS4.format(EGS4.AP[I - 1], 10)
					+ EGS4.format("", 6) + EGS4.format(EGS4.UE[I - 1], 10)
					+ EGS4.format("", 1) + EGS4.format(EGS4.UP[I - 1], 10);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		EGS4Geom.GEOMRZ_ISUMRY();

		EGS4SrcEns.SRCOUT();

		// "       PRINT A GRID OF THE ZONE DEPENDENT VARIABLES"
		// "       ============================================"

		// "Set defaults for non-used variables"
		// CDSTBL(1)='0';CTRTBL(1)='0';CABSRB(1)='0';
		// EGS4Grid.CDSTBL[0]=" ";
		// EGS4Grid.CTRTBL[0]="0";
		EGS4Grid.CABSRB[0] = "0";
		// "For identifying a cavity region"
		// REPLACE {$IRL} WITH {IZ+1+NZ*(IX-1)}
		/*
		 * for(int IRL=2;IRL<= EGS4Geom.NR*EGS4Geom.NZ+1;IRL++) { if
		 * (EGS4Geom.ntrack[IRL-1]==1) { EGS4Grid.CDSTBL[IRL-1]="C"; } else {
		 * EGS4Grid.CDSTBL[IRL-1]=" "; } }
		 */
		// "Make the material grid"
		EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, AMASS, 1, EGS4.ECUT,
				EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);

		if (EGS4Macro.IFULL == 2) {
			// OUTPUT (I,IPHR(I),I=1,NREG);
			// (///' PULSE HEIGHT DISTRIBUTION IS SCORED IN ',
			// 'THOSE REGIONS DENOTED WITH A 1'/(10(I3,'(',I1,'), ')));
			// WRITE(IOUT,311)(I,IPHR(I),I=1,NREG);
			s = " PULSE HEIGHT DISTRIBUTION IS SCORED!";// IN THOSE REGIONS
														// DENOTED WITH A 1";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			// for (int I=1;I<=EGS4Geom.nrg;I++)
			// {
			//
			// }

		} // "END IFULL=2 BLOCK"

	}

	@SuppressWarnings("unused")
	/**
	 * Unused here.
	 * @param P1 P1
	 * @param P2 P2
	 */
	private void EXCHANGE_STACK(int P1, int P2) {
		double FDUMMY = 0.0;
		int IDUMMY = 0;
		FDUMMY = EGS4.U[P2 - 1];
		EGS4.U[P2 - 1] = EGS4.U[P1 - 1];
		EGS4.U[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.V[P2 - 1];
		EGS4.V[P2 - 1] = EGS4.V[P1 - 1];
		EGS4.V[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.W[P2 - 1];
		EGS4.W[P2 - 1] = EGS4.W[P1 - 1];
		EGS4.W[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.E[P2 - 1];
		EGS4.E[P2 - 1] = EGS4.E[P1 - 1];
		EGS4.E[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.WT[P2 - 1];
		EGS4.WT[P2 - 1] = EGS4.WT[P1 - 1];
		EGS4.WT[P1 - 1] = FDUMMY;
		IDUMMY = EGS4.IQ[P2 - 1];
		EGS4.IQ[P2 - 1] = EGS4.IQ[P1 - 1];
		EGS4.IQ[P1 - 1] = IDUMMY;
		// "LATCH IS NOW STANDARD"
		IDUMMY = EGS4.LATCH[P2 - 1];
		EGS4.LATCH[P2 - 1] = EGS4.LATCH[P1 - 1];
		EGS4.LATCH[P1 - 1] = IDUMMY;
	}

	private void OSUMRY() {
		String s = "";
		int ll = 0;

		// $IMPLICIT-NONE;

		// COMIN/GEOM,IODAT1,IODAT2,PRINTC,SCORE,SOURCE,USER/;
		// COMIN/CH-Steps/;

	//	int I = 0;
		int IRL = 0;
		int IX = 0;
		int IZ = 0;
	//	double ASCT = 0.0;
	//	double ASCTUN = 0.0;
	//	double AATT = 0.0;
	//	double AATTUN = 0.0;
	//	double AWLL = 0.0;
	//	double AWLLUN = 0.0;
	//	double TDAW = 0.0;
	//	double TDAWUN = 0.0;
	//	double KSCT = 0.0;
	//	double KATT = 0.0;
	//	double KWLL = 0.0;
	//	double KSCTUN = 0.0;
	//	double KATTUN = 0.0;
	//	double KWLLUN = 0.0;

		// "SET UP THE PRINTER"
		// ICHPIN=12; "12 CHARACTERS/INCH"
		// ILPIN=6; "6 LINES/INCH"
		// IPAGE=0; "NO PAGE THROW"
		// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

		// "WRITE(IOUT,100)TITLE,DATEN,TIMEN; HEADER"

		// IF(ISOURC=21 | ISOURC=22)[
		// WRITE(IOUT,200) SCSTP,SCSTP2,
		// SCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		// NINCSRC),SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		// ]
		// ELSE[
		if (EGS4Macro.IFULL == 2) {
			// "PULSE HEIGHT DISTRIBUTION OUTPUT"
			// "IF(MAXBIN.LT.45)[CALL PRNTER(13,6,1,1);]ELSE[CALL PRNTER(13,8,1,1);]"
			// write(1,'(a)') '\f';
			// write(iout,101) title; call egs_fdate(iout);
			// write(iout,104);
			EGS4.seqStr = "=========================================================================";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "                   Summary of pulse height distribution";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "=========================================================================";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// " ********** IK: using wrapper for date and time routines.
			// " Was:
			// " WRITE(IOUT,101) TITLE, DATEN,TIMEN;
			PLOTPH();
			// CALL PLOTPH(TITLE,SCPDST,SCPDST2,SCPCUM,SCPCUM2,
			// SCPTOT,SCPTOT2,SCDFEP,SCDFEP2,MAXBIN,SLOTE,BINTOP,
			// IHSTRY,SCOMEG,SCOMEG2,SCPHEN,SCPHEN2);
		}// "END OF IFULL= 2 BLOCK"

		s = "                    # primary charged particle steps";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s + EGS4.format(SCSTP, 10, false) + " +/- "
				+ EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		s = "     # primary charged particle steps/initial history";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s
				+ EGS4.format(SCSTP / EGS4SrcEns.dble(IHSTRY), 10, false)
				+ " +/- " + EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		s = "# of presta-II steps/# primary charged particle steps";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s
				+ EGS4.format((EGS4.count_pII_steps + PIISTP) / SCSTP, 10,
						false) + " +/- " + EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,200) SCSTP,SCSTP2,SCSTP/dble(IHSTRY),
		// SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		// ]
		// 200 FORMAT(//' ' ,' # primary charged particle steps',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ',' # primary charged particle steps/initial history',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ','# of presta-II steps/# primary charged particle steps',
		// T58,F10.3,' +/- ',0PF6.3,'%');

		// "PRINT # CHARGED PARTICLE STEPS IN cavity REGION, etc"
		// IF(ISOURC=21 | ISOURC=22)[
		// WRITE(IOUT,210) SCCSTP,SCCSTP2,
		// SCCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		// NINCSRC),SCCSTP2;
		// ]
		// ELSE[
		// s="   # primary charged particle steps in cavity region";
		s = "   # primary charged particle steps in dose region";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		// EGS4.seqStr=s+EGS4.format(SCCSTP,10,false)+" +/- "+
		// EGS4.format(SCCSTP2,6)+"%";
		EGS4.seqStr = s + EGS4.format(SCDSTP, 10, false) + " +/- "
				+ EGS4.format(SCDSTP2, 6) + "%";

		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// s="    # primary steps in cavity region/initial history";
		s = "    # primary steps in dose region/initial history";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		// EGS4.seqStr=s+EGS4.format(SCCSTP/EGS4SrcEns.dble(IHSTRY),10,false)+" +/- "+
		// EGS4.format(SCCSTP2,6)+"%";
		EGS4.seqStr = s
				+ EGS4.format(SCDSTP / EGS4SrcEns.dble(IHSTRY), 10, false)
				+ " +/- " + EGS4.format(SCDSTP2, 6) + "%";

		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,210) SCCSTP,SCCSTP2,SCCSTP/dble(IHSTRY),SCCSTP2;
		// ]
		// 210 FORMAT(//' ',' # primary charged particle steps in cavity region'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ',' # primary steps in cavity region/initial history'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%');

		// if( EGS4SrcEns.ISOURC == 15 ) EGS4SrcEns.src15_out();//(iout);
		// "SCALE DOSE FRACTIONS"
		double TEMP = 0.0;
		// FRACS($MAXZREG,$MAXRADII,3:6);
		double[][][] FRACS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][6];
		if (EGS4Macro.IFULL == 1) {
			for (int IXD = 1; IXD <= NRDOSE; IXD++) {
				for (int IZD = 1; IZD <= NZDOSE; IZD++) {
					TEMP = SCDOSE[IZD - 1][IXD - 1][0];
					if (TEMP == 0) {
						for (int IT = 3; IT <= 6; IT++) {
							FRACS[IZD - 1][IXD - 1][IT - 1] = 0.0;
						}
					} else {
						for (int IT = 3; IT <= 6; IT++) {
							FRACS[IZD - 1][IXD - 1][IT - 1] = 100.
									* SCDOSE[IZD - 1][IXD - 1][IT - 1] / TEMP;
						}
					}
				}
			}
		}
		/*
		 * //"THE CAVITY SUMMARY" //"******************"
		 * 
		 * if(EGS4Geom.NSUMCV==1) {
		 * EGS4.seqStr="                   SUM OF RESULTS FOR THE CAVITY: 1 REGION"
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr="                   ***************************************"
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr); //WRITE(IOUT,300); }
		 * else {
		 * EGS4.seqStr="                   SUM OF RESULTS FOR THE CAVITY: " +
		 * EGS4Geom.NSUMCV+" REGION"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * EGS4.seqStr="                   ***************************************"
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,301)NSUMCV; }
		 * 
		 * //300 FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/
		 * //' ',T20,'***************************************'); //301 FORMAT(//
		 * ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
		 * //' ',T20,'*****************************************');
		 * //"******************" //"history by history" //"    EMH March 2002"
		 * //"******************" cav2_dose = cav2_dose/cav_dose; if
		 * (EGS4Macro.iifano == 1) {
		 * EGS4.seqStr=" This calculation was performed using regeneration ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" ================================================= ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * s=" D/Awall (grays/incident fluence): ";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(cav_dose,11,false)+" +/- "+
		 * EGS4.format(100*cav2_dose,6)+"%"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * //write(iout,'(a,t50,1PE11.4,a,0PF6.3,a)') // 'D/Awall
		 * (grays/incident fluence): ',cav_dose,' +/- ', // 100*cav2_dose,'%';
		 * return; }
		 * 
		 * if( EGS4Macro.use_enhance ) {
		 * EGS4.seqStr=" This calculation was performed using CS enhancement ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr="  enhancement factor was "
		 * +EGS4.format(EGS4Macro.cs_enhance,10,true); if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" ================================================= ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); } else if (
		 * EGS4Macro.n_split > 1 ) {
		 * EGS4.seqStr=" This calculation was performed using photon splitting "
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr="  splitting number was "
		 * +EGS4.format(EGS4Macro.n_split,6); if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" ================================================= ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr); } if (IFANO==1) {
		 * EGS4.seqStr=" This calculation was performed using regeneration ";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * } if (IFANO==2) { EGS4.seqStr=
		 * " This calculation was performed eliminating electrons originating in the cavity wall."
		 * ; if(EGS4.iprint>1) printSequence(EGS4.seqStr); } //"THE TOTAL DOSE"
		 * if
		 * (EGS4SrcEns.ISOURC==3||EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22)
		 * { s=" TOTAL DOSE (GRAYS/INCIDENT PARTICLE):";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * 
		 * EGS4.seqStr=s+
		 * EGS4.format(cav_dose,11,false)+" +/- "+EGS4.format(100*
		 * cav2_dose,5)+"%"; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,311)cav_dose,100*cav2_dose; } else {
		 * s=" TOTAL DOSE (GRAYS/INCIDENT FLUENCE):";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * 
		 * EGS4.seqStr=s+
		 * EGS4.format(cav_dose,11,false)+" +/- "+EGS4.format(100*
		 * cav2_dose,5)+"%"; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,310)cav_dose,100*cav2_dose; } //310 FORMAT(//'TOTAL DOSE
		 * (GRAYS/INCIDENT FLUENCE):', //T50,1PE11.4,' +/- ',0PF5.2,'%'); //311
		 * FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT PARTICLE):', //T50,1PE11.4,' +/-
		 * ',0PF5.2,'%');
		 * 
		 * if(EGS4Macro.IFULL==1) { //"CALCULATE Ascat, Aatt, Awall, DOSE/Awall"
		 * 
		 * cav2_dose1 = cav2_dose1/cav_dose1; cav_dosec =
		 * cav_dosec/cav_dose/cav_dose1;
		 * 
		 * //"total dose/Awall" TDAW=cav_dose1; TDAWUN=cav2_dose1;
		 * 
		 * s="  TOTAL DOSE/Awall:";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (TDAW,11,false)+" +/- "+EGS4.format(100*TDAWUN,5)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,350)TDAW,TDAWUN*100.; //350 FORMAT(' TOTAL DOSE/Awall:',
		 * //T50,1PE11.4,' +/- ',0PF5.2,'%');
		 * 
		 * cav_dosec = cav2_dose*cav2_dose + cav2_dose1*cav2_dose1 -
		 * 2*cav_dosec; if( cav_dosec > 0 ) cav_dosec = Math.sqrt(cav_dosec);
		 * //"Awall" AWLL=cav_dose/cav_dose1; AWLLUN=cav_dosec;
		 * //WRITE(IOUT,340)AWLL,AWLLUN*100.; //340 FORMAT(' Awall:',T50,F8.5,'
		 * +/- ',F5.2,'%'); s=" Awall:";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (AWLL,8,true)+"  +/- "+EGS4.format(100*AWLLUN,5,true)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1 ) { corr_02 =
		 * corr_02/(cav_dose0+cav_dose2);
		 * 
		 * cav2_dose0 = cav2_dose0/cav_dose0; cav2_dose2 = cav2_dose2/cav_dose2;
		 * 
		 * cav_dosec01 = cav_dosec01/cav_dose0/cav_dose1; cav_dosec02 =
		 * cav_dosec02/cav_dose0/cav_dose2;
		 * 
		 * cav_dosec01 = cav2_dose0*cav2_dose0 + cav2_dose1*cav2_dose1 -
		 * 2*cav_dosec01; cav_dosec02 = cav2_dose0*cav2_dose0 +
		 * cav2_dose2*cav2_dose2 - 2*cav_dosec02;
		 * 
		 * if( cav_dosec01 > 0 ) cav_dosec01 = Math.sqrt(cav_dosec01); if(
		 * cav_dosec02 > 0 ) cav_dosec02 = Math.sqrt(cav_dosec02);
		 * 
		 * //"... we had rel.error(x), get now rel.error(1+x)" cav_dosec02 =
		 * cav_dosec02/(1 + cav_dose0/cav_dose2);
		 * 
		 * //"Ascatt" ASCT=1+cav_dose2/cav_dose0; ASCTUN=cav_dosec02;
		 * //WRITE(IOUT,320)ASCT,ASCTUN*100.; /// 320 FORMAT('
		 * Ascat:',T50,F8.5,' +/- ',F5.2,'%'); s=" Ascat:";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (ASCT,8,true)+"  +/- "+EGS4.format(100*ASCTUN,5,true)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //"Aatt" AATT=cav_dose0/cav_dose1; AATTUN=cav_dosec01;
		 * //WRITE(IOUT,330)AATT,AATTUN*100.; //330 FORMAT(' Aatt :',T50,F8.5,'
		 * +/- ',F5.2,'%'); s=" Aatt :";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (AATT,8,true)+"  +/- "+EGS4.format(100*AATTUN,5,true)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * 
		 * //"Dpr + Dsec" //345 FORMAT(' Dprimary + Dsecondary:',T50,1PE11.4,'
		 * +/- ',0PF5.2,'%'); //WRITE(IOUT,345)cav_dose2+cav_dose0,corr_02*100.;
		 * s=" Dprimary + Dsecondary:";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(cav_dose2+cav_dose0,11,false)+"  +/- "+
		 * EGS4.format(corr_02*100.,5,true)+"%"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * } } //"END OF IFULL = 1"
		 * 
		 * //"OUTPUT Kscat,Katt,Kwall,Kpn,Kfl THE INVERSES OF THE Ai's"
		 * if(EGS4Macro.IFULL>0) { KWLL=1./AWLL;KWLLUN=AWLLUN/(AWLL*AWLL);
		 * s=" Kwall:"; ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr
		 * =s+EGS4.format(KWLL,8,true)+"  +/- "+EGS4.format(100*KWLLUN
		 * ,5,true)+"%"; if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,341)KWLL,KWLLUN*100.; // 341 FORMAT(' Kwall:',T50,F8.5,'
		 * +/- ',F5.2,'%'); if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1
		 * ) { KSCT=1./ASCT;KSCTUN=ASCTUN/(ASCT*ASCT); s=" Kscat:";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (KSCT,8,true)+"  +/- "+EGS4.format(100*KSCTUN,5,true)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,321)KSCT,KSCTUN*100.; /// 321 FORMAT('
		 * Kscat:',T50,F8.5,' +/- ',F5.2,'%');
		 * KATT=1./AATT;KATTUN=AATTUN/(AATT*AATT);
		 * //WRITE(IOUT,331)KATT,KATTUN*100.; //331 FORMAT(' Katt :',T50,F8.5,'
		 * +/- ',F5.2,'%'); s=" Katt :";
		 * ll=s.length();ll=50-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format
		 * (KATT,8,true)+"  +/- "+EGS4.format(100*KATTUN,5,true)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * } }
		 * 
		 * //"THE DETAILED OUTPUT FOR EACH CAVITY ZONE"
		 * //"****************************************" if(IOOPTN==1 &&
		 * EGS4Geom.NSUMCV>1) {
		 * //"ONLY IF REQUESTED AND MORE THAN ONE CAVITY ZONE"
		 * 
		 * //WRITE(IOUT,400) NSUMCV; // 400 FORMAT(// ,T15,'DETAILED RESULTS FOR
		 * EACH OF THE ',I4,' CAVITY REGIONS'/ ///
		 * ' ',T15,'****************************************************');
		 * EGS4.
		 * seqStr="   DETAILED RESULTS FOR EACH OF THE "+EGS4.format(EGS4Geom
		 * .NSUMCV,4)+ " CAVITY REGIONS"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * //"PRINT THE TABLE HEADER" if(EGS4Macro.IFULL==0) {
		 * //WRITE(IOUT,410); //410 FORMAT(//'Z# P# C# Total Dose '/ /// ' -- --
		 * -- -------------------');
		 * EGS4.seqStr=" Z# P# C#     Total Dose     "; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" -- -- -- -------------------"; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * } else if(EGS4Macro.IFULL>0) { //"SET UP THE PRINTER" //ICHPIN=16;
		 * "16 CHARACTERS/INCH" //ILPIN=6; "6 LINES/INCH" //IPAGE=0;
		 * "NO PAGE THROW" //"CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"
		 * 
		 * //WRITE(IOUT,420); /// 420 FORMAT(//'Z# P# C#', //' Total Dose ', //'
		 * Ascat ', //' Aatt ', //' Awall ', //' Total Dose/Awall'/ //' -- --
		 * --', //' ---------- ', //' ----- ', //' ---- ', //' ----- ', //'
		 * ----------------');
		 * EGS4.seqStr=" Z# P# C#"+"     Total Dose   "+"      Ascat      "+
		 * "       Aatt      "+"      Awall      "+"  Total Dose/Awall";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" -- -- --"+"     ----------   "+"      -----      "+
		 * "       ----      "+"      -----      "+"  ----------------";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * }
		 * 
		 * //"LOOP OVER THE CAVITY ZONES" for( I=1;I<=EGS4Geom.NSUMCV;I++) {
		 * IRL=EGS4Geom.ISUMCV[I-1]; //$GET-IX-IZ(IRL); IX=EGS4Geom.GET_IX(IRL);
		 * IZ=EGS4Geom.GET_IZC(IRL);
		 * 
		 * //"total dose" if(SCDOSE[IZ-1][IX-1][0]>0.)
		 * SCDOSE2[IZ-1][IX-1][0]=SCDOSE2[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][0];
		 * if(EGS4Macro.IFULL==1) { //"CALCULATE Ascat, Aatt, Awall, DOSE/Awall"
		 * 
		 * if(SCDOSE[IZ-1][IX-1][2]>0.) {
		 * SCDOSE2[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX-1][2]/SCDOSE[IZ-1][IX-1][2];
		 * SCDOSE_COV
		 * [IZ-1][IX-1][0]=SCDOSE_COV[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][0]/
		 * SCDOSE[IZ-1][IX-1][2]; }
		 * 
		 * //"TOTAL DOSE/Awall"
		 * TDAW=SCDOSE[IZ-1][IX-1][2];TDAWUN=SCDOSE2[IZ-1][IX-1][2];
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][0] =
		 * SCDOSE2[IZ-1][IX-1][0]*SCDOSE2[IZ-1][IX-1][0] +
		 * SCDOSE2[IZ-1][IX-1][2]*SCDOSE2[IZ-1][IX-1][2] -
		 * 2*SCDOSE_COV[IZ-1][IX-1][0];
		 * if(SCDOSE_COV[IZ-1][IX-1][0]>0.)SCDOSE_COV[IZ-1][IX-1][0]=
		 * Math.sqrt(SCDOSE_COV[IZ-1][IX-1][0]); //"Awall"
		 * AWLL=SCDOSE[IZ-1][IX-1][0]/SCDOSE[IZ-1][IX-1][2];
		 * AWLLUN=SCDOSE_COV[IZ-1][IX-1][0];
		 * 
		 * if( !EGS4Macro.use_enhance && EGS4Macro.n_split <= 1 ) {
		 * //"this condition not strictly necessary because"
		 * //"IOOPTN set = 0 if either of these options on"
		 * 
		 * if(SCDOSE[IZ-1][IX-1][1]>0.) {
		 * SCDOSE2[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][1];
		 * SCDOSE_COV
		 * [IZ-1][IX-1][1]=SCDOSE_COV[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][1]/
		 * SCDOSE[IZ-1][IX-1][2]; if(SCDOSE[IZ-1][IX-1][3]>0.) {
		 * SCDOSE_COV[IZ-1]
		 * [IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]/SCDOSE[IZ-1][IX-1][1]/
		 * SCDOSE[IZ-1][IX-1][3]; } }
		 * if(SCDOSE[IZ-1][IX-1][3]>0.)SCDOSE2[IZ-1][IX
		 * -1][3]=SCDOSE2[IZ-1][IX-1][3]/ SCDOSE[IZ-1][IX-1][3];
		 * 
		 * SCDOSE_COV[IZ-1][IX-1][1]=SCDOSE2[IZ-1][IX-1][1]*SCDOSE2[IZ-1][IX-1][1
		 * ]+ SCDOSE2[IZ-1][IX-1][2]*SCDOSE2[IZ-1][IX-1][2]-
		 * 2*SCDOSE_COV[IZ-1][IX-1][1];
		 * SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE2[IZ-1][IX
		 * -1][1]*SCDOSE2[IZ-1][IX-1][1]+
		 * SCDOSE2[IZ-1][IX-1][3]*SCDOSE2[IZ-1][IX-1][3]-
		 * 2*SCDOSE_COV[IZ-1][IX-1][2];
		 * 
		 * if(SCDOSE_COV[IZ-1][IX-1][1]>0)SCDOSE_COV[IZ-1][IX-1][1]=
		 * Math.sqrt(SCDOSE_COV[IZ-1][IX-1][1]);
		 * if(SCDOSE_COV[IZ-1][IX-1][2]>0)SCDOSE_COV[IZ-1][IX-1][2]=
		 * Math.sqrt(SCDOSE_COV[IZ-1][IX-1][2]);
		 * 
		 * //"... we had rel.error(x), get now rel.error(1+x)"
		 * if(SCDOSE[IZ-1][IX-1][3]>0.)
		 * SCDOSE_COV[IZ-1][IX-1][2]=SCDOSE_COV[IZ-1][IX-1][2]/
		 * (1+SCDOSE[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][3]);
		 * 
		 * //"Ascatt" ASCT=1+SCDOSE[IZ-1][IX-1][3]/SCDOSE[IZ-1][IX-1][1];
		 * ASCTUN=SCDOSE_COV[IZ-1][IX-1][2];
		 * 
		 * //"Aatt" AATT=SCDOSE[IZ-1][IX-1][1]/SCDOSE[IZ-1][IX-1][2];
		 * AATTUN=SCDOSE_COV[IZ-1][IX-1][1];
		 * 
		 * } } if(EGS4Macro.IFULL==0) {
		 * EGS4.seqStr=" "+EGS4.format(IRL,2)+" "+EGS4
		 * .format(IZ,2)+" "+EGS4.format(IX,2)+" "+
		 * EGS4.format(SCDOSE[IZ-1][IX-1
		 * ][0],11,false)+" +/-"+EGS4.format(100*SCDOSE2[IZ-1][IX-1][0],5)+"%";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,430) // IRL,IZ,IX, "POSITION" //
		 * SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100.; "TOTAL DOSE" }
		 * 
		 * //430 FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,' +/-',0PF5.2,'%'); //440
		 * FORMAT(' ',I2,2(1X,I2), //1X,1PE11.4,'(',0PF5.2,'%)',
		 * //3(F8.5,'(',F7.5,')'), //1X,1PE11.4,'(',0PF5.2,'%)');
		 * 
		 * else if(EGS4Macro.IFULL==1) {
		 * EGS4.seqStr=" "+EGS4.format(IRL,2)+" "+EGS4
		 * .format(IZ,2)+" "+EGS4.format(IX,2)+" "+
		 * EGS4.format(SCDOSE[IZ-1][IX-1
		 * ][0],11,false)+"("+EGS4.format(100*SCDOSE2[IZ-1][IX-1][0],5)+"%)"+
		 * EGS4.format(ASCT,8,true)+"("+EGS4.format(ASCTUN,7,true)+")"+
		 * EGS4.format(AATT,8,true)+"("+EGS4.format(AATTUN,7,true)+")"+
		 * EGS4.format(AWLL,8,true)+"("+EGS4.format(AWLLUN,7,true)+")"+
		 * EGS4.format(TDAW,11,false)+"("+EGS4.format(100*TDAWUN,5)+"%)";
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * //WRITE(IOUT,440) // IRL,IZ,IX, "POSITION" //
		 * SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100., "TOTAL DOSE" // ASCT,ASCTUN,
		 * "Ascat" // AATT,AATTUN, "Aatt" // AWLL,AWLLUN, "Awall" //
		 * TDAW,TDAWUN*100; "TOTAL DOSE/Awall" }
		 * 
		 * }// "END OF LOOP OVER CAVITY ZONES" }// "END OF DETAILED SUMMARY"
		 * 
		 * //"RESET UP THE PRINTER" //ICHPIN=12; "12 CHARACTERS/INCH" //ILPIN=6;
		 * "6 LINES/INCH" //IPAGE=0; "NO PAGE THROW"
		 * //"CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"
		 * 
		 * //RETURN; //%I0 //"FORMATS" //200 FORMAT(//' ' ,' # primary charged
		 * particle steps',T58, //1PE10.3,' +/- ',0PF6.3,'%'/ /// ' ',' #
		 * primary charged particle steps/initial history',T58, //1PE10.3,' +/-
		 * ',0PF6.3,'%'/ /// ' ','# of presta-II steps/# primary charged
		 * particle steps', ////T58,F10.3,' +/- ',0PF6.3,'%'); /// 202
		 * FORMAT(//' ' ,' # charged particle steps in run',T58, /// 1PE10.3,/
		 * /// ' ',' # charged particle steps in run/initial history',T58,
		 * //1PE10.3/ //' ','# of presta-II steps/# primary charged particle
		 * steps', /// T58,F10.3,' +/- ',0PF6.3,'%'); //210 FORMAT(//' ',' #
		 * primary charged particle steps in cavity region' //,T58,1PE10.3,' +/-
		 * ',0PF6.3,'%'/ //' ',' # primary steps in cavity region/initial
		 * history' //,T58,1PE10.3,' +/- ',0PF6.3,'%'); //220 FORMAT(//
		 * ,T8,'STEP COUNTING RESULTS FOR WALL MATERIAL IN THE CAVITY'/ /// //
		 * ,T8,'# primary charged particle steps',T51, /// I12,' +/-
		 * ',0PF5.2,'%'/ //' ',T8,'# OF Times mscat switched off',T51, ////I12,'
		 * +/- ',0PF5.2,'%'/ /// ' ',T8,'RATIO',T54,F7.3,' +/- ',0PF5.2,'%');
		 * //230 FORMAT(// ,T8,'# Primary charged particle steps in cavity
		 * region' /// ,T51,I12,' +/- ',0PF5.2,'%'/ //' ',T8,'# Times mscat
		 * switched off in cavity region.' //,T51,I12,' +/- ',0PF5.2,'%'/
		 * //' ',T8,'Ratio',T56,F7.3,' +/- ',0PF5.2,'%'); /// 300 FORMAT(//
		 * ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/ ///
		 * ' ',T20,'***************************************'); //301 FORMAT(//
		 * ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
		 * //' ',T20,'*****************************************'); //310
		 * FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT FLUENCE):', //T50,1PE11.4,' +/-
		 * ',0PF5.2,'%'); //311 FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT
		 * PARTICLE):', //T50,1PE11.4,' +/- ',0PF5.2,'%'); /// 320 FORMAT('
		 * Ascat:',T50,F8.5,' +/- ',F5.2,'%'); /// 321 FORMAT('
		 * Kscat:',T50,F8.5,' +/- ',F5.2,'%'); //330 FORMAT(' Aatt :',T50,F8.5,'
		 * +/- ',F5.2,'%'); //331 FORMAT(' Katt :',T50,F8.5,' +/- ',F5.2,'%');
		 * //340 FORMAT(' Awall:',T50,F8.5,' +/- ',F5.2,'%'); //341 FORMAT('
		 * Kwall:',T50,F8.5,' +/- ',F5.2,'%'); //345 FORMAT(' Dprimary +
		 * Dsecondary:',T50,1PE11.4,' +/- ',0PF5.2,'%'); //350 FORMAT(' TOTAL
		 * DOSE/Awall:', //T50,1PE11.4,' +/- ',0PF5.2,'%'); /// 360
		 * FORMAT(//'Apn :',T50,F8.5,' +/- ',F8.5); //361 FORMAT(' Kpn
		 * :',T50,F8.5,' +/- ',F8.5); //370 FORMAT(' TOTAL DOSE/([Apn]*Awall):',
		 * //T50,1PE11.4,' +/- ',0PF5.2,'%'); //380 FORMAT(//'Afl :',T50,F8.5,'
		 * +/- ',F8.5); //381 FORMAT(' Kfl :',T50,F8.5,' +/- ',F8.5); //390
		 * FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall):', //T50,1PE11.4,' +/-
		 * ',0PF5.2,'%'); //395 FORMAT(//'<s>g,w :',T50,F8.5,' +/- ',F8.5); ///
		 * 396 FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall*<s>g,w):', /// T50,1PE11.4,'
		 * +/- ',0PF5.2,'%'); //400 FORMAT(// ,T15,'DETAILED RESULTS FOR EACH OF
		 * THE ',I4,' CAVITY REGIONS'/ ///
		 * ' ',T15,'****************************************************');
		 * //410 FORMAT(//'Z# P# C# Total Dose '/ /// ' -- -- --
		 * -------------------'); /// 420 FORMAT(//'Z# P# C#', //' Total Dose ',
		 * //' Ascat ', //' Aatt ', //' Awall ', //' Total Dose/Awall'/ //' --
		 * -- --', //' ---------- ', //' ----- ', //' ---- ', //' ----- ', //'
		 * ----------------'); //430 FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,'
		 * +/-',0PF5.2,'%'); //440 FORMAT(' ',I2,2(1X,I2),
		 * //1X,1PE11.4,'(',0PF5.2,'%)', //3(F8.5,'(',F7.5,')'),
		 * //1X,1PE11.4,'(',0PF5.2,'%)'); //445 FORMAT(' ',I2,2(1X,I2),
		 * //1X,1PE11.4,'(',0PF5.2,'%)', //3(F8.5,'(',F7.5,')'),
		 * //1X,1PE11.4,'(',0PF5.2,'%)', //F8.5,'(',F7.5,')',
		 * //1X,1PE11.4,'(',0PF5.2,'%)'); /// 450 FORMAT(' ',I2,2(1X,I2),
		 * //1X,1PE11.4, //3(1X,F8.5,1X), //1X,1PE11.4, //1X,F8.5,1X,
		 * //1X,1PE11.4, //1X,F8.5,1X, //1X,1PE11.4, //1X,F8.5,1X, ///
		 * 1X,1PE11.4/ ////' ','ERRORS = ', //4X,'(',0PF5.2,'%)', ///
		 * 3(1X,'(',F7.5,')'), //4X,'(',0PF5.2,'%)', /// 1X,'(',F7.5,')',
		 * //4X,'(',0PF5.2,'%)', //1X,'(',F7.5,')', /// 4X,'(',0PF5.2,'%)', ///
		 * 1X,'(',F7.5,')', //4X,'(',0PF5.2,'%)'); // /// END;
		 * "END OF SUBROUTINE OSUMRY" //
		 */
		// "PRINT A SUMMARY OF THE DOSE REGION RESULTS"
		IRL = 0;
		if ((IOOPTN == 1) || (IOOPTN > 2)) {

			// "A COMPACT VERSION IF ONLY ONE DOSE REGION ZONE"
			if (NDOSE == 1) {
				IZ = NZDMIN;
				IX = NRDMAX;
				// $GET-IRL(IZ,IX);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				// if(EGS4SrcEns.ISOURC==3||EGS4SrcEns.ISOURC==21||EGS4SrcEns.ISOURC==22||EGS4SrcEns.ISOURC==23)
				// {
				// WRITE(IOUT,299)
				// IRL,IZ,IX,(SCDOSE(1,1,IT),SCDOSE2(1,1,IT),IT=1,2);
				// }
				// else
				// {
				// @@@@@@@@@@@@@@@WRITE(IOUT,300)
				// IRL,IZ,IX,(SCDOSE(1,1,IT),SCDOSE2(1,1,IT),IT=1,2);
				// 300 FORMAT(/' Geometrical zone number:',T53,I10/
				// ' Planar zone number:',T53,I10/
				// ' Cylndrical zone number:',T53,I10//
				// ' Total dose (Gray/(incident fluence)):',
				// T50,1PE11.4,' +/- ',0PF6.3,'%'/
				// ' Total dose minus stoppers:',
				// T50,1PE11.4,' +/- ',0PF6.3,'%');
				s = " Geometrical zone number:";
				ll = s.length();
				ll = 53 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(IRL, 3) + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " Planar zone number:";
				ll = s.length();
				ll = 53 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(IZ, 3) + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " Cylndrical zone number:";
				ll = s.length();
				ll = 53 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(IX, 3) + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " Total dose (Gray/(incident fluence)):";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDOSE[0][0][0], 11, false)
						+ " +/- " + EGS4.format(SCDOSE2[0][0][0], 6, true)
						+ "%" + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " Total dose minus stoppers:";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDOSE[0][0][1], 11, false)
						+ " +/- " + EGS4.format(SCDOSE2[0][0][1], 6, true)
						+ "%" + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				// }
				if (EGS4Macro.IFULL == 1) {
					// @@@@@@@@@@@@@@WRITE(IOUT,302)
					// (SCDOSE(1,1,IT),SCDOSE2(1,1,IT),FRACS(1,1,IT),IT=3,6);
					/*
					 * if(ISOURC==3) { if(IRL==NSRCRG) {
					 * //WRITE(IOUT,303)SCDOSE(1,1,7),SCDOSE2(1,1,7); } }
					 */
				}
			} else {// "OUTPUT FOR MORE THAN 1 SCORING ZONE"
			// @@@@@@@@@WRITE(IOUT,400);
			// 400 FORMAT(/' ',T20,'Z# : Geometrical zone number'/
			// ' ',T20,'P# : Planar zone number'/
			// ' ',T20,'C# : Cylndrical zone number');
				s = "	 Z# : Geometrical zone number:";
				// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = "	 P# : Planar zone number:";
				// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = "	 C# : Cylndrical zone number:";
				// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s + " \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				if (EGS4Macro.IFULL != 3) {
					// IF(EGS4SrcEns.ISOURC=3|EGS4SrcEns.ISOURC=21|ISOURC=22|ISOURC=23)[WRITE(IOUT,399);]
					// ELSE
					// {
					// WRITE(IOUT,401);
					// }
				}
				if (EGS4Macro.IFULL == 1) {// "dose per entrance region"
				// @@@@@WRITE(IOUT,410);
				// 410 FORMAT(' ',T20,'F : TOTAL DOSE FROM FRONT PLANAR WALL'/
				// ' ',T20,'O : TOTAL DOSE FROM OUTSIDE CURVED WALL'/
				// ' ',T20,'B : TOTAL DOSE FROM BACK PLANAR WALL'/
				// ' ',T20,'I : TOTAL DOSE FROM INNER CURVED WALL');
					s = "	 F  : TOTAL DOSE FROM FRONT PLANAR WALL";
					// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);
					s = "	 O  : TOTAL DOSE FROM OUTSIDE CURVED WALL";
					// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);
					s = "	 B  : TOTAL DOSE FROM BACK PLANAR WALL";
					// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);
					s = "	 I  : TOTAL DOSE FROM INNER CURVED WALL";
					// ll=s.length();ll=53-ll;s=s+EGS4.format("",ll);
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

					// ICHPIN=12;IPAGE=0;
					// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);  SET UP THE PRINTER"
					// @@@@@@@@@WRITE(IOUT,411);
					// 411 FORMAT(/' Z# P# C# ',
					// ' T T-S F O ',
					// ' B I'/
					// ' --- -- -- ',
					// '---------- ---------- ---------- ---------- ',
					// '---------- ----------');
					s = "  Z# P# C# "
							+ "     T         T-S         F          O     "
							+ "     B          I";
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);
					s = " --- -- -- "
							+ "---------- ---------- ---------- ---------- "
							+ "---------- ----------";
					EGS4.seqStr = s + " \n";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

					for (int IXD = 1; IXD <= NRDOSE; IXD++) {
						for (int IZD = 1; IZD <= NZDOSE; IZD++) {
							IZ = IZD - 1 + NZDMIN;
							IX = IXD + NRDMIN;
							// $GET-IRL(IZ,IX);
							IRL = EGS4Geom.GET_IRL(IZ, IX);
							// @@@@@@WRITE(IOUT,412) IRL,IZ,IX,
							// (SCDOSE(IZD,IXD,IT),IT=1,$MAXIT-1),
							// (SCDOSE2(IZ,IX,IT),IT=1,$MAXIT-1),
							// (FRACS(IZ,IX,IT),IT=3,6);
							// 412
							// FORMAT(' ',I3,2(1X,I2),1X,1PE11.4,5(1X,E10.3)/
							// ' %ERROR=',
							// 6(3X,0PF6.3,4X)/
							// ' %DOSE=',21X,4(2X,0PF7.3,2X));
							s = " "
									+ EGS4.format(IRL, 3)
									+ EGS4.format("", 1)
									+ EGS4.format(IZ, 2)
									+ EGS4.format("", 1)
									+ EGS4.format(IX, 2)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][0],
											11, false)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][1],
											10, false)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][2],
											10, false)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][3],
											10, false)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][4],
											10, false)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][5],
											10, false);
							EGS4.seqStr = s;// \n";
							if (EGS4.iprint > 0)
								printSequence(EGS4.seqStr);
							s = "   %ERROR="
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][0],
											6, true)
									+ EGS4.format("", 4)
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][1],
											6, true)
									+ EGS4.format("", 4)
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][2],
											6, true)
									+ EGS4.format("", 4)
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][3],
											6, true)
									+ EGS4.format("", 4)
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][4],
											6, true)
									+ EGS4.format("", 4)
									+ EGS4.format("", 3)
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][5],
											6, true);
							EGS4.seqStr = s;// \n";
							if (EGS4.iprint > 0)
								printSequence(EGS4.seqStr);
							s = "    %DOSE="
									+ EGS4.format("", 21)
									+ EGS4.format("", 2)
									+ EGS4.format(FRACS[IZD - 1][IXD - 1][2],
											6, true)
									+ EGS4.format("", 2)
									+ EGS4.format("", 2)
									+ EGS4.format(FRACS[IZD - 1][IXD - 1][3],
											6, true)
									+ EGS4.format("", 2)
									+ EGS4.format("", 2)
									+ EGS4.format(FRACS[IZD - 1][IXD - 1][4],
											6, true)
									+ EGS4.format("", 2)
									+ EGS4.format("", 2)
									+ EGS4.format(FRACS[IZD - 1][IXD - 1][5],
											6, true) + EGS4.format("", 2);
							EGS4.seqStr = s;// \n";
							if (EGS4.iprint > 0)
								printSequence(EGS4.seqStr);

						}
					}
					/*
					 * IF(ISOURC.EQ.3)[ IF(CDSTBL(NSRCRG).EQ.DCHAR)[
					 * $GET-IX-IZ(NSRCRG);IZD=IZ-NZDMIN+1;IXD=IX-NRDMIN;
					 * WRITE(IOUT,304)NSRCRG,IZ,IX,
					 * SCDOSE(IZD,IXD,7),SCDOSE2(IZD,IXD,7); ] ]
					 */
				} else {// "no dose per entrance region"
					if (EGS4Macro.IFULL == 3) {
						// IF(ISOURC=3|ISOURC=21|ISOURC=22|ISOURC=23)[WRITE(IOUT,:F399B:);]
						// ELSE[
						// @@@@@@@@@@@WRITE(IOUT,:F401B:);]
						// :F401B: FORMAT(' ',T20,'T : Total dose
						// (Gray/(incident fluence))'/
						// ' ',T20,
						// 'Sca: Scatter dose (after Compton or fluorecent
						// reabsorbed ');
						s = "	";
						ll = s.length();
						ll = 20 - ll;
						s = s + EGS4.format("", ll);
						EGS4.seqStr = s
								+ "T  : Total dose (Gray/(incident fluence))";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);
						s = "	";
						ll = s.length();
						ll = 20 - ll;
						s = s + EGS4.format("", ll);
						EGS4.seqStr = s
								+ "Sca: Scatter dose (after Compton or fluorecent  reabsorbed ";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);

						// @@@@@@@@@@@@@WRITE(IOUT,:F402B:);
						// :F402B: FORMAT(/' Z# P# C# T Sca'/
						// ' --- -- -- -------------------
						// -------------------');
						s = "  Z# P# C#          T                  Sca";
						EGS4.seqStr = s + " \n";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);
						s = " --- -- -- ------------------- -------------------";
						EGS4.seqStr = s + " \n";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);

					} else {
						// @@@@@@@WRITE(IOUT,402);
						// 402 FORMAT(/' Z# P# C# T T-S'/
						// ' --- -- -- -------------------
						// -------------------');
						s = "  Z# P# C#          T                  T-S";
						EGS4.seqStr = s + " \n";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);
						s = " --- -- -- ------------------- -------------------";
						EGS4.seqStr = s + " \n";
						if (EGS4.iprint > 0)
							printSequence(EGS4.seqStr);

					}
					for (int IXD = 1; IXD <= NRDOSE; IXD++) {
						for (int IZD = 1; IZD <= NZDOSE; IZD++) {
							IZ = IZD - 1 + NZDMIN;
							IX = IXD + NRDMIN;
							// $GET-IRL(IZ,IX);
							IRL = EGS4Geom.GET_IRL(IZ, IX);
							// @@@@@@@@@@@WRITE(IOUT,403)
							// IRL,IZ,IX,(SCDOSE(IZD,IXD,IT),SCDOSE2(IZD,IXD,IT),IT=1,2);
							// 403 FORMAT(' ',I3,2(1X,I2),1X,2(1PE11.4,'
							// +/-',0PF6.3,'% '));
							s = " "
									+ EGS4.format(IRL, 3)
									+ EGS4.format("", 1)
									+ EGS4.format(IZ, 2)
									+ EGS4.format("", 1)
									+ EGS4.format(IX, 2)
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][0],
											11, false)
									+ " +/- "
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][0],
											6, true)
									+ "% "
									+ EGS4.format("", 1)
									+ EGS4.format(SCDOSE[IZD - 1][IXD - 1][1],
											10, false)
									+ " +/- "
									+ EGS4.format(SCDOSE2[IZD - 1][IXD - 1][1],
											6, true) + "% ";
							EGS4.seqStr = s;// \n";
							if (EGS4.iprint > 0)
								printSequence(EGS4.seqStr);

						}
					}
					// "Here we output a brief summary to unit 10 for a database"
					// dose_unit = egs_open_file(10,0,1,'.egsdose');
					// WRITE(dose_unit,'(80A1)') TITLE;
					// WRITE(dose_unit,
					// '('' There are'',I5,'' radial zones, midpoints:'')')NR;
					// WRITE(dose_unit,'(
					// 8(F9.4,'',''))')((RCYL(I-1)+RCYL(I))/2., I=1,NR);
					// WRITE(dose_unit,'('' There are '',I5,'' depth
					// regions'')') NZ;
					// DO IZ = 1,NZ[
					// WRITE(dose_unit,'('' Depth centered at: '',F12.3)')
					// (ZPLANE(IZ)+ZPLANE(IZ+1))/2.;
					// WRITE(dose_unit,'( 4(1PE10.3,''+/-'',0PF5.2,''% ''))')
					// (SCDOSE(IZ,IR,1),SCDOSE2(IZ,IR,1), IR=1,NR);
					// ]
					// close(dose_unit);
				}
			}// "END OF DOSE SUMMARY"
		} // "END OF CONDITIONAL SUMMARY OF DOSE"
			// double[][][] RESULTS=new
			// double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		// double[][][] UNCRTY=new
		// double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];//($MAXZREG,
		// $MAXRADII, $MAXCMPTS),
		// double[] RADIAL_BINS=new double[EGS4Geom.$MAXRADII];
		// double[] DEPTH_BINS=new double[EGS4Geom.$MAXZPLANE];
		EGS4Grid.RESULTS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		EGS4Grid.UNCRT = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];// ($MAXZREG,
																						// $MAXRADII,
																						// $MAXCMPTS),
		EGS4Grid.RADIAL_BINS = new double[EGS4Geom.$MAXRADII];
		EGS4Grid.DEPTH_BINS = new double[EGS4Geom.$MAXZPLANE];

		// CHARACTER*60 EXPLANATIONS($MAXCMPTS);
		// CHARACTER*4 LABELS($MAXCMPTS);
		// String[] EXPLANATIONS=new String[$MAXCMPTS];
		// String[] LABELS=new String[$MAXCMPTS];
		// int NCOMP=0;

		if ((IOOPTN == 0) || (IOOPTN == 2) || (IOOPTN == 4)) {
			// " OUTPUT GRID ONLY IF REQUESTED"
			// "Grid routine by Aaron Merovitz, 1998"

			// "1) Set up the arrays"
			for (int IXD = 1; IXD <= NRDOSE; IXD++) {
				for (int IZD = 1; IZD <= NZDOSE; IZD++) {
					// System.out.println("@@@@@ "+SCDOSE[IZD-1][IXD-1][0]);
					EGS4Grid.RESULTS[IZD - 1][IXD - 1][0] = SCDOSE[IZD - 1][IXD - 1][0];
					// System.out.println("@@@@@ "+SCDOSE[IZD-1][IXD-1][1]);
					EGS4Grid.RESULTS[IZD - 1][IXD - 1][1] = SCDOSE[IZD - 1][IXD - 1][1];
					// System.out.println("@@@@@ "+SCKERMA[IZD-1][IXD-1][0]);
					EGS4Grid.RESULTS[IZD - 1][IXD - 1][2] = SCKERMA[IZD - 1][IXD - 1][0];
					// System.out.println("@@@@@ "+SCDOSE[IZD-1][IXD-1][0]/SCKERMA[IZD-1][IXD-1][0]);
					EGS4Grid.RESULTS[IZD - 1][IXD - 1][3] = SCDOSE[IZD - 1][IXD - 1][0]
							/ SCKERMA[IZD - 1][IXD - 1][0];
					EGS4Grid.UNCRT[IZD - 1][IXD - 1][0] = SCDOSE2[IZD - 1][IXD - 1][0];
					EGS4Grid.UNCRT[IZD - 1][IXD - 1][1] = SCDOSE2[IZD - 1][IXD - 1][1];
					EGS4Grid.UNCRT[IZD - 1][IXD - 1][2] = SCKERMA2[IZD - 1][IXD - 1][0];
					EGS4Grid.UNCRT[IZD - 1][IXD - 1][3] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][0];
					if ((EGS4Macro.IFULL == 1) && (IKERMA == 0)) {
						for (int IT = 3; IT <= 6; IT++) {
							EGS4Grid.RESULTS[IZD - 1][IXD - 1][IT - 2] = FRACS[IZD - 1][IXD - 1][IT - 1];
							EGS4Grid.UNCRT[IZD - 1][IXD - 1][IT - 2] = SCDOSE2[IZD - 1][IXD - 1][IT - 1];
						}
					}
					if ((EGS4Macro.IFULL == 1) && (IKERMA == 1)) {
						EGS4Grid.RESULTS[IZD - 1][IXD - 1][1] = SCKERMA[IZD - 1][IXD - 1][0];
						EGS4Grid.RESULTS[IZD - 1][IXD - 1][2] = SCDOSE[IZD - 1][IXD - 1][0]
								/ SCKERMA[IZD - 1][IXD - 1][0];
						EGS4Grid.UNCRT[IZD - 1][IXD - 1][1] = SCKERMA2[IZD - 1][IXD - 1][0];
						EGS4Grid.UNCRT[IZD - 1][IXD - 1][2] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][0];
						for (int IT = 4; IT <= 7; IT++) {
							EGS4Grid.RESULTS[IZD - 1][IXD - 1][IT - 1] = FRACS[IZD - 1][IXD - 1][IT - 2];
							EGS4Grid.UNCRT[IZD - 1][IXD - 1][IT - 1] = SCDOSE2[IZD - 1][IXD - 1][IT - 2];
						}
					}
					if (EGS4Macro.IFULL == 3) {
						EGS4Grid.RESULTS[IZD - 1][IXD - 1][3] = SCKERMA[IZD - 1][IXD - 1][1];
						EGS4Grid.RESULTS[IZD - 1][IXD - 1][4] = SCDOSE[IZD - 1][IXD - 1][0]
								/ SCKERMA[IZD - 1][IXD - 1][0];
						EGS4Grid.RESULTS[IZD - 1][IXD - 1][5] = SCDOSE[IZD - 1][IXD - 1][1]
								/ SCKERMA[IZD - 1][IXD - 1][1];
						EGS4Grid.UNCRT[IZD - 1][IXD - 1][3] = SCKERMA2[IZD - 1][IXD - 1][1];
						EGS4Grid.UNCRT[IZD - 1][IXD - 1][4] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][0];
						EGS4Grid.UNCRT[IZD - 1][IXD - 1][5] = SCDOSEtoKERMA2[IZD - 1][IXD - 1][1];
					}
				}
			}
			// "2) Determine the number of components"
			if (IKERMA == 1) {
				if (EGS4Macro.IFULL == 0) {
					EGS4Grid.NCOMP = 4;
				}
				if (EGS4Macro.IFULL == 1) {
					EGS4Grid.NCOMP = 7;
				}
				if (EGS4Macro.IFULL == 2) {
					EGS4Grid.NCOMP = 4;
				}
				if (EGS4Macro.IFULL == 3) {
					EGS4Grid.NCOMP = 6;
				}
			}
			if (IKERMA == 0) {
				if (EGS4Macro.IFULL == 1) {
					EGS4Grid.NCOMP = 5;
				} else {
					EGS4Grid.NCOMP = 2;
				}
			}

			// "3) Set up the bin indicators"
			for (IX = 1; IX <= NRDOSE + 1; IX++) {
				EGS4Grid.RADIAL_BINS[IX - 1] = EGS4Geom.RCYL[IX + NRDMIN - 1];
			}
			for (IZ = 1; IZ <= NZDOSE + 1; IZ++) {
				EGS4Grid.DEPTH_BINS[IZ - 1] = EGS4Geom.ZPLANE[IZ + NZDMIN - 2];
			}

			// "4) Set up the labels and the explanations"
			// IF(ISOURC=3|ISOURC=21|ISOURC=22|ISOURC=23)[
			// EXPLANATIONS(1)='Total dose (Gray/incident particle)';
			// ]
			// ELSE[
			EGS4Grid.EXPLANATIONS[0] = "Total dose (Gray/incident fluence)";
			// ]
			EGS4Grid.LABELS[0] = "T  :";
			if (EGS4Macro.IFULL == 3) {
				EGS4Grid.EXPLANATIONS[1] = "Scatter dose (after Compton or fluorecent reabsorbed)";
				EGS4Grid.LABELS[1] = "Sca:";
			} else {
				EGS4Grid.EXPLANATIONS[1] = "Total dose minus stoppers";
				EGS4Grid.LABELS[1] = "T-S:";
			}
			EGS4Grid.EXPLANATIONS[2] = "Kerma";
			EGS4Grid.LABELS[2] = "K  :";
			if (EGS4Macro.IFULL != 3) {
				EGS4Grid.EXPLANATIONS[3] = "Dose to kerma";
				EGS4Grid.LABELS[3] = "D/K:";
			} else {
				EGS4Grid.EXPLANATIONS[3] = "Kerma scatter";
				EGS4Grid.LABELS[3] = "Ksc:";
				EGS4Grid.EXPLANATIONS[4] = "Dose to kerma";
				EGS4Grid.LABELS[4] = "D/K:";
				EGS4Grid.EXPLANATIONS[5] = "Dose to kerma scatter";
				EGS4Grid.LABELS[5] = "DsKs";
			}
			if ((EGS4Macro.IFULL == 1) && (IKERMA == 0)) {
				EGS4Grid.EXPLANATIONS[1] = "% dose from front wall (error is % OF %)";
				EGS4Grid.LABELS[1] = "FT :";
				EGS4Grid.EXPLANATIONS[2] = "% dose from outer wall";
				EGS4Grid.LABELS[2] = "OUT:";
				EGS4Grid.EXPLANATIONS[3] = "% dose from back wall";
				EGS4Grid.LABELS[3] = "BK :";
				EGS4Grid.EXPLANATIONS[4] = "% dose from inner wall";
				EGS4Grid.LABELS[4] = "IN :";
			}
			if ((EGS4Macro.IFULL == 1) && (IKERMA == 1)) {
				EGS4Grid.EXPLANATIONS[1] = "Total kerma";
				EGS4Grid.LABELS[1] = "K  :";
				EGS4Grid.EXPLANATIONS[2] = "Total dose to kerma";
				EGS4Grid.LABELS[2] = "D/K:";
				EGS4Grid.EXPLANATIONS[3] = "% dose from front wall (error is % of %)";
				EGS4Grid.LABELS[3] = "FT :";
				EGS4Grid.EXPLANATIONS[4] = "% dose from outer wall";
				EGS4Grid.LABELS[4] = "OUT:";
				EGS4Grid.EXPLANATIONS[5] = "% dose from back wall";
				EGS4Grid.LABELS[5] = "BK :";
				EGS4Grid.EXPLANATIONS[6] = "% dose from inner wall";
				EGS4Grid.LABELS[6] = "IN :";
			}

			// "5) Make the grid"
			// EGS4Grid.ZONEGRID(NRDOSE, NZDOSE, NRDMIN, NZDMIN, EGS4Geom.NZ,
			// RESULTS,
			// UNCRTY, NCOMP, RADIAL_BINS, DEPTH_BINS);//, LABELS,
			// EXPLANATIONS);
			// System.out.println(NCOMP);
			EGS4Grid.ZONEGRID(NRDOSE, NZDOSE, NRDMIN, NZDMIN, EGS4Geom.NZ);// ,
																			// EGS4Grid.NCOMP);

		}// "end IF IOOPTN=0, 2 or 4"

		// IF(ILPIN.NE.6)[ write(iout,'(a)') '\f'; "CALL PRNTER(13,6,IOUT,1);"]
		// IF(IOPLOT.EQ.1)CALL PLOTEN;
		// "CALL PRNTER(13,6,IOUT,0);"

		return;

	}

	/**
	 * Return base 10 logarithm.
	 * @param x x
	 * @return the result
	 */
	public static double DLOG10(double x) {
		double result = Math.log(x);
		result = result / Math.log(10);
		return result;
	}

	// CALL PLOTPH(TITLE,SCPDST,SCPDST2,SCPCUM,SCPCUM2,
	// SCPTOT,SCPTOT2,SCDFEP,SCDFEP2,MAXBIN,SLOTE,BINTOP,
	// IHSTRY,SCOMEG,SCOMEG2,SCPHEN,SCPHEN2);

	// SUBROUTINE PLOTPH(TITLE,PDST,PDSTUN,PCUM,PCUMUN,PTOT,PTOTUN,DFEP,DFEPUN,
	// MAXBIN,SLOTE,BINTOP,IHSTRY,OMEG,OMEGUN,PHEN,PHENUN);
	/**
	 * Plot some useful information (pulse height distribution.
	 */
	public void PLOTPH() {
		String s = "";
		int ll = 0;
		// "GET PEAK IN SPECTRUM"
		double sfac = 0.0;
		double EB = 0.0;
		String SPACE = " ";// /,BAR/'|'/,SYMBOL/$S'*+$-#@'/;
		String BAR = "|";
		String[] SYMBOL = { "*", "+", "$", "-", "#", "@" };
		String[] LINE = new String[61];
		int NHIST = 61;
		double HIST = 61.;
		//int IOUT = 1;
		int ILEV = 0;

		for (int IB = 1; IB <= MAXBIN; IB++) {
			if (SCPDST[IB - 1] > sfac) {
				sfac = SCPDST[IB - 1];
			}
		}
		// WRITE(IOUT,4) IHSTRY;
		s = " " + EGS4.format(IHSTRY, 12) + " HISTORIES ANALYSED";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,5); "HEADER FOR PLOT"
		s = EGS4.format("", 65) + "EBIN";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = " 3" + EGS4.format("", 19) + "2" + EGS4.format("", 19) + "1"
				+ EGS4.format("", 7) + "-LOG" + EGS4.format("", 8) + "0    TOP"
				+ EGS4.format("", 9) + "PDST              CUMULATIVE";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		for (int IB = 1; IB <= MAXBIN; IB++) {
			if (SLOTE > 0.0) {
				EB = IB * SLOTE;
			} else {
				EB = BINTOP[IB - 1];
			}
			for (int J = 1; J <= NHIST; J++) {
				LINE[J - 1] = SPACE;
				// "BLANK ENTIRE LINE EACH TIME"
			}
			ILEV = 0;
			if (SCPDST[IB - 1] > 0.0) {
				Double dbl = new Double(HIST + DLOG10(SCPDST[IB - 1] / sfac)
						* 20. + 0.5);
				ILEV = dbl.intValue();
				// "RUNS FROM 1 TO 61 FOR 0.0 TO 1.0"
			}
			if ((ILEV > 0) && (ILEV <= NHIST)) {
				LINE[ILEV - 1] = SYMBOL[0];
			}
			// "ADD REFERENCE BARS"
			for (int J = 1; J <= 4; J++) {
				ILEV = 20 * J - 19;
				if (LINE[ILEV - 1].compareTo(SPACE) == 0)
					LINE[ILEV - 1] = BAR;
			}
			// WRITE(IOUT,10) LINE,EB,PDST(IB),PDSTUN(IB),PCUM(IB),PCUMUN(IB);
			// 10 FORMAT(1X,61A1,F8.4,1X,2(0PF11.4,' (',F6.3,'%)'));
			s = " ";// +EGS4.format(IHSTRY,12)+" HISTORIES ANALYSED";
			for (int i = 1; i <= 61; i++) {
				s = s + LINE[i - 1];
			}
			s = s + EGS4.format(EB, 8, true) + " "
					+ EGS4.format(SCPDST[IB - 1], 11, true) + " ("
					+ EGS4.format(SCPDST2[IB - 1], 6, true) + "%)"
					+ EGS4.format(SCPCUM[IB - 1], 11, true) + " ("
					+ EGS4.format(SCPCUM2[IB - 1], 6, true) + "%)";
			EGS4.seqStr = s;// \n";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}// "END OF IB LOOP"

		// 4 FORMAT(/I12,' HISTORIES ANALYSED'/);
		// 5 FORMAT(/65X,'EBIN'/' 3',19X,'2',19X,'1',7X,
		// '-LOG',8X,'0 TOP',9X,'PDST CUMULATIVE');
		// WRITE(IOUT,30);
		// 30 FORMAT(/' PLOT NORMALIZED TO PEAK OF ONE, ',
		// 'PDST IS NORMALIZED TO UNIT AREA'/);
		s = " PLOT NORMALIZED TO PEAK OF ONE, "
				+ "PDST IS NORMALIZED TO UNIT AREA";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		double SUM = 0.0;
		for (int IPK = 1; IPK <= 4; IPK++) {
			SUM = SUM + SCDFEP[IPK - 1];
		}
		if (SUM >= 0.00005) {
			// "THERE ARE SOME COUNTS IN THE PEAKS"
			// WRITE(IOUT,40) (DFEP(IPK),DFEPUN(IPK),IPK=1,4);
			// 40 FORMAT(' PEAK EFFICIENCIES PER COUNT IN SPECTRUM'//
			// ' FULL ENERGY PEAK',T30,F10.4,'(+-',F6.3,'%)'/
			// ' SINGLE ESCAPE PK',T30,F10.4,'(+-',F6.3,'%)'/
			// ' DOUBLE ESCAPE PK',T30,F10.4,'(+-',F6.3,'%)'/
			// ' 511 KEV PK' ,T30,F10.4,'(+-',F6.3,'%)'/);
			s = " PEAK EFFICIENCIES PER COUNT IN SPECTRUM";
			EGS4.seqStr = s;// \n";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " FULL ENERGY PEAK";
			ll = s.length();
			ll = 30 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(SCDFEP[0], 10, true) + "(+-"
					+ EGS4.format(SCDFEP2[0], 6, true) + "%)";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " SINGLE ESCAPE PK";
			ll = s.length();
			ll = 30 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(SCDFEP[1], 10, true) + "(+-"
					+ EGS4.format(SCDFEP2[1], 6, true) + "%)";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " DOUBLE ESCAPE PK";
			ll = s.length();
			ll = 30 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(SCDFEP[2], 10, true) + "(+-"
					+ EGS4.format(SCDFEP2[2], 6, true) + "%)";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " 511 KEV PK";
			ll = s.length();
			ll = 30 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(SCDFEP[3], 10, true) + "(+-"
					+ EGS4.format(SCDFEP2[3], 6, true) + "%)";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}
		double UNCERT = 0.0;
		if (EGS4SrcEns.SCOMEG == 0.0) {
			// "THIS IS A PARALLEL EXTERNAL BEAM OR INTERNAL SOURCE CASE"
			// WRITE(IOUT,50) PTOT,PTOTUN;
			// 50 FORMAT(/' FRACTION OF INITIAL PARTICLES (HITTING DETECTOR
			// HOUSING ',
			// '(if parallel beam) WHICH CAUSE PULSE = ',0PF14.4,'(',F6.3,'%)');
			s = "  FRACTION OF INITIAL PARTICLES (HITTING DETECTOR HOUSING "
					+ "(if parallel beam) WHICH CAUSE PULSE = "
					+ EGS4.format(SCPTOT, 14, true) + "("
					+ EGS4.format(SCPTOT2, 6, true) + "%)";
			// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		} else {
			// "EXTERNAL POINT SOURCE GEOMETRY"
			UNCERT = Math.sqrt(SCPTOT2 * SCPTOT2 + EGS4SrcEns.SCOMEG2
					* EGS4SrcEns.SCOMEG2); // "UNCERTAINTY ON RATIO"
			// WRITE(IOUT,60)
			// OMEG,OMEGUN,PTOT/(4.*3.141593),PTOTUN,PTOT/OMEG,UNCERT;
			// 60 FORMAT(/' SOLID ANGLE SUBTENDED BY DETECTOR HOUSING = ',
			// 1PE13.3,'(',0PF5.1,'%)'/
			// /' FRACTION OF PARTICLES INTO 4-PI WHICH CAUSE A PULSE = ',
			// 1PE13.3,'(',0PF5.1,'%)'/
			// /' FRACTION OF PARTICLES INCIDENT ON HOUSING WHICH CAUSE A PULSE
			// =',
			// 1PE13.3,'(',0PF5.1,'%)');
			s = " SOLID ANGLE SUBTENDED BY DETECTOR HOUSING =                    "
					+ EGS4.format(EGS4SrcEns.SCOMEG, 13, true)
					+ "("
					+ EGS4.format(EGS4SrcEns.SCOMEG2, 5, true) + "%)";
			// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			s = " FRACTION OF PARTICLES INTO 4-PI WHICH CAUSE A PULSE =          "
					+ EGS4.format(SCPTOT / (4. * 3.141593), 13, true)
					+ "("
					+ EGS4.format(SCPTOT2, 5, true) + "%)";
			// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			s = " FRACTION OF PARTICLES INCIDENT ON HOUSING WHICH CAUSE A PULSE ="
					+ EGS4.format(SCPTOT / EGS4SrcEns.SCOMEG, 13, true)
					+ "("
					+ EGS4.format(UNCERT, 5, true) + "%)";
			// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}
		// WRITE(IOUT,70) PHEN,PHENUN;//=>SCPHEN,SCPHEN2
		// 70 FORMAT(/' ENERGY DEPOSITED IN DETECTOR PER INITIAL PARTICLE:',
		// 0PF14.5,' MeV (+-',F6.3,'%)');
		s = "  ENERGY DEPOSITED IN DETECTOR PER INITIAL PARTICLE:"
				+ EGS4.format(SCPHEN, 14, true) + "(+-"
				+ EGS4.format(SCPHEN2, 6, true) + "%)";
		// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		/*
		 * CHARACTER SPACE,BAR,LINE(61),SYMBOL(6); REAL*8
		 * PDST($EBIN),PDSTUN($EBIN),PTOT,PTOTUN,DFEP(4),DFEPUN(4),
		 * PCUM($EBIN),PCUMUN($EBIN),omeg,omegun,phen,phenun; $REAL
		 * BINTOP($EBIN),xplot($EBIN),yplot($EBIN), errplot($EBIN),SLOTE;
		 * $LONG_INT ihstry; $INTEGER maxbin; COMIN/IODAT2/; character*1
		 * TITLE(80); character*60 xtitle,ytitle,subtitle,seriestitle,plottitle;
		 * 
		 * $INTEGER NHIST,IOUT,IB,j,ilev,ipk,iunit,icurvenum,iplttype,iaxistype;
		 * $REAL HIST,sfac,eb,sum,uncert,binw,histmin; $INTEGER egs_open_file;
		 * DATA SPACE/' '/,BAR/'|'/,SYMBOL/$S'*+$-#@'/; DATA
		 * NHIST/61/,HIST/61./,IOUT/1/;
		 * 
		 * "********* IK: on many compilers, arguments to the MAX intrinsic MUST
		 * be "              of the same type. "
		 * "DO IB=1,MAXBIN[SFAC = MAX(PDST(IB),SFAC);]" " " OUTPUT RESULTS TO
		 * DISK FILE ON UNIT 9 " "IUNIT=9;"
		 * "OPEN(UNIT=9,file='fort.9',STATUS='UNKNOWN'); "
		 * "needed by absoft compiler for some reason"
		 * 
		 * //"IK: open the plot file using a file name instead of fort.9 "
		 * //iunit = egs_open_file(9,0,1,'.egseff');
		 * 
		 * //WRITE(IUNIT,75) TITLE; //WRITE(IUNIT,80)MAXBIN;
		 * //"Feb 1992 corrected to output cts/MeV in bin and get file right"
		 * //"Jan 99 changed uncertainties to be absolute so plotxvgr can handle"
		 * //DO IB=1,MAXBIN[ //IF(SLOTE.GT.0.0)[EB=IB*SLOTE;BINW=SLOTE;] //ELSE[
		 * //EB=BINTOP(IB);IF(IB=1)[BINW=BINTOP(1);]ELSE[BINW=EB-BINTOP(IB-1);]
		 * //]"end block for SLOTE" //IF(BINW=0.0)[BINW = 1E-9;] //xplot(IB) =
		 * EB; //yplot(IB) = PDST(IB)/BINW; //IF (PDST(IB) ~= 0.0)
		 * [errplot(IB)=PDSTUN(IB)*PDST(IB)/(BINW*100.);] //ELSE [errplot(IB) =
		 * 0.0;] //WRITE(IUNIT,90) EB,PDST(IB)/BINW,errplot(IB); //]
		 * //WRITE(IUNIT,100) (DFEP(IPK),DFEPUN(IPK),IPK=1,4),PTOT,PTOTUN;
		 * //close(iunit);
		 * 
		 * //"Now setup to create a proper xvgr plot of the spectrum on unit 22"
		 * //iunit = egs_open_file(22,0,1,'.plotphd');
		 * 
		 * xtitle='energy/MeV'; ytitle='cts/MeV'; subtitle='pulse height
		 * distribution from dosrznrc'; icurvenum=0; iplttype=1;
		 * "histogram=1, XY-plot=0" histmin=0.0; iaxistype = 0; "no logs"
		 * seriestitle=' '; call
		 * xvgrplot(xplot,yplot,errplot,MAXBIN,icurvenum,seriestitle,
		 * xtitle,ytitle,title,subtitle,iunit,iplttype,histmin,iaxistype);
		 * "iaxistype is not set, is this intended??? " close(iunit);
		 * 
		 * "FORMATS" 4 FORMAT(/I12,' HISTORIES ANALYSED'/); 5
		 * FORMAT(/65X,'EBIN'/' 3',19X,'2',19X,'1',7X, '-LOG',8X,'0
		 * TOP',9X,'PDST CUMULATIVE'); 10 FORMAT(1X,61A1,F8.4,1X,2(0PF11.4,'
		 * (',F6.3,'%)')); 20 FORMAT(/62X,'TOTALS',2X,6(I4,3X)); 30 FORMAT(/'
		 * PLOT NORMALIZED TO PEAK OF ONE, ', 'PDST IS NORMALIZED TO UNIT
		 * AREA'/); 40 FORMAT(' PEAK EFFICIENCIES PER COUNT IN SPECTRUM'// '
		 * FULL ENERGY PEAK',T30,F10.4,'(+-',F6.3,'%)'/ ' SINGLE ESCAPE
		 * PK',T30,F10.4,'(+-',F6.3,'%)'/ ' DOUBLE ESCAPE
		 * PK',T30,F10.4,'(+-',F6.3,'%)'/ ' 511 KEV PK'
		 * ,T30,F10.4,'(+-',F6.3,'%)'/); 50 FORMAT(/' FRACTION OF INITIAL
		 * PARTICLES (HITTING DETECTOR HOUSING ', '(if parallel beam) WHICH
		 * CAUSE PULSE = ',0PF14.4,'(',F6.3,'%)'); 60 FORMAT(/' SOLID ANGLE
		 * SUBTENDED BY DETECTOR HOUSING = ', 1PE13.3,'(',0PF5.1,'%)'/ /'
		 * FRACTION OF PARTICLES INTO 4-PI WHICH CAUSE A PULSE = ',
		 * 1PE13.3,'(',0PF5.1,'%)'/ /' FRACTION OF PARTICLES INCIDENT ON HOUSING
		 * WHICH CAUSE A PULSE =', 1PE13.3,'(',0PF5.1,'%)'); 70 FORMAT(/' ENERGY
		 * DEPOSITED IN DETECTOR PER INITIAL PARTICLE:', 0PF14.5,' MeV
		 * (+-',F6.3,'%)'); 75 FORMAT(1X,80A1); 80 FORMAT(I5,T40,'# energy bins
		 * to follow, bin top,cts/MeV,uncertainty(abs)'); 90
		 * FORMAT(3(E12.4,',')); 100 FORMAT( 1X,'NOW INDIVIDUAL PEAK AND TOTAL
		 * EFFICIENCIES'/ (2E12.4));
		 * 
		 * RETURN;END;
		 */
	}
}
