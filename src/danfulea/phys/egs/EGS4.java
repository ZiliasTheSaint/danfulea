package danfulea.phys.egs;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

/**
 * Utility class for Monte-Carlo simulation of photon-electron transport. Based on 
 * EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada). 
 * @author Dan Fulea, 08 SEPT. 2005
 */

public class EGS4 {
	public static EgsQuestion eq;// sometimes may be necessary!!!
	// ------------------variables--------------------------------
	public static int $MXMED = 10; // "MAX. NO. OF DIFFERENT MEDIA (EXCL. VAC.)"
	public static int $MXREG = 2000; // "MAXIMUM NO. OF REGIONS ALLOCATED"
	public static int $MXSTACK = 40; // "STACK SIZE"
	// PARAMETER $MXVRT1=1000; //"NO. OF REPRESENTATIVE ANGLES IN VERT1"
	// PARAMETER $FMXVRT1=1000.; //"FLOATING $MXVRT1"
	public static int $MXPWR2I = 50; // "SIZE OF TABLE OF INVERSE POWERS OF TWO"
	// PARAMETER $MXJREFF=200; //"SIZE OF MULTIPLE SCATTERING JREFF MAP"
	// PARAMETER $MSSTEPS=16; //"NO. OF MULTIPLE SCATTERING STEP SIZES"
	// PARAMETER $MXVRT2=100; //"DISTRBTN OF NONOVERLAPPING PARTS OF VERT"
	// "FOLLOWING DEFINE MAX. NO. OF INTERVALS FOR FIT TO FUNCTIONS"
	public static int $MXSGE = 400; // "GAMMA SMALL ENERGY INTERVALS"
	public static int $MXGE = 200; // "GAMMA MAPPED ENERGY INTERVALS"
	// PARAMETER $MXSEKE=300; //"ELECTRON SMALL ENERGY INTERVALS"
	public static int $MXEKE = 150; // "ELECTRON MAPPED ENERGY INTERVALS"
	public static int $MXEL = 50; // "MAXIMUM # OF ELEMENTS IN A MEDIUM"
	// PARAMETER $MXLEKE=100; //"ELECTRON ENERGY INTERVALS BELOW EKELIM"
	// PARAMETER $MXCMFP=100; //"CUMULATIVE ELECTRON MEAN-FREE-PATH"
	public static int $MXRANGE = 100; // "ELECTRON RANGE"
	// PARAMETER $MXBLC=20; //"MOLIERE'S LOWER CASE B"
	// PARAMETER $MXRNTH=20; //"RANDOM NUMBER FOR SUBMEDIAN ANGLES"
	// PARAMETER $MXRNTHI=20; //"RANDOM NUMBER FOR SUPERMEDIAN ANGLES"
	public static int $MXRAYFF = 100; // "RAYLEIGH ATOMIC FORM FACTOR"
	public static int $MXSINC = 1002; // "ANGLE INTERVALS IN (0,5*PI/2) FOR SINE"
	// "THE FOLLOWING ARE PARAMETERS FOR AUSGAB CALLS:"
	public static int $MXAUS = 29; // "CHANGE IF MORE AUSGAB CALLS ARE ADDED"
	public static int $MXAUSM5 = 24; // "SET THIS TO $MXAUS VALUE LESS 5"
	// "FIRST FIVE AUSGAB CALLS BELOW TURNED ON IN BLOCK DATA (DEFAULT)"
	public static int $TRANAUSB = 0; // "BEFORE TRANSPORT"
	public static int $EGSCUTAUS = 1; // "ENERGY BELOW ECUT OR PCUT"
	public static int $PEGSCUTAUS = 2; // "ENERGY BELOW AE OR AP"
	public static int $USERDAUS = 3; // "USER REQUESTED DISCARD"
	public static int $PHOTXAUS = 4; // "FLUORESCENT PHOTON DISCARD"
	// "THE REMAINING 23 ARE TURNED OFF IN BLOCK DATA (DEFAULT)"
	public static int $TRANAUSA = 5; // "AFTER TRANSPORT"
	public static int $BREMAUSB = 6; // "BEFORE BREMS CALL"
	public static int $BREMAUSA = 7; // "AFTER BREMS CALL"
	public static int $MOLLAUSB = 8; // "BEFORE MOLLER CALL"
	public static int $MOLLAUSA = 9; // "AFTER MOLLER CALL"
	public static int $BHABAUSB = 10; // "BEFORE BHABHA CALL"
	public static int $BHABAUSA = 11; // "AFTER BHABHA CALL"
	public static int $ANNIHFAUSB = 12; // "BEFORE ANNIH CALL"
	public static int $ANNIHFAUSA = 13; // "AFTER ANNIH CALL"
	public static int $ANNIHRAUSB = 28; // "BEFORE POSITRON ANNIH AT REST"
	public static int $ANNIHRAUSA = 14; // "POSITRON ANNIHILATED AT REST"
	public static int $PAIRAUSB = 15; // "BEFORE PAIR CALL"
	public static int $PAIRAUSA = 16; // "AFTER PAIR CALL"
	public static int $COMPAUSB = 17; // "BEFORE COMPT CALL"
	public static int $COMPAUSA = 18; // "AFTER COMPT CALL"
	public static int $PHOTOAUSB = 19; // "BEFORE PHOTO CALL"
	public static int $PHOTOAUSA = 20; // "AFTER PHOTO CALL"
	public static int $UPHIAUSB = 21; // "ENTERED UPHI"
	public static int $UPHIAUSA = 22; // "LEFT UPHI"
	public static int $RAYLAUSB = 23; // "BEFORE RAYLEIGH EVENT"
	public static int $RAYLAUSA = 24; // "AFTER RAYLEIGH EVENT"
	public static int $FLUORTRA = 25; // "A fluorescent transition just occured"
	public static int $COSKROTRA = 26; // "A Coster-Kronig transition just occured"
	public static int $AUGERTRA = 27; // "An Auger transition just occured"

	// "------------------------------------------------------------------"
	// "*** BOUNDS--CUTOFF ENERGIES & VACUUM TRANSPORT DISTANCE           "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/BOUNDS/;} WITH
	// {
	// ;COMMON/BOUNDS/ECUT($MXREG),PCUT($MXREG),VACDST;
	public static double[] ECUT = new double[$MXREG];// "Minimum electron transport energy"
	public static double[] PCUT = new double[$MXREG];// "Minimum photon transport energy"
	public static double VACDST = 1.E8;// "Infinity (1E8)"
	// }

	// "------------------------------------------------------------------"
	// "*** STACK--INFORMATION KEPT ABOUT CURRENT PARTICLES               "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/STACK/;} WITH
	// {;
	// COMMON/STACK/
	// $LGN(E,X,Y,Z,U,V,W,DNEAR,WT,IQ,IR,LATCH($MXSTACK)),
	// LATCHI,NP,NPold;
	// $ENERGY PRECISION
	public static double[] E = new double[$MXSTACK];// "total particle energy"
	public static double[] X = new double[$MXSTACK];
	public static double[] Y = new double[$MXSTACK];
	public static double[] Z = new double[$MXSTACK];// "particle co-ordinates"
	public static double[] U = new double[$MXSTACK];
	public static double[] V = new double[$MXSTACK];
	public static double[] W = new double[$MXSTACK];// "particle direction cosines"
	public static double[] DNEAR = new double[$MXSTACK];// "perpendicular distance to nearest boundary"
	public static double[] WT = new double[$MXSTACK];// "particle weight"
	public static int[] IQ = new int[$MXSTACK];// "charge, -1 for electrons, 0 for photons, 1 for positrons"
	public static int[] IR = new int[$MXSTACK];// "current region"
	public static int[] LATCH = new int[$MXSTACK];// "extra phase space variable"
	public static int LATCHI = 0;// "needed because shower does not pass latch-BLOCK DATA sets 0"
	public static int NP = 0;// "stack pointer"
	public static int NPold = 0; // "stack pointer before an interaction"
	// }
	// REPLACE {;COMIN/MISC/;} WITH
	// {;
	public static double DUNIT = 0.0;// "unit scaling factor"
	public static double[] RHOR = new double[$MXREG];// "density of a given region"
	public static int[] MED = new int[$MXREG];// "medium number for a given region"
	public static int[] IRAYLR = new int[$MXREG];// "Rayleigh switch for a given region"
	// }
	// "------------------------------------------------------------------"
	// "*** MEDIA--NAMES OF MEDIA CURRENTLY BEING USED                    "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/MEDIA/;} WITH
	// {;
	// COMMON/MEDIA/
	// $LGN(RLC,RLDU,RHO,MSGE,MGE,MSEKE,
	// MEKE,MLEKE,MCMFP,MRANGE,//IRAYLM($MXMED)),
	// MEDIA(24,$MXMED),NMED;
	// public static String[][] MEDIA=new String[24][$MXMED];//"media names"
	public static String[] MEDIA = new String[$MXMED];// "media names"
	public static double[] RLC = new double[$MXMED];// "radiation length in centimeters for a given medium"
	public static double[] RLDU = new double[$MXMED];// "radiation length after user scaling over-ride"
	public static double[] RHO = new double[$MXMED];// "mass density of a given medium"
	public static int[] MSGE = new int[$MXMED];// "??? "
	public static int[] MGE = new int[$MXMED]; // "number of photon mapped energy intervals for a given medium"
	public static int[] MSEKE = new int[$MXMED];// "??? "
	public static int[] MEKE = new int[$MXMED]; // "number of e mapped energy intervals for a given medium"
	public static int[] MLEKE = new int[$MXMED];// "??? "
	public static int[] MCMFP = new int[$MXMED];// "??? "
	public static int[] MRANGE = new int[$MXMED];// "??? "
	public static int[] IRAYLM = new int[$MXMED];// "Rayleigh switch for a given medium"
	public static int NMED = 0; // "number of media"
	// }
	// "******************************************************************"
	// "                                                                  "
	// "       transport algorithm related stuff                          "
	// "                                                                  "
	// "******************************************************************"
	// "Macros to denote the various transport algorithms"
	// "These numbers just have to be distinct"
	// "Note that the distributed version of EGSnrc does not include the VMC option"
	public static int $PRESTA_II = 0;
	public static int $PRESTA__I = 1;
	public static int $VMC = 2;
	public static double $SKIN_DEPTH_FOR_BCA = 3.0;
	public static boolean $PRESTA_DEBUG = false;
	public static double $EXACT_BCA_XIMAX = 0.5;
	public static double $INEXACT_BCA_XIMAX = 0.5;// "this is not realy neccessary, "
	// "it remained from Alex's coding"
	public static double $MAX_ELOSS = 0.25;
	public static boolean $SUBSTEP_ELOSS_EVALUATION = false;
	public static double $MAX_SMAX = 1.e10;
	public static double $GLOBAL_ECUT = 0.;
	public static double $GLOBAL_PCUT = 0.;
	public static int $IBRDST_DEFAULT = 1;
	public static int $IBR_NIST_DEFAULT = 0;
	public static int $IPRDST_DEFAULT = 1;;
	public static int $IBCMP_DEFAULT = 1;
	public static int $IEDGFL_DEFAULT = 1;
	public static int $IPHTER_DEFAULT = 1;
	public static int $TRANSPORT_ALGORITHM_DEFAULT = $PRESTA_II;
	public static int $BCA_ALGORITHM_DEFAULT = 0;
	public static boolean $EXACT_BCA_DEFAULT = true;
	public static boolean $SPIN_EFFECTS_DEFAULT = true;
	public static int $IRAYLR_DEFAULT = 0;
	public static int $PAIR_NRC_DEFAULT = 0;
	public static int $TRIPLET_DEFAULT = 0;
	public static int pair_nrc = 0;
	public static int itriplet = 0;

	// "This macro sets the minimum step size for a condensed"
	// "history (CH) step. When the exact BCA is used, the minimum"
	// "CH step is determined by efficiency considerations only"
	// "At about 3 elastic MFP's single scattering becomes more"
	// "efficient than CH and so the algorithm switches off CH"
	// "If one of the various inexact BCA's is invoked, this macro"
	// "provides a simple way to include more sophisticated"
	// "decisions about the maximum acceptable approximated CH step"
	// -----------------------------------------------------------
	// "------------------------------------------------------------------------"
	// "*** COMPTON-DATA -- Incoherent scattering data                          "
	// "------------------------------------------------------------------------"
	public static int $MXTOTSH = 1538;// "Total number of shells for Z=1..100    "
	public static int $MXMDSH = 50;// "Max. number of shells per medium       "
	// REPLACE {;COMIN/COMPTON-DATA/;} WITH
	// {
	public static int[] iz_array = new int[$MXTOTSH];// "Atomic number for each shell"
	public static double[] be_array = new double[$MXTOTSH];// "Shell binding energies      "
	public static double[] Jo_array = new double[$MXTOTSH];// "Compton profile parameter   "
	public static double[] erfJo_array = new double[$MXTOTSH];// "needed for the calculation"
	// "of the incoherent scattering function
	public static int[] ne_array = new int[$MXTOTSH]; // "Occupation number           "
	public static int[] shn_array = new int[$MXTOTSH];// "shell type                  "
	// "(=1 for K;=2,3,4 for L1,L2,L3; =5 for M;=6 for N; =7     for all others      "
	public static int[][] shell_array = new int[$MXMDSH][$MXMED];
	public static double[][] eno_array = new double[$MXMDSH][$MXMED];
	public static int[] n_shell = new int[$MXMED];
	public static int[] ibcmp = new int[$MXREG];// "flag to turn on binding effects"
	// }
	// "------------------------------------------------------------------ "
	// "*** EDGE -- Containes binding energies for K,L1,L2,L3,             "
	// "             'average' M and 'average' N shells; photo-absorption  "
	// "             interaction probabilities with these shells;          "
	// "             + fluorescence, Auger, Coster-Kronig transition       "
	// "             probabilities                                         "
	// "             IEDGFL is a flag for turning on/off atomic relaxations"
	// "             IPHTER is a flag for turning on/off photo-lectron     "
	// "                    angular distribution                           "
	// "             both are left-overs from the previous coding          "
	// "             Have put now also data to calculate elemental PE      "
	// "             cross sections needed to sample the element the photon"
	// " is interacting with.
	// "------------------------------------------------------------------ "
	public static int $MXELEMENT = 100;// " Number of elements               "
	public static int $MXSHELL = 6;// " Number of shells treated         "
	public static int $MXINTER = 5;// " $MXSHELL-1                       "
	public static int $MXTRANS = 39;// " Number of possible transitions   "
	public static int $MXEDGE = 16;// " max. number of edges above 1 keV "
	// REPLACE {;COMIN/EDGE/;} WITH
	// {;
	// " K,L1,L2,L3,M,N binding energies  "
	public static double[][] binding_energies = new double[$MXSHELL][$MXELEMENT];
	// " prob. for interaction with one of the above shells (provided photon"
	// " energy is above be)              "
	public static double[][] interaction_prob = new double[$MXSHELL][$MXELEMENT];
	// " relaxation probabilities         "
	public static double[][] relaxation_prob = new double[$MXTRANS][$MXELEMENT];
	// " photo-absorption edge energies   "
	public static double[][] edge_energies = new double[$MXEDGE][$MXELEMENT];
	// " number of `edges' for each element"
	public static int[] edge_number = new int[$MXELEMENT];
	// " photo cross section fit parameters "
	public static double[][] edge_a = new double[$MXEDGE][$MXELEMENT];
	public static double[][] edge_b = new double[$MXEDGE][$MXELEMENT];
	public static double[][] edge_c = new double[$MXEDGE][$MXELEMENT];
	public static double[][] edge_d = new double[$MXEDGE][$MXELEMENT];
	// "flag for switching on fluorscent emission"
	public static int[] iedgfl = new int[$MXREG];
	// "flag for switching on photo-electron angular distr."
	public static int[] iphter = new int[$MXREG];
	// }
	// "***************************************************************************"
	// " ------------ common block for EII data -----------------                  "
	// " Added by Iwan Kawrakow, March 20 2004.
	// "****************************************************************************
	public static int $MAX_EII_SHELLS = 40;// "Maximum number of shells participating"
											// "in EII in a simulation                "
	public static int $N_EII_BINS = 250;// "Number of bins for EII x-section
										// interpolations
	public static int $MAX_EII_BINS = $N_EII_BINS * $MAX_EII_SHELLS;
	// "We store the EII x-section interpolation coefficients in 1D arrays  "
	// "The above is the dimension of these arrays required to hold the data"
	// REPLACE {;COMIN/EII-DATA/;} WITH {;
	// "EII x-section interpolation coeff."
	public static double[] eii_xsection_a = new double[$MAX_EII_BINS];
	public static double[] eii_xsection_b = new double[$MAX_EII_BINS];
	public static int[] eii_z = new int[$MAX_EII_SHELLS];// "Z of each shell                   "
	public static int[] eii_sh = new int[$MAX_EII_SHELLS];// "shell type (1=K, 2=LI, eyc.)      "
	// "energy grid coeff. for each shell "
	public static double[] eii_a = new double[$MAX_EII_SHELLS];
	public static double[] eii_b = new double[$MAX_EII_SHELLS];
	public static int[] eii_nshells = new int[$MXELEMENT];// "No. of EII shells for each element"
	public static int[] eii_nsh = new int[$MXMED];// "No. of EII shells for each medium "
	public static double[] eii_cons = new double[$MXMED];
	// "First EII shell in the list of shells"
	public static int[][] eii_first = new int[$MXMED][$MXEL];
	public static int[][] eii_no = new int[$MXMED][$MXEL];// "N. of EII shells     "
	public static int eii_flag = 0; // "EII flag    "
	// = 0 => no EII;= 1 => simple EII;> 1 => future use "
	// };
	// " The following common block is made available to the user so that  "
	// " he/she knows which shell was being relaxed when the call to ausgab"
	// " occured " Added by Iwan Kawrakow, March 22 2004. "
	// REPLACE {;COMIN/RELAX-USER/;} WITH {;
	public static double u_relax = 0.0;
	public static int ish_relax = 0;
	public static int iZ_relax = 0;
	// };
	// REPLACE {;COMIN/ET-Control/;} WITH "ET stands for Electron Transport"
	// {
	// "geom. step-size constrain for each region"
	public static double[] SMAXIR = new double[$MXREG];
	public static double estepe = 0.0;// "global energy loss constrain"
	public static double ximax = 0.0;// "max. first GS moment per step"
	// "ximin_for_bca,"
	// "distance from a boundary (in elastic MFP)to switch to one of the BCAs "
	public static double skindepth_for_bca = 0.0;
	public static int transport_algorithm = 0;// //"=$PRESTA-II or $PRESTA--I"
	public static int bca_algorithm = 0;// "will be used if other inexact BCAs"
										// "implemented in the future"
	public static boolean exact_bca = true;// "if .true. => BCA in single scattering mode"
	public static boolean spin_effects = true;// "if .true. electron/positron spin effects"
	// "are taken into account in the single and"
	// "multiple elasting scattering routines"
	// }
	// "***************************************************************************"
	// "         EGSnrc internal Variance Reduction Macros                         "
	// "***************************************************************************"
	// REPLACE {;COMIN/EGS-VARIANCE-REDUCTION/;} WITH {;
	// "max energy at which to do range rejection (RR)"
	public static double[] e_max_rr = new double[$MXREG];
	public static double prob_RR = 0.0;// "probability for survival in R. Roulette"
	public static int nbr_split = 0;// "do brems splitting if > 1"
	public static int i_play_RR = 0;// "0=>don't play Russian Roulette,1=>play Russian Roulette"
	public static int i_survived_RR = 0;// "0=> all particles survive RR, n=> n partilces were"
	// "eliminated by RR in this interaction"
	public static int n_RR_warning = 0;// "a counter for user errors"
	public static int[] i_do_rr = new int[$MXREG];// "0=>no RR, region by region, 1=>there is RR"
	// };
	public static int $MAX_RR_WARNING = 50;

	// "------------------------------------------------------------------"
	// "*** THRESH--THRESHOLD (AND OTHER) ENERGIES                        "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/THRESH/;} WITH
	// {;
	public static double RMT2 = 0.0;// "2*electron mass in MeV"
	public static double RMSQ = 0.0;// "electron mass squared in MeV**2"
	// $LGN(AP,AE,UP,UE,TE,THMOLL($MXMED));
	public static double[] AP = new double[$MXMED];// "photon creation threshold energy"
	public static double[] AE = new double[$MXMED];// "electron creation threshold energy (total)"
	public static double[] UP = new double[$MXMED]; // "upper photon energy in PEGS4 data set"
	public static double[] UE = new double[$MXMED]; // "upper electron energy in PEGS4 data set"
	public static double[] TE = new double[$MXMED]; // "electron creation threshold energy (kinetic)"
	public static double[] THMOLL = new double[$MXMED];// "Moller threshold = AE + TE"
	// }
	// "------------------------------------------------------------------"
	// "*** USEFUL--HEAVILY USED VARIABLES                                "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/USEFUL/;} WITH
	// {;
	public static double PZERO = 0.0;// "precise zero"
	public static double PRM = 0.0;// "precise electron mass in MeV"
	public static double PRMT2 = 0.0;// "2*PRM"
	public static double RM = 0.0;// "electron mass in MeV"
	public static int MEDIUM = 0;// "medium index of current region"
	public static int MEDOLD = 0;// "medium index of previous region"
	// }

	// "------------------------------------------------------------------"
	// "*** BREMPR--BREMSSTRAHLUNG AND PAIR PRODUCTION DATA               "
	// "------------------------------------------------------------------"
	public static int $MXBREN = 57;
	public static int $MXBRXX = 30;
	public static int $MXBREL = 100;
	public static int $MXGAUSS = 64;
	public static int $MXBRES = 100;
	public static int $MXBRXS = 50;
	// REPLACE {;COMIN/BREMPR/;} WITH
	// {
	// $LGN(DL(8,$MXMED)/1,2,3,4,5,6/),
	// "Parameter for the fit of the screening rejection function, eq. (2.7.14 and 15)"
	public static double[][] DL1 = new double[8][$MXMED];// 8???see PEGS4
	public static double[][] DL2 = new double[8][$MXMED];
	public static double[][] DL3 = new double[8][$MXMED];
	public static double[][] DL4 = new double[8][$MXMED];
	public static double[][] DL5 = new double[8][$MXMED];
	public static double[][] DL6 = new double[8][$MXMED];
	// "Prob. for the (1-BR)/BR part in BREMS, eq. (2.7.64)"
	public static double[][] ALPHI = new double[2][$MXMED];
	// "Prob. for the 12*(BR-1/2)**2 part in PAIR, eq. (2.7.105)"
	public static double[][] BPAR = new double[2][$MXMED];
	// "maximum delta, eq. (2.7.31)"
	public static double[][] DELPOS = new double[2][$MXMED];
	public static double[][] WA = new double[$MXMED][$MXEL];// "atomic weight"
	// "atomic fraction of an element in a compound"
	public static double[][] PZ = new double[$MXMED][$MXEL];
	public static double[][] ZELEM = new double[$MXMED][$MXEL];// "Z for a given component"
	// "density of an element in a compound"
	public static double[][] RHOZ = new double[$MXMED][$MXEL];
	// "powers of 1/2 (used for sampling (1-BR)/BR"
	public static double[] PWR2I = new double[$MXPWR2I];
	public static double[] DELCM = new double[$MXMED];// "136*m*exp(Zg), eq. (2.7.51)"
	// "composite factor for angular distributions"
	public static double[] ZBRANG = new double[$MXMED];
	public static double[] LZBRANG = new double[$MXMED];// "-Log(ZBRANG)"
	public static int[] NNE = new int[$MXMED];// "number of elements/compound"
	// "flag to switch on bremsstrahlung angular distributions"
	public static int ibrdst = 0;
	// "flag to switch on pair angular distributions"
	public static int iprdst = 0;
	// "use the NIST bremsstrahlung cross sections"
	public static int ibr_nist = 0;
	// ASYM($MXMED,$MXEL,2);->??????????????????
	public static String[][] ASYM = new String[$MXMED][$MXEL];
	// };
	// REPLACE {;COMIN/NIST-BREMS/;} WITH {;
	public static double[][][] nb_fdata = new double[$MXBRXS + 1][$MXBRES][$MXMED];
	public static double[][][] nb_xdata = new double[$MXBRXS + 1][$MXBRES][$MXMED];
	public static double[][][] nb_wdata = new double[$MXBRXS][$MXBRES][$MXMED];
	public static int[][][] nb_idata = new int[$MXBRXS][$MXBRES][$MXMED];
	public static double[] nb_emin = new double[$MXMED];
	public static double[] nb_emax = new double[$MXMED];
	public static double[] nb_lemin = new double[$MXMED];
	public static double[] nb_lemax = new double[$MXMED];
	public static double[] nb_dle = new double[$MXMED];
	public static double[] nb_dlei = new double[$MXMED];
	public static double[] log_ap = new double[$MXMED];
	// common/nist_brems/ nb_fdata(0:$MXBRXS,$MXBRES,$MXMED),//???????????????
	// nb_xdata(0:$MXBRXS,$MXBRES,$MXMED),
	// };
	// "------------------------------------------------------------------"
	// "*** EPCONT--ELECTRON-PHOTON CONTROL VARIABLES                     "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/EPCONT/;} WITH
	// {;
	public static double EDEP = 0.0;// "energy deposition in MeV"
	public static double TSTEP = 0.0;// "distance to a discrete interaction"
	public static double TUSTEP = 0.0;// "intended step length, befor check with geometry"
	public static double USTEP = 0.0;// "transport distance calculated from TUSTEP"
	public static double VSTEP = 0.0;// "transport distance after truncation by HOWFAR"
	public static double TVSTEP = 0.0;// "curved path-length calculated from TVSTEP"
	public static double RHOF = 0.0;// "mass density ratio"
	public static double EOLD = 0.0;// "energy before deduction of energy loss"
	public static double ENEW = 0.0;// "energy after  deduction of energy loss"
	public static double EKE = 0.0;// "kinetic energy"
	public static double ELKE = 0.0;// "Log(EKE)"
	public static double GLE = 0.0;// "Log(energy) in PHOTON"
	public static double E_RANGE = 0.0;// "range of electron before an iarg=0 ausgab call"
	public static int[] iausfl = new int[$MXAUS];// "flags for AUSGAB calls"
	public static int IDISC = 0;// "flag indicating user discard"
	public static int IROLD = 0;// "region before transport"
	public static int IRNEW = 0;// "region after transport"

	// }

	// " ======================== multiple scattering commons ================= "
	// " Screened Rutherford MS data "

	public static int $MAXL_MS = 63;
	public static int $MAXQ_MS = 7;
	public static int $MAXU_MS = 31;
	// REPLACE {$0_MAXL_MS} WITH {0:63};
	// REPLACE {$0-MAXQ_MS} WITH {0:7}
	// REPLACE {$0-MAXU_MS} WITH {0:31}
	public static double $LAMBMIN_MS = 1.;
	public static double $LAMBMAX_MS = 1.e5;
	public static double $QMIN_MS = 1.e-3;
	public static double $QMAX_MS = 0.5;

	public static double[][][] ums_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
	public static double[][][] fms_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
	public static double[][][] wms_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
	public static int[][][] ims_array = new int[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
	// REPLACE {COMIN/MS-Data/;} WITH {
	// common/ms_data/
	// ums_array($0-MAXL_MS,$0-MAXQ_MS,$0-MAXU_MS),
	// fms_array($0-MAXL_MS,$0-MAXQ_MS,$0-MAXU_MS),
	// wms_array($0-MAXL_MS,$0-MAXQ_MS,$0-MAXU_MS),
	// ims_array($0-MAXL_MS,$0-MAXQ_MS,$0-MAXU_MS),
	public static double llammin = 0.0;
	public static double llammax = 0.0;
	public static double dllamb = 0.0;
	public static double dllambi = 0.0;
	public static double dqms = 0.0;
	public static double dqmsi = 0.0;
	// }

	// " spin effect data used in an additional rejection loop "
	public static int $MAXE_SPIN = 15;
	public static int $MAXE_SPI1 = 2 * $MAXE_SPIN + 1;// $MAXE_SPIN=15
	public static int $MAXQ_SPIN = 15;
	public static int $MAXU_SPIN = 31;
	// REPLACE {$0-MAXE_SPI1} WITH {0:$MAXE_SPI1}
	// REPLACE {$0-MAXQ_SPIN} WITH {0:$MAXQ_SPIN}
	// REPLACE {$0-MAXU_SPIN} WITH {0:$MAXU_SPIN}

	// REPLACE {COMIN/Spin-Data/;} WITH {
	// common/spin_data/
	// spin_rej($MXMED,0:1,$0-MAXE_SPI1,$0-MAXQ_SPIN,$0-MAXU_SPIN),
	public static double[][][][][] spin_rej = new double[$MXMED][2][$MAXE_SPI1 + 1][$MAXQ_SPIN + 1][$MAXU_SPIN + 1];
	public static double espin_min = 0.0;
	public static double espin_max = 0.0;
	public static double espml = 0.0;
	public static double b2spin_min = 0.0;
	public static double b2spin_max = 0.0;
	public static double dbeta2 = 0.0;
	public static double dbeta2i = 0.0;
	public static double dlener = 0.0;
	public static double dleneri = 0.0;
	public static double dqq1 = 0.0;
	public static double dqq1i = 0.0;
	public static boolean fool_intel_optimizer = false;
	// }

	// REPLACE {COMIN/CH-Steps/;} WITH
	// {
	public static double count_pII_steps = 0.0;
	public static double count_all_steps = 0.0;
	public static boolean is_ch_step = false;
	// }

	// " default numer of media. "
	public static int $default_nmed = 1;
	// "------------------------------------------------------------------"
	// "*** UPHIIN--SINE TABLES FOR UPHI                                  "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/UPHIIN/;} WITH
	// {;
	// COMMON/UPHIIN/SINC0,SINC1,$LGN(SIN($MXSINC)/0,1/);
	public static double SINC0 = 0.0;
	public static double SINC1 = 0.0;
	public static double[] SIN0 = new double[$MXSINC];
	public static double[] SIN1 = new double[$MXSINC];
	// }

	// "------------------------------------------------------------------"
	// "*** UPHIOT--UPHI'S INPUT/OUTPUT WITH ITS USERS                    "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/UPHIOT/;} WITH
	// {;
	public static double THETA = 0.0;// "polar scattering angle"
	public static double SINTHE = 0.0;// "sin(THETA)"
	public static double COSTHE = 0.0;// "cos(THETA)"
	public static double SINPHI = 0.0;// "sine of the azimuthal scattering angle"
	public static double COSPHI = 0.0;// "cosine of the azimuthal scattering angle"
	public static double PI = 0.0;//
	public static double TWOPI = 0.0;//
	public static double PI5D2 = 0.0;//
	// }
	// "------------------------------------------------------------------"
	// "*** ELECIN--ELECTRON TRANSPORT INPUT                              "
	// "        MODIFIED 1989/12/19 TO INCLUDE IUNRST,EPSTFL AND IAPRIM   "
	// "        NRC DWOR                                                  "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/ELECIN/;} WITH
	// {
	// "maximum electron cross section per energy loss for each medium"
	public static double[] esig_e = new double[$MXMED];
	// "maximum positron cross section per energy loss for each medium"
	public static double[] psig_e = new double[$MXMED];
	// "maximum electron cross section per energy loss"
	public static double esige_max = 0.0;
	// "maximum positron cross section per energy loss"
	public static double psige_max = 0.0;
	// range_ep(0:1,$MXEKE,$MXMED),//"electron (0) or positron (1) range"
	public static double[][][] range_ep = new double[2][$MXEKE][$MXMED];
	public static double[][] E_array = new double[$MXEKE][$MXMED];// "table energies"
	// "for interpolation of screening parameter (e-)"
	public static double[][] etae_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] etae_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of screening parameter (e+)"
	public static double[][] etap_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] etap_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of q1 correction due to spin (e-)"
	public static double[][] q1ce_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] q1ce_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of q1 correction due to spin (e+)"
	public static double[][] q1cp_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] q1cp_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of q2 correction due to spin (e-)"
	public static double[][] q2ce_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] q2ce_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of q2 correction due to spin (e+)"
	public static double[][] q2cp_ms0 = new double[$MXEKE][$MXMED];
	public static double[][] q2cp_ms1 = new double[$MXEKE][$MXMED];
	// "for interpolation of scattering power correction necessary to account
	// for
	// scattering already taken into account in discrete Moller/Bhabha
	public static double[][] blcce0 = new double[$MXEKE][$MXMED];
	public static double[][] blcce1 = new double[$MXEKE][$MXMED];
	public static double[] expeke1 = new double[$MXMED]; // "Exp(1/eke1)-1"
	// "table for kinetic energy indexing"
	public static double[] EKE0 = new double[$MXMED];
	public static double[] EKE1 = new double[$MXMED];
	public static double[] XR0 = new double[$MXMED];// "unused, but read in HATCH"
	public static double[] TEFF0 = new double[$MXMED];// "unused, but read in HATCH"
	public static double[] BLCC = new double[$MXMED];// "b lower case sub c"
	public static double[] XCC = new double[$MXMED];// "chi sub-c-c"
	// "used for electron cross section interpolation"
	public static double[][] ESIG0 = new double[$MXEKE][$MXMED];
	public static double[][] ESIG1 = new double[$MXEKE][$MXMED];
	// "used for positron cross section interpolation"
	public static double[][] PSIG0 = new double[$MXEKE][$MXMED];
	public static double[][] PSIG1 = new double[$MXEKE][$MXMED];
	// "used for electron dE/dx interpolation"
	public static double[][] EDEDX0 = new double[$MXEKE][$MXMED];
	public static double[][] EDEDX1 = new double[$MXEKE][$MXMED];
	// "used for positron dE/dx interpolation"
	public static double[][] PDEDX0 = new double[$MXEKE][$MXMED];
	public static double[][] PDEDX1 = new double[$MXEKE][$MXMED];
	// "used for e- branching into brems interpolation"
	public static double[][] EBR10 = new double[$MXEKE][$MXMED];
	public static double[][] EBR11 = new double[$MXEKE][$MXMED];
	// "used for e+ branching into brems interpolation"
	public static double[][] PBR10 = new double[$MXEKE][$MXMED];
	public static double[][] PBR11 = new double[$MXEKE][$MXMED];
	// "used for e+ branching into Bhabha interpolation"
	public static double[][] PBR20 = new double[$MXEKE][$MXMED];
	public static double[][] PBR21 = new double[$MXEKE][$MXMED];
	// "used for maximum step-size interpolation"
	public static double[][] TMXS0 = new double[$MXEKE][$MXMED];
	public static double[][] TMXS1 = new double[$MXEKE][$MXMED];
	// "flag for type of stopping power (see PEGS4)"
	public static int[] IUNRST = new int[$MXMED];
	// "flag for ICRU37 collision stopping powers"
	public static int[] EPSTFL = new int[$MXMED];
	// "flag for ICRU37 radiative stopping powers"
	public static int[] IAPRIM = new int[$MXMED];
	// }
	// "------------------------------------------------------------------"
	// "*** PHOTIN--PHOTON TRANSPORT DATA                                 "
	// "------------------------------------------------------------------"
	// REPLACE {;COMIN/PHOTIN/;} WITH
	// {;
	// COMMON/PHOTIN/
	// "energy of the K-edge for a given medium"
	public static double[] EBINDA = new double[$MXMED];
	// "used for indexing in logarithmic interpolations"
	public static double[] GE0 = new double[$MXMED];
	public static double[] GE1 = new double[$MXMED];
	// "used for gamma MFP interpolation"
	public static double[][] GMFP0 = new double[$MXGE][$MXMED];
	public static double[][] GMFP1 = new double[$MXGE][$MXMED];
	// "used for branching into pair interpolation"
	public static double[][] GBR10 = new double[$MXGE][$MXMED];
	public static double[][] GBR11 = new double[$MXGE][$MXMED];
	// "used for branching into Compton interpolation"
	public static double[][] GBR20 = new double[$MXGE][$MXMED];
	public static double[][] GBR21 = new double[$MXGE][$MXMED];
	// "used for indexing in momentum trans. sampling in Rayleigh"
	public static double[] RCO0 = new double[$MXMED];
	public static double[] RCO1 = new double[$MXMED];
	// "used for interpolation of momentum trans. func. in R"
	public static double[][] RSCT0 = new double[$MXRAYFF][$MXMED];
	public static double[][] RSCT1 = new double[$MXRAYFF][$MXMED];
	// "used for Rayleigh modification interpolation"
	public static double[][] COHE0 = new double[$MXGE][$MXMED];
	public static double[][] COHE1 = new double[$MXGE][$MXMED];
	// "number of MFP's to go to the next interaction"
	public static double DPMFP = 0.0;
	public static int[][] MPGEM = new int[$MXSGE][$MXMED];
	// "array size for Rayleigh scattering data"
	public static int[] NGR = new int[$MXMED];
	// }

	// -------------------
	public static int IRAYL = 0;// "Rayleigh switch read in from PEGS"
	public static int NEKE = 0;// "array size input from PEGS."
	public static int I1ST = 1;// "flag = 0 on first pass"
	public static boolean got_data = false;
	// --------------------
	private static final String datas = "Data";
	private static final String egsDatas = "egs";
	// private static final String userdirs = System.getProperty("user.dir");
	private static final String file_sep = System.getProperty("file.separator");
	private static final String intDir = "interactiv";
	private static final String pegs4datext = ".pegs4dat";
	private static final String defaultext = ".data";
	private static final String msfile = "msnew";
	private static final String spinDir = "V3_spinms";
	private static final String zfile = "z";
	private final static String z00file = "00";
	private static final String z0file = "0";
	private static final String phrelaxfile = "photo_relax";
	private static final String phcsfile = "photo_cs";
	private static final String incohfile = "incoh";
	private static final String nistbremsfile = "nist_brems";
	public static String eiifile = "eii_ik";// default;
	private static final String radcompton = "rad_compton1";
	private static final String tripletfile = "triplet";
	private static final String spinfile = "spinms.data";
	private static final String pair_nrc_file = "pair_nrc1.data";

	public static double tperp = 0.0;// USED in HOWNEAR!!!!!!!!!!!!!!!!!!!!!!
	public static int irl = 0;// USED in HOWNEAR!!!!!!!!!!!!!!!!!!!!!!
	public static boolean STOPPROGRAM = false;
	public static int RandomUse = 0;// if 0- old Java random
	// if 1- EGS based RANLUX random generator, if 2-a pure java Math random
	public static boolean inEMF = false;// if is an electromagnetic field
										// present
	public static int i_EMF = 0;// default em field configuration-> cylinder of
								// 0.5 cm radius
	public static long startSimulationTime = 0;// a java used for time elapsed!!
	// ***************RANDOM
	// data****************************************************************
	public static int rng_seed = 999999;// "current pointer for rng_array            "
	public static double[] rng_array = new double[24];// "containes 24 random numbers              "
	public static int[] seeds = new int[24];// "for storing the rng state                "
	public static int seedin = 0;
	public static int luxury_level = 0;
	public static int[] state = new int[25];
	public static int carry = 0;// " The state of the generator "
	public static int i24 = 0;
	public static int j24 = 0;// " The rng seeds "
	public static int[] next = new int[24];// " for convinience "
	public static int jseed_dflt = 314159265;// 0;
	public static int nskipRnd = 0;
	public static int icon = 2147483563;// 0;
	public static int status = 0;
	public static int jseed = 0;
	public static int[] nskipll = { 0, 24, 73, 199, 365 };// new int[5];//(0:4),
	public static int icarry = 0;
	public static int kRnd = 0;
	public static int jRnd = 0;
	public static boolean not_initialized = true;
	public static int uni = 0;
	public static double twom24 = 0.0;
	public static double twop24 = 0.0;
	public static int $DEFAULT_LL = 1;// REPLACE {$DEFAULT-LL} WITH {1}
	public static int $MXRNGDIM = 5;// REPLACE {$MXRNGDIM} WITH {5}
	// ranmar
	public static int[] urndm = new int[97];
	public static int crndm = 0;
	public static int cdrndm = 0;
	public static int cmrndm = 0;
	public static int i4opt = 0;
	public static int ixx = 0;
	public static int jxx = 0;
	public static int fool_optimizer = 0;
	// *********************************************************************************************
	// ******************AUX
	// data****************************************************************
	public static int $MXDATA = 10;// 10 arrays of date to handle
	public static int $STAT = 200;// maximum data available per array
	// ex Data[n][i] means n siruri de cate i date maxim fiecare
	public static double[][] DATA = new double[$MXDATA][$STAT];
	// at output data will contain for each array, [n][0] for avg and and [n][1]
	// for errors
	// private static double[] wws_array;
	// private static double[] ffs_array;
	// private static int[] iibin_array;
	public static int IERRsgm = 0;
	public static int ICOUNT = 0;
	public static int JHSTRY = 1;
	public static double KE = 0.0;
	public static NumberFormat nf = NumberFormat.getInstance(Locale.US);// NumberFormat.getInstance();
	public static String pattern = "0.###E0";
	// ====================
	public static DecimalFormatSymbols dfs = new DecimalFormatSymbols(Locale.US);
	// ======================
	public static DecimalFormat nff = new DecimalFormat(pattern, dfs);// new
																		// DecimalFormat(pattern);
	public static int idigits = 3;
	// rad comptom
	public static int radc_flag = 0;
	public static int $RADC_NE = 128; // "Number of initialized energies"
	// REPLACE {$0-RADC_NE} WITH {0:128}
	public static int $RADC_NU = 32;// "Number of angular bins"
	// REPLACE {$0-RADC_NU} WITH {0:32}
	public static int $RADC_NBOX = 13917;// "Number of cells (boxes) for all energies"
	public static int $RADC_NX = 8929; // "Number of cell boundaries for all energies"
	/*****************************************************************************
	 * Common block with pre-calculated data read from rad_compton1.data
	 ****************************************************************************/
	// REPLACE {;COMIN/RAD_COMPTON/;} WITH {;
	// "single Compton total cross section in units of KN"
	public static double[] radc_sigs = new double[$RADC_NE + 1];// ($0-RADC_NE)
	// "double Compton total cross section in units of KN"
	public static double[] radc_sigd = new double[$RADC_NE + 1];// ($0-RADC_NE),
	// "rejection function for single Compton scattering"
	public static double[][] radc_frej = new double[$RADC_NE + 1][$RADC_NU + 1];// ($0-RADC_NE,$0-RADC_NU),
	// "box boundaries"
	public static double[] radc_x = new double[$RADC_NX];
	// "probabilities for double Compton scattering"
	public static double[] radc_fdat = new double[$RADC_NBOX];
	// "maximums for double Compton scattering"
	public static double[] radc_Smax = new double[$RADC_NBOX];
	public static double radc_emin = 0.0;// "minimum initialized enrgy"
	public static double radc_emax = 0.0;// "maximum initialized enrgy"
	public static double radc_lemin = 0.0;// "ln(radc_emin)"
	public static double radc_dw = 0.0;// "soft photon threshold"
	public static double radc_dle = 0.0;// "energy grid bin width"
	public static double radc_dlei = 0.0;// "inverse energy grid bin width"
	public static double radc_le1 = 0.0;

	public static int[] radc_bins = new int[$RADC_NBOX];
	// "indexes of box boundaries"
	public static int[] radc_ixmin1 = new int[$RADC_NBOX];
	public static int[] radc_ixmax1 = new int[$RADC_NBOX];
	public static int[] radc_ixmin2 = new int[$RADC_NBOX];
	public static int[] radc_ixmax2 = new int[$RADC_NBOX];
	public static int[] radc_ixmin3 = new int[$RADC_NBOX];
	public static int[] radc_ixmax3 = new int[$RADC_NBOX];
	public static int[] radc_ixmin4 = new int[$RADC_NBOX];
	public static int[] radc_ixmax4 = new int[$RADC_NBOX];
	public static int[] radc_startx = new int[$RADC_NE + 1];// $0-RADC_NE);
	public static int[] radc_startb = new int[$RADC_NE + 1];// $0-RADC_NE);
	// };
	public static int $N_TRIPLET_DATA = 55;
	public static int $N_ELEMENT = 100;
	public static int $MAX_TRIPLET = 250;
	public static double[][] a_triplet = new double[$MAX_TRIPLET][$MXMED];
	public static double[][] b_triplet = new double[$MAX_TRIPLET][$MXMED];
	public static double dl_triplet = 0.0;
	public static double dli_triplet = 0.0;
	public static double bli_triplet = 0.0;
	public static double log_4rm = 0.0;
	// pair nrc
	public static int $NRC_PAIR_NXX = 65;
	public static int $NRC_PAIR_NEE = 84;
	public static int $NRC_PAIR_NX_1 = 64;
	public static int $NRC_PAIR_NE_1 = 83;
	// //($NRC-PAIR-NXX,$NRC-PAIR-NEE,$MXMED)
	public static int[][][] nrcp_idata = new int[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
	public static double[] nrcp_xdata = new double[$NRC_PAIR_NXX];
	public static double[][][] nrcp_fdata = new double[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
	public static double[][][] nrcp_wdata = new double[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
	public static double nrcp_emin = 0.0;
	public static double nrcp_emax = 0.0;
	public static double nrcp_dle = 0.0;
	public static double nrcp_dlei = 0.0;
	// -----------------PAssed var for macro----------------------
	public static int $NRANMAR = 128;
	public static double[] rng_array1 = new double[$NRANMAR];
	public static boolean ranluxB = false;// use ranlux if true, ranmar
											// otherwise
	public static int hatchindex = 0;// default
	public static int USER_CONTROLS_TSTEP_RECURSION = 0;// default
	public static boolean nextB = false;

	public static final int iTut7 = -1;
	public static final int iCavity = 1;
	public static final int iCavitySPH = 2;
	public static final int iDose = 3;
	public static final int iEdk = 4;
	public static final int iFlur = 5;
	public static final int iSpr = 6;
	public static final int iGe = 7;
	public static int iraycorr = 0;// default;1=>cav
	public static int isemfp = 0;// default;//for $SELECT-ELECTRON-MFP;1=>cav
	public static int iurd = 0;// default;//for user_range_discard;1=>cav
	public static int ispmfp = 0;// default;//for $SELECT-PHOTON-MFP;1=>cav
	public static int iGeom = 0;// default no geometry!
	// -----------------------------------------------------------
	// what to print 0=nothing,1=results only(eg from watch), 2=summary, 3=all
	// shit
	public static int iprint = 0;
	public static String seqStr = "";
	// ----------------------------------------------------------
	// sig_ismonotone(0,medium),sig_ismonotone(0:1,$MXMED);
	public static boolean[][] sig_ismonotone = new boolean[2][$MXMED];

	// ******************************************************************************************
	// set the number of medium->default 10:must be
	// >=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	/**
	 * Set the maximum number of mediums. Default is 10. Must be at least 1.
	 * @param n n
	 */
	public static void setMXMED(int n) {
		$MXMED = n;
	}

	// set the number of geometric regions->default 2000:must be
	// >=1!!!!!!!!!!!!!!!!!!!!!!!!!!
	/**
	 * Set the maximum number of geometric regions. Default is 2000. Must be at least 1.
	 * @param n n
	 */
	public static void setMXREG(int n) {
		$MXREG = n;
	}

	// set the number of particles in stack at once->default 40:must be
	// >=1!!!!!!!!!!!!!!!!!!!!!
	/**
	 * Set the maximum number of particles in stack at once. Default is 40. Must be at least 1.
	 * @param n n
	 */
	public static void setMXSTACK(int n) {
		$MXSTACK = n;
	}

	/**
	 * Set up the electron maximum range (for range rejection). Default is 100.
	 * @param n n
	 */
	public static void setMXRANGE(int n) {
		$MXRANGE = n;
	}

	// ---------------------------------
	// TO BE CONSISTENT WITH OTHER IMPLEMENTATION OF EQ, THE PASS OF EgsQuestion
	// MUST BE
	// MADE BY STATEMENT:"EGS4.eq=this;" FROM APPZ AND THEN SOMETIME:
	// EGS4.Hatch!!!
	// public EGS4(EgsQuestion eq)
	// {
	// this.eq=eq;
	// HATCH();//only purpose of this class
	// }

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * "Hatch" the media information.
	 * Setup which the user is expected to do before calling HATCH is: <br>
	 * 1. SET 'NMED' TO THE NUMBER OF MEDIA TO BE USED. <br>
	 * 2. SET THE ARRAY 'MEDIA', WHICH CONTAINS THE NAMES OF THE 
	 *  MEDIA THAT ARE DESIRED.  THE CHARACTER FORMAT IS A1, SO THAT MEDIA(IB,IM) CONTAINS THE 
	 *  IB'TH BYTE OF THE NAME OF THE IM'TH MEDIUM IN A1 FORMAT.<br>
	 * 3. SET 'DUNIT', THE DISTANCE UNIT TO BE USED. 
	 * DUNIT.GT.0 MEANS VALUE OF DUNIT IS LENGTH OF DISTANCE UNIT. <br>
	 * 4. FILL THE ARRAY 'MED' WITH THE MEDIUM INDICES FOR THE REGIONS. <br>
	 * 5. FILL ARRAYS 'ECUT' AND 'PCUT' WITH THE ELECTRON AND PHOTON 
	 * CUT-OFF ENERGIES FOR EACH REGION RESPECTIVELY. SETUP WILL 
	 * RAISE THESE IF NECESSARY TO MAKE THEM AT LEAST AS LARGE AS 
	 * THE REGION'S MEDIUM'S AE AND AP RESPECTIVELY. <br>
	 * 6. FILL 'MED' ARRAY.  MED(IR) IS THE MEDIUM INDEX FOR REGION 
	 * IR.  A ZERO MEDIUM INDEX MEANS THE REGION IS IN A VACUUM. <br>
	 * 7. FILL THE ARRAY 'IRAYLR' WITH 1 FOR EACH REGION IN WHICH 
	 * RAYLEIGH (COHERENT) SCATTERING IS TO BE INCLUDED.
	 */
	public static void HATCH()// private void HATCH()
	{
		// "******************************************************************"
		// "   Setup which the user is expected to do before calling HATCH is:"
		// "     1. SET 'NMED' TO THE NUMBER OF MEDIA TO BE USED."
		// "     2. SET THE ARRAY 'MEDIA', WHICH CONTAINS THE NAMES OF THE"
		// "        MEDIA THAT ARE DESIRED.  THE CHARACTER FORMAT IS A1, SO"
		// "        THAT MEDIA(IB,IM) CONTAINS THE IB'TH BYTE OF THE NAME OF"
		// "        THE IM'TH MEDIUM IN A1 FORMAT."
		// "     3. SET 'DUNIT', THE DISTANCE UNIT TO BE USED."
		// "        DUNIT.GT.0 MEANS VALUE OF DUNIT IS LENGTH OF DISTANCE UNIT"
		// "        CENTIMETERS.  DUNIT.LT.0 MEANS USE THE RADIATION LENGTH OF"
		// "        THE ABS(DUNIT)'TH MEDIUM FOR THE DISTANCE UNIT."
		// "     4. FILL THE ARRAY 'MED' WITH THE MEDIUM INDICES FOR THE"
		// "        REGIONS."
		// "     5. FILL ARRAYS 'ECUT' AND 'PCUT' WITH THE ELECTRON AND PHOTON"
		// "        CUT-OFF ENERGIES FOR EACH REGION RESPECTIVELY.  SETUP WILL"
		// "        RAISE THESE IF NECESSARY TO MAKE THEM AT LEAST AS LARGE AS"
		// "        THE REGION'S MEDIUM'S AE AND AP RESPECTIVELY."
		// "     6. FILL 'MED' ARRAY.  MED(IR) IS THE MEDIUM INDEX FOR REGION"
		// "        IR.  A ZERO MEDIUM INDEX MEANS THE REGION IS IN A VACUUM."
		// "     7. FILL THE ARRAY 'IRAYLR' WITH 1 FOR EACH REGION IN WHICH"
		// "        RAYLEIGH (COHERENT) SCATTERING IS TO BE INCLUDED."
		// "
		// "   KMPO = 8 and KMPI = 12 are set in BLOCK DATA"
		// "   The echo to unit 8 has been removed since it is sent"
		// "                                          to /dev/null anyway!"
		// "           To put it back search $UOUTPUT  and $ECHO and uncomment"
		// "******************************************************************"
		// ; Copyright NRC;
		// $DEFINE-LOCAL-VARIABLES-HATCH;
		// ------------------------------local
		// var------------------------------------
		// String[] MBUF = new String[72];
		// String[] MDLABL = new String[8];
		double ACD = 0.0;// , "used to test goodness of sine-table look-up"
		double ADEV = 0.0;// "absolute deviation in sine-table look-up"
		double ASD = 0.0;// "used to test goodness of sine-table look-up"
		double COST = 0.0;// "cos(theta) from instrinsic library function"
		double CTHET = 0.0;// "use to calculate cos(theta) according to look-up tables"
		double DEL = 0.0;// "leat squares delta for sine-table look-up"
		double DFACT = 0.0;// "converts rl to dunits"
		double DFACTI = 0.0;// "converts rl**-1 to dunits**-1"
		double DUNITO = 0.0;// "units scaling varable"
		double DUNITR = 0.0;// "saved value of dunit"
		double FNSSS = 0.0;// "real form of integer nsinss"
		double P = 0.0;// "counter used in the pwr2i(i) = 1/2**(i - 1) construction"
		double PZNORM = 0.0;// "used in $INITIALIZE-BREMS-ANGLE"
		double RDEV = 0.0;// "relative deviation in sine-table look-up"
		double S2C2 = 0.0;// "sinthe**2 + costhe**2, used to test look-up table"
		double S2C2MN = 0.0;// "min(s2c2)"
		double S2C2MX = 0.0;// "max(s2c2)"
		double SINT = 0.0;// "sin(theta) from instrinsic library function"
		double SX = 0.0;// "sum of angles for least squared analysis of look-up table errors"
		double SXX = 0.0;// "sum**2 of angles for least square analysis of look-up table errors"
		double SXY = 0.0;// "sum of angle*sin(angle) for least squared analysis of look-up"
							// "table errors"
		double SY = 0.0;// "sum of sin(angle) for least squared analysis of look-up table "
						// "errors"
		double WID = 0.0;// "width of sine-table mesh sub-interval (sine-table algorithm)"
		double XS = 0.0;// "angle value in a sub-sub-interval (sine-table algorithm)"
		double XS0 = 0.0;// "lower limit of a sub-sub-interval (sine-table algorithm)"
		double XS1 = 0.0;// "upwer limit of a sub-sub-interval (sine-table algorithm)"
		double XSI = 0.0;// "beginning angle of a sun-interval (sine-table algorithm)"
		double WSS = 0.0;// "width of a sub-sub-interval (sine-table algorithm)"
		double YS = 0.0;// "sin(angle) for least squared analysis of look-up table errors"
		double[] ZEROS = new double[3];// "zeros of sine, 0,pi,twopi"

		int I = 0;// , "generic do-loop variable"
		// int I1ST=0;//"flag = 0 on first pass"
		// int IB = 0;// "do-loop variable used for reading the medium type"
		int ID = 0;// "integer value of -dunit, when dunit is negative"
		int IE = 0;// "do-loop variable for reading over elements in a compound/mixture"
		// int IL = 0;// "do-loop variable used for reading the medium type"
		int IM = 0;// "do-loop variable looping over nmed, number of media"
		// int IRAYL=0;//"Rayleigh switch read in from PEGS"
		int IRN = 0;// "do-loop variable over random set of sine-table look-ups"
		int ISTEST = 0;// "flag that switches on test of sine function fit"
		int ISUB = 0;// "do-loop variable over sub-intervals of the sine look-up table"
		int ISS = 0;// "do-loop variable over sub-sub-intervals of the sine look-up table"
		int IZ = 0;// "used to locate an exact zero of a sun-interval mesh point in the"
					// "sine-table look-up"
		int IZZ = 0;// "do-loop variable over the exact zeros of the sine-table look-up"
		int J = 0;// "do-loop variable looping over nmed, number of media"
		int JR = 0;// "do-loop variable looping over number of regions"
		int LCTHET = 0;// "$SET INTERVAL index for cos(theta) from look-up table"
		// int LMDL = 0;// "character width of medium header ' MEDIUM='"
		// int LMDN = 0;// "character width of medium description"
		int LTHETA = 0;// "$SET INTERVAL index for sin(theta) from look-up table"
		int MD = 0;// "temporary storage for the medium number"
		int MXSINC = 0;// "number of intervals approximating the sine function"
		// int NCMFP = 0;//
		// "array size input from PEGS. Currently 0, but probably intended"
		// "to be cumulative electron mean free path. Presently unused."
		// int NEKE=0;//"array size input from PEGS."
		// "Number of electron mapped energy intervals."
		// int NGE = 0;// "array size input from PEGS."
		// "Number of photon mapped energy intervals."
		// int NGRIM = 0;// "Rayleigh cross section array size."
		int NISUB = 0;// "mxsinc - 2. Size of array with endpoints removed."
		// int NLEKE = 0;//
		// "array size input from PEGS. Currently 0, but probably intended"
		// "to be number of electron energy intervals below threshold."
		// "Presently unused."
		int NM = 0;// "number of media found in the "
		// int NRANGE = 0;//
		// "array size input from PEGS. Currently 0, but probably intended"
		// "to be number of intervals in an array giving the electron range."
		// "Presently unused."
		int NRNA = 0;// "number of random angles testing sine function fit"
		// int NSEKE = 0;//
		// "array size input from PEGS. Currently 0, but probably intended"
		// "to be number of electron small energy intervals. Presently unused."
		// int NSGE = 0;//
		// "array size input from PEGS. Currently 0, but probably intended"
		// "to be number of gamma small energy intervals. Presently unused."
		int NSINSS = 0;// "number of sub-intervals for each sine function interval"
		int[] LOK = new int[$MXMED]; // "flag indicating that medium has been found in the PEGS "
		// "datafile"

		// --------------------------------------------------------------------------
		STOPPROGRAM = false;
		ISTEST = 0;
		// ISTEST=1;
		// LMDL = 8;
		// LMDN = 24;
		DUNITO = 1.;
		// I1ST=1;//default
		NSINSS = 37;
		MXSINC = $MXSINC;
		NRNA = 1000;
		if (I1ST != 0) {
			I1ST = 0;
			// "RESET FIRST TIME FLAG"
			// "   DO FIRST TIME INITIALIZATION"
			EGS4Macro.HATCH_USER_INPUT_INIT();
			// "MACRO TO ALLOW USER TO INITIALIZE IN SUBROUTINE HATCH"
			// "   SET UP ENERGY PRECISION VARIABLES"
			PRM = RM; // "PRECISE REST MASS"
			PRMT2 = 2.0 * PRM;// "TWICE THE PRECISION REST MASS"
			PZERO = 0.0; // "PRECISE ZERO"
			// "   NOW CONSTRUCT PIECEWISE LINEAR FIT TO SINE FUNCTION OVER THE"
			// "   INTERVAL (0,5*PI/2).  DIVIDE THIS INTERVAL INTO MXSINC SUB-"
			// "   INTERVALS.  EACH OF THESE SUBINTERVALS IS THEN SUBDIVIDED INTO"
			// "   NSINSS SUB-SUB-INTERVALS.  THE ANGLES AT THE BOUNDARIES OF"
			// "   THESE SUB-SUB-INTERVALS AND THEIR SINES ARE USED TO COMPUTE"
			// "   LEAST SQUARES COEFFICIENTS FOR THE SUBINTERVAL.  AN EXTRA"
			// "   SUBINTERVAL ON EACH SIDE OF THE INTERVAL (0,5*PI/2) IS INCLUDED"
			// "   FOR GOOD MEASURE."
			NISUB = MXSINC - 2;
			FNSSS = (new Integer(NSINSS)).doubleValue();
			Integer in = new Integer(NISUB);
			// WID=PI5D2/FLOAT(NISUB);
			WID = PI5D2 / in.doubleValue();
			WSS = WID / (FNSSS - 1.0);
			// ZEROS(1)=0.;ZEROS(2)=PI; ZEROS(3)=TWOPI;
			ZEROS[0] = 0.;
			ZEROS[1] = PI;
			ZEROS[2] = TWOPI;
			for (ISUB = 1; ISUB <= MXSINC; ISUB++) {// "LOOP OVER SUBINTERVALS"
				SX = 0.;
				SY = 0.;
				SXX = 0.;
				SXY = 0.;// "ZERO SUMS"
				Integer intg = new Integer(ISUB - 2);
				// XS0=WID*FLOAT(ISUB-2);XS1=XS0+WID;"LOWER & UPPER LIMITS"
				XS0 = WID * intg.doubleValue();
				XS1 = XS0 + WID;// "LOWER & UPPER LIMITS"
				// "   NOW CHECK TO SEE IF ANY ZEROS ARE IN THE INTERVAL"
				IZ = 0;
				// DO IZZ=1,3 [
				for (IZZ = 1; IZZ <= 3; IZZ++) {
					if ((XS0 <= ZEROS[IZZ - 1]) && (ZEROS[IZZ - 1] <= XS1)) {
						IZ = IZZ;
						break;
					}
				}// "END OF LOOP OVER ZEROS"
				if (IZ == 0) {
					XSI = XS0;
				} else {
					XSI = ZEROS[IZ - 1];
				}
				for (ISS = 1; ISS <= NSINSS; ISS++) {// "LOOP OVER SUB-SUBINTERVALS"
					Integer i1 = new Integer(ISUB - 2);
					Integer i2 = new Integer(ISS - 1);
					XS = WID * i1.doubleValue() + WSS * i2.doubleValue() - XSI; // "ANGLE VALUE"
					YS = Math.sin(XS + XSI); // "SINE OF ANGLE"
					SX = SX + XS; // "ACCUMULATE SUMS"
					SY = SY + YS;
					SXX = SXX + XS * XS;
					SXY = SXY + XS * YS;
				}// "END SUB-SUBINTERVAL LOOP"
					// "   NOW COMPUTE LEAST SQUARES COEFFICIENTS"
				if (IZ != 0) {// "FORCE FIT THROUGH SINES' ZEROS,"
								// "             FOR SMALL REL.ERR.&GOOD"
								// "   VALUES OF SINTHE/THETA NEAR ZERO"
					SIN1[ISUB - 1] = SXY / SXX;// SIN1(ISUB)=SXY/SXX;
					SIN0[ISUB - 1] = -SIN1[ISUB - 1] * XSI;// SIN0(ISUB)=-SIN1(ISUB)*XSI;
				} else {// "DO FULL LEAST SQUARES"
					DEL = FNSSS * SXX - SX * SX;
					SIN1[ISUB - 1] = (FNSSS * SXY - SY * SX) / DEL;// SIN1(ISUB)=(FNSSS*SXY-SY*SX)/DEL;
					SIN0[ISUB - 1] = (SY * SXX - SX * SXY) / DEL
							- SIN1[ISUB - 1] * XSI;// SIN0(ISUB)=(SY*SXX-SX*SXY)/DEL
													// - SIN1(ISUB)*XSI ;
				}
			} // "END SUB-INTERVAL LOOP"

			SINC0 = 2.0;// "SET COEFFICIENTS WHICH DETERMINE INTERVAL"
			SINC1 = 1.0 / WID;
			// "   NOW TEST FIT, IF REQUESTED"
			if (ISTEST != 0) {
				// "   FIRST TEST AT POINTS PREVIOUSLY COMPUTED, EXCLUDING"
				// "   END SUBINTERVALS"
				ADEV = 0.;
				RDEV = 0.;
				S2C2MN = 10.;
				S2C2MX = 0.;
				for (ISUB = 1; ISUB <= NISUB; ISUB++) {
					for (ISS = 1; ISS <= NSINSS; ISS++) {
						Integer i1 = new Integer(ISUB - 1);
						Integer i2 = new Integer(ISS - 1);
						THETA = WID * i1.doubleValue() + WSS * i2.doubleValue();
						CTHET = PI5D2 - THETA;
						// $SET INTERVAL THETA,SINC;$SET INTERVAL CTHET,SINC;
						// SPECIFY SNAME AS
						// ['sinc'|'blc'|'rthr'|'rthri'|'SINC'|'BLC'|'RTHR'|'RTHRI'];
						// SPECIFY SNAME1 AS ['sin'|'SIN'];

						// REPLACE {$SETINTERVAL#,#;} WITH {
						// [IF] '{P2}'=SNAME [L{P1}={P2}1*{P1}+{P2}0;]
						// [ELSE] [L{P1}={P2}1(MEDIUM)*{P1}+{P2}0(MEDIUM);]}
						// "TWO ARGUMENT SET INTERVAL CALL (ABOVE) CURRENTLY MEANS"
						// "INTERVAL INDEX IS LINEAR FUNCTION OF THE LINEAR VARIABLE,"
						// "WITH NO MAPPING.  {P1} IS THE LINEAR VARIABLE AND {P2} IS"
						// "THE NAME OF THE INTERVAL(WHICH IS USED TO CONSTRUCT THE"
						// "COEFFICIENTS USED IN COMPUTING THE INTERVAL INDEX)."
						// "BUT, IF {P2} IS SINC OR BLC OR RTHR OR RTHRI, IT DOES"
						// "NOT DEPEND ON MEDIUM"
						// SPECIFY SNAME AS
						// ['sinc'|'blc'|'rthr'|'rthri'|'SINC'|'BLC'|'RTHR'|'RTHRI'];
						// SPECIFY SNAME1 AS ['sin'|'SIN'];
						Double dth = new Double(SINC1 * THETA + SINC0);// System.out.println(dth);
						LTHETA = dth.intValue();
						Double dcth = new Double(SINC1 * CTHET + SINC0);
						LCTHET = dcth.intValue();
						// REPLACE {$EVALUATE#USING#(#);} WITH {
						// [IF] '{P2}'=SNAME1
						// [{P1}={P2}1(L{P3})*{P3}+{P2}0(L{P3});] [ELSE]
						// [{P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM);]}
						// "{P1} IS VARIABLE TO BE ASSIGNED VALUE."
						// "{P2} IS THE FUNCTION BEING APPROXIMATED."
						// "{P3} IS THE ARGUMENT OF THE FUNCTION. IN THE CURRENT"
						// "PWLF METHOD, THE ARGUMENT DETERMINES AN INTERVAL USING THE"
						// "$SET INTERVAL MACROS.   WITH IN THIS INTERVAL THE"
						// "FUNCTION IS APPROXIMATED AS A LINEAR FUNCTION OF"
						// "THE ARGUMENT. BUT"
						// "IF {P2}=SIN IT DOES NOT DEPEND ON MEDIUM"
						// REPLACE {$EVALUATE#USING SIN(#);} WITH
						// {{P1}=sin({P2});}
						// $EVALUATE SINTHE USING SIN(THETA);
						SINTHE = SIN1[LTHETA - 1] * THETA + SIN0[LTHETA - 1];
						// $EVALUATE COSTHE USING SIN(CTHET);
						COSTHE = SIN1[LCTHET - 1] * CTHET + SIN0[LCTHET - 1];

						SINT = Math.sin(THETA);
						COST = Math.cos(THETA);
						ASD = Math.abs(SINTHE - SINT);
						ACD = Math.abs(COSTHE - COST);
						ADEV = max(ADEV, ASD, ACD);
						if (SINT != 0.0)
							RDEV = Math.max(RDEV, ASD / Math.abs(SINT));
						if (COST != 0.0)
							RDEV = Math.max(RDEV, ACD / Math.abs(COST));
						S2C2 = SINTHE * SINTHE + COSTHE * COSTHE;
						S2C2MN = Math.min(S2C2MN, S2C2);
						S2C2MX = Math.max(S2C2MX, S2C2);
						if (ISUB < 11) {
							seqStr = "th=" + THETA + " isin=" + SINTHE
									+ " rsin=" + SINT + " icos=" + COSTHE
									+ " rcos=" + COST;// +" \n";
							if (iprint > 1)
								eq.printSequence(seqStr);
						}
					}
				} // "END OF FIXED INTERVAL TEST-OUTPUT RESULTS"
					// OUTPUT MXSINC,NSINSS;(' SINE TESTS,MXSINC,NSINSS=',2I5);
					// OUTPUT ADEV,RDEV,S2C2MN,S2C2MX;
					// (' ADEV,RDEV,S2C2(MN,MX) =',1PE16.8,3E16.8);
				seqStr = "SINE TESTS, MXSINC:" + MXSINC + " NSINSS= " + NSINSS;// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "AbsoluteDEV=" + ADEV + " RDEV= " + RDEV + " S2C2MN= "
						+ S2C2MN + " S2C2MX= " + S2C2MX;// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "**************************RANDOM TEST*******************************";
				// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
				// "   NOW DO RANDOM TEST"
				ADEV = 0.;
				RDEV = 0.;
				S2C2MN = 10.;
				S2C2MX = 0.;
				for (IRN = 1; IRN <= NRNA; IRN++) {
					// $RANDOMSET THETA;THETA=THETA*PI5D2;
					THETA = random01();
					THETA = THETA * PI5D2;
					CTHET = PI5D2 - THETA;
					// $SET INTERVAL THETA,SINC;$SET INTERVAL CTHET,SINC;
					Double dth = new Double(SINC1 * THETA + SINC0);// System.out.println(dth);
					LTHETA = dth.intValue();
					Double dcth = new Double(SINC1 * CTHET + SINC0);
					LCTHET = dcth.intValue();
					// $EVALUATE SINTHE USING SIN(THETA);
					SINTHE = SIN1[LTHETA - 1] * THETA + SIN0[LTHETA - 1];
					// $EVALUATE COSTHE USING SIN(CTHET);
					COSTHE = SIN1[LCTHET - 1] * CTHET + SIN0[LCTHET - 1];

					SINT = Math.sin(THETA);
					COST = Math.cos(THETA);
					ASD = Math.abs(SINTHE - SINT);
					ACD = Math.abs(COSTHE - COST);
					ADEV = max(ADEV, ASD, ACD);
					if (SINT != 0.0)
						RDEV = Math.max(RDEV, ASD / Math.abs(SINT));
					if (COST != 0.0)
						RDEV = Math.max(RDEV, ACD / Math.abs(COST));
					S2C2 = SINTHE * SINTHE + COSTHE * COSTHE;
					S2C2MN = Math.min(S2C2MN, S2C2);
					S2C2MX = Math.max(S2C2MX, S2C2);

					seqStr = "thR=" + THETA + " isinR=" + SINTHE + " rsinR="
							+ SINT + " icosR=" + COSTHE + " rcosR=" + COST;
					// +" \n";
					if (iprint > 1)
						eq.printSequence(seqStr);
				} // "END RANDOM ANGLE LOOP"
				seqStr = "TEST AT " + NRNA + " RANDOM ANGLES IN (0,5*PI/2)";
				// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "AbsoluteDEV=" + ADEV + " RDEV= " + RDEV + " S2C2MN= "
						+ S2C2MN + " S2C2MX= " + S2C2MX;
				// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
			}// "END OF SINE TABLE TEST"
				// "   NOW FILL IN POWER OF TWO TABLE.  PWR2I(I)=1/2**(I-1)"
			P = 1.;
			// DO I=1,$MXPWR2I [PWR2I(I)=P; P=P/2.;]
			for (I = 0; I < $MXPWR2I; I++) {
				PWR2I[I] = P;
				P = P / 2.;
			}

		}// "END OF FIRST TIME INITIALIZATION"
			// "FILL IRAYLM ARRAY BASED ON IRAYLR INPUTS"
			// $need_rayleigh_data;
			// REPLACE {$need_rayleigh_data;} WITH {;
			// for(J=1,NMED [
		for (J = 1; J <= NMED; J++) {
			LOOP_OVER_REGIONS:
			// DO I=1,$MXREG [
			for (I = 0; I < $MXREG; I++) {
				if ((IRAYLR[I] == 1) && (MED[I] == J))// /because default
														// 0==VACUUM!!
				{
					// "REGION I = MEDIUM J AND WE WANT RAYLEIGH SCATTERING, SO"
					// "SET FLAG TO PICK UP DATA FOR MEDIUM J AND TRY NEXT MEDIUM."
					IRAYLM[J - 1] = 1;// iraylm[J-1]=1;
					break LOOP_OVER_REGIONS;
				}
			}// "END OF REGION-LOOP"
		}// "END OF MEDIA-LOOP"
			// };
			// "   NOW SEARCH FILE FOR DATA FOR REQUESTED MATERIALS"
		NM = 0; // "NUMBER OF MEDIA FOUND"
		for (IM = 1; IM <= NMED; IM++) {
			LOK[IM - 1] = 0;// "SET FLAG TELLING WHICH MEDIA ARE OK"
							// "NOW TELL USER IF RAYLEIGH OPTION HAS BEEN REQUESTED"
			if (IRAYLM[IM - 1] == 1) {
				seqStr = " RAYLEIGH OPTION REQUESTED FOR MEDIUM NUMBER " + IM;
				// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
			}
		}

		// :MEDIUM:
		// LOOP["MEDIUM SEARCH LOOP"
		// UNTIL NM.GE.NMED; "LOOP UNTIL WE HAVE ENOUGH.  END :MEDIUM: LOOP"
		while (NM < NMED) {
			NM = NM + 1;// "SET FOUND FLAG AND STEP MEDIUM COUNTER"
			IM = NM;// "'IM' IS THE INDEX OF THE MEDIUM READY TO BE READ"<->FORCED
					// ALWAYS EXISTS!!!
			// "   NOW READY TO READ IN DATA FOR THIS MEDIUM"
			seqStr = "Data for medium " + IM + " which is:" + MEDIA[IM - 1];
			// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
			readPegs4Dat(MEDIA[IM - 1], IM);
			TE[IM - 1] = AE[IM - 1] - RM;
			THMOLL[IM - 1] = TE[IM - 1] * 2. + RM;
			// NSGE = MSGE[IM - 1];
			// NGE = MGE[IM - 1];
			// NSEKE = MSEKE[IM - 1];
			NEKE = MEKE[IM - 1];
			// NLEKE = MLEKE[IM - 1];
			// NCMFP = MCMFP[IM - 1];
			// NRANGE = MRANGE[IM - 1];
			// NGRIM = NGR[IM - 1];

			seqStr = "****************************************************************************";
			// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
			seqStr = "RHO " + " = " + format(RHO[IM - 1], 8) + " NNE " + " = "
					+ NNE[IM - 1] + " IUNRST " + " = " + IUNRST[IM - 1]
					+ " EPSTFL " + " = " + EPSTFL[IM - 1] + " IAPRIM " + " = "
					+ IAPRIM[IM - 1];// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			for (IE = 0; IE < NNE[IM - 1]; IE++) {
				int ie = IE + 1;
				seqStr = "index:" + format(ie, 2) + " ASYM " + " = "
						+ format(ASYM[IM - 1][IE], 4) + " Z " + " = "
						+ format(ZELEM[IM - 1][IE], 8) + " WA " + " = "
						+ format(WA[IM - 1][IE], 8) + " PZ " + " = "
						+ format(PZ[IM - 1][IE], 8) + " RHOZ " + " = "
						+ format(RHOZ[IM - 1][IE], 8);// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
			}
			seqStr = "RLC " + " = " + format(RLC[IM - 1], 8) + " AE " + " = "
					+ format(AE[IM - 1], 8) + " AP " + " = "
					+ format(AP[IM - 1], 8) + " UE " + " = "
					+ format(UE[IM - 1], 8) + " UP " + " = "
					+ format(UP[IM - 1], 8);// +" TE "+" = "+format(TE[IM-1],8)+
			// " THMOLL "+" = "+format(THMOLL[IM-1],8);//+" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
			seqStr = "TE " + " = " + format(TE[IM - 1], 8) + " THMOLL " + " = "
					+ format(THMOLL[IM - 1], 8);// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			// "   THAT'S ALL FOR THIS MEDIUM"
		}// while(NM<NMED)

		// "   WE NOW HAVE DATA FOR ALL MEDIA REQUESTED.  NOW DO DISTANCE UNIT"
		// "   CHANGE.  DATA FROM PEGS IS IN UNITS OF RADIATION LENGTHS."
		// "   EGS IS RUN IN UNITS OF 'DUNIT' CENTIMETERS, IF DUNIT.GT.0 "
		// "   OR IN UNITS OF RLC(-DUNIT) CENTIMETERS IF DUNIT.LT.0."
		// "   THAT IS, A NEGATIVE DUNIT MEANS UNIT IS TO BE THE RADIATION"
		// "   LENGTH OF THE MEDIUM WHOSE INDEX IS -DUNIT"
		DUNITR = DUNIT; // "SAVE REQUESTED"
		Double cbl = new Double(-DUNIT);
		int duni = cbl.intValue();
		if (DUNIT < 0.0) {
			ID = Math.max(1, Math.min($MXMED, duni));
			DUNIT = RLC[ID - 1];
		}
		if (DUNIT != 1.0) {
			seqStr = " DUNIT REQUESTED&USED ARE:" + DUNITR + " & " + DUNIT
					+ " (cm.) ";// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
		}

		seqStr = " DUNIT REQUESTED ARE:" + DUNIT;// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		for (IM = 1; IM <= NMED; IM++) {
			DFACT = RLC[IM - 1] / DUNIT; // "CONVERTS RL TO DUNITS"
			DFACTI = 1.0 / DFACT; // "CONVERT RL**-1 TO DUNITS**-1"
			// REPLACE {$SCALE# BY #;} WITH {{P1}={P1}*{P2};}
			for (I = 0; I < MEKE[IM - 1]; I++) {
				// $SCALE $LGN(ESIG,PSIG,EDEDX,PDEDX(I,IM)/0,1/) BY DFACTI;
				// 1/rlc=>cm!!!!
				ESIG0[I][IM - 1] = ESIG0[I][IM - 1] * DFACTI;
				ESIG1[I][IM - 1] = ESIG1[I][IM - 1] * DFACTI;
				PSIG0[I][IM - 1] = PSIG0[I][IM - 1] * DFACTI;
				PSIG1[I][IM - 1] = PSIG1[I][IM - 1] * DFACTI;
				EDEDX0[I][IM - 1] = EDEDX0[I][IM - 1] * DFACTI;
				EDEDX1[I][IM - 1] = EDEDX1[I][IM - 1] * DFACTI;
				PDEDX0[I][IM - 1] = PDEDX0[I][IM - 1] * DFACTI;
				PDEDX1[I][IM - 1] = PDEDX1[I][IM - 1] * DFACTI;
				// $SCALE $LGN(TMXS(I,IM)/0,1/) BY DFACT;
				TMXS0[I][IM - 1] = TMXS0[I][IM - 1] * DFACT;
				TMXS1[I][IM - 1] = TMXS1[I][IM - 1] * DFACT;

			}
			// $SCALE TEFF0(IM) BY DFACT;
			TEFF0[IM - 1] = TEFF0[IM - 1] * DFACT;
			// $SCALE BLCC(IM) BY DFACTI;
			BLCC[IM - 1] = BLCC[IM - 1] * DFACTI;
			// $SCALE XCC(IM) BY SQRT(DFACTI);
			XCC[IM - 1] = XCC[IM - 1] * Math.sqrt(DFACTI);
			RLDU[IM - 1] = RLC[IM - 1] / DUNIT;
			for (I = 0; I < MGE[IM - 1]; I++) {
				// $SCALE $LGN(GMFP(I,IM)/0,1/) BY DFACT;
				GMFP0[I][IM - 1] = GMFP0[I][IM - 1] * DFACT;
				GMFP1[I][IM - 1] = GMFP1[I][IM - 1] * DFACT;
			}

		}// "END IM DO"
			// "   SCALE VACDST.  UNDO PREVIOUS SCALE, THEN DO NEW."
		VACDST = VACDST * DUNITO / DUNIT;
		DUNITO = DUNIT; // "SAVE OLD DUNIT"
		// "   NOW MAKE SURE ECUT AND PCUT ARE NOT LOWER THAN ANY AE OR AP"
		// "   ALSO SET DEFAULT DENSITIES"
		// $adjust_rhor_ecut_pcut;
		// REPLACE {$adjust_rhor_ecut_pcut;} WITH {;
		// DO JR=1,$MXREG [
		for (JR = 0; JR < $MXREG; JR++) {
			MD = MED[JR];// medium 0=vacuum else medium 1 or 2 etc
			if ((MD >= 1) && (MD <= NMED)) {// "IT IS LEGAL NON-VACUUM MEDIUM."
				ECUT[JR] = Math.max(ECUT[JR], AE[MD - 1]);
				PCUT[JR] = Math.max(PCUT[JR], AP[MD - 1]);
				// "   USE STANDARD DENSITY FOR REGIONS NOT SPECIALLY SET UP"
				if (RHOR[JR] == 0.0) {
					RHOR[JR] = RHO[MD - 1];
				}
			}
		}
		// };
		// "BREMSSTRAHLUNG ANGULAR DISTRIBUTION INITIALIZATION - DEFAULT IS NULL"
		// "NEXT LINE ADDED AFB 88/05/31"
		// $INITIALIZE-BREMS-ANGLE;
		// "Macro to initialize data for bremsstrahlung production               "
		// "The quantity ZBRANG is ( (1/111)*Zeff**(1/3) )**2                    "
		// "where Zeff is defined in equation (7) OF PIRS0203                    "
		// "This macro goes in SUBROUTINE HATCH                                  "
		// REPLACE {$INITIALIZE-BREMS-ANGLE;} WITH {
		if (ibrdst == 1) {
			// DO IM=1,NMED[
			for (IM = 1; IM <= NMED; IM++) {
				ZBRANG[IM - 1] = 0.0;
				PZNORM = 0.0;
				// DO IE=1,NNE(IM)[
				for (IE = 0; IE < NNE[IM - 1]; IE++) {
					ZBRANG[IM - 1] = ZBRANG[IM - 1] + PZ[IM - 1][IE]
							* ZELEM[IM - 1][IE] * (ZELEM[IM - 1][IE] + 1.0);
					PZNORM = PZNORM + PZ[IM - 1][IE];
				}
				ZBRANG[IM - 1] = (8.116224E-05) * Math.pow(ZBRANG[IM - 1]
						/ PZNORM, 1. / 3.);
				// ok,Z^2->Z(Z+1.0)
				// ZBRANG[IM-1]=(8.116224E-05)*Math.pow(ZBRANG[IM-1]/PZNORM,2./3.);
				LZBRANG[IM - 1] = -Math.log(ZBRANG[IM - 1]);
			}
		}
		// }
		// "PAIR ANGULAR DISTRIBUTION INITIALIZATION - DEFAULT IS NULL"
		// "NEXT LINE ADDED AFB 91/05/29"
		// $INITIALIZE-PAIR-ANGLE;
		// "Macro to initialize data for PAIR PRODUCTION                         "
		// "THE QUANTITY ZBRANG IS ( (1/111)*Zeff**(1/3) )**2                    "
		// "WHERE Zeff IS DEFINED IN EQUATION (7) OF PIRS0287                    "
		// "THIS MACRO GOES IN SUBROUTINE HATCH                                  "
		// "THIS MACRO IS IDENTICAL TO THE $INITIALIZE-BREMS-ANGLE DEFINED ABOVE "
		// "                                                                     "
		// REPLACE {$INITIALIZE-PAIR-ANGLE;} WITH {
		if (iprdst > 0) {
			// DO IM=1,NMED[
			for (IM = 1; IM <= NMED; IM++) {
				ZBRANG[IM - 1] = 0.0;
				PZNORM = 0.0;
				// DO IE=1,NNE(IM)[
				for (IE = 0; IE < NNE[IM - 1]; IE++) {
					ZBRANG[IM - 1] = ZBRANG[IM - 1] + PZ[IM - 1][IE]
							* ZELEM[IM - 1][IE] * (ZELEM[IM - 1][IE] + 1.0);
					PZNORM = PZNORM + PZ[IM - 1][IE];
				}
				ZBRANG[IM - 1] = (8.116224E-05) * Math.pow(ZBRANG[IM - 1]
						/ PZNORM, 1. / 3.);
				// ZBRANG[IM-1]=(8.116224E-05)*Math.pow(ZBRANG[IM-1]/PZNORM,2./3.);
				// ADDED------------------------->
				LZBRANG[IM - 1] = -Math.log(ZBRANG[IM - 1]);
			}
		}
		// }
		// @@@@@@@@@@@@@@@@@@@@@@
		// " See if the user has requested to use a different set of photon "
		// " cross section data "
		// IF( photon_xsections(1:1) ~= ' ' ) [
		// call egs_init_user_photon(photon_xsections,1);
		// ]
		if (photon_xsections.compareTo("") != 0) {
			egs_init_user_photon(photon_xsections, 1);
		}

		// =-====================
		// "Initialize new MS, step-sizes, etc, IK Oct 97"
		mscati();
		// "Initialize relaxations and photo-absorption data if requested
		EDGSET(1, 1);
		// "Initialize bound compton scattering, IK, Jan 99 if requested
		init_compton();
		// "Re-calculate dl1,... for the different technique"
		// "employed in BREMS. Note that the old EGS sampling"
		// "technique for BREMS had a bug that shows up only"
		// "if AP is not much smaller than electron kinetic energy"
		fix_brems();
		// "initializes the sampling tables and modifies the total"
		// "brems cross sections if the NIST brems data base is to"
		// "be used                                             "
		if (ibr_nist == 1) {
			init_nist_brems();
		}

		if (pair_nrc == 1) {
			init_nrc_pair();
		}
		// " Load and initialize EII data if needed. "
		eii_init();

		// " Load and initialize the triplet data if needed "
		init_triplet();

		// "   SETUP IS NOW COMPLETE"
		if (NMED == 1) {
			seqStr = " SUCCESSFULLY ''HATCHED'' FOR ONE MEDIUM.";// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
		} else {
			seqStr = " SUCCESSFULLY ''HATCHED'' FOR " + NMED + " MEDIA.";// +
																			// " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);
		}

		// "END OF SUBROUTINE HATCH"
	}

	// "**************************************************************************"
	// "Init EII. This subroutine is called from HATCH after all media are known, "
	// "threshold energies and interpolation data have been initialized.          "
	// "**************************************************************************"
	/**
	 * Init EII (electron impact ionization). This subroutine is called from HATCH 
	 * after all media are known, threshold energies and interpolation data have been initialized. 
	 */
	public static void eii_init() {
		// "**************************************************************************"
		int iZ = 0;
		int nsh = 0;
		// int nsh_tot = 0;
		int itmp = 0;
		int nbin = 0;
		int ii = 0;
		int iii = 0;
		int[] tmp_array = new int[$MXELEMENT];
		double emax = 0.0;
		double fmax = 0.0;
		double[] aux_array = new double[$N_EII_BINS];
		double sigo = 0.0;
		double loge = 0.0;
		double tau = 0.0;
		double beta2 = 0.0;
		double p2 = 0.0;
		double uwm = 0.0;
		double Wmax = 0.0;
		double ss_0 = 0.0;
		double ss_1 = 0.0;
		double sh_0 = 0.0;
		double sh_1 = 0.0;
		double aux = 0.0;
		// double av_e = 0.0;
		double con_med = 0.0;
		double dedx_old = 0.0;
		// double sigm_old = 0.0;
		double dedx = 0.0;
		double e = 0.0;
		double sig = 0.0;
		double sigm = 0.0;
		double wbrem = 0.0;
		double sum_a = 0.0;
		double sum_z = 0.0;
		double sum_pz = 0.0;
		double sum_wa = 0.0;
		double Ec = 0.0;
		double Ecc = 0.0;
		double sum_sh = 0.0;
		double sum_occn = 0.0;
		double U = 0.0;
		double sum_sigma = 0.0;
		double sum_dedx = 0.0;
		double sigma = 0.0;
		double sigma_old = 0.0;
		double wbrem_old = 0.0;
		double sig_j = 0.0;
		double de = 0.0;
		int lloge = 0;
		double cons = 0.153536; // " 2*Pi*Re^2*rm/u "
		int[] occn_numbers = { 2, 2, 2, 4 };
		int i = 0;
		int jj = 0;
		int jjj = 0;
		double av_E = 0.0;
		double sigma_max = 0.0;

		for (int j = 1; j <= $MXELEMENT; j++) {
			eii_nshells[j - 1] = 0;
		}
		for (int j = 1; j <= $MXMED; j++) {
			eii_nsh[j - 1] = 0;
		}
		if (eii_flag == 0) {
			return;
		}

		// find minimum of all threshold energies
		double e_eii_min = 1.e30;

		for (int imed = 1; imed <= NMED; imed++) {
			if (AE[imed - 1] - RM < e_eii_min)
				e_eii_min = AE[imed - 1] - RM;
			if (AP[imed - 1] < e_eii_min)
				e_eii_min = AP[imed - 1];
		}
		seqStr = "eii_init: minimum threshold energy found: " + e_eii_min;// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// determine elements that need to load EII data
		for (int imed = 1; imed <= NMED; imed++) {
			for (int iele = 1; iele <= NNE[imed - 1]; iele++) {
				Double dbl = new Double(ZELEM[imed - 1][iele - 1] + 0.5);
				iZ = dbl.intValue(); // int(zelem(imed,iele)+0.5);
				if (eii_nshells[iZ - 1] == 0)// if( eii_nshells[iZ] == 0 ),
												// ZELEM=min 1!!!
				{
					nsh = 0;
					for (int ish = 1; ish <= 4; ish++) {
						if (binding_energies[ish - 1][iZ - 1] > e_eii_min)
							nsh = nsh + 1;
					}
					eii_nshells[iZ - 1] = nsh;
				}
			}
		}

		// total number of shells that need to be loaded
		nsh = 0;
		for (iZ = 1; iZ <= $MXELEMENT; iZ++) {
			nsh = nsh + eii_nshells[iZ - 1];
		}
		if (nsh == 0) {
			seqStr = "*** EII requested but no shells with binding energies "
					+ "above the specified threshold found"
					+ "    => turning off EII";// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			eii_flag = 0;
		}
		if (nsh > $MAX_EII_SHELLS) {
			STOPPROGRAM = true;
			seqStr = "*** Number of shells with binding energies greater than "
					+ "the specified thresholds is " + nsh
					+ "    This is more than the allocated arrays can hold"
					+ "    Increase the macro $MAX_EII_SHELLS and retry";
			eq.printSequence(seqStr);

			return;
		}
		seqStr = "eii_init: number of shells to simulate EII: " + nsh;// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// nsh_tot = nsh;
		tmp_array[0] = 0;// tmp_array(1) = 0;
		// DO j=2,$MXELEMENT [ tmp_array(j) = tmp_array(j-1) + eii_nshells(j-1);
		// ]
		for (int j = 2; j <= $MXELEMENT; j++) {
			tmp_array[j - 1] = tmp_array[j - 2] + eii_nshells[j - 2];
		}
		// Get relaxation data if necessary so that the binding energies
		// are available
		itmp = iedgfl[0];// iedgfl(1);
		iedgfl[0] = 1;// iedgfl(1) = 1;
		EDGSET(1, 1);
		iedgfl[0] = itmp;// iedgfl(1) = itmp;
		// set EII active shells per medium and for each element
		for (int imed = 1; imed <= NMED; imed++) {
			nsh = 0;
			for (int iele = 1; iele <= NNE[imed - 1]; iele++) {
				Double dbl = new Double(ZELEM[imed - 1][iele - 1] + 0.5);
				iZ = dbl.intValue(); // int(zelem(imed,iele)+0.5);
				eii_no[imed - 1][iele - 1] = eii_nshells[iZ - 1];
				nsh = nsh + eii_nshells[iZ - 1];
				if (eii_nshells[iZ - 1] > 0) {
					eii_first[imed - 1][iele - 1] = tmp_array[iZ - 1] + 1;
				} else {
					eii_first[imed - 1][iele - 1] = 0;
				}
			}
			eii_nsh[imed - 1] = nsh;
		}

		// read EII data
		String filename = datas + file_sep + egsDatas + file_sep + eiifile
				+ defaultext;
		int iread = 0;
		int lnr = 0;// data number
		int lnrr = 0;// line number
		int indx = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals = "=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		int nskip = 0;
		boolean emaxb = false;
		boolean nbinb = false;
		int maini = 0;
		boolean mainb = false;

		seqStr = " eii_init: reading EII data ... ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// read(eii_unit,*,err=:eii-reading-error:,end=:eii-reading-error:)
		// nskip;
		FileInputStream in =null;
		try {
			//FileInputStream in = new FileInputStream(filename);
			in = new FileInputStream(filename);
			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if (lnr == 1)// read nskip
						{
							String s = desc.toString();
							nskip = stringToInt(s);// System.out.println(nskip);
						} else if (lnrr == nskip + 1)// =>have data for
														// lnrr>nskip!!
						{
							// read emax,nbin;
							if (!emaxb) {
								String s = desc.toString();
								emax = stringToDouble(s);
								// System.out.println("emax: "+emax);
								emaxb = true;
							} else if (!nbinb) {
								String s = desc.toString();
								nbin = stringToInt(s);
								// System.out.println("nbin: "+nbin);
								nbinb = true;
								if (nbin != $N_EII_BINS) {
									STOPPROGRAM = true;
									seqStr = "Inconsistent EII data file";
									// if(iprint>1)
									eq.printSequence(seqStr);

									return;
								}
								ii = 0;
							}
						}
						// else if(lnrr==nskip+2)//=>have data for lnrr>nskip!!
						else if (lnrr > nskip + 1)// =>have data for
													// lnrr>nskip!!
						{
							if (!mainb) {
								// -----------------
								// read iZ,nsh;
								String s = desc.toString();
								iZ = stringToInt(s);
								// System.out.println("iZ="+iZ);

								desc = new StringBuffer();// reset
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											appendB = false;// read data so
															// ready to
															// loop!!!!!
											String s1 = desc.toString();
											nsh = stringToInt(s1);
											// System.out.println("nsh: "+nsh);
											b = false;// exit the loop
										}

									}
								}

								if (nsh < eii_nshells[iZ - 1]) {
									STOPPROGRAM = true;
									seqStr = "EII data file has data for "
											+ nsh
											+ " shells for element "
											+ iZ
											+ ", but according"
											+ "to binding energies and thresholds "
											+ eii_nshells[iZ - 1]
											+ "shells are required. This is a fatal error.";
									// if(iprint>1)
									eq.printSequence(seqStr);

									return;
								}

								// ish=1;
								if ((nsh > 0))// &&(ish<=nsh))
								{
									int ish = 1;
									while (ish <= nsh) {
										desc = new StringBuffer();// reset
										// read fmax
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													appendB = false;// read data
																	// so ready
																	// to
																	// loop!!!!!
													String s1 = desc.toString();
													fmax = stringToDouble(s1);
													// System.out.println("fmax: "+fmax);
													b = false;// exit the loop
												}

											}
										}
										// read aux_array
										int kindx = 0;
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													appendB = false;// read data
																	// so ready
																	// to
																	// loop!!!!!
													String s1 = desc.toString();
													aux_array[kindx] = stringToDouble(s1);
													// System.out.println("i= "+ish+" k= "+kindx+" aux[k] "+aux_array[kindx]);
													kindx++;
													if (kindx == $N_EII_BINS) {
														b = false;// exit the
																	// loop
														if (ish <= eii_nshells[iZ - 1]) {
															ii = ii + 1;
															eii_z[ii - 1] = iZ;
															eii_sh[ii - 1] = ish;
															eii_a[ii - 1] = (new Integer(
																	nbin))
																	.doubleValue();
															eii_a[ii - 1] = eii_a[ii - 1]
																	/ Math.log(emax
																			/ binding_energies[ish - 1][iZ - 1]);
															eii_b[ii - 1] = 1.
																	- eii_a[ii - 1]
																	* Math.log(binding_energies[ish - 1][iZ - 1]);
															for (int k = 1; k <= nbin; k++) {
																if (k > 1) {
																	sigo = fmax
																			* aux_array[k - 2];
																} else {
																	sigo = 0.;
																}
																loge = (k - eii_b[ii - 1])
																		/ eii_a[ii - 1];
																iii = nbin
																		* (ii - 1)
																		+ k;// ok
																eii_xsection_a[iii - 1] = (fmax
																		* aux_array[k - 1] - sigo)
																		* eii_a[ii - 1];
																eii_xsection_b[iii - 1] = sigo
																		- eii_xsection_a[iii - 1]
																		* loge;
															}
														}

													}
												}
												desc = new StringBuffer();// reset
											}
										}
										// -----------
										ish++;
									}// while ish<=nsh
										// if( ii == nsh_tot ) [ break; ]

								}// nsh>0
									// -------------------
								maini++;
								if (maini == $MXELEMENT) {
									mainb = true;// exit
								}

							}
						}

					}// have data
					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			//in.close();
		}// try
		catch (Exception exc) {
			
		}
		finally {
	        if( null != in ) {
	            try 
	            {
	                in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
		}

		seqStr = " OK !";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// Adjust restricted stopping power and discrete inelastic cross
		// sections

		// Discrete interaction cross sections for all shells that come from
		// PEGS4
		// are calculated using the Moller cross section
		// => we must subtract the Moller cross section for the shells that will
		// have EII and then add the EII cross sections for these shells
		// The restricted stopping power that comes from PEGS4 is calculated
		// assuming
		// that all shells will be producing secondaries according to Moller
		// => we must add the energy lost in Moller events in EII shells to the
		// restricted stopping power and then subtract the average energy lost
		// in EII events.
		int medium = 0;
		for (int imed = 1; imed <= NMED; imed++) {
			// -------------------------mine
			medium = imed;
			// ------------------------------
			Ec = AE[imed - 1] - RM;
			Ecc = Math.min(Ec, AP[medium - 1]);
			sum_z = 0.0;
			sum_pz = 0.0;
			sum_a = 0.0;
			sum_wa = 0.0;
			for (int iele = 1; iele <= NNE[imed - 1]; iele++) {
				sum_z = sum_z + PZ[imed - 1][iele - 1]
						* ZELEM[imed - 1][iele - 1];
				sum_pz = sum_pz + PZ[imed - 1][iele - 1];
				sum_wa = sum_wa + RHOZ[imed - 1][iele - 1];
				sum_a = sum_a + PZ[imed - 1][iele - 1] * WA[imed - 1][iele - 1];
			}
			con_med = RHO[imed - 1] / 1.6605655 / sum_a;
			eii_cons[imed - 1] = con_med;
			if (eii_nsh[imed - 1] > 0) {
				sigma_max = 0.0;
				for (int j = 1; j <= MEKE[imed - 1]; j++) {
					loge = (j - EKE0[imed - 1]) / EKE1[imed - 1];
					e = Math.exp(loge);
					tau = e / RM;
					beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
					p2 = 2. * RM * tau * (tau + 2.);
					Wmax = (e + U) / 2.0;
					uwm = U / Wmax;
					lloge = j;
					medium = imed;
					// eg $EVALUATE sigt USING esig(eil);
					// //" sigt is the total cross section "
					// eg
					// sigt=ESIG1[leil-1][medium-1]*eil+ESIG0[leil-1][medium-1];

					// $EVALUATE dedx USING ededx(loge);
					dedx = EDEDX1[lloge - 1][medium - 1] * loge
							+ EDEDX0[lloge - 1][medium - 1];
					if ((e > AP[medium - 1]) || (e > 2 * Ec)) {
						// $EVALUATE sig USING esig(loge);
						sig = ESIG1[lloge - 1][medium - 1] * loge
								+ ESIG0[lloge - 1][medium - 1];
					} else {
						sig = 0.0;
					}
					if (e > 2. * Ec) {
						// $EVALUATE wbrem USING ebr1(loge);
						wbrem = EBR11[lloge - 1][medium - 1] * loge
								+ EBR10[lloge - 1][medium - 1];
						sigm = sig * (1.0 - wbrem);
					} else {
						sigm = 0.0;
						wbrem = 1.0;
					}
					sum_occn = 0.0;
					sum_sigma = 0.0;
					sum_dedx = 0.0;
					for (int iele = 1; iele <= NNE[imed - 1]; iele++) {
						Double dbl = new Double(ZELEM[imed - 1][iele - 1] + 0.5);
						iZ = dbl.intValue(); // int(zelem(imed,iele)+0.5);
						sum_sh = 0.0;
						for (int ish = 1; ish <= eii_no[imed - 1][iele - 1]; ish++) {
							// "jj is the shell index in the list of EII shells "
							jj = eii_first[imed - 1][iele - 1] + ish - 1;
							// "jjj is shell type (1 = K, 2 = LI, 3 = LII, etc.)
							jjj = eii_sh[jj - 1];
							U = binding_energies[jjj - 1][iZ - 1];// (jjj,iZ);
							Wmax = (e + U) / 2.0;
							uwm = U / Wmax;// @@@@@@@@@@@@@@@@@@@@@@@@MUST BE
											// SET AGAIN HERE!!!
							// "IF( Uj >= Ecc ) sum_sh = sum_sh +
							// occn_numbers(jjj);
							if ((U < e) && (U > Ecc)) {
								// " At this energy interactions with this shell will "
								// " be done using the EII differential x-section "
								sum_sh = sum_sh + occn_numbers[jjj - 1];
								ss_0 = 2.0
										* (Math.log(p2 / U) - uwm * uwm * uwm
												* Math.log(p2 / Wmax) - (beta2 + 0.833333)
												* (1.0 - uwm * uwm * uwm))
										/ 3.0 / U;
								sh_0 = ((1.0 - uwm) * (1.0 + uwm / (2.0 - uwm))
										+ U * (Wmax - U)
										/ ((e + RM) * (e + RM)) - (2.0 * tau + 1.0)
										/ ((tau + 1.0) * (tau + 1.0))
										* uwm
										/ 2.0 * Math.log((2.0 - uwm) / uwm))
										/ U;
								ss_1 = Math.log(p2 / U) - uwm * uwm
										* Math.log(p2 / Wmax) - (beta2 + 1.0)
										* (1.0 - uwm * uwm);
								sh_1 = Math.log(Wmax / U / (2.0 - uwm)) + 2.0
										* (Wmax - U) / (2.0 * Wmax - U)
										+ (Wmax * Wmax - U * U)
										/ ((e + RM) * (e + RM)) / 2.0
										- (2.0 * tau + 1.0)
										/ ((tau + 1.0) * (tau + 1.0))
										* Math.log((2.0 * Wmax - U) / Wmax);
								av_E = (ss_1 + sh_1) / (ss_0 + sh_0);
								// "av_E is the average energy lost in a collision"
								// "with this shell"
								Double dbl1 = new Double(eii_a[jjj - 1] * loge
										+ eii_b[jjj - 1]);
								i = dbl1.intValue();// eii_a[jjj-1]*loge +
													// eii_b[jjj-1];
								i = (jj - 1) * $N_EII_BINS + i;
								sig_j = eii_xsection_a[i - 1] * loge
										+ eii_xsection_b[i - 1];
								sig_j = sig_j * PZ[imed - 1][iele - 1]
										* con_med;
								sum_sigma = sum_sigma + sig_j;
								sum_dedx = sum_dedx + sig_j * av_E;
							}
						}
						sum_occn = sum_occn + sum_sh * PZ[imed - 1][iele - 1];
					}
					sigm = sigm + sum_sigma;
					dedx = dedx - sum_dedx;
					aux = Ec / e;
					if (e > 2.0 * Ec) {
						sigo = cons
								* sum_occn
								* RHO[imed - 1]
								/ (beta2 * Ec)
								* ((1.0 - 2.0 * aux)
										* (1.0 + aux / (1.0 - aux) + (tau / (tau + 1.0))
												* (tau / (tau + 1.0))
												* aux
												/ 2.0) - (2.0 * tau + 1.0)
										/ ((tau + 1.0) * (tau + 1.0)) * aux
										* Math.log((1.0 - aux) / aux)) / sum_a;
						de = cons
								* sum_occn
								* RHO[imed - 1]
								/ beta2
								* (Math.log(0.25 / aux / (1.0 - aux))
										+ (1.0 - 2.0 * aux) / (1.0 - aux)
										+ (tau / (tau + 1.0))
										* (tau / (tau + 1.0))
										* (1.0 - 4.0 * aux * aux) / 8.0 - (2.0 * tau + 1.0)
										/ ((tau + 1.0) * (tau + 1.0))
										* Math.log(2.0 * (1.0 - aux))) / sum_a;
						sigm = sigm - sigo;
						// "sigm = sig*(1-wbrem)*(1-sum_occn/sum_z);
						dedx = dedx + de;
					}
					sigma = sigm + wbrem * sig;
					if (sigma / dedx > sigma_max) {
						sigma_max = sigma / dedx;
					}
					if (sigma > 0.0) {
						wbrem = wbrem * sig / sigma;
					} else {
						wbrem = 1.0;
					}
					if (j > 1) {
						EDEDX1[j - 2][imed - 1] = (dedx - dedx_old)
								* EKE1[imed - 1];
						EDEDX0[j - 2][imed - 1] = dedx
								- EDEDX1[j - 2][imed - 1] * loge;
						ESIG1[j - 2][imed - 1] = (sigma - sigma_old)
								* EKE1[imed - 1];
						ESIG0[j - 2][imed - 1] = sigma - ESIG1[j - 2][imed - 1]
								* loge;
						EBR11[j - 2][imed - 1] = (wbrem - wbrem_old)
								* EKE1[imed - 1];
						EBR10[j - 2][imed - 1] = wbrem - EBR11[j - 2][imed - 1]
								* loge;
					}
					dedx_old = dedx;
					// sigm_old = sigm;
					sigma_old = sigma;
					wbrem_old = wbrem;
				}
				EDEDX1[MEKE[imed - 1] - 1][imed - 1] = EDEDX1[MEKE[imed - 1] - 2][imed - 1];
				EDEDX0[MEKE[imed - 1] - 1][imed - 1] = EDEDX0[MEKE[imed - 1] - 2][imed - 1];
				ESIG1[MEKE[imed - 1] - 1][imed - 1] = ESIG1[MEKE[imed - 1] - 2][imed - 1];
				ESIG0[MEKE[imed - 1] - 1][imed - 1] = ESIG0[MEKE[imed - 1] - 2][imed - 1];
				EBR11[MEKE[imed - 1] - 1][imed - 1] = EBR11[MEKE[imed - 1] - 2][imed - 1];
				EBR10[MEKE[imed - 1] - 1][imed - 1] = EBR10[MEKE[imed - 1] - 2][imed - 1];

				esig_e[imed - 1] = sigma_max;// (imed) = sigma_max;
			}
		}

	}

	// "*************************************************************************
	// " The following is a set of macros and subroutines that implement
	// " bremsstrahlung sampling from the S. Seltzer (NIST) cross sections,
	// " which are the basis for ICRU radiative stopping powers,
	// " into the EGSnrc environment.
	// " In order to use it, you have to `turn on' this option by
	// " setting ibr_nist = 1 (which is in COMON/BREMPR/)
	// " If this option is turned on, subroutine HATCH will call
	// " subroutine init_nist_brems.
	// " In init_nist_brems the NIST cross sections are read in,
	// " total bremsstrahlung cross sections are calculated using
	// " 64 point Gauss-Legendre quadrature, the interpolation arrays
	// " that are used for total cross sections and brems fraction
	// interpolations
	// " (esig0, esig1, ebr10, ebr11 for electrons)
	// " (psig0, psig1, pbr10, pbr11, pbr20, pbr21 for positrons)
	// " are updated and alias sampling tables for rapid sampling of brems
	// " energies are created. These alias sampling tables are then used
	// " during the simulation in subroutine BREMS.
	// " Be aware that there is a slight inconsistency when using this option
	// " as resttricted radiative stopping powers used are the ones coming
	// " from PEGS and so, they are calculated using Bethe-Heitler.
	// " This will not matter at all if
	// " - AP is much smaller than the electron energy
	// " and/or
	// " - the restricted radiative stopping power is much smaller
	// " then the restricted collision stopping power
	// " Both conditions are usually satisfied.
	// " I. Kawrakow, NRC, January 2000.
	// "****************************************************************************
	public static double $NIST_ENERGY_SCALE = 1.0;// scale factor

	/**
	 * The following is a set of macros and subroutines that implement 
	 * bremsstrahlung sampling from the S. Seltzer (NIST) cross sections, 
	 * which are the basis for ICRU radiative stopping powers. 
	 * In order to use it, you have to `turn on' this option by setting ibr_nist = 1. 
	 * If this option is turned on, subroutine HATCH will call subroutine init_nist_brems. 
	 * In init_nist_brems the NIST cross sections are read in, total bremsstrahlung cross sections are calculated using 
	 * 64 point Gauss-Legendre quadrature, the interpolation arrays that are used for total cross sections and brems fraction 
	 * interpolations (esig0, esig1, ebr10, ebr11 for electrons) (psig0, psig1, pbr10, pbr11, pbr20, pbr21 for positrons) 
	 * are updated and alias sampling tables for rapid sampling of brems energies are created. These alias sampling tables are then used 
	 * during the simulation in subroutine BREMS. Be aware that there is a slight inconsistency when using this option 
	 * as resttricted radiative stopping powers used are the ones coming from PEGS and so, they are calculated using Bethe-Heitler. 
	 * This will not matter at all if AP is much smaller than the electron energy and/or 
	 * the restricted radiative stopping power is much smaller then the restricted collision stopping power. Both conditions are usually satisfied. 
	 * 
	 */
	public static void init_nist_brems() {
		// "**************************"
		// ; Copyright NRC;
		// int $MXBREN=57; $MXBRXX=30; $MXBREL=100; $MXGAUSS=64;

		double[] energy_array = new double[$MXBREN];
		double[] x_array = new double[$MXBRXX];
		double[][][] cs_array = new double[$MXBREN][$MXBRXX][$MXBREL];
		double[] xi_array = new double[$MXBRXX];
		double[] x_gauss = new double[$MXGAUSS];
		double[] w_gauss = new double[$MXGAUSS];
		// -------------
		int nmix = 0;
		int kmix = 0;
		int ngauss = 0;
		double emin = 0.0;
		int ifirst = 0;
		int i = 0;
		int ilast = 0;
		int ii = 0;
		double sumA = 0.0;
		double Z = 0.0;
		int iz = 0;
		double xi = 0.0;
		double res = 0.0;
		// double spline = 0.0;
		double eil = 0.0;
		double ei = 0.0;
		double beta2 = 0.0;
		double aux = 0.0;
		double sigb = 0.0;
		double sigt = 0.0;
		double ebr1 = 0.0;
		double ebr2 = 0.0;
		int i_gauss = 0;
		int nener = 0;
		double sigee = 0.0;
		double sigep = 0.0;
		double sige = 0.0;
		double si_esig = 0.0;
		double si1_esig = 0.0;
		double si_ebr1 = 0.0;
		double si1_ebr1 = 0.0;
		double ededx = 0.0;
		double sig_bhabha = 0.0;
		double si_psig = 0.0;
		double si1_psig = 0.0;
		double si_pbr1 = 0.0;
		double si1_pbr1 = 0.0;
		double si_pbr2 = 0.0;
		double si1_pbr2 = 0.0;
		int leil = 0;
		double ple = 0.0;
		double qle = 0.0;
		double x = 0.0;
		double f = 0.0;
		double error = 0.0;
		double max_error = 0.0;
		double x_max_error = 0.0;
		double f_max_error = 0.0;
		int ndat = 0;
		int k_max_error = 0;
		// ---------------------
		double[][] cs = new double[$MXBREN][$MXBRXX];
		double[] ee = new double[$MXBREN];
		double[] ele = new double[$MXBREN];
		double[] csx = new double[$MXBRXX];
		double[] afx = new double[$MXBRXX];
		double[] bfx = new double[$MXBRXX];
		double[] cfx = new double[$MXBRXX];
		double[] dfx = new double[$MXBRXX];
		double[] cse = new double[$MXBREN];
		double[] afe = new double[$MXBREN];
		double[] bfe = new double[$MXBREN];
		double[] cfe = new double[$MXBREN];
		double[] dfe = new double[$MXBREN];

		double amu = 1660.5655;// "converts the cross sections from mB/per atom to cm^2/g"
		// $declare_max_medium;//--->nothing
		// "Get the S. Seltzer brems cross sections"->nistbremsfile
		String filename = datas + file_sep + egsDatas + file_sep
				+ nistbremsfile + defaultext;
		int iread = 0;
		int lnr = 0;// data number
		int lnrr = 0;// line number
		// int indx = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals = "=";
		// char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		int eneri = 0;
		boolean enerb = false;
		int xarai = 0;
		boolean xarab = false;
		int csarai = 0;
		boolean csarab = false;
		int n = 0;
		int k = 0;
		FileInputStream in = null;//new FileInputStream(filename);
		try {
			//FileInputStream in = new FileInputStream(filename);
			in = new FileInputStream(filename);
			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if ((lnrr != 0) && ((lnrr != 17)))// skip those lines
						{
							lnr++;
							if (lnr == 1) {
								String s = desc.toString();// System.out.println(s);
								nmix = stringToInt(s);
								// System.out.println("nmix= "+nmix);
							}
							if (lnr == 2) {
								String s = desc.toString();// System.out.println(s);
								kmix = stringToInt(s);
								// System.out.println("kmix= "+kmix);
								if ((kmix != $MXBRXX) || (nmix != $MXBREN)) {
									STOPPROGRAM = true;
									seqStr = " init_nist_brems: wrong data file!";
									// if(iprint>1)
									eq.printSequence(seqStr);

									return;
								}

							}
							if (lnr > 2) {
								if (!enerb) {
									String s = desc.toString();
									energy_array[eneri] = stringToDouble(s);
									// System.out.println("ener "+energy_array[eneri]);
									eneri++;
									if (eneri == $MXBREN) {
										enerb = true;
										// =====================================
										// DO n=1,nmix [ energy_array(n) =
										// $NIST-ENERGY-SCALE*energy_array(n); ]
										for (int nn = 1; nn <= nmix; nn++) {
											energy_array[nn - 1] = $NIST_ENERGY_SCALE
													* energy_array[nn - 1];
										}
										// ======================================
									}
								} else if (!xarab) {
									String s = desc.toString();
									x_array[xarai] = stringToDouble(s);
									// System.out.println("xara "+x_array[xarai]);
									xarai++;
									if (xarai == $MXBRXX)
										xarab = true;

								} else if (!csarab) {
									// cs_array(n,k,i),n=1,nmix),k=1,kmix);
									String s = desc.toString();
									cs_array[n][k][csarai] = stringToDouble(s);
									// System.out.println("cs "+cs_array[n][k][csarai]+" n "+n+" k "+k+" i "+csarai);
									n++;
									if (n == nmix) {
										k++;
										n = 0;
									}
									if (k == kmix) {
										csarai++;
										k = 0;
									}
									if (csarai == $MXBREL)
										csarab = true;

								}

							}// lbr>2

						}

					}// have data
					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			//in.close();
		}// try
		catch (Exception exc) {
			//try 
            //{
             //   in.close();
            //} catch(IOException ex) {
                // log or fail if you like
            //}
		}
		finally {
	        if( null != in ) {
	            try 
	            {
	                in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
		}
		for (k = 1; k <= kmix; k++) {
			xi_array[k - 1] = Math.log(1.0 - x_array[k - 1] + 1.e-6);
			if (fool_intel_optimizer) {
				seqStr = "xi_array(k): " + xi_array[k - 1];// +" \n";
				if (iprint > 2)
					eq.printSequence(seqStr);
			}
		}

		// "Get abscissas and weights for Gauss-Legendre quadrature"

		ngauss = $MXGAUSS;
		gauss_legendre(0.0, 1.0, x_gauss, w_gauss, ngauss);

		// "Calculate total brems cross sections and sampling tables"
		// "for all media                                           "
		seqStr = " Using NIST brems cross sections! ";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		for (int medium = 1; medium <= NMED; medium++) {

			log_ap[medium - 1] = Math.log(AP[medium - 1]);
			seqStr = " Initializing brems data for medium " + medium + "...";// +
																				// " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			emin = Math.max(AE[medium - 1] - RM, AP[medium - 1]);
			for (i = 1; i <= nmix; i++) {
				if (energy_array[i - 1] >= emin)
					break;
			}
			ifirst = i;
			for (i = nmix; i >= 1; i--) {
				if (energy_array[i - 1] < UE[medium - 1] - RM)
					break;
			}
			ilast = i + 1;
			if ((ifirst < 1) || (ilast > nmix)) {
				seqStr = " init_nist_brems: data available only for "
						+ energy_array[i - 1] + " <= E <= "
						+ energy_array[nmix - 1]
						+ " will use spline interpolations to get cross"
						+ " sections beyond the available data but this may"
						+ " produce nonsense!";// + " \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

				// 1-biased!!
				if (ifirst < 1)
					ifirst = 1;
				if (ilast > nmix)
					ilast = nmix;
			}
			for (i = ifirst; i <= ilast; i++) {
				ii = i + 1 - ifirst;
				ee[ii - 1] = energy_array[i - 1];
				ele[ii - 1] = Math.log(ee[ii - 1]);
				sumA = 0;
				for (int j = 1; j <= NNE[medium - 1]; j++) {
					sumA = sumA + PZ[medium - 1][j - 1] * WA[medium - 1][j - 1];
				}
				sumA = sumA * amu;
				for (k = 1; k <= kmix; k++) {
					cs[ii - 1][k - 1] = 0;
					for (int j = 1; j <= NNE[medium - 1]; j++) {
						Z = ZELEM[medium - 1][j - 1];
						Double dbl = new Double(Z + 0.1);
						iz = dbl.intValue();
						Z = Z * Z / sumA;
						cs[ii - 1][k - 1] = cs[ii - 1][k - 1]
								+ PZ[medium - 1][j - 1] * Z
								* cs_array[i - 1][k - 1][iz - 1];
					}
					csx[k - 1] = Math.log(cs[ii - 1][k - 1]);
				}
				set_spline(xi_array, csx, afx, bfx, cfx, dfx, kmix);

				// " Integrate

				cse[ii - 1] = 0;
				aux = Math.log(ee[ii - 1] / AP[medium - 1]);
				for (i_gauss = 1; i_gauss <= ngauss; i_gauss++) {
					xi = Math.log(1.0 - AP[medium - 1] / ee[ii - 1]
							* Math.exp(x_gauss[i_gauss - 1] * aux) + 1.e-6);
					res = spline(xi, xi_array, afx, bfx, cfx, dfx, kmix);
					// res = spline(xi,kmix);
					cse[ii - 1] = cse[ii - 1] + w_gauss[i_gauss - 1]
							* Math.exp(res);
				}
			}
			nener = ilast - ifirst + 1;
			set_spline(ele, cse, afe, bfe, cfe, dfe, nener);

			// " Now replace the PEGS cross sections "

			NEKE = MEKE[medium - 1];
			sigee = 1.E-15;
			sigep = 1.E-15;
			for (i = 1; i <= NEKE; i++) {
				Integer intg = new Integer(i);
				eil = (intg.doubleValue() - EKE0[medium - 1])
						/ EKE1[medium - 1];
				ei = Math.exp(eil);
				leil = i;
				beta2 = ei * (ei + 2.0 * RM) / ((ei + RM) * (ei + RM));
				if (ei <= AP[medium - 1]) {
					sigb = 1.e-30;
				} else {
					sigb = spline(eil, ele, afe, bfe, cfe, dfe, nener);
					// sigb = spline(eil,nener);
					sigb = sigb * Math.log(ei / AP[medium - 1]) / beta2
							* RHO[medium - 1];
				}
				// $EVALUATE sigt USING esig(eil);
				// //" sigt is the total cross section "
				sigt = ESIG1[leil - 1][medium - 1] * eil
						+ ESIG0[leil - 1][medium - 1];
				// $EVALUATE ededx USING ededx(eil);
				// ededx =
				// EDEDX1[leil-1][medium-1]*eil+EDEDX0[leil-1][medium-1];

				// $EVALUATE ebr1 USING ebr1(eil);
				// //" coming from PEGS, ebr1*sigt is "
				// " then the brems cross section "
				ebr1 = EBR11[leil - 1][medium - 1] * eil
						+ EBR10[leil - 1][medium - 1];

				if (sigt < 0.0)
					sigt = 0.0;
				if (ebr1 > 1.0)
					ebr1 = 1.0;
				if (ebr1 < 0.0)
					ebr1 = 0.0;
				if (i > 1) {
					si_esig = si1_esig;
					si_ebr1 = si1_ebr1;
					si1_esig = sigt * (1.0 - ebr1) + sigb;
					si1_ebr1 = sigb / si1_esig;
					ESIG1[i - 2][medium - 1] = (si1_esig - si_esig)
							* EKE1[medium - 1];
					ESIG0[i - 2][medium - 1] = si1_esig
							- ESIG1[i - 2][medium - 1] * eil;
					EBR11[i - 2][medium - 1] = (si1_ebr1 - si_ebr1)
							* EKE1[medium - 1];
					EBR10[i - 2][medium - 1] = si1_ebr1
							- EBR11[i - 2][medium - 1] * eil;
				} else {
					si1_esig = sigt * (1.0 - ebr1) + sigb;
					si1_ebr1 = sigb / si1_esig;
				}

				// " Positrons "

				// $EVALUATE sigt USING psig(eil);
				sigt = PSIG1[leil - 1][medium - 1] * eil
						+ PSIG0[leil - 1][medium - 1];
				// $EVALUATE ebr1 USING pbr1(eil);
				ebr1 = PBR11[leil - 1][medium - 1] * eil
						+ PBR10[leil - 1][medium - 1];
				// $EVALUATE ebr2 USING pbr2(eil);
				ebr2 = PBR21[leil - 1][medium - 1] * eil
						+ PBR20[leil - 1][medium - 1];
				if (sigt < 0.0)
					sigt = 0.0;
				if (ebr1 > 1.0)
					ebr1 = 1.0;
				if (ebr1 < 0.0)
					ebr1 = 0.0;
				if (ebr2 > 1.0)
					ebr2 = 1.0;
				if (ebr2 < 0.0)
					ebr2 = 0.0;
				sig_bhabha = sigt * (ebr2 - ebr1);
				if (sig_bhabha < 0.0)
					sig_bhabha = 0.0;
				if (i > 1) {
					si_psig = si1_psig;
					si_pbr1 = si1_pbr1;
					si_pbr2 = si1_pbr2;
					si1_psig = sigt * (1.0 - ebr1) + sigb;
					si1_pbr1 = sigb / si1_psig;
					si1_pbr2 = (sigb + sig_bhabha) / si1_psig;
					PSIG1[i - 2][medium - 1] = (si1_psig - si_psig)
							* EKE1[medium - 1];
					PSIG0[i - 2][medium - 1] = si1_psig
							- PSIG1[i - 2][medium - 1] * eil;
					PBR11[i - 2][medium - 1] = (si1_pbr1 - si_pbr1)
							* EKE1[medium - 1];
					PBR10[i - 2][medium - 1] = si1_pbr1
							- PBR11[i - 2][medium - 1] * eil;
					PBR21[i - 2][medium - 1] = (si1_pbr2 - si_pbr2)
							* EKE1[medium - 1];
					PBR20[i - 2][medium - 1] = si1_pbr2
							- PBR21[i - 2][medium - 1] * eil;
				} else {
					si1_psig = sigt * (1.0 - ebr1) + sigb;
					si1_pbr1 = sigb / si1_psig;
					si1_pbr2 = (sigb + sig_bhabha) / si1_psig;
				}
				// $EVALUATE ededx USING ededx(eil);
				ededx = EDEDX1[leil - 1][medium - 1] * eil
						+ EDEDX0[leil - 1][medium - 1];
				sige = si1_esig / ededx;
				if (sige > sigee)
					sigee = sige;
				// $EVALUATE ededx USING pdedx(eil);
				ededx = PDEDX1[leil - 1][medium - 1] * eil
						+ PDEDX0[leil - 1][medium - 1];
				sige = si1_psig / ededx;
				if (sige > sigep)
					sigep = sige;
				// "write(6,*) ei,leil,si1_esig,si1_psig;
			}
			ESIG1[NEKE - 1][medium - 1] = ESIG1[NEKE - 2][medium - 1];
			ESIG0[NEKE - 1][medium - 1] = ESIG0[NEKE - 2][medium - 1];
			EBR11[NEKE - 1][medium - 1] = EBR11[NEKE - 2][medium - 1];
			EBR10[NEKE - 1][medium - 1] = EBR10[NEKE - 2][medium - 1];
			PSIG1[NEKE - 1][medium - 1] = PSIG1[NEKE - 2][medium - 1];
			PSIG0[NEKE - 1][medium - 1] = PSIG0[NEKE - 2][medium - 1];
			PBR11[NEKE - 1][medium - 1] = PBR11[NEKE - 2][medium - 1];
			PBR10[NEKE - 1][medium - 1] = PBR10[NEKE - 2][medium - 1];
			PBR21[NEKE - 1][medium - 1] = PBR21[NEKE - 2][medium - 1];
			PBR20[NEKE - 1][medium - 1] = PBR20[NEKE - 2][medium - 1];

			seqStr = " Max. new cross sections per energy loss: " + sigee
					+ " , " + sigep;// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			esig_e[medium - 1] = sigee;
			psig_e[medium - 1] = sigep;
			if (sigee > esige_max)
				esige_max = sigee;
			if (sigep > psige_max)
				psige_max = sigep;

			// " Now prepare the arrays for brems sampling

			nb_emin[medium - 1] = energy_array[ifirst - 1];
			if (nb_emin[medium - 1] <= AP[medium - 1]) {
				nb_emin[medium - 1] = energy_array[ifirst];
			}
			nb_emax[medium - 1] = energy_array[ilast - 1];
			nb_lemin[medium - 1] = Math.log(nb_emin[medium - 1]);
			nb_lemax[medium - 1] = Math.log(nb_emax[medium - 1]);
			nb_dle[medium - 1] = (nb_lemax[medium - 1] - nb_lemin[medium - 1])
					/ ($MXBRES - 1.0);
			nb_dlei[medium - 1] = 1.0 / nb_dle[medium - 1];
			// "
			eil = nb_lemin[medium - 1] - nb_dle[medium - 1];
			for (i = 1; i <= $MXBRES; i++) {
				eil = eil + nb_dle[medium - 1];
				ei = Math.exp(eil);
				for (ii = 1; ii <= nener; ii++) {
					if (ei < ee[ii - 1])
						break;
				}
				ii = ii - 1;
				if (ii > nener - 1)
					ii = nener - 1;
				// -----------
				if (ii == 0)
					ii = 1;
				// ----------
				// " ple and qle are energy interpolation coefficients

				ple = (eil - ele[ii - 1]) / (ele[ii] - ele[ii - 1]);
				qle = 1.0 - ple;
				for (k = 1; k <= $MXBRXX; k++) {
					csx[k - 1] = Math.log(qle * cs[ii - 1][k - 1] + ple
							* cs[ii][k - 1]);
				}
				set_spline(xi_array, csx, afx, bfx, cfx, dfx, kmix);

				// " fill the abscissas for this energy

				x = AP[medium - 1] / ei;
				aux = -Math.log(x);
				xi = Math.log(1.0 - x + 1.e-6);
				res = spline(xi, xi_array, afx, bfx, cfx, dfx, kmix);
				// res = spline(xi,kmix);
				nb_xdata[0][i - 1][medium - 1] = 0.0;
				nb_fdata[0][i - 1][medium - 1] = Math.exp(res);

				for (k = 1; k <= $MXBRXX; k++) {
					if (x_array[k - 1] > x)
						break;
				}
				if (k > $MXBRXX)
					k = $MXBRXX;
				ndat = 0;
				for (int j = k + 1; j <= $MXBRXX - 1; j++) {
					ndat = ndat + 1;
					nb_xdata[ndat][i - 1][medium - 1] = Math.log(x_array[j - 1]
							/ x)
							/ aux;
					nb_fdata[ndat][i - 1][medium - 1] = Math.exp(csx[j - 1]);
					if (fool_intel_optimizer) {

						seqStr = "nb_xdata(ndat,i,medium): "
								+ nb_xdata[ndat][i - 1][medium - 1];// + " \n";
						if (iprint > 2)
							eq.printSequence(seqStr);

					}
				}
				ndat = ndat + 1;
				nb_xdata[ndat][i - 1][medium - 1] = 1.0;
				nb_fdata[ndat][i - 1][medium - 1] = Math.exp(csx[$MXBRXX - 1]);

				// " Now expand the arrays by filling intermediate points
				// " at the positions that show the maxium relative error
				// " when using linear interpolation in x.
				// " If arrays were allocated dynamically one could use
				// " a certain condition to stop the iteration but in our case
				// " memory is allocated anyway and so we use the maximum
				// " space provided

				// if( ndat = $MXBRXS ) //goto :SKIP-LOOP:;
				if (ndat != $MXBRXS) {
					while (ndat < $MXBRXS) {
						x_max_error = 0.0;
						f_max_error = 0.0;
						k_max_error = 0;
						max_error = 0.0;
						for (k = 0; k <= ndat - 1; k++) {
							x = 0.5 * (nb_xdata[k][i - 1][medium - 1] + nb_xdata[k + 1][i - 1][medium - 1]);
							f = 0.5 * (nb_fdata[k][i - 1][medium - 1] + nb_fdata[k + 1][i - 1][medium - 1]);
							xi = Math.log(1.0 - AP[medium - 1] / ei
									* Math.exp(x * aux) + 1.e-6);
							res = spline(xi, xi_array, afx, bfx, cfx, dfx, kmix);
							// res = spline(xi,kmix);
							res = Math.exp(res);
							error = Math.abs(1.0 - f / res);
							if (error > max_error) {
								x_max_error = x;
								f_max_error = res;
								max_error = error;
								k_max_error = k;
							}
						}
						ndat = ndat + 1;
						for (k = ndat; k >= k_max_error + 2; k--) {
							nb_xdata[k][i - 1][medium - 1] = nb_xdata[k - 1][i - 1][medium - 1];
							nb_fdata[k][i - 1][medium - 1] = nb_fdata[k - 1][i - 1][medium - 1];
						}
						nb_xdata[k_max_error + 1][i - 1][medium - 1] = x_max_error;
						nb_fdata[k_max_error + 1][i - 1][medium - 1] = f_max_error;
					}// UNTIL (ndat = $MXBRXS);
				}

				// :SKIP-LOOP:

				// " Now generate the alias tables for rapid brems sampling
				// " during run time

				// prepare_alias_table($MXBRXS,nb_xdata[0][i-1][medium-1],
				// nb_fdata[0][i-1][medium-1],nb_wdata[1][i-1][medium-1],nb_idata[1][i-1][medium-1]);
				double[] patx = new double[$MXBRXS + 1];
				double[] patf = new double[$MXBRXS + 1];
				double[] patw = new double[$MXBRXS];
				int[] pati = new int[$MXBRXS];
				for (int ll = 0; ll <= $MXBRXS; ll++) {
					patx[ll] = nb_xdata[ll][i - 1][medium - 1];
					patf[ll] = nb_fdata[ll][i - 1][medium - 1];
				}
				prepare_alias_table($MXBRXS, patx, patf, patw, pati);

				for (int ll = 0; ll < $MXBRXS; ll++) {
					nb_wdata[ll][i - 1][medium - 1] = patw[ll];
					nb_idata[ll][i - 1][medium - 1] = pati[ll];
					// System.out.println(nb_wdata[ll][i-1][medium-1]);
				}

			}
		}

	}

	/**
	 * Calculates the parameter for the rejection function used in 
	 * the current implementation of bremsstrahlung sampling. It is called by HATCH.
	 */
	private static void fix_brems() {
		// "******************************************************************"
		// " Calculates the parameter for the rejection function used in
		// " the current implementation of bremsstrahlung sampling
		// " I Kawrakow, January 2000
		// "*******************************************************************"
		// ; Copyright NRC;
		double Zt = 0.0;
		double Zb = 0.0;
		double Zf = 0.0;
		double Zg = 0.0;
		double Zv = 0.0;
		double fmax1 = 0.0;
		double fmax2 = 0.0;
		double Zi = 0.0;
		double pi = 0.0;
		double fc = 0.0;
		double xi = 0.0;
		double aux = 0.0;

		// $declare_max_medium;
		for (int medium = 1; medium <= NMED; medium++) {

			log_ap[medium - 1] = Math.log(AP[medium - 1]);
			Zt = 0;
			Zb = 0;
			Zf = 0;
			for (int i = 1; i <= NNE[medium - 1]; i++) {
				Zi = ZELEM[medium - 1][i - 1];
				pi = PZ[medium - 1][i - 1];
				fc = FCOULC(Zi);
				xi = XSIF(Zi);
				aux = pi * Zi * (Zi + xi);// ZZX[I] = PZ[I]*Z[I]*(Z[I]+XSI[I]);
				Zt = Zt + aux;// ZT = ZT + ZZX[I];
				Zb = Zb - aux * Math.log(Zi) / 3.0;// ZB = ZB +
													// ZZX[I]*Math.log(Math.pow(Z[I],-1./3.));
				Zf = Zf + aux * fc;// ZF = ZF + ZZX[I]*FCOUL[I];
			}
			Zv = (Zb - Zf) / Zt;// ZV = (ZB-ZF)/ZT;
			Zg = Zb / Zt;// ZG = ZB/ZT;
			fmax1 = 2.0 * (20.863 + 4.0 * Zg) - 2.0 * (20.029 + 4.0 * Zg) / 3.0;
			fmax2 = 2.0 * (20.863 + 4.0 * Zv) - 2.0 * (20.029 + 4.0 * Zv) / 3.0;

			// "These are used in BREMS"
			DL1[0][medium - 1] = (20.863 + 4.0 * Zg) / fmax1;
			DL2[0][medium - 1] = -3.242 / fmax1;
			DL3[0][medium - 1] = 0.625 / fmax1;
			DL4[0][medium - 1] = (21.12 + 4.0 * Zg) / fmax1;
			DL5[0][medium - 1] = -4.184 / fmax1;
			DL6[0][medium - 1] = 0.952;
			DL1[1][medium - 1] = (20.029 + 4.0 * Zg) / fmax1;
			DL2[1][medium - 1] = -1.93 / fmax1;
			DL3[1][medium - 1] = -0.086 / fmax1;
			DL4[1][medium - 1] = (21.12 + 4.0 * Zg) / fmax1;
			DL5[1][medium - 1] = -4.184 / fmax1;
			DL6[1][medium - 1] = 0.952;
			DL1[2][medium - 1] = (20.863 + 4.0 * Zv) / fmax2;
			DL2[2][medium - 1] = -3.242 / fmax2;
			DL3[2][medium - 1] = 0.625 / fmax2;
			DL4[2][medium - 1] = (21.12 + 4.0 * Zv) / fmax2;
			DL5[2][medium - 1] = -4.184 / fmax2;
			DL6[2][medium - 1] = 0.952;
			DL1[3][medium - 1] = (20.029 + 4.0 * Zv) / fmax2;
			DL2[3][medium - 1] = -1.93 / fmax2;
			DL3[3][medium - 1] = -0.086 / fmax2;
			DL4[3][medium - 1] = (21.12 + 4.0 * Zv) / fmax2;
			DL5[3][medium - 1] = -4.184 / fmax2;
			DL6[3][medium - 1] = 0.952;

			// "and these in PAIR"
			DL1[4][medium - 1] = (3.0 * (20.863 + 4.0 * Zg) - (20.029 + 4.0 * Zg));
			DL2[4][medium - 1] = (3.0 * (-3.242) - (-1.930));
			DL3[4][medium - 1] = (3.0 * (0.625) - (-0.086));
			DL4[4][medium - 1] = (2.0 * 21.12 + 8.0 * Zg);
			DL5[4][medium - 1] = (2.0 * (-4.184));
			DL6[4][medium - 1] = 0.952;
			DL1[5][medium - 1] = (3.0 * (20.863 + 4.0 * Zg) + (20.029 + 4.0 * Zg));
			DL2[5][medium - 1] = (3.0 * (-3.242) + (-1.930));
			DL3[5][medium - 1] = (3.0 * 0.625 + (-0.086));
			DL4[5][medium - 1] = (4.0 * 21.12 + 16.0 * Zg);
			DL5[5][medium - 1] = (4.0 * (-4.184));
			DL6[5][medium - 1] = 0.952;
			DL1[6][medium - 1] = (3.0 * (20.863 + 4.0 * Zv) - (20.029 + 4.0 * Zv));
			DL2[6][medium - 1] = (3.0 * (-3.242) - (-1.930));
			DL3[6][medium - 1] = (3.0 * (0.625) - (-0.086));
			DL4[6][medium - 1] = (2.0 * 21.12 + 8.0 * Zv);
			DL5[6][medium - 1] = (2.0 * (-4.184));
			DL6[6][medium - 1] = 0.952;
			DL1[7][medium - 1] = (3.0 * (20.863 + 4.0 * Zv) + (20.029 + 4.0 * Zv));
			DL2[7][medium - 1] = (3.0 * (-3.242) + (-1.930));
			DL3[7][medium - 1] = (3.0 * 0.625 + (-0.086));
			DL4[7][medium - 1] = (4.0 * 21.12 + 16.0 * Zv);
			DL5[7][medium - 1] = (4.0 * (-4.184));
			DL6[7][medium - 1] = 0.952;

			BPAR[1][medium - 1] = DL1[6][medium - 1]
					/ (3.0 * DL1[7][medium - 1] + DL1[6][medium - 1]);
			BPAR[0][medium - 1] = 12.0 * DL1[7][medium - 1]
					/ (3.0 * DL1[7][medium - 1] + DL1[6][medium - 1]);

		}
	}

	/**
	 * Used by fix_brems. See EGSnrc documentation (PIRS docs).
	 * @param Z Z
	 * @return the result
	 */
	public static double FCOULC(double Z) {
		// "************************"
		// ; Copyright NRC;
		double FCOULC = 0.0;
		double fine = 137.03604;
		double asq = Z / fine;
		asq = asq * asq;
		FCOULC = asq
				* (1.0 / (1.0 + asq) + 0.20206 + asq
						* (-0.0369 + asq * (0.0083 + asq * (-0.002))));

		return FCOULC;
	}

	/**
	 * Used by fix_brems. See EGSnrc documentation (PIRS docs).
	 * @param Z Z
	 * @return the result
	 */
	public static double XSIF(double Z) {
		// "**********************"
		// ; Copyright NRC;
		double xsif = 0.0;

		double[] alrad = { 5.31, 4.79, 4.74, 4.71 };
		double[] alradp = { 6.144, 5.621, 5.805, 5.924 };
		double a1440 = 1194.0;
		double a183 = 184.15;
		if (Z <= 4.0) {
			Double dbl = new Double(Z);
			int iZ = dbl.intValue();
			xsif = alradp[iZ - 1] / (alrad[iZ - 1] - FCOULC(Z));
		} else {
			xsif = Math.log(a1440 * Math.pow(Z, -0.666667))
					/ (Math.log(a183 * Math.pow(Z, -0.33333)) - FCOULC(Z));
		}
		return xsif;
	}

	// "******************************************************************"
	/**
	 * Reads in bound Compton scattering data and performs necessary initializations.
	 * Called by HATCH.
	 */
	public static void init_compton() {
		// " Reads in bound comptin scattering data from unit $INCOHUNIT
		// " and performs necessary initializations
		// " See definitions of variables in egsnrc.macros with definition
		// " of COMIN/COMPTON-DATA/
		// " I.Kawrakow, January 99
		// "******************************************************************"
		// ; Copyright NRC;
		int medium = 0;
		double rm = 0.5110034;
		double aux_erf = 0.0;
		double aux = 0.0;
		double pztot = 0.0;
		int nsh = 0;
		int iz = 0;
		boolean getd = false;

		// " Initialize radiative Compton corrections, if needed "
		RADC_HATCH();
		// $need_bound_compton_data(getd);
		// REPLACE {$need_bound_compton_data(#);} WITH {
		getd = false;
		// {P1} = .false.;
		for (int j = 0; j < $MXREG; j++) {
			medium = MED[j];
			if ((medium > 0) && (medium <= NMED)) {
				if (ibcmp[j] > 0) {
					getd = true;// {P1} = .true.;
					break;
				}
			}
		}
		// };

		if (!getd) {
			seqStr = " Bound Compton scattering not requested! ";// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			return;
		}
		seqStr = " Bound Compton scattering requested, reading data ......";// +
																			// " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		String filename = datas + file_sep + egsDatas + file_sep + incohfile
				+ defaultext;
		int iread = 0;
		int lnr = 0;// data number
		int indx = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals = "=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		int eni = 0;
		boolean enb = false;

		try {
			FileInputStream in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;

				} else {

					if (haveData)// we have data
					{
						haveData = false;// reset

						if (lnr >= 18)// first 18 lines are skipped
						{
							// String s=desc.toString();
							// System.out.println(s);
							if ((eni < $MXTOTSH) && (!enb)) {
								// first data is already readed
								String ss1 = desc.toString();
								iz_array[eni] = stringToInt(ss1);
								// System.out.println("iz "+iz_array[eni]);

								// int kindx=0;
								int j = 0;
								desc = new StringBuffer();
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											appendB = false;// read data so
															// ready to
															// loop!!!!!
											String s1 = desc.toString();// System.out.println(" ## "+s1);
											if (j == 0) {
												shn_array[eni] = stringToInt(s1);
											} else if (j == 1) {
												ne_array[eni] = stringToInt(s1);
											} else if (j == 2) {
												Jo_array[eni] = stringToDouble(s1);
												Jo_array[eni] = Jo_array[eni] * 137.;
												aux_erf = 0.70710678119 * (1.0 + 0.3 * Jo_array[eni]);
												erfJo_array[eni] = 0.82436063535 * (ERF1(aux_erf) - 1.0);
												// "0.82436063535 is exp(0.5)/2"
											} else if (j == 3) {
												be_array[eni] = stringToDouble(s1);
												be_array[eni] = be_array[eni]
														* 1.e-6 / rm;
											}

											j++;

											if (j == 4) {
												// System.out.println("  i: "+eni+"  shn  "+shn_array[eni]+" ne "+ne_array[eni]+
												// " jo "+Jo_array[eni]+" be "+be_array[eni]+" errf "+erfJo_array[eni]);

												j = 0;
												b = false;// exit the loop
											}

										}
										desc = new StringBuffer();// reset
									}
								}

								// next set
								eni++;
								if (eni == $MXTOTSH)
									enb = true;
							}

						}

					}// have data
					if ((char) iread == lineSep)
						lnr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
		}// try
		catch (Exception exc) {

		}

		seqStr = " Done";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		seqStr = " Initializing Bound Compton scattering ......";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		for (medium = 1; medium <= NMED; medium++) {
			pztot = 0;
			nsh = 0;
			for (int i = 1; i <= NNE[medium - 1]; i++) {
				Double dbl = new Double(ZELEM[medium - 1][i - 1]);
				iz = dbl.intValue();
				for (int j = 1; j <= $MXTOTSH; j++) {
					if (iz == iz_array[j - 1]) {
						nsh = nsh + 1;
						if (nsh > $MXMDSH) {
							STOPPROGRAM = true;
							seqStr = " For medium " + medium
									+ " the number of shells is > " + $MXMDSH
									+ "!" + " \n"
									+ " Increase the parameter $MXMDSH! ";
							// if(iprint>1)
							eq.printSequence(seqStr);

							return;
						}
						shell_array[nsh - 1][medium - 1] = j;
						aux = PZ[medium - 1][i - 1] * ne_array[j - 1];
						eno_array[nsh - 1][medium - 1] = aux;
						pztot = pztot + aux;
					}
				}
			}
			if (nsh == 0) {
				seqStr = " Medium " + medium + " has zero shells! ";// + " \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

				return;
			}
			n_shell[medium - 1] = nsh;

			seqStr = " Medium " + medium + " has " + nsh + " shells: ";// +
																		// " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			for (int i = 1; i <= nsh; i++) {
				int j = shell_array[i - 1][medium - 1];
				eno_array[i - 1][medium - 1] = eno_array[i - 1][medium - 1]
						/ pztot;
				// ----------
				aux = be_array[j - 1] * rm * 1000.;
				// -----------
				// OUTPUT i,j,shn_array(j),eno_array(i,medium),
				// Jo_array(j),be_array(j)*rm*1000.;
				// (i3,i4,i3,f9.5,e10.3,f10.3);
				seqStr = "i:" + format(i, 3) + " j:" + format(j, 4) + " shn:"
						+ format(shn_array[j - 1], 3) + " eno:"
						+ format(eno_array[i - 1][medium - 1], 9, true)
						+ " jo:" + format(Jo_array[j - 1], 10, false) + " be:"
						+ format(aux, 10, true);// + " \n";
				if (iprint > 1)
					eq.printSequence(seqStr);
			}
		}

		seqStr = "...... Done.";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// $need_relaxation_data(getd);
		// REPLACE {$need_relaxation_data(#);} WITH {
		// {P1} = .false.;
		getd = false;
		for (int j = 0; j < $MXREG; j++) {
			if ((iedgfl[j] > 0) && (iedgfl[j] <= $MXELEMENT)) {
				// {P1} = .true.;
				getd = true;
				break;
			}
		}
		// };
		if (getd)
			return;

		seqStr = " In subroutine init_compton: fluorescence not set but relaxation"
				+ "data are required for bound Compton scattering.-->   calling EDGSET. ";// +
																							// " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		iedgfl[0] = 1;// iedgfl[1] = 1; //"This was (2) originally  DR"
		EDGSET(1, 1);
		iedgfl[0] = 0;// iedgfl[1] = 0; //"This was (2) originally  DR"

	}

	// "******************************************************************"
	/**
	 * Set up parameters for atomic relaxation and proper handling of photo-electric absorption. 
	 * Input:  NREGLO and NREGHI, not needed but left there for compatibility with older user codes. 
	 * This routine is called from HATCH, it checks whether one of the 
	 * elements of IEDGFL has been set to an integer number between 
	 * 1 and 100 and if so reads in photo-absorption and relaxation 
	 * data. Note that the array IEDGFL, which used to be the 
	 * `effective' atomic number of a medium is not used for any 
	 * purpose other than to indicate that relaxations are requested (if non-zero).
	 * @param NREGLO NREGLO
	 * @param NREGHI NREGHI
	 */
	public static void EDGSET(int NREGLO, int NREGHI) {
		// "******************************************************************"
		// " SUBPROGRAM TO SET UP PARAMETERS FOR ATOMIC RELAXATIONS           "
		// " and proper handling of photo-electric absorption
		// "******************************************************************"
		// " Programmer:   I. Kawrakow, (NRC)                                 "
		// "******************************************************************"
		// "  Input:  NREGLO and NREGHI, not needed but left there for        "
		// "                             compatibility with older user codes  "
		// "  This routine is called from HATCH, it checks whether one of the "
		// "  elements of IEDGFL has been set to an integer number between    "
		// "  1 and 100 and if so reads in photo-absorption and relaxation    "
		// "  data. Note that the array IEDGFL, which used to be the          "
		// "  `effecvtive' atomic number of a medium is not used for any      "
		// "  purpose other than to indicate that relaxations are requested   "
		// " (if non-zero).
		// "******************************************************************"
		// ; Copyright NRC;
		// ;COMIN/EDGE,X-OPTIONS/;
		// "Input variables"
		// integer NREGLO,NREGHI;

		// $INTEGER i,j,k,jj,iz;
		// logical do_relax;
		// logical got_data;
		// save got_data;
		// data got_data/.false./;
		String filename = datas + file_sep + egsDatas + file_sep + phrelaxfile
				+ defaultext;
		int iread = 0;
		@SuppressWarnings("unused")
		int lnr = 0;// data number
		int indx = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals = "=";
		char comma = ',';
		int beni = 0;
		boolean benb = false;
		int ipri = 0;
		boolean iprb = false;
		int rki = 0;
		boolean rkb = false;
		int rl1i = 0;
		boolean rl1b = false;
		int rl2i = 0;
		boolean rl2b = false;
		int rl3i = 0;
		boolean rl3b = false;
		int rmi = 0;
		boolean rmb = false;

		// boolean got_data=false;
		if (got_data)
			return;
		// "EDGSET is now called from HATCH. In older user codes it was called"
		// "from within the user code. If this happens, and the data is already"
		// "available, we don't need to read it again. That's why the above"
		// "statement."
		seqStr = " Output from subroutine EDGSET: ==============================";// +
																					// " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// $need_relaxation_data(do_relax);
		// REPLACE {$need_relaxation_data(#);} WITH {
		// {P1} = .false.;
		boolean do_relax = false;
		for (int j = 0; j < $MXREG; j++) {
			if ((iedgfl[j] > 0) && (iedgfl[j] <= $MXELEMENT)) {
				// {P1} = .true.;
				do_relax = true;
				break;
			}
		}
		// };

		if (!do_relax) {
			seqStr = " Atomic relaxations not requested! ";// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			return;
		}
		seqStr = "  Atomic relaxations requested! ";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		seqStr = " Reading photo-absorption data .....";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		got_data = true;

		try {
			FileInputStream in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if ((beni < $MXELEMENT) && (!benb)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->binding_energies(k,i),k=1,$MXSHELL);//$MXSHELL=6
							int kindx = 0;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										binding_energies[kindx][beni] = stringToDouble(s1);
										// System.out.println("ben  "+binding_energies[kindx][beni]+"  i: "+beni+"  k: "+kindx);
										kindx++;
										if (kindx == $MXSHELL) {
											for (int k = 0; k < $MXSHELL; k++) {
												binding_energies[k][beni] = binding_energies[k][beni] * 1.e-6; // "Convert to MeV"
												// System.out.println("ben  "+binding_energies[k][beni]+"  i: "+beni+"  k: "+k);
											}

											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							beni++;
							if (beni == $MXELEMENT)
								benb = true;
						} else if ((ipri < $MXELEMENT) && (!iprb)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->(interaction_prob(k,i),k=1,$MXINTER-1);//$MXINTER=5
							int kindx = 0;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										interaction_prob[kindx][ipri] = stringToDouble(s1);
										// System.out.println("ipr  "+interaction_prob[kindx][ipri]+"  i: "+ipri+"  k: "+kindx);
										kindx++;
										if (kindx == $MXINTER) {
											b = false;// exit the loop
											// interaction_prob($MXSHELL,i)=1.01;
											interaction_prob[kindx][ipri] = 1.01;
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							ipri++;
							if (ipri == $MXELEMENT) {
								iprb = true;

								seqStr = "Done";// + " \n";
								if (iprint > 1)
									eq.printSequence(seqStr);

								seqStr = " Reading relaxation data .... ";// +
																			// " \n";
								if (iprint > 1)
									eq.printSequence(seqStr);

							}
						} else if ((rki < $MXELEMENT) && (!rkb)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->(relaxation_prob(k,i),k=1,19); //"K-shell"
							int kindx = 0;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										relaxation_prob[kindx][rki] = stringToDouble(s1);
										// System.out.println("rki  "+relaxation_prob[kindx][rki]+"  i: "+rki+"  k: "+kindx);
										kindx++;
										if (kindx == 19) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							rki++;
							if (rki == $MXELEMENT) {
								rkb = true;
							}
						} else if ((rl1i < $MXELEMENT) && (!rl1b)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->(relaxation_prob(k,i),k=20,26); "L1-shell"
							int kindx = 19;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										relaxation_prob[kindx][rl1i] = stringToDouble(s1);
										// System.out.println("rl1i  "+relaxation_prob[kindx][rl1i]+"  i: "+rl1i+"  k: "+kindx);
										kindx++;
										if (kindx == 26) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							rl1i++;
							if (rl1i == $MXELEMENT) {
								rl1b = true;
							}
						} else if ((rl2i < $MXELEMENT) && (!rl2b)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->(relaxation_prob(k,i),k=27,32); "L2-shell"
							int kindx = 26;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										relaxation_prob[kindx][rl2i] = stringToDouble(s1);
										// System.out.println("rl2i  "+relaxation_prob[kindx][rl2i]+"  i: "+rl2i+"  k: "+kindx);
										kindx++;
										if (kindx == 32) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							rl2i++;
							if (rl2i == $MXELEMENT) {
								rl2b = true;
							}
						} else if ((rl3i < $MXELEMENT) && (!rl3b)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->(relaxation_prob(k,i),k=33,37); "L3-shell"
							int kindx = 32;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										relaxation_prob[kindx][rl3i] = stringToDouble(s1);
										// System.out.println("rl3i  "+relaxation_prob[kindx][rl3i]+"  i: "+rl3i+"  k: "+kindx);
										kindx++;
										if (kindx == 37) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							rl3i++;
							if (rl3i == $MXELEMENT) {
								rl3b = true;
							}
						} else if ((rmi < $MXELEMENT) && (!rmb)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->relaxation_prob(38,i); "M-shell"
							int kindx = 37;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										relaxation_prob[kindx][rmi] = stringToDouble(s1);
										// System.out.println("rmi  "+relaxation_prob[kindx][rmi]+"  i: "+rmi+"  k: "+kindx);
										kindx++;
										if (kindx == 38) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							rmi++;
							if (rmi == $MXELEMENT) {
								rmb = true;
							}
						}

					}// have data
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
		}// try
		catch (Exception exc) {

		}

		seqStr = " Done";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		seqStr = " Reading photo cross section data .... ";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		filename = datas + file_sep + egsDatas + file_sep + phcsfile
				+ defaultext;
		iread = 0;
		lnr = 0;// data number
		indx = 0;
		desc = new StringBuffer();
		haveData = false;
		// equals="=";comma=',';
		// edge_a=new double[$MXEDGE][$MXELEMENT];$MXEDGE=16
		int eni = 0;
		boolean enb = false;

		try {
			FileInputStream in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if ((eni < $MXELEMENT) && (!enb)) {
							// first index is already readed
							// String ss1=desc.toString();
							// System.out.println("test "+ss1);

							// read->edge_number(i);
							int kindx = 0;
							desc = new StringBuffer();
							boolean b = true;
							boolean appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										edge_number[eni] = stringToInt(s1);
										// System.out.println("eni  "+edge_number[eni]+"  i: "+eni);
										kindx++;
										if (kindx == 1) {
											b = false;// exit the loop
										}

									}
									desc = new StringBuffer();// reset
								}
							}
							// read($PHOCSUNIT,*)
							// edge_a(j,i),edge_b(j,i),edge_c(j,i),
							// edge_d(j,i),edge_energies(j,i);
							kindx = 0;
							int j = 0;
							desc = new StringBuffer();
							b = true;
							appendB = false;
							while (b) {
								indx = in.read();
								if (!Character.isWhitespace((char) indx)
										&& ((char) indx != comma)) {
									desc.append((char) indx);
									appendB = true;
								} else {
									if ((indx == -1) || (appendB))// EOF or
																	// normal
																	// exit
									{
										appendB = false;// read data so ready to
														// loop!!!!!
										String s1 = desc.toString();// System.out.println(" ## "+s1);
										if (j == 0) {
											edge_a[kindx][eni] = stringToDouble(s1);
										} else if (j == 1) {
											edge_b[kindx][eni] = stringToDouble(s1);
										} else if (j == 2) {
											edge_c[kindx][eni] = stringToDouble(s1);
										} else if (j == 3) {
											edge_d[kindx][eni] = stringToDouble(s1);
										} else if (j == 4) {
											edge_energies[kindx][eni] = stringToDouble(s1);
										}

										j++;

										if (j == 5) {
											// System.out.println("  i: "+eni+" k  "+kindx+"  a  "+edge_a[kindx][eni]+" b "+edge_b[kindx][eni]+
											// " c "+edge_c[kindx][eni]+" d "+edge_d[kindx][eni]+" en "+edge_energies[kindx][eni]);
											kindx++;
											j = 0;
										}
										if (kindx == edge_number[eni]) {
											b = false;// exit the loop

										}

									}
									desc = new StringBuffer();// reset
								}
							}

							// next set
							eni++;
							if (eni == $MXELEMENT)
								enb = true;
						}

					}// have data
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
		}// try
		catch (Exception exc) {

		}

		seqStr = " Done";// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

	}

	/**
	 * Subroutine to read the pre-calculated q^(2+)-surface, prepare data 
	 * required by the mscat and msdist subroutines, initialize spin effect. Called by HATCH.
	 */
	public static void mscati() {
		// "**********************************************************************"
		// " Subroutine to read the pre-calculated q^(2+)-surface, prepare data
		// " required by the mscat and msdist subroutines, initialize spin
		// effect
		// " corrections
		// " I.Kawrakow, NRC
		// "**********************************************************************"
		// ; Copyright NRC;

		// " moved the following from prestaII_inputs, "
		// " if transport_algorithm = presta-I, exact_bca = .false. and
		// " skin_depth_for_bca <= 1  ==> calculate default presta-I tmin for bca"
		double ecutmn = 0.0;
		double tstbmn = 0.0;
		double tstbm = 0.0;
		double aux = 0.0;
		double ei = 0.0;
		double eil = 0.0;
		int leil = 0;
		double ededx = 0.0;
		double sig = 0.0;
		double sigee = 0.0;
		double sigep = 0.0;
		double rm = 0.0;
		double eip1 = 0.0;
		double eke = 0.0;
		double elke = 0.0;
		int lelke = 0;
		double p2 = 0.0;
		double beta2 = 0.0;
		double chi_a2 = 0.0;
		double dedx0 = 0.0;
		double estepx = 0.0;
		double si = 0.0;
		double ekef = 0.0;
		double sip1 = 0.0;
		double elkef = 0.0;
		int lelkef = 0;
		int leip1l = 0;
		int lelktmp = 0;
		double eip1l = 0.0;
		double ektmp = 0.0;
		double elktmp = 0.0;

		// $LOGICAL ise_monoton, isp_monoton;
		boolean ise_monoton = false;
		boolean isp_monoton = false;
		double sige_old = 0.0;
		double sigp_old = 0.0;

		if (bca_algorithm == 0) {
			exact_bca = true;
		} else {
			exact_bca = false;
		}
		if ((estepe <= 0.0) || (estepe >= 1.0)) {
			estepe = $MAX_ELOSS;
		}
		if ((ximax <= 0.0) || (ximax >= 1.0)) {
			if (exact_bca) {
				ximax = $EXACT_BCA_XIMAX;
			} else {
				ximax = $INEXACT_BCA_XIMAX;
			}
		}

		if (((transport_algorithm != $PRESTA_II) && (transport_algorithm != $PRESTA__I))
				&& (transport_algorithm != $VMC)) {// prevent errors
			transport_algorithm = $PRESTA_II;
		}
		if (skindepth_for_bca <= 1.e-4) {
			// "IF( transport_algorithm = $PRESTA--I & ~exact_bca ) ["
			if (!exact_bca) {

				seqStr = " old PRESTA calculates default min. step-size for BCA: ";// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

				// $set_ecutmn;
				// REPLACE {$set_ecutmn;} WITH {
				ecutmn = 1.e30;
				for (int i = 0; i < $MXREG; i++) {
					if ((MED[i] > 0) && (MED[i] <= NMED)) {
						ecutmn = Math.min(ecutmn, ECUT[i]);
					}
				}
				// };

				seqStr = "minimum ECUT found: " + ecutmn;// +" \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

				tstbmn = 1.e30;
				for (int medium = 0; medium < NMED; medium++) {
					tstbm = (ecutmn - 0.5110034) * (ecutmn + 0.5110034)
							/ (ecutmn * ecutmn);
					tstbm = BLCC[medium] * tstbm * (ecutmn / XCC[medium])
							* (ecutmn / XCC[medium]);
					aux = Math.log(tstbm);
					if (aux > 300.0) {
						seqStr = "aux > 300 ? " + aux;// + " \n";
						if (iprint > 1)
							eq.printSequence(seqStr);
					}
					tstbm = Math.log(tstbm / aux);
					// "Changed the following to the above so that the Intel compiler"
					// " does not vectorize the loop with -xK. Vectorizing this loop"
					// " on an Athlon CPU results in segmentation fault."
					// " IK, Jan 29 2004."
					// "tstbm = Log(tstbm/Log(tstbm));"
					tstbmn = Math.min(tstbmn, tstbm);
				}

				seqStr = "  default BLCMIN is: " + tstbmn;// + " \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

				skindepth_for_bca = Math.exp(tstbmn);

				seqStr = "  this corresponds to " + skindepth_for_bca
						+ " elastic MFPs ";// + " \n";
				if (iprint > 1)
					eq.printSequence(seqStr);

			} else {
				skindepth_for_bca = $SKIN_DEPTH_FOR_BCA;
			}
		}
		// " read MS data for screened Rutherford scattering "
		init_ms_SR();

		for (int medium = 0; medium < NMED; medium++) {
			// "Absorb Euler constant into the multiple scattering parameter"
			// "1.16699413758864573 = Exp[2 EulerGamma - 1]"//see eq 2.14.30
			// SLAC
			BLCC[medium] = 1.16699413758864573 * BLCC[medium];
			// "Take its square as this is employed throughout"
			XCC[medium] = XCC[medium] * XCC[medium];
		}

		if (spin_effects) {
			init_spin();// read binary data spinms.data.
			// init_spin_old();
		}
		// "Determine maximum cross section per energy loss for every medium"
		// write(6,*);
		esige_max = 0.0;
		psige_max = 0.0;
		for (int medium = 1; medium <= NMED; medium++) {

			sigee = 1.E-15;
			sigep = 1.E-15;
			NEKE = MEKE[medium - 1]; // "Number of elements in storage array"
			// ise_monoton = .true.; isp_monoton = .true.;
			ise_monoton = true;
			isp_monoton = true;
			sige_old = -1.0;
			sigp_old = -1.0;
			for (int i = 1; i <= NEKE; i++) {
				Integer intg = new Integer(i);
				double idbl = intg.doubleValue();// exp((float(i) -
													// eke0(medium))/eke1(medium));
				ei = Math.exp((idbl - EKE0[medium - 1]) / EKE1[medium - 1]);
				eil = Math.log(ei);
				leil = i;
				// [{P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM);]}
				// ex: sig =
				// ESIG1[leil-1][medium-1]*eil+ESIG0[leil-1][medium-1];
				// $EVALUATE ededx USING ededx(eil);
				ededx = EDEDX1[leil - 1][medium - 1] * eil
						+ EDEDX0[leil - 1][medium - 1];
				// $EVALUATE sig USING esig(eil);
				sig = ESIG1[leil - 1][medium - 1] * eil
						+ ESIG0[leil - 1][medium - 1];
				sig = sig / ededx;
				if (sig > sigee)
					sigee = sig;
				if (sig < sige_old)
					ise_monoton = false;
				sige_old = sig;

				// $EVALUATE ededx USING pdedx(eil);
				ededx = PDEDX1[leil - 1][medium - 1] * eil
						+ PDEDX0[leil - 1][medium - 1];
				// $EVALUATE sig USING psig(eil);
				sig = PSIG1[leil - 1][medium - 1] * eil
						+ PSIG0[leil - 1][medium - 1];
				sig = sig / ededx;
				if (sig > sigep)
					sigep = sig;
				if (sig < sigp_old)
					isp_monoton = false;
				sigp_old = sig;

			}

			seqStr = " Medium " + medium + " sige = " + sigee + "    " + sigep;// +
																				// " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			sig_ismonotone[0][medium - 1] = ise_monoton;
			sig_ismonotone[1][medium - 1] = isp_monoton;

			esig_e[medium - 1] = sigee;
			psig_e[medium - 1] = sigep;
			if (sigee > esige_max)
				esige_max = sigee;
			if (sigep > psige_max)
				psige_max = sigep;
		}

		seqStr = " Initializing tmxs for estepe = " + estepe + " and ximax = "
				+ ximax;// + " \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// "Determine upper limit in step size for multiple scattering
		rm = 0.5110034;
		for (int medium = 1; medium <= NMED; medium++) {
			// " Calculate range array first "
			// " =========================== "
			ei = Math.exp((1.0 - EKE0[medium - 1]) / EKE1[medium - 1]); // "Energy of first table entry"
			eil = Math.log(ei);
			leil = 1;
			E_array[0][medium - 1] = ei;// E_array[1][medium-1] = ei;
			expeke1[medium - 1] = Math.exp(1. / EKE1[medium - 1]) - 1.0;
			// range_ep[0][1,medium) = 0; range_ep(1,1,medium) = 0;
			range_ep[0][0][medium - 1] = 0.0;
			range_ep[1][0][medium - 1] = 0.0;
			NEKE = MEKE[medium - 1]; // "Number of elements in storage array"
			for (int i = 1; i <= NEKE - 1; i++) {
				Integer intg = new Integer(i + 1);
				double idbl = intg.doubleValue();
				eip1 = Math.exp((idbl - EKE0[medium - 1]) / EKE1[medium - 1]); // "Energy at i+1"
				E_array[i][medium - 1] = eip1;// E_array(i+1,medium) = eip1;
				// " Calculate range. The following expressions result from the"
				// " logarithmic interpolation for the (restricted) stopping power "
				// " and a power power series expansion of the integral "
				eke = 0.5 * (eip1 + ei);
				elke = Math.log(eke);
				// $SET INTERVAL elke,eke;
				// else-> L{P1}={P2}1(MEDIUM)*{P1}+{P2}0(MEDIUM);]}
				// SPECIFY SNAME AS
				// ['sinc'|'blc'|'rthr'|'rthri'|'SINC'|'BLC'|'RTHR'|'RTHRI'];
				// SPECIFY SNAME1 AS ['sin'|'SIN'];
				// eg: $SET INTERVAL THETA,SINC;
				// eg Double dth=new
				// Double(SINC1*THETA+SINC0);//System.out.println(dth);
				// eg LTHETA=dth.intValue();
				Double dbl = new Double(EKE1[medium - 1] * elke
						+ EKE0[medium - 1]);
				lelke = dbl.intValue();
				// ------------------------------
				// {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM);]}
				// eg $EVALUATE ededx USING pdedx(eil);
				// eg ededx =
				// PDEDX1[leil-1][medium-1]*eil+PDEDX0[leil-1][medium-1];
				// $EVALUATE ededx USING pdedx(elke);
				ededx = PDEDX1[lelke - 1][medium - 1] * elke
						+ PDEDX0[lelke - 1][medium - 1];
				aux = PDEDX1[i - 1][medium - 1] / ededx;// aux =
														// PDEDX1(i,medium)/ededx;
				// range_ep(1,i+1,medium) = range_ep(1,i,medium) +
				range_ep[1][i][medium - 1] = range_ep[1][i - 1][medium - 1]
						+ (eip1 - ei)
						/ ededx
						* (1.0 + aux * (1.0 + 2.0 * aux) * ((eip1 - ei) / eke)
								* ((eip1 - ei) / eke) / 24.0);
				// $EVALUATE ededx USING ededx(elke);
				ededx = EDEDX1[lelke - 1][medium - 1] * elke
						+ EDEDX0[lelke - 1][medium - 1];
				aux = EDEDX1[i - 1][medium - 1] / ededx;// aux =
														// ededx1(i,medium)/ededx;
				// range_ep(0,i+1,medium) = range_ep(0,i,medium) +
				range_ep[0][i][medium - 1] = range_ep[0][i - 1][medium - 1]
						+ (eip1 - ei)
						/ ededx
						* (1.0 + aux * (1.0 + 2.0 * aux) * ((eip1 - ei) / eke)
								* ((eip1 - ei) / eke) / 24.0);
				ei = eip1;
			}

			// " Now tmxs "
			// " ======== "
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			ei = Math.exp(eil);
			leil = 1; // " As in $SET INTERVAL but avoids roundoff"
			p2 = ei * (ei + 2.0 * rm);
			beta2 = p2 / (p2 + rm * rm);
			chi_a2 = XCC[medium - 1] / (4.0 * p2 * BLCC[medium - 1]);
			// $EVALUATE dedx0 USING ededx(eil);
			dedx0 = EDEDX1[leil - 1][medium - 1] * eil
					+ EDEDX0[leil - 1][medium - 1];
			estepx = 2.0 * p2 * beta2 * dedx0 / ei / XCC[medium - 1]
					/ (Math.log(1.0 + 1. / chi_a2) * (1.0 + chi_a2) - 1.0);
			estepx = estepx * ximax;
			if (estepx > estepe) {
				estepx = estepe;
			}
			si = estepx * ei / dedx0;
			// "write(1,'(3g15.6)') ei,si,range_ep(0,1,medium);"
			for (int i = 1; i <= NEKE - 1; i++)// DO i = 1,neke - 1
			{

				elke = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				eke = Math.exp(elke);
				lelke = i + 1;
				p2 = eke * (eke + 2.0 * rm);
				beta2 = p2 / (p2 + rm * rm);
				chi_a2 = XCC[medium - 1] / (4.0 * p2 * BLCC[medium - 1]);
				// $EVALUATE ededx USING ededx(elke);
				ededx = EDEDX1[lelke - 1][medium - 1] * elke
						+ EDEDX0[lelke - 1][medium - 1];
				estepx = 2.0 * p2 * beta2 * ededx / eke / XCC[medium - 1]
						/ (Math.log(1.0 + 1. / chi_a2) * (1.0 + chi_a2) - 1.0);
				estepx = estepx * ximax;
				if (estepx > estepe) {
					estepx = estepe;
				}
				ekef = (1.0 - estepx) * eke;
				// if( ekef <= E_array(1,medium) )
				if (ekef <= E_array[0][medium - 1]) {
					// sip1 = (E_array(1,medium) - ekef)/dedx0;
					sip1 = (E_array[0][medium - 1] - ekef) / dedx0;
					ekef = E_array[0][medium - 1];// ekef = E_array(1,medium);
					elkef = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
					lelkef = 1;
				} else {
					elkef = Math.log(ekef);
					// $SET INTERVAL elkef,eke;
					Double dbl = new Double(EKE1[medium - 1] * elkef
							+ EKE0[medium - 1]);
					lelkef = dbl.intValue();

					leip1l = lelkef + 1;
					eip1l = (leip1l - EKE0[medium - 1]) / EKE1[medium - 1];
					eip1 = E_array[leip1l - 1][medium - 1];// eip1 =
															// E_array(leip1l,medium);
					aux = (eip1 - ekef) / eip1;
					elktmp = 0.5 * (elkef + eip1l + 0.25 * aux * aux
							* (1.0 + aux * (1.0 + 0.875 * aux)));
					ektmp = 0.5 * (ekef + eip1);
					lelktmp = lelkef;
					// $EVALUATE ededx USING ededx(elktmp);
					ededx = EDEDX1[lelktmp - 1][medium - 1] * elktmp
							+ EDEDX0[lelktmp - 1][medium - 1];
					aux = EDEDX1[lelktmp - 1][medium - 1] / ededx;// aux =
																	// ededx1(lelktmp,medium)/ededx;
					sip1 = (eip1 - ekef)
							/ ededx
							* (1.0 + aux * (1.0 + 2.0 * aux)
									* ((eip1 - ekef) / ektmp)
									* ((eip1 - ekef) / ektmp) / 24.0);
				}
				// sip1 = sip1 + range_ep(0,lelke,medium) -
				// range_ep(0,lelkef+1,medium);
				sip1 = sip1 + range_ep[0][lelke - 1][medium - 1]
						- range_ep[0][lelkef][medium - 1];

				// "Now solve these equations"
				// "  si   = tmxs1 * eil   + tmxs0"
				// "  sip1 = tmxs1 * eip1l + tmxs0"

				// TMXS1(i,medium) = (sip1 - si)*EKE1[medium-1];
				// TMXS0(i,medium) = sip1 - TMXS1(i,medium)*elke;
				TMXS1[i - 1][medium - 1] = (sip1 - si) * EKE1[medium - 1];
				TMXS0[i - 1][medium - 1] = sip1 - TMXS1[i - 1][medium - 1]
						* elke;

				// "write(1,'(3g15.6)') eke,sip1,range_ep(0,i+1,medium);"
				si = sip1;
			}
			// "Now pick up last table entry which applies only to last energy"
			TMXS0[NEKE - 1][medium - 1] = TMXS0[NEKE - 2][medium - 1];
			TMXS1[NEKE - 1][medium - 1] = TMXS1[NEKE - 2][medium - 1];
		}

	}

	@Deprecated
	public static void init_spin_old() {
		// "======================================================================="
		// " Reads in spin rejection data for multiple elastic scattering and      "
		// " initializes interpolation arrays for the screening parameter,         "
		// " elastic cross section, first and second MS moments                    "
		// " This version uses the old EGSnrc ASCII spin data format.              "
		// " If the installation failed to figure out the byte order of your       "
		// " machine, rename this subroutine to init_spin, the current init_spin   "
		// " to init_spin_new (or whatever) and download the old spinms data       "
		// " directory.
		// " I. Kawrakow, NRC                                                      "
		// "======================================================================="
		// ; Copyright NRC;
		double[][] eta_array = new double[2][$MAXE_SPI1 + 1];
		double[][] c_array = new double[2][$MAXE_SPI1 + 1];
		double[][] g_array = new double[2][$MAXE_SPI1 + 1];
		double[] earray = new double[$MAXE_SPI1 + 1];
		double[] tmp_array = new double[$MAXE_SPI1 + 1];
		double sum_Z2 = 0.0;
		double sum_Z = 0.0;
		double sum_A = 0.0;
		double sum_pz = 0.0;
		double Z = 0.0;
		double tmp = 0.0;
		double Z23 = 0.0;
		double g_m = 0.0;
		double g_r = 0.0;
		double sig = 0.0;
		double dedx = 0.0;
		double dum1 = 0.0;
		double dum2 = 0.0;
		double dum3 = 0.0;
		double aux_o = 0.0;
		double tau = 0.0;
		double tauc = 0.0;
		double beta2 = 0.0;
		double eta = 0.0;
		double gamma = 0.0;
		double fmax = 0.0;
		double eil = 0.0;
		double e = 0.0;
		double si1e = 0.0;
		double si2e = 0.0;
		double si1p = 0.0;
		double si2p = 0.0;
		double aae = 0.0;
		double etap = 0.0;
		double[] elarray = new double[$MAXE_SPI1 + 1];
		double[] farray = new double[$MAXE_SPI1 + 1];
		double[] af = new double[$MAXE_SPI1 + 1];
		double[] bf = new double[$MAXE_SPI1 + 1];
		double[] cf = new double[$MAXE_SPI1 + 1];
		double[] df = new double[$MAXE_SPI1 + 1];
		// double spline = 0.0;
		double fine = 137.03604;
		double TF_constant = 0.88534138;
		String endzs = "";
		int Zi = 0;
		// public static double[][][][][] spin_rej=
		// new double[$MXMED][2][$MAXE_SPI1+1][$MAXQ_SPIN+1][$MAXU_SPIN+1];
		int n_ener = 0;
		int n_q = 0;
		int n_point = 0;
		int je = 0;
		int ndata = 0;
		int leil = 0;
		// " First construct the path to the spin dbase directory "
		String spin_file = datas + file_sep + egsDatas + file_sep + spinDir
				+ file_sep + zfile;
		int iread = 0;
		int lnr = 0;// data number
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		String equals = "=";
		char comma = ',';

		int iqq = 0;

		for (int medium = 1; medium <= NMED; medium++) {

			seqStr = "  Initializing spin data for medium " + medium;// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			// spin_file = datas+file_sep+spinDir+file_sep+zfile;
			for (int iq = 0; iq <= 1; iq++) {
				for (int i = 0; i <= $MAXE_SPI1; i++) {
					eta_array[iq][i] = 0.0;
					c_array[iq][i] = 0.0;
					g_array[iq][i] = 0.0;
					for (int j = 0; j <= $MAXQ_SPIN; j++) {
						for (int k = 0; k <= $MAXU_SPIN; k++) {
							spin_rej[medium - 1][iq][i][j][k] = 0.0;
						}// k

					}// j
				}// i
			}// iq
			sum_Z2 = 0.0;
			sum_A = 0.0;
			sum_pz = 0.0;
			sum_Z = 0.0;
			for (int i_ele = 1; i_ele <= NNE[medium - 1]; i_ele++) {
				spin_file = datas + file_sep + egsDatas + file_sep + spinDir
						+ file_sep + zfile;
				Z = ZELEM[medium - 1][i_ele - 1];
				Double Zdbl = new Double(Z);
				Zi = Zdbl.intValue();
				tmp = PZ[medium - 1][i_ele - 1] * Z * (Z + 1);
				// "For now, we take into account the contribution of atomic"
				// "electrons to elastic scattering by replacing Z**2 with  "
				// "Z*(Z+1). The part of the scattering power that is taken "
				// "into account by discrete Moller/Bhabha events is        "
				// "substracted below => bc is energy dependent. We will    "
				// "worry about better approaches in the future (a realistic"
				// "inelastic scattering model is needed first)             "
				if (Zi < 10)
					endzs = z00file + intToString(Zi);
				else if (Zi < 100)
					endzs = z0file + intToString(Zi);
				else
					endzs = intToString(Zi);
				spin_file = spin_file + endzs;
				iread = 0;
				lnr = 0;// reset
				iqq = 0;
				boolean rectb = true;
				boolean startb = true;// reset
				int indx = 0;
				int SPI1i = 0;
				int QSPINi = 0;
				int USPINi = 0;
				int tmpi = 0;
				haveData = false;
				desc = new StringBuffer();// reset
				// System.out.println(spin_file+"    "+$MAXE_SPI1);
				// open(61,file=spin_file,status='old',err=:SPIN-DBASE-ERROR:);
				// **********************READ and set DATA
				try {
					FileInputStream in = new FileInputStream(spin_file);

					while ((iread = in.read()) != -1) {
						if (!Character.isWhitespace((char) iread)) {
							desc.append((char) iread);
							haveData = true;
						} else {
							if (haveData)// we have data
							{
								haveData = false;// reset
								lnr++;
								if (lnr == 1) {
									// read(61,*)
									// espin_min,espin_max,b2spin_min,b2spin_max;
									String s = desc.toString();
									espin_min = stringToDouble(s);
									// System.out.println(espin_min);
								}
								if (lnr == 2) {
									// read(61,*)
									// espin_min,espin_max,b2spin_min,b2spin_max;
									String s = desc.toString();
									espin_max = stringToDouble(s);
									// System.out.println(espin_max);
								}
								if (lnr == 3) {
									// read(61,*)
									// espin_min,espin_max,b2spin_min,b2spin_max;
									String s = desc.toString();
									b2spin_min = stringToDouble(s);
									// System.out.println(b2spin_min);
								}
								if (lnr == 4) {
									// read(61,*)
									// espin_min,espin_max,b2spin_min,b2spin_max;
									String s = desc.toString();
									b2spin_max = stringToDouble(s);
									// System.out.println(b2spin_max);
								}
								if (lnr == 5) {
									// read(61,*) n_ener,n_q,n_point;
									String s = desc.toString();
									n_ener = stringToInt(s);
									// System.out.println(n_ener);
								}
								if (lnr == 6) {
									// read(61,*) n_ener,n_q,n_point;
									String s = desc.toString();
									n_q = stringToInt(s);
									// System.out.println(n_q);
								}
								if (lnr == 7) {
									// read(61,*) n_ener,n_q,n_point;
									String s = desc.toString();
									n_point = stringToInt(s);
									// --------------------
									if (((n_ener != $MAXE_SPIN) || (n_q != $MAXQ_SPIN))
											|| (n_point != $MAXU_SPIN)) {
										STOPPROGRAM = true;
										seqStr = " Wrong spin file for Z = "
												+ Zi;
										// if(iprint>1)
										eq.printSequence(seqStr);

										return;
									}
									sum_Z2 = sum_Z2 + tmp;
									sum_Z = sum_Z + PZ[medium - 1][i_ele - 1]
											* Z;
									sum_A = sum_A + PZ[medium - 1][i_ele - 1]
											* WA[medium - 1][i_ele - 1];
									sum_pz = sum_pz + PZ[medium - 1][i_ele - 1];
									Z23 = Math.pow(Z, 0.6666667);
									// System.out.println(n_point);
								}
								if (lnr > 7) {
									if (startb)// electron or positron is next
									{
										iqq++;// 1 for e and 2 for positron
										if (iqq == 2) {
											earray = new double[$MAXE_SPI1 + 1];// in
																				// fact
																				// is
																				// the
																				// same!!!!
											SPI1i = 0;
										}
										startb = false;
										// auto second # is readed so
										rectb = false;
										SPI1i++;
										// String s=desc.toString();
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													// System.out.println("test"+s1);
													if (s1.compareTo(equals) == 0)
														b = false;// exit the
																	// loop
												}
												desc = new StringBuffer();// reset
											}
										}
										// read e array
										desc = new StringBuffer();
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													earray[SPI1i - 1] = stringToDouble(s1);
													// System.out.println("e  "+earray[SPI1i-1]+"  i: "+SPI1i);
													b = false;// exit the loop
												}
											}
										}
										// read dum1
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum1 = stringToDouble(s1);
													// System.out.println("d1  "+dum1);
													b = false;// exit the loop
												}
											}
										}
										// read dum2
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum2 = stringToDouble(s1);
													// System.out.println("d2  "+dum2);
													b = false;// exit the loop
												}
											}
										}
										// read dum3
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum3 = stringToDouble(s1);
													// System.out.println("d3  "+dum3);
													b = false;// exit the loop
												}
											}
										}
										// read aux_o
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													aux_o = stringToDouble(s1);
													// System.out.println("aux_o  "+aux_o);
													b = false;// exit the loop
												}
											}
										}
										// ********************
										eta_array[iqq - 1][SPI1i - 1] = eta_array[iqq - 1][SPI1i - 1]
												+ tmp * Math.log(Z23 * aux_o);
										tau = earray[SPI1i - 1] / 511.0034; // "energy in the file is in keV"
										beta2 = tau * (tau + 2.0)
												/ ((tau + 1.0) * (tau + 1.0));
										eta = Z23
												/ ((fine * TF_constant) * (fine * TF_constant))
												* aux_o / 4.0 / tau
												/ (tau + 2.0);
										c_array[iqq - 1][SPI1i - 1] = c_array[iqq - 1][SPI1i - 1]
												+ tmp
												* (Math.log(1.0 + 1.0 / eta) - 1.0 / (1.0 + eta))
												* dum1 * dum3;
										g_array[iqq - 1][SPI1i - 1] = g_array[iqq - 1][SPI1i - 1]
												+ tmp * dum2;

									}// startb
									else if (rectb)// read line
									{
										rectb = false;
										SPI1i++;
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													// System.out.println("test"+s1);
													if (s1.compareTo(equals) == 0)
														b = false;// exit the
																	// loop
												}
												desc = new StringBuffer();// reset
											}
										}
										// read e array
										desc = new StringBuffer();
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													earray[SPI1i - 1] = stringToDouble(s1);
													// System.out.println("e  "+earray[SPI1i-1]+"  i: "+SPI1i);
													b = false;// exit the loop
												}
											}
										}
										// read dum1
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum1 = stringToDouble(s1);
													// System.out.println("d1  "+dum1);
													b = false;// exit the loop
												}
											}
										}
										// read dum2
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum2 = stringToDouble(s1);
													// System.out.println("d2  "+dum2);
													b = false;// exit the loop
												}
											}
										}
										// read dum3
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													dum3 = stringToDouble(s1);
													// System.out.println("d3  "+dum3);
													b = false;// exit the loop
												}
											}
										}
										// read aux_o
										desc = new StringBuffer();// reset
										b = true;
										appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();
													aux_o = stringToDouble(s1);
													// System.out.println("aux_o  "+aux_o);
													b = false;// exit the loop
												}
											}
										}

										// ********************
										eta_array[iqq - 1][SPI1i - 1] = eta_array[iqq - 1][SPI1i - 1]
												+ tmp * Math.log(Z23 * aux_o);
										tau = earray[SPI1i - 1] / 511.0034; // "energy in the file is in keV"
										beta2 = tau * (tau + 2.0)
												/ ((tau + 1.0) * (tau + 1.0));
										eta = Z23
												/ ((fine * TF_constant) * (fine * TF_constant))
												* aux_o / 4.0 / tau
												/ (tau + 2.0);
										c_array[iqq - 1][SPI1i - 1] = c_array[iqq - 1][SPI1i - 1]
												+ tmp
												* (Math.log(1.0 + 1.0 / eta) - 1.0 / (1.0 + eta))
												* dum1 * dum3;
										g_array[iqq - 1][SPI1i - 1] = g_array[iqq - 1][SPI1i - 1]
												+ tmp * dum2;

									} else// read data
									{
										// System.out.println("enter");
										tmpi++;
										USPINi++;
										// QSPINi=maxim 15
										if (QSPINi <= $MAXQ_SPIN) {
											String s = desc.toString();
											tmp_array[tmpi - 1] = stringToDouble(s);
											if (tmpi == $MAXE_SPI1 + 1) {
												// ----------------
												for (int k = 0; k < USPINi; k++) {
													spin_rej[medium - 1][iqq - 1][SPI1i - 1][QSPINi][k] = spin_rej[medium - 1][iqq - 1][SPI1i - 1][QSPINi][k]
															+ tmp
															* tmp_array[k];
													// System.out.println("med1 "+medium+" iq1 "+iqq+" i1 "+SPI1i+" j "+QSPINi+" k "+k+" s= "+tmp_array[k]);
												}
												QSPINi++;
												// ------------------
												tmpi = 0;
												USPINi = 0;
												tmp_array = new double[$MAXE_SPI1 + 1];
											}

										}
										// reset
										if (QSPINi == $MAXQ_SPIN + 1) {
											QSPINi = 0;
											rectb = true;
											if (SPI1i == $MAXE_SPI1 + 1)// 32
											{
												startb = true;
											}
										}
									}

								}// lnr>7

							}// have data
							desc = new StringBuffer();
						}// else
					}// main while
					in.close();
				}// try
				catch (Exception exc) {

				}

			}// for(int i_ele=1;i_ele<=NNE[medium-1];i_ele++)
				// " spin_rej will be used as a rejection function in MS sampling, "
				// " so scale maximum to unity"
			for (int iq = 0; iq <= 1; iq++) {
				for (int i = 0; i <= $MAXE_SPI1; i++) {
					for (int j = 0; j <= $MAXQ_SPIN; j++) {
						fmax = 0.0;
						for (int k = 0; k <= $MAXU_SPIN; k++) {
							if (spin_rej[medium - 1][iq][i][j][k] > fmax) {
								fmax = spin_rej[medium - 1][iq][i][j][k];
							}
						}
						for (int k = 0; k <= $MAXU_SPIN; k++) {
							spin_rej[medium - 1][iq][i][j][k] = spin_rej[medium - 1][iq][i][j][k]
									/ fmax;
						}

					}

				}

			}
			// " Process eta_array, c_array and g_array to their final form "
			for (int i = 0; i <= $MAXE_SPI1; i++) {
				tau = earray[i] / 511.0034;
				beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
				for (int iq = 0; iq <= 1; iq++) {
					aux_o = Math.exp(eta_array[iq][i] / sum_Z2)
							/ ((fine * TF_constant) * (fine * TF_constant));
					eta_array[iq][i] = 0.26112447 * aux_o * BLCC[medium - 1]
							/ XCC[medium - 1];
					eta = aux_o / 4.0 / tau / (tau + 2.0);
					gamma = 3.0
							* (1.0 + eta)
							* (Math.log(1.0 + 1.0 / eta) * (1.0 + 2.0 * eta) - 2.0)
							/ (Math.log(1.0 + 1.0 / eta) * (1.0 + eta) - 1.0);
					g_array[iq][i] = g_array[iq][i] / sum_Z2 / gamma;
					c_array[iq][i] = c_array[iq][i] / sum_Z2
							/ (Math.log(1.0 + 1.0 / eta) - 1.0 / (1.0 + eta));
				}

			}
			// " Convert to MeV and set interpolation interavals"
			espin_min = espin_min / 1000.0;
			espin_max = espin_max / 1000.0;
			dlener = Math.log(espin_max / espin_min) / $MAXE_SPIN;
			dleneri = 1.0 / dlener;
			espml = Math.log(espin_min);
			dbeta2 = (b2spin_max - b2spin_min) / $MAXE_SPIN;
			dbeta2i = 1.0 / dbeta2;
			dqq1 = 0.5 / $MAXQ_SPIN;
			dqq1i = 1.0 / dqq1;

			// " Prepare interpolation table for the screening parameter "
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			e = Math.exp(eil);
			if (e <= espin_min) {
				// System.out.println("here");
				si1e = eta_array[0][0];
				si1p = eta_array[1][0];
			} else {
				if (e <= espin_max) {
					aae = (eil - espml) * dleneri;
					Double dbl = new Double(aae);
					je = dbl.intValue();
					// System.out.println("aae "+aae+" je "+je);
					aae = aae - je;
				} else {
					tau = e / 0.5110034;
					beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
					aae = (beta2 - b2spin_min) * dbeta2i;
					Double dbl = new Double(aae);
					je = dbl.intValue();
					// System.out.println("aae "+aae+" je "+je);
					aae = aae - je;
					je = je + $MAXE_SPIN + 1;
				}
				si1e = (1 - aae) * eta_array[0][je] + aae
						* eta_array[0][je + 1];
				si1p = (1 - aae) * eta_array[1][je] + aae
						* eta_array[1][je + 1];
			}
			NEKE = MEKE[medium - 1];
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				e = Math.exp(eil);
				if (e <= espin_min) {
					si2e = eta_array[0][0];
					si2p = eta_array[1][0];
				} else {
					if (e <= espin_max) {
						aae = (eil - espml) * dleneri;
						Double dbl = new Double(aae);
						je = dbl.intValue();
						aae = aae - je;
					} else {
						tau = e / 0.5110034;
						beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
						aae = (beta2 - b2spin_min) * dbeta2i;
						Double dbl = new Double(aae);
						je = dbl.intValue();
						aae = aae - je;
						je = je + $MAXE_SPIN + 1;
					}
					si2e = (1 - aae) * eta_array[0][je] + aae
							* eta_array[0][je + 1];
					si2p = (1 - aae) * eta_array[1][je] + aae
							* eta_array[1][je + 1];
				}
				etae_ms1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				etae_ms0[i - 1][medium - 1] = si2e
						- etae_ms1[i - 1][medium - 1] * eil;
				etap_ms1[i - 1][medium - 1] = (si2p - si1p) * EKE1[medium - 1];
				etap_ms0[i - 1][medium - 1] = si2p
						- etap_ms1[i - 1][medium - 1] * eil;
				// "write(6,'(6g13.5)') e,si2e,si2p,"
				// "  etae_ms1(i,medium),etae_ms0(i,medium),"
				// "  etap_ms1(i,medium),etap_ms0(i,medium);"
				si1e = si2e;
				si1p = si2p;
			}
			etae_ms1[NEKE - 1][medium - 1] = etae_ms1[NEKE - 2][medium - 1];
			etae_ms0[NEKE - 1][medium - 1] = etae_ms0[NEKE - 2][medium - 1];
			etap_ms1[NEKE - 1][medium - 1] = etap_ms1[NEKE - 2][medium - 1];
			etap_ms0[NEKE - 1][medium - 1] = etap_ms0[NEKE - 2][medium - 1];

			// "Prepare correction to the first MS moment due to spin effects"
			// "first electrons"
			for (int i = 0; i <= $MAXE_SPIN; i++) {
				elarray[i] = Math.log(earray[i] / 1000.0);
				farray[i] = c_array[0][i];
			}
			for (int i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
				elarray[i] = Math.log(earray[i + 1] / 1000.0);
				farray[i] = c_array[0][i + 1];
			}
			ndata = $MAXE_SPI1 + 1;
			if (UE[medium - 1] > 1.e5) {
				elarray[ndata - 1] = Math.log(UE[medium - 1]);
			} else {
				elarray[ndata - 1] = Math.log(1.e5);
			}
			farray[ndata - 1] = 1.0;
			set_spline(elarray, farray, af, bf, cf, df, ndata);
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			si1e = spline(eil, elarray, af, bf, cf, df, ndata);
			// si1e = spline(eil,ndata);
			// "write(6,*) ' Interpolation table for q1 correction (e-) ';"
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				si2e = spline(eil, elarray, af, bf, cf, df, ndata);
				// si2e = spline(eil,ndata);
				q1ce_ms1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				q1ce_ms0[i - 1][medium - 1] = si2e
						- q1ce_ms1[i - 1][medium - 1] * eil;
				// "write(6,'(4g13.5)')"
				// "  Exp(eil),si2e,q1ce_ms1(i,medium),q1ce_ms0(i,medium);"
				si1e = si2e;
			}
			// "write(6,*);"
			q1ce_ms1[NEKE - 1][medium - 1] = q1ce_ms1[NEKE - 2][medium - 1];
			q1ce_ms0[NEKE - 1][medium - 1] = q1ce_ms0[NEKE - 2][medium - 1];
			// "now positrons"
			for (int i = 0; i <= $MAXE_SPIN; i++) {
				farray[i] = c_array[1][i];
			}
			for (int i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
				farray[i] = c_array[1][i + 1];
			}
			set_spline(elarray, farray, af, bf, cf, df, ndata);
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			si1e = spline(eil, elarray, af, bf, cf, df, ndata);
			// si1e = spline(eil,ndata);
			// "write(6,*) ' Interpolation table for q1 correction (e+) ';"
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				si2e = spline(eil, elarray, af, bf, cf, df, ndata);
				// si2e = spline(eil,ndata);
				q1cp_ms1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				q1cp_ms0[i - 1][medium - 1] = si2e
						- q1cp_ms1[i - 1][medium - 1] * eil;
				// "write(6,'(4g13.5)')"
				// "  Exp(eil),si2e,q1cp_ms1(i,medium),q1cp_ms0(i,medium);"
				si1e = si2e;
			}
			// "write(6,*);"
			q1cp_ms1[NEKE - 1][medium - 1] = q1cp_ms1[NEKE - 2][medium - 1];
			q1cp_ms0[NEKE - 1][medium - 1] = q1cp_ms0[NEKE - 2][medium - 1];
			// "prepare interpolation table for the second MS moment correction"
			// "e-"
			for (int i = 0; i <= $MAXE_SPIN; i++) {
				farray[i] = g_array[0][i];
			}
			for (int i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
				farray[i] = g_array[0][i + 1];
			}
			set_spline(elarray, farray, af, bf, cf, df, ndata);
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			si1e = spline(eil, elarray, af, bf, cf, df, ndata);
			// si1e = spline(eil,ndata);
			// "write(6,*) ' Interpolation table for q2 correction (e-) ';"
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				si2e = spline(eil, elarray, af, bf, cf, df, ndata);
				// si2e = spline(eil,ndata);
				q2ce_ms1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				q2ce_ms0[i - 1][medium - 1] = si2e
						- q2ce_ms1[i - 1][medium - 1] * eil;
				// "write(6,'(4g13.5)')"
				// "  Exp(eil),si2e,q2ce_ms1(i,medium),q2ce_ms0(i,medium);"
				si1e = si2e;
			}
			// "write(6,*);"
			q2ce_ms1[NEKE - 1][medium - 1] = q2ce_ms1[NEKE - 2][medium - 1];
			q2ce_ms0[NEKE - 1][medium - 1] = q2ce_ms0[NEKE - 2][medium - 1];
			// "e+"
			for (int i = 0; i <= $MAXE_SPIN; i++) {
				farray[i] = g_array[1][i];
			}
			for (int i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
				farray[i] = g_array[1][i + 1];
			}
			set_spline(elarray, farray, af, bf, cf, df, ndata);
			eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
			si1e = spline(eil, elarray, af, bf, cf, df, ndata);
			// si1e = spline(eil,ndata);
			// "write(6,*) ' Interpolation table for q2 correction (e+) ';"
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				si2e = spline(eil, elarray, af, bf, cf, df, ndata);
				// si2e = spline(eil,ndata);
				q2cp_ms1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				q2cp_ms0[i - 1][medium - 1] = si2e
						- q2cp_ms1[i - 1][medium - 1] * eil;
				// "write(6,'(4g13.5)')"
				// "  Exp(eil),si2e,q2cp_ms1(i,medium),q2cp_ms0(i,medium);"
				si1e = si2e;
			}
			q2cp_ms1[NEKE - 1][medium - 1] = q2cp_ms1[NEKE - 2][medium - 1];
			q2cp_ms0[NEKE - 1][medium - 1] = q2cp_ms0[NEKE - 2][medium - 1];

			// "Now substract scattering power that is already taken into account in"
			// "discrete Moller/Bhabha events"
			tauc = TE[medium - 1] / 0.5110034;
			si1e = 1;
			for (int i = 1; i <= NEKE - 1; i++) {
				eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
				e = Math.exp(eil);
				leil = i + 1;
				tau = e / 0.5110034;
				if (tau > 2.0 * tauc) {
					// $EVALUATE sig USING esig(eil);
					// [{P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM);]}
					sig = ESIG1[leil - 1][medium - 1] * eil
							+ ESIG0[leil - 1][medium - 1];
					// $EVALUATE dedx USING ededx(eil);
					dedx = EDEDX1[leil - 1][medium - 1] * eil
							+ EDEDX0[leil - 1][medium - 1];
					sig = sig / dedx;
					if (sig > 1.e-6) {// "To be sure that this is not a CSDA calc."
										// $EVALUATE etap USING etae_ms(eil);
						etap = etae_ms1[leil - 1][medium - 1] * eil
								+ etae_ms0[leil - 1][medium - 1];
						eta = 0.25 * etap * XCC[medium - 1] / BLCC[medium - 1]
								/ tau / (tau + 2.0);
						g_r = (1.0 + 2.0 * eta) * Math.log(1.0 + 1.0 / eta)
								- 2.0;
						g_m = Math.log(0.5 * tau / tauc)
								+ (1.0 + ((tau + 2.0) / (tau + 1.0))
										* ((tau + 2.0) / (tau + 1.0)))
								* Math.log(2.0 * (tau - tauc + 2.0)
										/ (tau + 4.0))
								- 0.25
								* (tau + 2.0)
								* (tau + 2.0 + 2.0 * (2.0 * tau + 1.0)
										/ ((tau + 1.0) * (tau + 1.0)))
								* Math.log((tau + 4.0) * (tau - tauc) / tau
										/ (tau - tauc + 2.0))
								+ 0.5
								* (tau - 2.0 * tauc)
								* (tau + 2.0)
								* (1.0 / (tau - tauc) - 1.0 / ((tau + 1.0) * (tau + 1.0)));
						if (g_m < g_r) {
							g_m = g_m / g_r;
						} else {
							g_m = 1.0;
						}
						si2e = 1.0 - g_m * sum_Z / sum_Z2;
					} else {
						si2e = 1.0;
					}
				} else {
					si2e = 1.0;
				}
				// "write(6,'(2g14.6)') e,si2e;"
				// "si2e=1;"
				blcce1[i - 1][medium - 1] = (si2e - si1e) * EKE1[medium - 1];
				blcce0[i - 1][medium - 1] = si2e - blcce1[i - 1][medium - 1]
						* eil;
				si1e = si2e;
			}
			blcce1[NEKE - 1][medium - 1] = blcce1[NEKE - 2][medium - 1];
			blcce0[NEKE - 1][medium - 1] = blcce0[NEKE - 2][medium - 1];
			// System.out.println(blcce0[NEKE-1][medium-1]);
			// "We will not bother to do the same for positrons at this time"

			seqStr = " done";// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

		}// for (int medium=1;medium<=NMED;medium++)

	}

	/**
	 * Called by mscati. Reads in pre-calculated screened Rutherford multiple elastic 
	 * scattering data.
	 */
	private static void init_ms_SR() {
		// "================================================================"
		// " Reads in pre-calculated screened Rutherford multiple elastic   "
		// " scattering data                                                "
		// " I. Kawrakow, NRC                                               "
		// "================================================================"
		// ; Copyright NRC;
		// public static double [][][] ums_array=new
		// double[$MAXL_MS+1][$MAXQ_MS+1][$MAXU_MS+1];
		// $MAXU_MS=31; $MAXL_MS=63;$MAXQ_MS=7;
		seqStr = "  Reading screened Rutherford MS data ............... ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		String filename = datas + file_sep + egsDatas + file_sep + msfile
				+ defaultext;
		int i = 0;
		StringBuffer desc = new StringBuffer();
		// String descs = null;
		boolean haveData = false;

		int uki = 0;
		int uii = 0;
		int uji = 0;
		boolean ukb = false;
		int fki = 0;
		int fii = 0;
		int fji = 0;
		boolean fkb = false;
		int wki = 0;
		int wii = 0;
		int wji = 0;
		boolean wkb = false;
		int iki = 0;
		int iii = 0;
		int iji = 0;
		boolean ikb = false;
		try {
			FileInputStream in = new FileInputStream(filename);

			while ((i = in.read()) != -1) {
				if (!Character.isWhitespace((char) i)) {
					desc.append((char) i);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if (((uki < $MAXU_MS + 1) && (!ukb))
								&& (uii < $MAXL_MS + 1)) {
							String s = desc.toString();
							// System.out.println("i "+uii+" j "+uji+" k "+uki+" value: "+s);
							ums_array[uii][uji][uki] = stringToDouble(s);
							uki++;
							if (uki == $MAXU_MS + 1) {
								ukb = true;// next
								if (uji < $MAXQ_MS)
									uji++;
								else {
									uji = 0;
									uii++;
								}
								uki = 0;
							}
						} else if (((fki < $MAXU_MS + 1) && (!fkb))
								&& (fii < $MAXL_MS + 1)) {
							String s = desc.toString();
							// System.out.println("i "+fii+" j "+fji+" k "+fki+" value: "+s);
							fms_array[fii][fji][fki] = stringToDouble(s);
							fki++;
							if (fki == $MAXU_MS + 1) {
								fkb = true;// next
								if (fji < $MAXQ_MS)
									fji++;
								else {
									fji = 0;
									fii++;
								}
								fki = 0;
							}
						} else if (((wki < $MAXU_MS) && (!wkb))
								&& (wii < $MAXL_MS + 1)) {
							String s = desc.toString();
							// System.out.println("i "+wii+" j "+wji+" k "+wki+" value: "+s);
							wms_array[wii][wji][wki] = stringToDouble(s);
							wki++;
							if (wki == $MAXU_MS) {
								wkb = true;// next
								if (wji < $MAXQ_MS)
									wji++;
								else {
									wji = 0;
									wii++;
								}
								wki = 0;
							}
						} else if (((iki < $MAXU_MS) && (!ikb))
								&& (iii < $MAXL_MS + 1)) {
							String s = desc.toString();
							// System.out.println("i "+iii+" j "+iji+" k "+iki+" value: "+s);
							ims_array[iii][iji][iki] = stringToInt(s);
							iki++;
							if (iki == $MAXU_MS) {
								ikb = true;// next
								if (iji < $MAXQ_MS)
									iji++;
								else {
									iji = 0;
									iii++;
								}
								iki = 0;
								// reconstruct:
								for (int k = 0; k <= $MAXU_MS - 1; k++) {
									fms_array[fii][fji][k] = fms_array[fii][fji][k + 1]
											/ fms_array[fii][fji][k] - 1;
									ims_array[iii][iji][k] = ims_array[iii][fji][k] - 1;
								}
								fms_array[fii][fji][$MAXU_MS] = fms_array[fii][fji][$MAXU_MS - 1];
								// reset and next read
								ukb = false;
								fkb = false;
								wkb = false;
								ikb = false;
							}
						}

					}
					desc = new StringBuffer();
				}
			}
			in.close();
		} catch (Exception e) {

		}
		seqStr = " done ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		llammin = Math.log($LAMBMIN_MS);
		llammax = Math.log($LAMBMAX_MS);
		dllamb = (llammax - llammin) / $MAXL_MS;
		dllambi = 1. / dllamb;
		dqms = $QMAX_MS / $MAXQ_MS;
		dqmsi = 1. / dqms;

	}

	// ---------------------------------
	// read a .pegs4.dat file generated by PEGS4 SLAC if true or generated by
	// java PEGS4A if false
	// title= filename(without extension) and index of medium ->index -1 for 0
	// biased var!!!
	// ALWAYS medium name is the same with filename!!!
	/**
	 * Called by HATCH. Read a .pegs4.dat file generated by 
	 * java PEGS4A. title= filename(without extension) and index of medium [index -1] 
	 * ALWAYS medium name is the same with filename.	
	 * @param title title
	 * @param index index
	 */
	public static void readPegs4Dat(String title, int index) {
		String filename = datas + file_sep + egsDatas + file_sep + intDir
				+ file_sep + title + pegs4datext;
		char comma = ',';
		String nulls = "";
		String equalSeq = "=";
		boolean haveData = false;
		int i = 0;
		int indx = 0;// for special cases
		int lnr = 0;// data number
		int asymi = 0;// at least 1 element!!
		boolean asymb = false;

		int zi = 0;// at least 1 element!!
		boolean zb = false;
		int ai = 0;// at least 1 element!!
		boolean ab = false;
		int pzi = 0;// at least 1 element!!
		boolean pzb = false;
		int rhozi = 0;// at least 1 element!!
		boolean rhozb = false;
		boolean rlcb = false;
		boolean aeb = false;
		boolean apb = false;
		boolean ueb = false;
		boolean upb = false;
		boolean inrstb = false;
		boolean iraylb = false;
		boolean mrangeb = false;
		boolean mcmfpb = false;
		boolean mlekeb = false;
		boolean mekeb = false;
		boolean msekeb = false;
		boolean mgeb = false;
		boolean msgeb = false;
		int dl1i = 0;
		int dl2i = 0;
		int dl3i = 0;
		int dl4i = 0;
		int dl5i = 0;
		int dl6i = 0;
		boolean dl1b = false;
		boolean dl2b = false;
		boolean dl3b = false;
		boolean dl4b = false;
		boolean dl5b = false;
		boolean dl6b = false;
		boolean delcmb = false;
		int alphii = 0;
		int bpari = 0;
		int delposi = 0;
		boolean alphib = false;
		boolean bparb = false;
		boolean delposb = false;
		boolean xr0b = false;
		boolean teff0b = false;
		boolean blccb = false;
		boolean xccb = false;
		boolean eke0b = false;
		boolean eke1b = false;
		int esig0i = 0;
		int psig0i = 0;
		int ededx0i = 0;
		int pdedx0i = 0;
		int ebr10i = 0;
		int pbr10i = 0;
		int pbr20i = 0;
		int tmxs0i = 0;
		int esig1i = 0;
		int psig1i = 0;
		int ededx1i = 0;
		int pdedx1i = 0;
		int ebr11i = 0;
		int pbr11i = 0;
		int pbr21i = 0;
		int tmxs1i = 0;
		boolean esig0b = false;
		boolean psig0b = false;
		boolean ededx0b = false;
		boolean pdedx0b = false;
		boolean ebr10b = false;
		boolean pbr10b = false;
		boolean pbr20b = false;
		boolean tmxs0b = false;
		boolean esig1b = false;
		boolean psig1b = false;
		boolean ededx1b = false;
		boolean pdedx1b = false;
		boolean ebr11b = false;
		boolean pbr11b = false;
		boolean pbr21b = false;
		boolean tmxs1b = false;
		boolean ebindab = false;
		boolean ge0b = false;
		boolean ge1b = false;
		int gmfp0i = 0;
		int gmfp1i = 0;
		int gbr10i = 0;
		int gbr11i = 0;
		int gbr20i = 0;
		int gbr21i = 0;
		boolean gmfp0b = false;
		boolean gmfp1b = false;
		boolean gbr10b = false;
		boolean gbr11b = false;
		boolean gbr20b = false;
		boolean gbr21b = false;
		boolean ngrb = false;
		boolean rco0b = false;
		boolean rco1b = false;
		int rsct0i = 0;
		int rsct1i = 0;
		int cohe0i = 0;
		int cohe1i = 0;
		boolean rsct0b = false;
		boolean rsct1b = false;
		boolean cohe0b = false;
		boolean cohe1b = false;

		StringBuffer desc = new StringBuffer();

		try {
			FileInputStream in = new FileInputStream(filename);

			while ((i = in.read()) != -1) {
				if (!Character.isWhitespace((char) i) && ((char) i != comma)) {
					desc.append((char) i);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if (lnr == 1)// MEDIUM=XXXX
						{
							// String s = desc.toString();//
							// System.out.println(s);
							// String med = stripAfterSeq(s, equalSeq);//
							// System.out.println(med);
							// NAPRZ=Convertor.stringToInt(NAPRZs);//System.out.println(NAPRZ);
							// NAPRZ=stringToInt(NAPRZs);//System.out.println(NAPRZ);
							// ZPRIM=new double[NAPRZ];
						}
						if (lnr == 2)// STERNCID=XXXX
						{
							// String s = desc.toString();//
							// System.out.println(s);
							// String med = stripAfterSeq(s, equalSeq);//
							// System.out.println(med);
						}
						if (lnr == 3)// ELEM or MIXT
						{
							// String s = desc.toString();//
							// System.out.println(s);
						}
						if (lnr == 4)// RHO=XXXX
						{
							String s = desc.toString();// System.out.println(s);
							String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
							if (med.compareTo(nulls) == 0) {
								// do nothing just set a flag
								// indx++;//eg 1
								// System.out.println("este ba");

								desc = new StringBuffer();
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);// System.out.println(desc.toString());
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											String s1 = desc.toString();// System.out.println("test"+s1);
											RHO[index - 1] = stringToDouble(s1);
											b = false;// exit the loop
										}
									}
								}
							} else {
								RHO[index - 1] = stringToDouble(med);
								// System.out.println(med);
							}
						}
						if (lnr == 5)// NE=XXXX
						{
							String s = desc.toString();// System.out.println("dupa 4"+s);
							String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
							if (med.compareTo(nulls) == 0) {
								desc = new StringBuffer();
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											String s1 = desc.toString();// System.out.println("test"+s1);
											NNE[index - 1] = stringToInt(s1);
											b = false;// exit the loop
										}
									}
								}
							} else {
								NNE[index - 1] = stringToInt(med);
								// System.out.println(med);
							}
						}

						if (lnr == 6)// IUNRST=XXXX but could be GASP
						{
							String s = desc.toString();// System.out.println(s);
							String bef = stripBeforeSeq(s, equalSeq);
							if (bef.compareTo("GASP") == 0) {
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0)// pegs4 based
																// GASP= xxx
								{
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											appendB = true;
										} else {
											// reach to the end of GASP=xxx
											// string so read next IUNRST
											b = false;// exit the loop
										}
									}
									b = true;
									appendB = false;
									while (b) {
										indx = in.read();
										if (((char) indx != comma))// no space
																	// serch for
																	// comma!!
										{
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												String med1 = stripAfterSeq(s1,
														equalSeq);
												IUNRST[index - 1] = stringToInt(med1);// System.out.println("pegs4-read");
												b = false;// exit the loop
											}
										}
									}

								} else// end of GASP=xxx string so read next
										// IUNRST=XXX in JAVA case
								{
									desc = new StringBuffer();// read next
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												String med1 = stripAfterSeq(s1,
														equalSeq);
												IUNRST[index - 1] = stringToInt(med1);// System.out.println("java-read");
												b = false;// exit the loop
											}
										}
									}
								}

							} else {
								// regular no GASP case
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												IUNRST[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									IUNRST[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
							}
						}

						if (lnr == 7)// EPSTFL=XXXX
						{
							String s = desc.toString();// System.out.println("dupa 4"+s);
							String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
							if (med.compareTo(nulls) == 0) {
								desc = new StringBuffer();
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											String s1 = desc.toString();// System.out.println("test"+s1);
											EPSTFL[index - 1] = stringToInt(s1);
											b = false;// exit the loop
										}
									}
								}
							} else {
								EPSTFL[index - 1] = stringToInt(med);
								// System.out.println(med);
							}
						}

						if (lnr == 8)// IAPRIM=XXXX
						{
							String s = desc.toString();// System.out.println("dupa 4"+s);
							String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
							if (med.compareTo(nulls) == 0) {
								desc = new StringBuffer();
								boolean b = true;
								boolean appendB = false;
								while (b) {
									indx = in.read();
									if (!Character.isWhitespace((char) indx)
											&& ((char) indx != comma)) {
										desc.append((char) indx);
										appendB = true;
									} else {
										if ((indx == -1) || (appendB))// EOF or
																		// normal
																		// exit
										{
											String s1 = desc.toString();// System.out.println("test"+s1);
											IAPRIM[index - 1] = stringToInt(s1);
											b = false;// exit the loop
										}
									}
								}
							} else {
								IAPRIM[index - 1] = stringToInt(med);
								// System.out.println(med);
							}
						}

						if (lnr > 8) {
							// *******************************************medium
							// composition
							if ((asymi < NNE[index - 1]) && (!asymb)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												ASYM[index - 1][asymi] = s1;
												b = false;// exit the loop
											}
										}
									}
								} else {
									ASYM[index - 1][asymi] = med;
									// System.out.println(med);
								}
								asymb = true;// set flag for next read
								asymi++;
							} else if ((zi < NNE[index - 1]) && (!zb))// //medium
																		// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												ZELEM[index - 1][zi] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									ZELEM[index - 1][zi] = stringToDouble(med);
									// System.out.println(med);
								}
								zb = true;// set flag for next read
								zi++;
							} else if ((ai < NNE[index - 1]) && (!ab))// //medium
																		// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												WA[index - 1][ai] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									WA[index - 1][ai] = stringToDouble(med);
									// System.out.println(med);
								}
								ab = true;// set flag for next read
								ai++;
							}// medium comp
							else if ((pzi < NNE[index - 1]) && (!pzb))// //medium
																		// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PZ[index - 1][pzi] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PZ[index - 1][pzi] = stringToDouble(med);
									// System.out.println(med);
								}
								pzb = true;// set flag for next read
								pzi++;
							}// medium comp
							else if ((rhozi < NNE[index - 1]) && (!rhozb))// //medium
																			// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												RHOZ[index - 1][rhozi] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									RHOZ[index - 1][rhozi] = stringToDouble(med);
									// System.out.println(med);
								}
								rhozb = true;// set flag for next read
								rhozi++;
								// repeat loop on ASYM
								asymb = false;
								rhozb = false;
								zb = false;
								ab = false;
								pzb = false;
							}// medium comp
								// *********************** "   MEDIA AND THRESH"
							else if ((!rlcb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												RLC[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									RLC[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								rlcb = true;
							} else if ((!aeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												AE[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									AE[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								aeb = true;
							} else if ((!apb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												AP[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									AP[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								apb = true;
							} else if ((!ueb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												UE[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									UE[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ueb = true;
							} else if ((!upb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												UP[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									UP[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								upb = true;
							}
							// **********************
							// "   ACTUAL ARRAY SIZES FROM PEGS"
							else if ((!msgeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MSGE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MSGE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								msgeb = true;
							} else if ((!mgeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MGE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MGE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								mgeb = true;
							} else if ((!msekeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MSEKE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MSEKE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								msekeb = true;
							} else if ((!mekeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MEKE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MEKE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								mekeb = true;
							} else if ((!mlekeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MLEKE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MLEKE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								mlekeb = true;
							} else if ((!mcmfpb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MCMFP[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MCMFP[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								mcmfpb = true;
							} else if ((!mrangeb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												MRANGE[index - 1] = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									MRANGE[index - 1] = stringToInt(med);
									// System.out.println(med);
								}
								mrangeb = true;
							} else if ((!iraylb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												IRAYL = stringToInt(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									IRAYL = stringToInt(med);
									// System.out.println(med);
								}
								iraylb = true;
							} else if ((!inrstb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												// String s1 =
												// desc.toString();//
												// System.out.println("test"+s1);
												// read nothing
												b = false;// exit the loop
											}
										}
									}
								} else {
									// read nothing//IRAYL=stringToInt(med);
									// System.out.println(med);
								}
								inrstb = true;
							}
							// ***********************
							// "   BREMPR"*************************************
							else if ((dl1i < 6) && (!dl1b))// //medium
															// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL1[dl1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL1[dl1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl1b = true;// set flag for next read
								dl1i++;
							} else if ((dl2i < 6) && (!dl2b))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL2[dl2i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL2[dl2i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl2b = true;// set flag for next read
								dl2i++;
							} else if ((dl3i < 6) && (!dl3b))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL3[dl3i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL3[dl3i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl3b = true;// set flag for next read
								dl3i++;
							} else if ((dl4i < 6) && (!dl4b))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL4[dl4i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL4[dl4i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl4b = true;// set flag for next read
								dl4i++;
							} else if ((dl5i < 6) && (!dl5b))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL5[dl5i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL5[dl5i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl5b = true;// set flag for next read
								dl5i++;
							} else if ((dl6i < 6) && (!dl6b))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DL6[dl6i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DL6[dl6i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								dl6b = true;// set flag for next read
								dl6i++;
								dl1b = false;
								dl2b = false;
								dl3b = false;
								dl4b = false;
								dl5b = false;
								dl6b = false;
							} else if ((!delcmb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DELCM[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DELCM[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								delcmb = true;
							} else if ((alphii < 2) && (!alphib))// //medium
																	// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												ALPHI[alphii][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									ALPHI[alphii][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								alphib = true;// set flag for next read
								alphii++;
							} else if ((bpari < 2) && (!bparb))// //medium
																// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												BPAR[bpari][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									BPAR[bpari][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								bparb = true;// set flag for next read
								bpari++;
							} else if ((delposi < 2) && (!delposb))// //medium
																	// composition
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												DELPOS[delposi][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									DELPOS[delposi][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								delposb = true;// set flag for next read
								delposi++;
								alphib = false;
								bparb = false;
								delposb = false;
							} else if ((!xr0b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												XR0[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									XR0[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								xr0b = true;
							} else if ((!teff0b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												TEFF0[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									TEFF0[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								teff0b = true;
							} else if ((!blccb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												BLCC[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									BLCC[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								blccb = true;
							} else if ((!xccb))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												XCC[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									XCC[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								xccb = true;
							} else if ((!eke0b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EKE0[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EKE0[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								eke0b = true;
							} else if ((!eke1b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EKE1[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EKE1[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								eke1b = true;
							} else if ((esig0i < MEKE[index - 1]) && (!esig0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												ESIG0[esig0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									ESIG0[esig0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								esig0b = true;// set flag for next read
								esig0i++;
							} else if ((esig1i < MEKE[index - 1]) && (!esig1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												ESIG1[esig1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									ESIG1[esig1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								esig1b = true;// set flag for next read
								esig1i++;
							} else if ((psig0i < MEKE[index - 1]) && (!psig0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PSIG0[psig0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PSIG0[psig0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								psig0b = true;// set flag for next read
								psig0i++;
							} else if ((psig1i < MEKE[index - 1]) && (!psig1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PSIG1[psig1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PSIG1[psig1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								psig1b = true;// set flag for next read
								psig1i++;
							} else if ((ededx0i < MEKE[index - 1])
									&& (!ededx0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EDEDX0[ededx0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EDEDX0[ededx0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ededx0b = true;// set flag for next read
								ededx0i++;
							} else if ((ededx1i < MEKE[index - 1])
									&& (!ededx1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EDEDX1[ededx1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EDEDX1[ededx1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ededx1b = true;// set flag for next read
								ededx1i++;
							} else if ((pdedx0i < MEKE[index - 1])
									&& (!pdedx0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PDEDX0[pdedx0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PDEDX0[pdedx0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pdedx0b = true;// set flag for next read
								pdedx0i++;
							} else if ((pdedx1i < MEKE[index - 1])
									&& (!pdedx1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PDEDX1[pdedx1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PDEDX1[pdedx1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pdedx1b = true;// set flag for next read
								pdedx1i++;
							} else if ((ebr10i < MEKE[index - 1]) && (!ebr10b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EBR10[ebr10i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EBR10[ebr10i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ebr10b = true;// set flag for next read
								ebr10i++;
							} else if ((ebr11i < MEKE[index - 1]) && (!ebr11b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EBR11[ebr11i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EBR11[ebr11i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ebr11b = true;// set flag for next read
								ebr11i++;
							} else if ((pbr10i < MEKE[index - 1]) && (!pbr10b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PBR10[pbr10i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PBR10[pbr10i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pbr10b = true;// set flag for next read
								pbr10i++;
							} else if ((pbr11i < MEKE[index - 1]) && (!pbr11b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PBR11[pbr11i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PBR11[pbr11i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pbr11b = true;// set flag for next read
								pbr11i++;
							} else if ((pbr20i < MEKE[index - 1]) && (!pbr20b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PBR20[pbr20i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PBR20[pbr20i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pbr20b = true;// set flag for next read
								pbr20i++;
							} else if ((pbr21i < MEKE[index - 1]) && (!pbr21b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												PBR21[pbr21i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									PBR21[pbr21i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								pbr21b = true;// set flag for next read
								pbr21i++;
							} else if ((tmxs0i < MEKE[index - 1]) && (!tmxs0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												TMXS0[tmxs0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									TMXS0[tmxs0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								tmxs0b = true;// set flag for next read
								tmxs0i++;
							} else if ((tmxs1i < MEKE[index - 1]) && (!tmxs1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												TMXS1[tmxs1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									TMXS1[tmxs1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								tmxs1b = true;// set flag for next read
								tmxs1i++;
								// null and re-read next set
								esig0b = false;
								esig1b = false;
								psig0b = false;
								psig1b = false;
								ededx0b = false;
								ededx1b = false;
								pdedx0b = false;
								pdedx1b = false;
								ebr10b = false;
								ebr11b = false;
								pbr10b = false;
								pbr11b = false;
								pbr20b = false;
								pbr21b = false;
								tmxs0b = false;
								tmxs1b = false;
							}
							// **************
							// PHOTIN*******************************************************************
							else if ((!ebindab))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												EBINDA[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									EBINDA[index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								ebindab = true;
							} else if ((!ge0b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();
												// System.out.println("ge0 test  "+s1);
												GE0[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GE0[index - 1] = stringToDouble(med);
									// System.out.println("ge0  "+med);
								}
								ge0b = true;
							} else if ((!ge1b))// flag
							{
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();
												// System.out.println("ge1.test  "+s1);
												GE1[index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GE1[index - 1] = stringToDouble(med);
									// System.out.println("ge1   "+med);
								}
								ge1b = true;
							} else if ((gmfp0i < MGE[index - 1]) && (!gmfp0b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GMFP0[gmfp0i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GMFP0[gmfp0i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gmfp0b = true;// set flag for next read
								gmfp0i++;
							} else if ((gmfp1i < MGE[index - 1]) && (!gmfp1b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GMFP1[gmfp1i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GMFP1[gmfp1i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gmfp1b = true;// set flag for next read
								gmfp1i++;
							} else if ((gbr10i < MGE[index - 1]) && (!gbr10b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GBR10[gbr10i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GBR10[gbr10i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gbr10b = true;// set flag for next read
								gbr10i++;
							} else if ((gbr11i < MGE[index - 1]) && (!gbr11b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GBR11[gbr11i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GBR11[gbr11i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gbr11b = true;// set flag for next read
								gbr11i++;
							} else if ((gbr20i < MGE[index - 1]) && (!gbr20b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GBR20[gbr20i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GBR20[gbr20i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gbr20b = true;// set flag for next read
								gbr20i++;
							} else if ((gbr21i < MGE[index - 1]) && (!gbr21b)) {
								String s = desc.toString();// System.out.println("dupa 4"+s);
								String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
								if (med.compareTo(nulls) == 0) {
									desc = new StringBuffer();
									boolean b = true;
									boolean appendB = false;
									while (b) {
										indx = in.read();
										if (!Character
												.isWhitespace((char) indx)
												&& ((char) indx != comma)) {
											desc.append((char) indx);
											appendB = true;
										} else {
											if ((indx == -1) || (appendB))// EOF
																			// or
																			// normal
																			// exit
											{
												String s1 = desc.toString();// System.out.println("test"+s1);
												GBR21[gbr21i][index - 1] = stringToDouble(s1);
												b = false;// exit the loop
											}
										}
									}
								} else {
									GBR21[gbr21i][index - 1] = stringToDouble(med);
									// System.out.println(med);
								}
								gbr21b = true;// set flag for next read
								gbr21i++;
								// reset
								gmfp0b = false;
								gmfp1b = false;
								gbr10b = false;
								gbr11b = false;
								gbr20b = false;
								gbr21b = false;
							}
							// ******************"   PHOTIN (CONTINUED)---OPTIONAL RAYLEIGH SCATTERING INPUT"
							else if (IRAYL == 1) {
								if ((!ngrb))// flag
								{
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													NGR[index - 1] = stringToInt(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										NGR[index - 1] = stringToInt(med);
										// System.out.println(med);
									}
									ngrb = true;
								} else if ((!rco0b))// flag
								{
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													RCO0[index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										RCO0[index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									rco0b = true;
								} else if ((!rco1b))// flag
								{
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													RCO1[index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										RCO1[index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									rco1b = true;
								} else if ((rsct0i < NGR[index - 1])
										&& (!rsct0b)) {
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													RSCT0[rsct0i][index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										RSCT0[rsct0i][index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									rsct0b = true;// set flag for next read
									rsct0i++;
								} else if ((rsct1i < NGR[index - 1])
										&& (!rsct1b)) {
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													RSCT1[rsct1i][index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										RSCT1[rsct1i][index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									rsct1b = true;// set flag for next read
									rsct1i++;
									// reset
									rsct0b = false;
									rsct1b = false;
								} else if ((cohe0i < MGE[index - 1])
										&& (!cohe0b)) {
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													COHE0[cohe0i][index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										COHE0[cohe0i][index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									cohe0b = true;// set flag for next read
									cohe0i++;
								} else if ((cohe1i < MGE[index - 1])
										&& (!cohe1b)) {
									String s = desc.toString();// System.out.println("dupa 4"+s);
									String med = stripAfterSeq(s, equalSeq);// System.out.println(med);
									if (med.compareTo(nulls) == 0) {
										desc = new StringBuffer();
										boolean b = true;
										boolean appendB = false;
										while (b) {
											indx = in.read();
											if (!Character
													.isWhitespace((char) indx)
													&& ((char) indx != comma)) {
												desc.append((char) indx);
												appendB = true;
											} else {
												if ((indx == -1) || (appendB))// EOF
																				// or
																				// normal
																				// exit
												{
													String s1 = desc.toString();// System.out.println("test"+s1);
													COHE1[cohe1i][index - 1] = stringToDouble(s1);
													b = false;// exit the loop
												}
											}
										}
									} else {
										COHE1[cohe1i][index - 1] = stringToDouble(med);
										// System.out.println(med);
									}
									cohe1b = true;// set flag for next read
									cohe1i++;
									// reset
									cohe0b = false;
									cohe1b = false;
								}

							}// if (IRAYL==1)
						}// if(lnr>8)
					}// we have data

					desc = new StringBuffer();
				}// else
			}// while
			in.close();
		}// try
		catch (Exception e) {

		}

	}

	/**
	 * Called by readPEGS4Dat.
	 * @param toBeStripped toBeStripped
	 * @param seqToLook seqToLook
	 * @return the result
	 */
	public static String stripAfterSeq(String toBeStripped, String seqToLook) {
		String ss = "";
		int idx = toBeStripped.lastIndexOf(seqToLook);
		if (idx < 1)
			ss = toBeStripped;
		else
			ss = toBeStripped.substring(idx + 1);
		return ss;
	}

	/**
	 * Called by readPEGS4Dat.
	 * @param toBeStripped toBeStripped
	 * @param seqToLook seqToLook
	 * @return the result
	 */
	public static String stripBeforeSeq(String toBeStripped, String seqToLook) {
		String ss = "";
		int idx = toBeStripped.lastIndexOf(seqToLook);
		if (idx < 1)
			ss = toBeStripped;
		else
			ss = toBeStripped.substring(0, idx);
		return ss;
	}

	// max from 3 var
	/**
	 * Return the maximum from 3 doubles.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double max(double a, double b, double c) {
		double r = Math.max(a, b);
		r = Math.max(r, c);
		return r;
	}

	/**
	 * Return the minimum from 3 doubles.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double min(double a, double b, double c) {
		double r = Math.min(a, b);
		r = Math.min(r, c);
		return r;
	}

	// return [0,1] random number
	/**
	 * Return a random number in [0,1] interval.
	 * @return the result.
	 */
	public static double random01() {
		// To get <0,1> instead of <0,1), we substract the random number from 1
		// in half cases
		double result = 0.0;
		if (RandomUse == 0) {
			double randomNumber = Math.random();
			if (Math.random() < 0.5)
				randomNumber = 1.0 - randomNumber;
			return randomNumber;
		} else if (RandomUse == 1) {
			result = RANDOMSET();
		} else if (RandomUse == 2) {
			result = Math.random();
		}

		return result;
	}

	/**
	 * Sets cubic spline interpolation coefficients for the data contained 
	 * in the array f(n) at the abscissas x(n). 
	 * @param x x
	 * @param f f
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param n n
	 */
	public static void set_spline(double[] x, double[] f, double[] a,
			double[] b, double[] c, double[] d, int n) {
		// "======================================================================"
		// "
		// " Sets cubic spline interpolation coefficients for the data contained  "
		// " in the array f(n) at the abscissas x(n)                              "
		// " original fortran by I.Kawrakow, NRC,                                 "
		// "======================================================================"
		// ; Copyright NRC;

		// aa=new double[n];bb=new double[n];cc=new double[n];dd=new
		// double[n];xx=new double[n];

		int m1 = 0;
		int m2 = 0;
		int m = 0;
		int mr = 0;
		double s = 0.0;
		double r = 0.0;
		m1 = 2;
		m2 = n - 1;
		s = 0.0;
		// DO m=1,m2 [
		for (m = 0; m < m2; m++)// 0 - n-2
		{
			d[m] = x[m + 1] - x[m];
			r = (f[m + 1] - f[m]) / d[m];
			c[m] = r - s;
			s = r;
		}
		// /s,r,c(1),c(n)/=0;
		s = 0.0;
		r = 0.0;
		c[0] = 0.0;// reset to null
		c[n - 1] = 0.0;// reset to null
		// DO m=m1,m2 [
		for (m = m1 - 1; m < m2; m++)// 1 - n-2
		{
			c[m] = c[m] + r * c[m - 1];
			b[m] = 2.0 * (x[m - 1] - x[m + 1]) - r * s;
			s = d[m];
			r = s / b[m];
		}
		// mr = m2;
		mr = m2 - 1;// n-2
		// DO m=m1,m2 [
		for (m = m1 - 1; m < m2; m++)// 1 - n-2
		{
			c[mr] = (d[mr] * c[mr + 1] - c[mr]) / b[mr];
			mr = mr - 1;
		}
		// DO m=1,m2 [
		for (m = 0; m < m2; m++)// 0 - n-2
		{
			s = d[m];
			r = c[m + 1] - c[m];
			d[m] = r / s;
			c[m] = 3.0 * c[m];
			b[m] = (f[m + 1] - f[m]) / s - (c[m] + r) * s;
			a[m] = f[m];
		}
	}

	/**
	 * Returns the value of the function at s using the spline coefficients 
	 * a,b,c,d, which must have been set using set_spline.
	 * @param s s
	 * @param x x
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param n n
	 * @return the result
	 */
	public static double spline(double s, double[] x, double[] a, double[] b,
			double[] c, double[] d, int n) {
		// "======================================================================"
		// "                                                                      "
		// " Returns the value of the function at s using the spline coefficients "
		// " a,b,c,d, which must have been set using set_spline                   "
		// "                                                                      "
		// " original fortran I.Kawrakow, NRC                                     "
		// "======================================================================"

		// ; Copyright NRC;

		// $INTEGER n;
		// $REAL s,x(n),a(n),b(n),c(n),d(n);

		int m_lower = 0;
		int m_upper = 0;
		int direction = 0;
		int m = 0;
		int ml = 0;
		int mu = 0;
		int mav = 0;
		double q = 0.0;

		// if( x(1) > x(n) ) [ direction = 1; m_lower = n; m_upper = 0; ]
		if (x[0] > x[n - 1]) {
			direction = 1;
			m_lower = n;
			m_upper = 0;
		} else {
			direction = 0;
			m_lower = 0;
			m_upper = n;
		}
		// if ( s >= x(m_upper + direction) ) [
		if (s >= x[m_upper + direction - 1]) {
			m = m_upper + 2 * direction - 1;
		}
		// else if ( s <= x(m_lower+1-direction) ) [
		else if (s <= x[m_lower - direction]) {
			m = m_lower - 2 * direction + 1;
		} else {// " Perform a binary search to find the interval s is in "
			ml = m_lower;
			mu = m_upper;
			// while ( iabs(mu-ml) > 1 ) [
			while (Math.abs(mu - ml) > 1) {
				mav = (ml + mu) / 2;
				// if( s < x(mav) )
				if (s < x[mav - 1]) {
					mu = mav;
				} else {
					ml = mav;
				}
			}
			m = mu + direction - 1;
		}

		q = s - x[m - 1];
		// double spline = a(m) + q*(b(m) + q*(c(m) + q*d(m)));
		double spline = a[m - 1] + q
				* (b[m - 1] + q * (c[m - 1] + q * d[m - 1]));
		return spline;
	}

	// "****************************************************************************
	/**
	 * Prepare an alias sampling table, given the histogram probabilities xs_array,fs_array.
	 * @param nsbin nsbin
	 * @param xs_array xs_array
	 * @param fs_array fs_array
	 * @param ws_array ws_array
	 * @param ibin_array ibin_array
	 */
	public static void prepare_alias_table(int nsbin, double[] xs_array,
			double[] fs_array, double[] ws_array, int[] ibin_array) {

		// " Prepare an alias sampling table, given the histogram probabilities
		// " xs_array,fs_array.
		// "*****************************************************************************
		// ; Copyright NRC;
		// integer nsbin;
		// $INTEGER ibin_array(nsbin);
		// $REAL xs_array(0:nsbin),fs_array(0:nsbin),ws_array(nsbin);
		// wss_array=new double[nsbin];
		// ibinn_array=new int[nsbin];
		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;
		boolean e1 = false;
		boolean e2 = false;
		sum = 0.0;
		// DO i=1,nsbin [
		for (i = 1; i <= nsbin; i++)// !!!!!!!!!!!!!!!see above definition
		// xs,fs-array of nsbin+1!!!!!
		{
			aux = 0.5 * (fs_array[i] + fs_array[i - 1])
					* (xs_array[i] - xs_array[i - 1]);
			if (aux < 1.e-30)
				aux = 1.e-30;
			ws_array[i - 1] = -aux;// ws_array[i] = -aux;
			ibin_array[i - 1] = 1;

			sum = sum + aux;
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			e1 = false;
			e2 = false;
			// DO j_h = 1,nsbin [
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (Math.abs(ws_array[j_h - 1]) > sum) {// GOTO :AT_EXIT_1:;
						e1 = true;
						break;// exit the loop
					}
				}
			}
			// j_h = nsbin;
			if (!e1)
				j_h = nsbin;
			// :AT_EXIT_1:
			// DO j_l = 1,nsbin [
			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (Math.abs(ws_array[j_l - 1]) < sum) {// GOTO :AT_EXIT_2:;
						e2 = true;
						break;// exit the loop
					}
				}
			}
			// j_l = nsbin;
			if (!e2)
				j_l = nsbin;
			// :AT_EXIT_2:
			aux = sum - Math.abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;
			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1.0;
			}
		}

	}

	// "******************************************************************************
	/**
	 * Sample a random variable from the alias table xs_array,fs_array,ws_array,ibin_array 
	 * which must have been prepared with prepare_alias_table.
	 * @param nsbin nsbin
	 * @param xs_array xs_array
	 * @param fs_array fs_array
	 * @param ws_array ws_array
	 * @param ibin_array ibin_array
	 * @return the result.
	 */
	public static double alias_sample1(int nsbin, double[] xs_array,
			double[] fs_array, double[] ws_array, int[] ibin_array) {
		// " Sample a random variable from the alias table
		// " xs_array,fs_array,ws_array,ibin_array
		// " which must have been prepared with prepare_alias_table
		// "******************************************************************************
		// ; Copyright NRC;

		// integer nsbin;
		// $INTEGER ibin_array(nsbin);
		// $REAL xs_array(0:nsbin),fs_array(0:nsbin),ws_array(nsbin);

		// ;COMIN/RANDOM/;
		double alias_sample1 = 0.0;
		int j = 0;
		double r1 = 0.0;
		double r2 = 0.0;
		double aj = 0.0;
		double x = 0.0;
		double dx = 0.0;
		double a = 0.0;
		double rnno1 = 0.0;

		r1 = random01();
		r2 = random01();
		aj = 1 + r1 * nsbin;
		Double de = new Double(aj);
		j = de.intValue();// j = aj;->minim 1 si maxim nsbin
		aj = aj - j;
		// -------------------------------
		// if( aj > ws_array(j) ) j = ibin_array(j);
		if (aj > ws_array[j - 1])
			j = ibin_array[j - 1];

		x = xs_array[j - 1];// 0 biased OK!!
		dx = xs_array[j] - x;

		if (fs_array[j - 1] > 0.0) {
			a = fs_array[j] / fs_array[j - 1] - 1.0;
			if (Math.abs(a) < 0.2) {
				rnno1 = 0.5 * (1.0 - r2) * a;
				alias_sample1 = x + r2 * dx * (1.0 + rnno1 * (1.0 - r2 * a));
			} else {
				alias_sample1 = x - dx / a
						* (1.0 - Math.sqrt(1.0 + r2 * a * (2.0 + a)));
			}
		} else {
			alias_sample1 = x + dx * Math.sqrt(r2);
		}

		return alias_sample1;
	}

	/**
	 * Given the lower and upper limit of integration, x1 and x2, 
	 * and given n, this routine returns arrays x and w, containing the abscissas and weights of the Gauss-Legendre 
	 * n - point quadrature formula. Called by init_nist_brems.
	 * @param x1 x1
	 * @param x2 x2
	 * @param x x
	 * @param w w
	 * @param n n
	 */
	public static void gauss_legendre(double x1, double x2, double[] x,
			double[] w, int n) {

		// " Given the lower and upper limit of integration, x1 and x2,
		// " and given n, this routine returns arrays x and w,
		// " containing the abscissas and weights of the Gauss-Legendre
		// " n - point quadrature formula
		// "******************************************************************************
		// ; Copyright NRC;

		double eps = 3.E-14;
		double Pi = Math.PI;
		// parameter (eps = 3.D-14, Pi = 3.141592654D0);

		int i = 0;
		int m = 0;
		int j = 0;
		double xm = 0.0;
		double xl = 0.0;
		double z = 0.0;
		double z1 = 0.0;
		double p1 = 0.0;
		double p2 = 0.0;
		double p3 = 0.0;
		double pp = 0.0;

		m = (n + 1) / 2;
		// xm=0.5d0*(x2+x1); xl=0.5d0*(x2-x1);
		xm = 0.5 * (x2 + x1);
		xl = 0.5 * (x2 - x1);
		for (i = 1; i <= m; i++) {
			// z=cos(Pi*(i-.25d0)/(n+.5d0));
			z = Math.cos(Pi * (i - 0.25) / (n + 0.5));
			while (Math.abs(z - z1) >= eps)
			// LOOP
			{
				// p1=1.d0;
				// p2=0.d0;
				p1 = 1.0;
				p2 = 0.0;
				for (j = 1; j <= n; j++) {
					p3 = p2;
					p2 = p1;
					// p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j;
					p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
				}
				// pp=n*(z*p1-p2)/(z*z-1.d0);
				pp = n * (z * p1 - p2) / (z * z - 1.0);
				z1 = z;
				z = z1 - p1 / pp;
			}// UNTIL (abs(z-z1) < eps);
				// x(i)=xm-xl*z; x(n+1-i)=xm+xl*z;
				// w(i)=2.d0*xl/((1.d0-z*z)*pp*pp); w(n+1-i)=w(i);
			x[i - 1] = xm - xl * z;
			x[n - i] = xm + xl * z;
			// w(i-1)=2.d0*xl/((1.d0-z*z)*pp*pp);
			w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
			w[n - i] = w[i];
		}
	}

	// "************************************************************************
	// " an error function routine which is needed since some of
	// " the compiler don't have it as an intrinsic
	// " Originally came from some library somewhere (Harwell I think)
	// " recoded in mortran
	// "************************************************************************
	/**
	 * An error function routine which is needed since some of the compiler don't have it as an intrinsic 
	 * Originally came from some library somewhere (Harwell I think) recoded in mortran. Called by init_compton.
	 * @param X X
	 * @return the result
	 */
	public static double ERF1(double X) {
		double erf1 = 0.0;
		double x = X;

		// double precision A(0:22,2); //" Coefficients in expansion for erf(x)
		// if x<3
		// " (K=1) and for erfc(x) x>3 (K=2)
		// double precision
		// CONST, //" 2/sqrt(pi)
		// BN,BN1,BN2, //" Recursion coefficients B(n),B(n+1),B(n+2)
		// Y,FAC; //" y=x/3 or 3/x and FAC = 2(2y**2-1)
		// $INTEGER N, //" recursion index n
		// K, //" K=1,2 for x <= 3 or x > 3
		// NLIM(2); //" Maximum value of n in sum for K=1,2
		double y = 0.0;
		double FAC = 0.0;
		double BN = 0.0;
		double BN1 = 0.0;
		double BN2 = 0.0;
		int k = 0;
		int n = 0;
		double[][] A = { { 1.0954712997776232, 0.9750834237085559 },
				{ -0.2891754011269890, -0.0240493938504146 },
				{ 0.1104563986337951, 0.0008204522408804 },
				{ -0.0412531882278565, -0.0000434293081303 },
				{ 0.0140828380706516, 0.0000030184470340 },
				{ -0.0043292954474314, -0.0000002544733193 },
				{ 0.0011982719015923, 0.0000000248583530 },
				{ -0.0002999729623532, -0.0000000027317201 },
				{ 0.0000683258603789, 0.0000000003308472 },
				{ -0.0000142469884549, 0.0000000000001464 },
				{ 0.0000027354087728, -0.0000000000000244 },
				{ -0.0000004861912872, 0.0000000000000042 },
				{ 0.0000000803872762, -0.0000000000000008 },
				{ -0.0000000124184183, 0.0000000000000001 },
				{ 0.0000000017995326, 0.0 }, { -0.0000000002454795, 0.0 },
				{ 0.0000000000316251, 0.0 }, { -0.0000000000038590, 0.0 },
				{ 0.0000000000004472, 0.0 }, { -0.0000000000000493, 0.0 },
				{ 0.0000000000000052, 0.0 }, { -0.0000000000000005, 0.0 },
				{ 0.0000000000000001, 0.0 } };
		int[] NLIM = { 22, 16 };
		double CONST = 2 / Math.sqrt(Math.PI);// 1.128379167095513 /;

		if (x > 3) {
			y = 3.0 / x;
			k = 2;
		} else {
			y = x / 3.0;
			k = 1;
		}
		double Y = y;
		int K = k;
		// int N = n;
		// " Calculate sum of Chebyshev polynomials by backwards recursion
		// "
		// " sum { A(n)*T(2n+1;y) : n=0,N } = y * ( B(0) - B(1) )
		// " sum { A(n)*T(2n;y) : n=0,N } = ( B(0) - (2*y**2-1) * B(1) ) / 2
		// " = ( B(0) - B(2) + A(0) ) / 2
		// "
		// " where B(N+2) = B(N+1) = 0
		// " and B(n) = 2*(2*y**2-1)*B(n+1) - B(n+2) + A(n) for
		// n=N,(N-1),...,1,0
		// "
		FAC = 2.0 * (2.0 * Y * Y - 1.0);
		BN1 = 0.0; // " Initialise B(N+2) = 0
		BN = 0.0; // " Initialise B(N+1) = 0

		// for(n = NLIM(K),0,-1 [
		for (n = NLIM[K - 1]; n >= 0; n--) {
			BN2 = BN1;
			BN1 = BN;
			// BN = FAC * BN1 - BN2 + A(N,K)
			BN = FAC * BN1 - BN2 + A[n][K - 1];
		}

		if (k == 1) {
			erf1 = CONST * Y * (BN - BN1);
		} else {
			// erf1 = 1 - CONST * Math.exp(-X*X) * ( BN - BN2 + A(0,K) )/(4.0 *
			// X);
			erf1 = 1.0 - CONST * Math.exp(-X * X) * (BN - BN2 + A[0][K - 1])
					/ (4.0 * X);
		}

		return erf1;
	}

	// ****************************************************************************
	// erf from c
	/**
	 * Complementary error function.
	 * @param x x
	 * @return the result
	 */
	public static double erfc1(double x)
	// ****************************************************************************
	//
	// Purpose:
	//
	// ERF1 evaluates the error function.
	//
	// Parameters:
	//
	// Input, double *X, the argument.
	//
	// Output, double ERF1, the value of the error function at X.
	//
	{
		double c = .564189583547756e0;
		double[] a = { .771058495001320e-04, -.133733772997339e-02,
				.323076579225834e-01, .479137145607681e-01,
				.128379167095513e+00 };
		double[] b = { .301048631703895e-02, .538971687740286e-01,
				.375795757275549e+00 };
		double[] p = { -1.36864857382717e-07, 5.64195517478974e-01,
				7.21175825088309e+00, 4.31622272220567e+01,
				1.52989285046940e+02, 3.39320816734344e+02,
				4.51918953711873e+02, 3.00459261020162e+02 };
		double[] q = { 1.00000000000000e+00, 1.27827273196294e+01,
				7.70001529352295e+01, 2.77585444743988e+02,
				6.38980264465631e+02, 9.31354094850610e+02,
				7.90950925327898e+02, 3.00459260956983e+02 };
		double[] r = { 2.10144126479064e+00, 2.62370141675169e+01,
				2.13688200555087e+01, 4.65807828718470e+00,
				2.82094791773523e-01 };
		double[] s = { 9.41537750555460e+01, 1.87114811799590e+02,
				9.90191814623914e+01, 1.80124575948747e+01 };
		double erf1, ax, bot, t, top, x2;

		ax = Math.abs(x);
		// if(ax > 0.5e0)
		// {
		// goto S10;
		// }
		if (ax <= 0.5e0) {
			t = x * x;
			top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]
					+ 1.0e0;
			bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
			erf1 = x * (top / bot);
			return erf1;
		}
		// S10:
		// if(ax > 4.0e0) goto S20;
		else if (ax <= 4.0e0) {
			top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4])
					* ax + p[5])
					* ax + p[6])
					* ax + p[7];
			bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4])
					* ax + q[5])
					* ax + q[6])
					* ax + q[7];
			erf1 = 0.5e0 + (0.5e0 - Math.exp(-(x * x)) * top / bot);
			if (x < 0.0e0)
				erf1 = -erf1;
			return erf1;
		}
		// S20:
		// if(ax >= 5.8e0) goto S30;
		else if (ax < 5.8e0) {
			x2 = x * x;
			t = 1.0e0 / x2;
			top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
			bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
			erf1 = (c - top / (x2 * bot)) / ax;
			erf1 = 0.5e0 + (0.5e0 - Math.exp(-x2) * erf1);
			if (x < 0.0e0)
				erf1 = -erf1;
			return erf1;
		}
		// S30:
		// erf1 = fifdsign(1.0e0,*x);
		else {
			erf1 = 1.0e0;
			if (x < 0.0e0)
				erf1 = -erf1;
		}
		return erf1;
	}

	// --------------convertors to be independent with other class
	// file--------------------------
	/**
	 * Convert an int to other base (than 10).
	 * @param value value
	 * @param base base
	 * @return the result
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertIntToOtherBase(String value, int base)
			throws NumberFormatException// 3a3f15f7=format Hexazecimal ex.
	{
		BigInteger bi = new BigInteger(value);
		return bi.toString(base);
	}

	/**
	 * Convert a number from other (than 10) base to an int.
	 * @param value value
	 * @param base base
	 * @return the result
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertOtherBaseToInt(String value, int base)
			throws NumberFormatException {
		BigInteger bi = new BigInteger(value, base);
		return bi.toString(10);
	}

	/**
	 * Converts a string to int value.
	 * 
	 * @param value
	 *            the string
	 * @return the int value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static int stringToInt(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Integer.parseInt(value);
		}
	}

	/**
	 * Converts an int number to string.
	 * 
	 * @param i
	 *            the int value
	 * @return the string representation
	 */
	public static String intToString(int i) {
		return Integer.toString(i);
	}

	/**
	 * Converts a string to float value.
	 * 
	 * @param value
	 *            the string
	 * @return the float value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static float stringToFloat(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Float.parseFloat(value);
		}
	}

	/**
	 * Converts a float number to string.
	 * 
	 * @param f
	 *            the float value
	 * @return the string representation
	 */
	public static String floatToString(float f) {
		return Float.toString(f);
	}

	/**
	 * Converts a string to double value.
	 * 
	 * @param value
	 *            the string
	 * @return the double value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static double stringToDouble(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Double.parseDouble(value);
		}
	}

	/**
	 * Converts a double number to string.
	 * 
	 * @param d
	 *            the double value
	 * @return the string representation
	 */
	public static String doubleToString(double d) {
		return Double.toString(d);
	}

	// BIT OPS OVER LATCH!!!
	// Function to set bit "index" in LATCH to 1. index start with 1, 2,...
	// ilatch is index of int value from LATCH array starting with 1.
	/**
	 * Function to set bit "index" in LATCH to 1. index start with 1, 2,... 
	 * ilatch is index of int value from LATCH array starting with 1.
	 * @param ilatch ilatch
	 * @param index index
	 * @return the result
	 */
	public static int IBSET_LATCH(int ilatch, int index) {
		int result = 0;
		int il = ilatch - 1;// real index starting from 0!!
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCH[il]);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		BigInteger bi2 = bi.setBit(ib);// set the bit to 1.
		result = bi2.intValue();// the result

		return result;
	}

	/**
	 * Similar with above, LATCHI is considered here.
	 * @param index index
	 * @return the result
	 */
	public static int IBSET_LATCHI(int index) {
		int result = 0;
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCHI);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		BigInteger bi2 = bi.setBit(ib);// set the bit to 1.
		result = bi2.intValue();// the result

		return result;
	}

	// Function to set bit "index" in LATCH to 0. index start with 1, 2,...
	// ilatch is index of int value from LATCH array starting with 1.
	/**
	 * Inverse of IBSET. Clear bit operation is performed here.
	 * @param ilatch ilatch
	 * @param index index
	 * @return the result
	 */
	public static int IBCLR_LATCH(int ilatch, int index) {
		int result = 0;
		int il = ilatch - 1;// real index starting from 0!!
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCH[il]);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		BigInteger bi2 = bi.clearBit(ib);// set the bit to 0.
		result = bi2.intValue();// the result

		return result;
	}

	/**
	 * Inverse of IBSET. Clear bit operation is performed here.
	 * @param index index
	 * @return the result
	 */
	public static int IBCLR_LATCHI(int index) {
		int result = 0;
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCHI);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		BigInteger bi2 = bi.clearBit(ib);// set the bit to 0.
		result = bi2.intValue();// the result

		return result;
	}

	// Function to test bit "index" in LATCH
	// ilatch is index of int value from LATCH array starting with 1.
	/**
	 * Function to test bit "index" in LATCH.
	 * ilatch is index of int value from LATCH array starting with 1.
	 * @param ilatch ilatch
	 * @param index index
	 * @return the result
	 */
	public static boolean BTEST_LATCH(int ilatch, int index) {
		boolean result = false;
		int il = ilatch - 1;// real index starting from 0!!
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCH[il]);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		result = bi.testBit(ib);// test the bit if 1.

		return result;
	}

	/**
	 * Similar with above, LATCHI is considered here.
	 * @param index index
	 * @return the result
	 */
	public static boolean BTEST_LATCHI(int index) {
		boolean result = false;
		int ib = index - 1;// bit index starting from 0!!

		String s = intToString(LATCHI);// string representation of int value
		BigInteger bi = new BigInteger(s);// in order to construct the
											// BigInteger object
		result = bi.testBit(ib);// test the bit if 1.

		return result;
	}

	// --------------END convertors to be independent with other class
	// file--------------------------

	/**
	 * Set EGS4 variables to default values.
	 */
	public static void egs_set_defaults() {
		VACDST = 1.e8;
		// $set-region-by-region-defaults;
		// REPLACE {$set-region-by-region-defaults;} WITH {;
		for (int i = 0; i < $MXREG; i++) {
			ECUT[i] = $GLOBAL_ECUT;
			PCUT[i] = $GLOBAL_PCUT; // "cut-off energies"
			ibcmp[i] = $IBCMP_DEFAULT; // "Compton "
			iedgfl[i] = $IEDGFL_DEFAULT;// "Relaxations"
			iphter[i] = $IPHTER_DEFAULT;// "photo-electron angular distribution"
			SMAXIR[i] = $MAX_SMAX; // "maximum step size"
			i_do_rr[i] = 0; // "range rejection flag"
			e_max_rr[i] = 0.0; // "`save' energy for range rejection"
			MED[i] = 1; // "default medium"
			RHOR[i] = 0.0; // "default mass density"
			IRAYLR[i] = $IRAYLR_DEFAULT;// "Rayleigh flag"
		}
		eii_flag = 0; // "No EII by default. "
		for (int i = 0; i < $MXMED; i++) {
			IRAYLM[i] = 0; // "Rayleigh data available?"
			// " set all thresholds to zero "
			AE[i] = 0.0;
			AP[i] = 0.0;
			UE[i] = 0.0;
			UP[i] = 0.0;
			TE[i] = 0.0;
			THMOLL[i] = 0.0;
		}
		for (int i = 0; i < $MXSHELL; i++) {
			for (int j = 0; j < $MXELEMENT; j++) {
				binding_energies[i][j] = 0.0;
			}
		}
		ibrdst = $IBRDST_DEFAULT; // " brems angular sampling"
		ibr_nist = $IBR_NIST_DEFAULT; // " flag for brems from NIST data base "
		pair_nrc = $PAIR_NRC_DEFAULT; // " flag for pair from the NRC data base "
		itriplet = $TRIPLET_DEFAULT; // " flag for triplet production "
		iprdst = $IPRDST_DEFAULT; // " pair angular sampling "

		RHOF = 1.0;
		for (int i = 0; i < 5; i++) {
			iausfl[i] = 1;
		}
		for (int i = 5; i < $MXAUS; i++) {
			iausfl[i] = 0;
		}
		ximax = $EXACT_BCA_XIMAX;
		estepe = $MAX_ELOSS;
		skindepth_for_bca = $SKIN_DEPTH_FOR_BCA;
		transport_algorithm = $TRANSPORT_ALGORITHM_DEFAULT;
		bca_algorithm = $BCA_ALGORITHM_DEFAULT;
		exact_bca = $EXACT_BCA_DEFAULT;
		spin_effects = $SPIN_EFFECTS_DEFAULT;
		count_pII_steps = 0.;
		count_all_steps = 0.;

		NMED = $default_nmed;
		DUNIT = 1.0;
		rng_seed = 999999;
		LATCHI = 0;

		RM = 0.5110034;
		RMT2 = 2 * RM;
		RMSQ = RM * RM;
		// pi = 4*datan(1d0); twopi = 2*pi; pi5d2 = 2.5*pi;

		PI = Math.PI;
		TWOPI = 2.0 * PI;
		PI5D2 = 2.5 * PI;
		nbr_split = 1;
		i_play_RR = 0;
		i_survived_RR = 0;
		prob_RR = -1;
		n_RR_warning = 0;

		RandomUse = 1;// turn on default EGS4 random generator numbers!!!
		ranluxB = true;// ranlux use!!!

		inEMF = false;// NO EM field
		i_EMF = 0;// if it is and no futher config-> the default one
		// createOutputFile=false;//default no output file
		// nf.setMinimumFractionDigits(3);//default is 2!!
		// nf.setMaximumFractionDigits(3);//default is 2!!
		idigits = 3;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
		radc_flag = 0;// default no radiative Compton corr
		// "  For now we exclude such corrections by default. They can be        "
		// "  included by adding the file rad_compton.mortran to the list of     "
		// "  files used to build EGSnrc just before egsnrc.mortran              "
		// "  The reason is that there is a fairly large amount of data needed   "
		// "  and this would be wasteful if the effect turns out to be small     "
		iprint = 2;// 0=no printing,1=results only(eg from WATCH), 2=summary,
					// 3=print all
		seqStr = "";
		hatchindex = 0;// default no user input in hatch
		nextB = false;// no resample by default
		USER_CONTROLS_TSTEP_RECURSION = 0;// default no user tstep control
											// override
		iraycorr = 0;// default no raycorr override
		isemfp = 0;// default no $SELECT-ELECTRON-MFP override
		iurd = 0;// default no user-range-discard override
		ispmfp = 0;// default no $SELECT-PHOTON-MFP override
		iGeom = 0;
	}

	// ******************RANDOM from
	// EGS**************************************************
	/**
	 * The EGSnrc implementation for RNG (RandomNumberGenerator).
	 * @return the random number
	 */
	public static double RANDOMSET() {
		double result = 0.0;
		if (ranluxB) {
			if (rng_seed > 24) {
				ranlux(rng_array);
				rng_seed = 1;
			}

			result = rng_array[rng_seed - 1];
			rng_seed = rng_seed + 1;
		} else// ranmar
		{
			if (rng_seed > $NRANMAR)
				ranmar_get();
			result = rng_array1[rng_seed - 1] * twom24;
			rng_seed = rng_seed + 1;
		}

		return result;
		// " i.e. take the rng_seed'th random number from rng_array, "
		// " if all 24 numbers used, generate a new set of 24        "
	}

	/**
	 * Internal use by RANDOMSET.
	 */
	public static void ranmar_get() {
		int i = 0;
		int iopt = 0;
		if (rng_seed == 999999)
			init_ranmar();
		for (i = 1; i <= $NRANMAR; i++) {
			iopt = urndm[ixx - 1] - urndm[jxx - 1];
			if (iopt < 0)
				iopt = iopt + 16777216;
			urndm[ixx - 1] = iopt;
			ixx = ixx - 1;
			jxx = jxx - 1;
			if (ixx == 0) {
				ixx = 97;
			} else if (jxx == 0) {
				jxx = 97;
			}
			crndm = crndm - cdrndm;
			if (crndm < 0)
				crndm = crndm + cmrndm;
			iopt = iopt - crndm;
			if (iopt < 0)
				iopt = iopt + 16777216;
			rng_array1[i - 1] = iopt;
		}
		rng_seed = 1;
		return;
	}

	/**
	 * Internal use by ranmar_get.
	 */
	public static void init_ranmar() {
		int s = 0;
		int t = 0;
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int m = 0;
		int ii = 0;
		int jj = 0;

		if (ixx <= 0 || ixx > 31328)
			ixx = 1802; // "Sets Marsaglia default"
		if (jxx <= 0 || jxx > 30081)
			jxx = 9373; // "sets Marsaglia default"

		i = ixx / 177;
		i = i % 177;
		i = i + 2;// i = mod(ixx/177,177) + 2;
		j = ixx % 177;
		j = j + 2;// j = mod(ixx, 177) + 2;
		k = jxx / 169;
		k = k % 178;
		k = k + 1;// k = mod(jxx/169,178) + 1;
		l = jxx % 169;// l = mod(jxx, 169) ;

		for (ii = 1; ii <= 97; ii++) {

			s = 0;
			t = 8388608;// "t is 2**23 i.e. half of the maximum allowed"
						// "(note that only 24 bits are used)          "

			for (jj = 1; jj <= 24; jj++) {

				// "The if( fool_optimizer ...) statements below are"
				// "to prevent re-arangement of statements for high "
				// "level optimizations and thus different sequences"
				// "on different architectures                      "
				m = i * j;
				m = m % 179;
				m = m * k;
				m = m % 179;
				// m = mod(mod(i*j,179)*k,179);
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				i = j;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				j = k;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				k = m;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				l = 53 * l + 1;
				l = l % 169;// l = mod(53*l+1,169);
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				int lm = l * m;
				lm = lm % 64;
				if (lm >= 32)
					s = s + t;// IF(mod(l*m,64) >= 32) s = s + t;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				t = t / 2;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
			}
			urndm[ii - 1] = s;
		}

		crndm = 362436;
		cdrndm = 7654321;
		cmrndm = 16777213;

		twom24 = 1. / 16777216.;

		ixx = 97;
		jxx = 33;

		rng_seed = $NRANMAR + 1;

		return;
	}

	/**
	 * Internal use by RANDOMSET
	 * @param rng_array rng_array
	 */
	public static void ranlux(double[] rng_array) {

		if (not_initialized) {
			not_initialized = false;
			nskipRnd = nskipll[$DEFAULT_LL];
			twom24 = 1.0;
			twop24 = 1.0;
			jseed = jseed_dflt;
			for (jRnd = 1; jRnd <= 24; jRnd++) {
				twom24 = twom24 * 0.5;
				twop24 = twop24 * 2.0;
				kRnd = jseed / 53668;
				jseed = 40014 * (jseed - kRnd * 53668) - kRnd * 12211;
				if (jseed < 0) {
					jseed = jseed + icon;
				}
				seeds[jRnd - 1] = jseed % 16777216;// mod(jseed,16777216);
				next[jRnd - 1] = jRnd - 1;
			}
			next[0] = 24;// next(1) = 24;
			i24 = 24;
			j24 = 10;
			carry = 0;
			// if( seeds(24) = 0 )
			if (seeds[23] == 0) {
				carry = 1;
			}
		}

		for (jRnd = 1; jRnd <= 24; jRnd++) // j=1,24
		{
			uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;// seeds(j24) -
															// seeds(i24) -
															// carry;
			if (uni < 0) {
				uni = uni + 16777216;
				carry = 1;
			} else {
				carry = 0;
			}
			seeds[i24 - 1] = uni;// seeds(i24) = uni;
			// "IF( uni = 0 ) [ uni = twom24*twom24; ]"
			i24 = next[i24 - 1];// i24 = next(i24);
			j24 = next[j24 - 1];// j24 = next(j24);
			if (uni >= 4096) {
				rng_array[jRnd - 1] = uni * twom24;// rng_array(j) = uni*twom24;
			} else {
				// rng_array(j) = uni*twom24 + seeds(j24)*twom24*twom24;
				rng_array[jRnd - 1] = uni * twom24 + seeds[j24 - 1] * twom24
						* twom24;
			}
		}

		if (nskipRnd > 0) {
			for (jRnd = 1; jRnd <= nskipRnd; jRnd++)// DO jRnd=1,nskipRnd
			{
				// uni = seeds(j24) - seeds(i24) - carry;
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0) {
					uni = uni + 16777216;
					carry = 1;
				} else {
					carry = 0;
				}
				seeds[i24 - 1] = uni;// seeds(i24) = uni;
				i24 = next[i24 - 1];// i24 = next[i24];
				j24 = next[j24 - 1];// j24 = next[j24];
			}
		}
		return;
	}

	/**
	 * Initialize ranlux.
	 * @param luxury_level1 luxury_level1
	 * @param seedin1 seedin1
	 */
	public static void init_ranlux(int luxury_level1, int seedin1) {
		seedin = seedin1;// $$$$$$$$$$$
		luxury_level = luxury_level1;// $$$$$$$$$$$

		jseed = seedin;
		if (jseed <= 0)
			jseed = jseed_dflt;
		if ((luxury_level < 0) || (luxury_level > 4)) {
			luxury_level = $DEFAULT_LL;
		}

		nskipRnd = nskipll[luxury_level];

		// OUTPUT luxury_level,jseed;
		// System.out.println(" ** RANLUX initialization **"+
		// " luxury level: "+luxury_level+
		// " initial seed: "+jseed+
		// "*******");

		not_initialized = false;
		twom24 = 1;
		twop24 = 1;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			twom24 = twom24 * 0.5;
			twop24 = twop24 * 2.0;
			kRnd = jseed / 53668;
			jseed = 40014 * (jseed - kRnd * 53668) - kRnd * 12211;
			if (jseed < 0) {
				jseed = jseed + icon;
			}
			seeds[jRnd - 1] = jseed % 16777216;// mod(jseed,16777216);
			next[jRnd - 1] = jRnd - 1;
		}
		next[0] = 24;
		i24 = 24;
		j24 = 10;
		carry = 0;
		if (seeds[23] == 0) {
			carry = 1;
		}

		return;
	}

	/**
	 * Get ranlux state based on seeds. It is stored in input state variable.
	 * @param state state
	 */
	public static void get_ranlux_state(int[] state) {

		for (jRnd = 1; jRnd <= 24; jRnd++) {
			state[jRnd - 1] = seeds[jRnd - 1];
		}
		state[24] = i24 + 100 * (j24 + 100 * nskipRnd);
		if (carry > 0)
			state[24] = -state[24];
		return;
	}

	/**
	 * Set ranlux status.
	 * @param state state
	 */
	public void set_ranlux_state(int[] state)// public static void
												// set_ranlux_state(int[] state)
	{
		twom24 = 1;
		twop24 = 1;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			twom24 = twom24 * 0.5;
			twop24 = twop24 * 2;
			next[jRnd - 1] = jRnd - 1;
		}
		next[0] = 24;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			seeds[jRnd - 1] = state[jRnd - 1];
		}
		if (state[24] <= 0) {
			status = -state[24];
			carry = 1;
		} else {
			status = state[24];
			carry = 0;
		}
		nskipRnd = status / 10000;
		status = status - nskipRnd * 10000;
		j24 = status / 100;
		i24 = status - 100 * j24;
		if (((j24 < 1) || (j24 > 24)) || ((i24 < 1) || (i24 > 24))) {
			// stop;
			STOPPROGRAM = true;

			seqStr = "*** Error in set_ranlux_state: seeds outside of allowed range!"
					+ " \n"
					+ "   status = "
					+ state[24]
					+ " \n"
					+ "   nskip = "
					+ nskipRnd
					+ " \n"
					+ "   i24 = "
					+ i24
					+ " \n" + "   j24 = " + j24;// +" \n";
			// if(iprint>1)
			eq.printSequence(seqStr);

			return;
		}
		not_initialized = false;
		return;
	}

	/**
	 * Show summary of Ranlux seeds
	 */
	public static void show_ranlux_seeds() {
		if (carry > 0) {
			icarry = 1;
		} else {
			icarry = 0;
		}
		seqStr = "' skip = " + format(nskipRnd, 4) + ",  ix jx = "
				+ format(i24, 3) + " ," + format(j24, 3) + " carry = "
				+ format(icarry, 2);
		if (iprint > 1)
			eq.printSequence(seqStr);

	}

	/**
	 * Show summary of RNG state.
	 */
	public static void SHOW_RNG_STATE() {
		if (ranluxB) {
			show_ranlux_seeds();
		} else// ranmar
		{
			seqStr = "  ixx jxx = " + format(ixx, 4) + " ," + format(jxx, 4);
			if (iprint > 1)
				eq.printSequence(seqStr);
		}
	}

	// ;" ************************* end of ranlux.mortran  *************************"
	// input System.currentTimeMillis() at start of event->e.g. simulation
	// output : time elapsed in String format
	/**
	 * Compute and print the time elapsed from start simulation time.
	 * @return the result
	 */
	public static String timeElapsed()// public static String timeElapsed()
	{
		String times = "";
		long currentTime = System.currentTimeMillis();
		// tipul int spre deosebire de long--->aprox 8-9 zile suporta ca numar
		// de sec!!
		int delta = (new Long(currentTime - startSimulationTime)).intValue();// time
																				// elapsed
																				// in
																				// milliseconds
		int sec = delta / 1000;// impartire intreaga->catul!!!
		int milis = delta % 1000;// restul impartirii intregi!!
		if (sec > 60) {
			int min = sec / 60;
			sec = sec % 60;
			if (min > 60) {
				int h = min / 60;
				min = min % 60;
				if (h > 24) {
					int z = h / 24;
					h = h % 24;
					times = z + " days " + h + " h, " + min + " min, " + sec
							+ " sec, " + milis + " milis";
				} else {
					times = h + " h, " + min + " min, " + sec + " sec, "
							+ milis + " milis";
				}
			} else {
				times = min + " min, " + sec + " sec, " + milis + " milis";
			}
		} else {
			times = sec + " sec, " + milis + " milis";
		}

		seqStr = "******************************************" + " \n"
				+ "Simulation time elapsed: " + times;// + " \n";
		if (iprint > 0)// show sim time
			eq.printSequence(seqStr);

		return times;
	}

	/**
	 * Compute and print the time elapsed from given start time.
	 * @param startTime startTime
	 * @return the result
	 */
	public static String timeElapsed(long startTime) {
		String times = "";
		long currentTime = System.currentTimeMillis();
		// tipul int spre deosebire de long--->aprox 8-9 zile suporta ca numar
		// de sec!!
		int delta = (new Long(currentTime - startTime)).intValue();// time
																	// elapsed
																	// in
																	// milliseconds
		int sec = delta / 1000;// impartire intreaga->catul!!!
		int milis = delta % 1000;// restul impartirii intregi!!
		if (sec > 60) {
			int min = sec / 60;
			sec = sec % 60;
			if (min > 60) {
				int h = min / 60;
				min = min % 60;
				if (h > 24) {
					int z = h / 24;
					h = h % 24;
					times = z + " days " + h + " h, " + min + " min, " + sec
							+ " sec, " + milis + " milis";
				} else {
					times = h + " h, " + min + " min, " + sec + " sec, "
							+ milis + " milis";
				}
			} else {
				times = min + " min, " + sec + " sec, " + milis + " milis";
			}
		} else {
			times = sec + " sec, " + milis + " milis";
		}

		return times;
	}

	/**
	 * Compute the time elapsed from given start time.
	 * @param startTime startTime
	 * @return the result
	 */
	public static String timeElapsedShort(long startTime) {
		String times = "";
		long currentTime = System.currentTimeMillis();
		// tipul int spre deosebire de long--->aprox 8-9 zile suporta ca numar
		// de sec!!
		int delta = (new Long(currentTime - startTime)).intValue();// time
																	// elapsed
																	// in
																	// milliseconds
		int sec = delta / 1000;// impartire intreaga->catul!!!
		int milis = delta % 1000;// restul impartirii intregi!!
		if (sec > 60) {
			int min = sec / 60;
			sec = sec % 60;
			if (min > 60) {
				int h = min / 60;
				min = min % 60;

				times = h + ":" + min + ":" + sec + "," + milis;
			} else {
				int h = min / 60;
				times = h + ":" + min + ":" + sec + "," + milis;
			}
		} else {
			int min = sec / 60;
			int h = min / 60;
			times = h + ":" + min + ":" + sec + "," + milis;
		}

		return times;
	}

	// *********************************AUX=WATCH+SIGMA***************************************
	// "**********************************************************************"
	// "                                                                      "
	// "  These are auxilliary routines used in many NRC user-codes.          "
	// "                                                                      "
	// "                                                                      "
	// "  If you include them via your .configuration file, your user code    "
	// "  will need to define $MXDATA and $STAT.                              "
	// "  eg  REPLACE {$MXDATA} WITH {1}; REPLACE{$STAT} WITH {1};            "
	// "      if you are not using the statistical analysis package           "
	// "                                                                      "
	// "                                                                      "
	// "**********************************************************************"
	/**
	 * Set the number of significant digits.
	 * @param n n
	 */
	public static void setDigits(int n) {
		idigits = n;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
	}

	// returns a string representing the number format (true) or scientific
	// format (false)
	// of a number d. If numbers of letters are smaller than offset, the
	// differences are
	// filled with blank characters.
	/**
	 * Returns a string representing the number format (true) or scientific format (false) 
	 * of a number d. If numbers of letters are smaller than offset, the differences are 
	 * filled with blank characters.
	 * @param d d
	 * @param offset offset
	 * @param nrformat nrformat
	 * @return the result
	 */
	public static String format(double d, int offset, boolean nrformat) {
		String result = "";
		/*
		 * String f = ""; if (nrformat) { f = nf.format(d); } else f =
		 * nff.format(d); int k = f.length(); if (k <= offset) { int i = offset
		 * - k; char[] chs = new char[i]; for (int j = 0; j < i; j++) chs[j] =
		 * ' ';// blank // f=f+new String(chs);//after f = new String(chs) +
		 * f;// before } else {// 1 character blank by default //
		 * f=f+" ";//after f = " " + f;// before }
		 * 
		 * result = f;
		 */

		String fmt = "";// +offset;
		if (nrformat) {
			fmt = "%" + offset + "." + idigits + "f";// "%20.5f"
		} else {
			fmt = "%" + offset + "." + idigits + "E";// "%20.5f"
		}
		result = String.format(fmt, d);

		return result;
	}

	/**
	 * Format using integer d format.
	 * @param d d
	 * @param offset offset
	 * @return the result
	 */
	public static String format(int d, int offset) {
		String result = "";
		/*
		 * String f = "" + d; int k = f.length(); if (k <= offset) { int i =
		 * offset - k; char[] chs = new char[i]; for (int j = 0; j < i; j++)
		 * chs[j] = ' ';// blank // f=f+new String(chs);//after f = new
		 * String(chs) + f;// before } else {// 1 character blank by default //
		 * f=f+" ";//after f = " " + f;// before }
		 * 
		 * result = f;
		 */

		String fmt = "%" + offset + "d";
		result = String.format(fmt, d);

		return result;
	}

	/**
	 * Format using float f format.
	 * @param d d
	 * @param offset offset
	 * @return the result
	 */
	public static String format(double d, int offset) {

		String result = "";
		int k = 0;
		String f = "";
		f = nf.format(d);

		double dbl = stringToDouble(f);
		/*
		 * if (dbl == 0.0) { f = nff.format(d); k = f.length(); if (k <= offset)
		 * { int i = offset - k; char[] chs = new char[i]; for (int j = 0; j <
		 * i; j++) chs[j] = ' ';// blank // f=f+new String(chs);//after f = new
		 * String(chs) + f;// before return f; } else { // f=f+" ";//after f =
		 * " " + f;// before return f; } }
		 * 
		 * k = f.length(); if (k <= offset) { int i = offset - k; char[] chs =
		 * new char[i]; for (int j = 0; j < i; j++) chs[j] = ' ';// blank //
		 * f=f+new String(chs);//after f = new String(chs) + f;// before return
		 * f; } else { f = nff.format(d); k = f.length(); if (k <= offset) { int
		 * i = offset - k; char[] chs = new char[i]; for (int j = 0; j < i; j++)
		 * chs[j] = ' ';// blank // f=f+new String(chs);//after f = new
		 * String(chs) + f;// before return f; } else { // f=f+" ";//after f =
		 * " " + f;// before } }
		 * 
		 * result = f;
		 */

		String fmt = "";// +offset;
		if (dbl == 0.0) {
			fmt = "%" + offset + "." + idigits + "E";// "%20.5f"
		} else {
			k = f.length();
			if (k <= offset) {
				fmt = "%" + offset + "." + idigits + "f";// "%20.5f"
			} else {
				fmt = "%" + offset + "." + idigits + "E";// "%20.5f"
			}
		}

		// result=new Format(fmt).form(d);
		result = String.format(fmt, d);

		return result;
	}

	/**
	 * Same as format but with 1 significant digit.
	 * @param d d
	 * @param offset offset
	 * @param nrformat nrformat
	 * @return the result
	 */
	public static String format1(double d, int offset, boolean nrformat) {
		String result = "";
		/*
		 * String f = ""; NumberFormat nf1 =
		 * NumberFormat.getInstance(Locale.US);
		 * nf1.setMinimumFractionDigits(1);// 3);//default is 2!!
		 * nf1.setMaximumFractionDigits(1);// (3);//default is 2!!
		 * nf1.setGroupingUsed(false);// no 4,568.02 but 4568.02
		 * 
		 * if (nrformat) { f = nf1.format(d); } else f = nff.format(d); int k =
		 * f.length(); if (k <= offset) { int i = offset - k; char[] chs = new
		 * char[i]; for (int j = 0; j < i; j++) chs[j] = ' ';// blank // f=f+new
		 * String(chs);//after f = new String(chs) + f;// before } // else //
		 * {//1 character blank by default // f=f+" ";//after //
		 * f=" "+f;//before // }
		 * 
		 * result = f;
		 */

		String fmt = "";// +offset;
		if (nrformat) {
			fmt = "%" + offset + "." + 1 + "f";// "%20.5f"
		} else {
			fmt = "%" + offset + "." + 1 + "E";// "%20.5f"
		}
		// result=new Format(fmt).form(d);
		result = String.format(fmt, d);// new Format(fmt).form(d);

		return result;
	}

	/**
	 * Format using String s format.
	 * @param d d
	 * @param offset offset
	 * @return the result
	 */
	public static String format(String d, int offset) {
		String result = "";
		/*
		 * String f = "" + d; int k = f.length(); if (k <= offset) { int i =
		 * offset - k; char[] chs = new char[i]; for (int j = 0; j < i; j++)
		 * chs[j] = ' ';// blank // f=f+new String(chs);//after f = new
		 * String(chs) + f;// before
		 * 
		 * } else { // f=f+" ";//after f = " " + f;// before
		 * 
		 * }
		 * 
		 * result = f;
		 */

		String fmt = "%" + offset + "s";
		// result=new Format(fmt).form(d);
		result = String.format(fmt, d);

		return result;
	}

	// private static void CNTOUT(int N, String s)
	/**
	 * Called by WATCH
	 * @param N N
	 * @param s s
	 */
	private static void CNTOUT(int N, String s) {
		// OUTPUT
		// {P1},KE,IQ({P1}),IR({P1}),X({P1}),Y({P1}),Z({P1}),U({P1}),V({P1}),
		// W({P1}),LATCH({P1}),WT({P1});
		// ({P2},I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
		ICOUNT = ICOUNT + 1;

		seqStr = s + format(N, 5) + format(KE, 9, true) + format(IQ[N - 1], 4)
				+ format(IR[N - 1], 4) + format(X[N - 1], 8, true)
				+ format(Y[N - 1], 8, true) + format(Z[N - 1], 8, true)
				+ format(U[N - 1], 7, true) + format(V[N - 1], 7, true)
				+ format(W[N - 1], 7, true) + format(LATCH[N - 1], 10)
				+ format(WT[N - 1], 10, false);// + " \n";
		if (iprint > 0)
			eq.printSequence(seqStr);

	}

	// "*****************************************************************************"
	// "                                                                             "
	// "                        WATCH                                                "
	// "                                                                             "
	/**
	 * A general purpose auxiliary routine. It prints out information about the particle transport: <br>
	 * For IWATCH = 1 it prints information about each discrete interaction <br>
	 * For IWATCH = 2 or 3 it prints information about each step as well <br>
	 * For IWATCH = 4 it prints graphing data<br>
	 * 
	 * Routine is used via two mandatory and 1 optional call from the user's code. <br>
	 * 1)The routine must be initialized by a call with IARG=-99 before the first 
	 * call to SHOWER. It should be after all inputs are in place.<br>
	 * 2)The routine must be called near the beginning of the AUSGAB subroutine 
	 * IF (IWATCH greatere than 0 ) CALL WATCH(IARG,IWATCH); <br>
	 * 3)The routine may be called at the end of each history with IARG = - 1 so 
	 * a message will get printed stated history is complete.<br>
	 * 
	 *
	 * 
	 * @param IARG IARG
	 * @param IWATCH IWATCH
	 */
	public static void WATCH(int IARG, int IWATCH)// public static void
													// WATCH(int IARG,int
													// IWATCH)
	{
		// "============================================================================="
		// "                                                                             "
		// " A general purpose auxiliary routine for use with the EGSnrc system
		// "
		// "     It prints out information about the particle transport                  "
		// "                                                                             "
		// "       For IWATCH = 1 it prints information about each discrete interaction  "
		// "       For IWATCH = 2 or 3 it prints information about each step as well     "
		// "       For IWATCH = 4 it prints graphing data for use with EGS_Windows       "
		// "        see http://www.irs.inms.nrc.ca/inms/irs/EGS_Windows/distribution.html"
		// "                                                                             "
		// "                                                                             "
		// "    Routine is used via two mandatory and 1 optional call from the user's    "
		// "          code                                                               "
		// "                                                                             "
		// "   1)The routine must be initialized by a call with IARG=-99 before the first"
		// "          call to SHOWER. It should be after all inputs are in place.        "
		// "   2)The routine must be called near the beginning of the AUSGAB subroutine  "
		// "          IF (IWATCH > 0 ) CALL WATCH(IARG,IWATCH);                          "
		// "   3)The routine may be called at the end of each history with IARG = - 1 so "
		// "          a message will get printed stated history is complete              "
		// "                                                                             "
		// "    Since WATCH cannot output values related to the initial values in a      "
		// "    shower call, it is useful to also put something like the following       "
		// "    immediately prior to the CALL SHOWER stmt                                "
		// "           IF((IWATCH ~= 0) & (IWATCH ~= 4))[                                "
		// "              OUTPUT 1,EIN,IQI,IRI,XI,YI,ZI,UI,VI,WI,LATCHI,WTI;              "
		// "               (/' INITIAL SHOWER VALUES',T36,':',                           "
		// "               I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);                         "
		// "           ]                                                                 "
		// "    Note EIN is the kinetic energy of the incident particle                  "
		// "                                                                             "
		// "                                                                             "
		// "   The routine uses up to 132 columns for output.                            "
		// "                                                                             "
		// "     JAN 1984  GENERALIZED VERSION WITH INITIALIZATION                       "
		// "                              DAVE ROGERS NRCC                               "
		// "     JUN 1987  PUT IN IWATCH = 4 OPTION     AFB                              "
		// "     JUL 1988  COMPATIBLE WITH X-RAY FLUORESCENCE  DWOR                      "
		// "     SEP 1990  ADDED ENERGY OUTPUT TO IWATCH = 4 OPTION     AFB              "
		// "     OCT 1990  UNIX compatible carriage control   DWOR                       "
		// "     JAN 2000  Rewritten to output relaxation particles and also             "
		// "               so some of the output makes more sense BW                     "
		// "                                                                             "
		// "*****************************************************************************"
		// ;Copyright NRC;
		// "define a local macro"

		// $IMPLICIT-NONE;
		// $INTEGER iarg,iwatch,IP,ICOUNT,JHSTRY,J,N;
		// $REAL KE;
		// save ICOUNT,JHSTRY; "if we ever decide not to use static variables"

		// ;COMIN/BOUNDS, STACK,EPCONT,EGS-VARIANCE-REDUCTION,USEFUL/;

		// ICOUNT=0;JHSTRY=1;
		String s = "";
		int ll = 0;
		int N = 0;

		if (IARG == -99) {// "Initialize flags so we will get calls thru AUSGAB"
			for (int J = 1; J <= 29; J++) {
				iausfl[J - 1] = 1;
			}
			// IAUSFL(22),IAUSFL(23),IAUSFL(24)/=0;
			iausfl[21] = 0;
			iausfl[22] = 0;
			iausfl[23] = 0;// uphib,uphia,rayb turned off
		}

		if (IARG == -1) {// "main is assumed to call AUSGAB with IARG=-1 at end of history"
			if (IWATCH == 4) {
				// WRITE(13,:GRAPHICS_FORMAT:) 0,0,0,0.0,0.0,0.0,0.0,JHSTRY;
				// :GRAPHICS_FORMAT:FORMAT(2I4,1X,I6,4G15.8,I12);
				seqStr = format(0, 4) + format(0, 4) + format("", 1)
						+ format(0, 6) + format(0.0, 15) + format(0.0, 15)
						+ format(0.0, 15) + format(0.0, 15)
						+ format(JHSTRY, 12);// +"\n";
				if (iprint > 0)
					eq.printSequence(seqStr);

				JHSTRY = JHSTRY + 1;
			} else {
				// OUTPUT JHSTRY;(' END OF HISTORY',I8,3X,40('*')/);
				seqStr = " END OF HISTORY" + format(JHSTRY, 8) + format("", 3)
						+ "****************************************";// +"\n";
				if (iprint > 0)
					eq.printSequence(seqStr);

				JHSTRY = JHSTRY + 1;
				ICOUNT = ICOUNT + 2;
				return;
			}
		}

		if ((IWATCH != 4)
				&& ((ICOUNT >= 50) || ((ICOUNT == 0) || (IARG == -99)))) {
			// "PRINT HEADER"
			ICOUNT = 1;
			// OUTPUT;(//T39,' NP',3X,'ENERGY Q REGION X',7X,
			// 'Y',7X,'Z',6X,'U',6X,'V',6X,'W',6X,'LATCH',2X,'WEIGHT'/);
			seqStr = format("", 39) + " NP" + format("", 3)
					+ "ENERGY  Q REGION    X" + format("", 7) + "Y"
					+ format("", 7) + "Z" + format("", 6) + "U" + format("", 6)
					+ "V" + format("", 6) + "W" + format("", 6) + "LATCH"
					+ format("", 2) + "WEIGHT";// +"\n";
			if (iprint > 0)
				eq.printSequence(seqStr);

		}

		if ((IWATCH == 4) && (IARG >= 0) && (IARG != 5)) {// "GRAPHICS OUTPUT"
															// WRITE(13,:GRAPHICS_FORMAT:)
															// NP,IQ(NP),IR(NP),X(NP),Y(NP),Z(NP),E(NP);
															// :GRAPHICS_FORMAT:FORMAT(2I4,1X,I6,4G15.8,I12);
			seqStr = format(NP, 4) + format(IQ[NP - 1], 4) + format("", 1)
					+ format(IR[NP - 1], 6) + format(X[NP - 1], 15)
					+ format(Y[NP - 1], 15) + format(Z[NP - 1], 15)
					+ format(E[NP - 1], 15);// +"\n";
			if (iprint > 0)
				eq.printSequence(seqStr);

		}

		if ((IARG == 5) || (IARG < 0))
			return;
		if (IWATCH == 4)
			return; // "NONE OF THE REST NEEDED FOR GRAPHICS OUTPUT"

		KE = E[NP - 1];
		if (IQ[NP - 1] != 0) {
			KE = E[NP - 1] - 0.511;
		}

		if (IARG == 0 && IWATCH == 2) {
			s = format("", 11) + "STEP ABOUT TO OCCUR";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (T11,'STEP ABOUT TO OCCUR', T36,':');
		} else if (IARG == 0) {
			return;
		}

		if (IARG == 1) {
			s = " Discard  AE,AP<E<ECUT";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Discard AE,AP<E<ECUT',T36,':');
		} else if (IARG == 2) {
			s = " Discard  E<AE,AP";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Discard E<AE,AP',T36,':');
		} else if (IARG == 3) {
			s = " Discard -user request";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Discard -user request',T36,':');
		} else if (IARG == 4) {
			s = "Local energy deposition";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);

			seqStr = format("", 10) + s + format(EDEP, 12, true)
					+ " MeV in region " + format(IR[NP - 1], 6);// +"\n";
			if (iprint > 0)
				eq.printSequence(seqStr);

			// OUTPUT EDEP,IR(NP);
			// (T10,'Local energy deposition',T36,':',F12.5,' MeV in region
			// ',I6);
		}

		else if (IARG == 6) {
			s = " bremsstrahlung  about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' bremsstrahlung about to occur',T36,':');
		} else if (IARG == 7) {
			if (nbr_split == 1) {// "no splitting or SBS is on in BEAMnrc"
				for (int IP = NPold; IP <= NP; IP++) {
					if (IQ[IP - 1] == -1) {
						KE = E[IP - 1] - RM;
						s = format("", 10) + "Resulting electron";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T10,'Resulting electron',T36,':');
					} else {
						KE = E[IP - 1];
						s = format("", 10) + "Resulting photon";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T10,'Resulting photon',T36,':');
					}
				}
			} else {// "splitting case--e- is always at NPold"
				KE = E[NPold - 1] - RM;
				s = format("", 10) + "Resulting electron";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NPold, s);// (T10,'Resulting electron',T36,':');
				for (int IP = NPold + 1; IP <= NP; IP++) {
					KE = E[IP - 1];
					if (IP == NPold + 1) {// "print info for first one only"
						s = format("", 10) + "Split photons";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T10,'Split photons',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
			}// " end of splitting block"
		}

		else if (IARG == 8) {
			s = " Moller   about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Moller about to occur',T36,':');
		} else if (IARG == 9) {
			// "surely this logic not needed?"
			if (NP == NPold) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else {
				for (int IP = NPold; IP <= NP; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[NP - 1]) * RM;
					if (IP == NPold) {
						s = format("", 11) + "Resulting electrons";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'Resulting electrons',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
			}
		}

		else if (IARG == 10) {
			s = " Bhabba   about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Bhabba about to occur',T36,':');
		} else if (IARG == 11) {
			if (NP == NPold) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else {
				for (int IP = NPold; IP <= NP; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[IP - 1]) * RM;
					if (IP == NPold) {
						s = format("", 11) + "Resulting e- or e+";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'Resulting e- or e+',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
			}
		}

		else if (IARG == 12) {
			s = " Positron about to decay in flight";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Positron about to decay in flight',T36,':');
		} else if (IARG == 13) {
			if (NP == NPold) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else {
				for (int IP = NPold; IP <= NP; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[IP - 1]) * RM;
					if (IP == NPold) {
						s = format("", 11) + "Resulting photons";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'Resulting photons',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
			}
		}

		else if (IARG == 28) {
			s = " Positron will annihilate at rest";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Positron will annihilate at rest',T36,':');
		} else if (IARG == 14) {
			if (NP == NPold) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else {
				for (int IP = NPold; IP <= NP; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[IP - 1]) * RM;
					if (IP == NPold) {
						s = " Positron annihilates at rest";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (' Positron annihilates at
										// rest',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
			}
		}

		else if (IARG == 15) {
			s = " Pair production about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Pair production about to occur',T36,':');
		} else if (IARG == 16) {// "after pair production"
			if ((NP == NPold) && (i_survived_RR == 0)) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else if ((NP == NPold) && (i_survived_RR > 0)) {// "we have cleared the stack"
																// OUTPUT
																// i_survived_rr,prob_rr;
																// (T10,'Russian
																// Roulette
																// eliminated
																// ',I2,
																// ' particle(s)
																// with
																// probability
																// ',F8.5)
				seqStr = format("", 10) + "Russian Roulette eliminated "
						+ format(i_survived_RR, 2)
						+ " particle(s) with probability "
						+ format(prob_RR, 8, true);// +"\n";
				if (iprint > 0)
					eq.printSequence(seqStr);

				s = format("", 10) + "Now on top of stack";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T10,'Now on top of stack',T36,':');
			} else {
				for (int IP = NPold; IP <= NP; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[IP - 1]) * RM;
					if (IP == NPold) {
						s = format("", 11) + "Resulting pair";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'Resulting pair',T36,':');
					} else {
						s = format(":", 36);
						CNTOUT(IP, s);// (T36,':');
					}
				}
				if (i_survived_RR > 0) {
					// OUTPUT i_survived_rr,prob_rr;
					// (T10,'Russian Roulette eliminated ',I2,'
					// particle(s) with probability ',F8.5);
					seqStr = format("", 10) + "Russian Roulette eliminated "
							+ format(i_survived_RR, 2)
							+ " particle(s) with probability "
							+ format(prob_RR, 8, true);// +"\n";
					if (iprint > 0)
						eq.printSequence(seqStr);

					s = format("", 10) + "Now on top of stack";
					ll = s.length();
					ll = 36 - ll;
					s = s + format(":", ll);
					CNTOUT(NP, s);// (T10,'Now on top of stack',T36,':');
				}
			}
		}

		else if (IARG == 17) {
			s = " Compton  about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Compton about to occur',T36,':');
		} else if (IARG == 18) {// "after call to COMPT"
			if ((NP == NPold) && (i_survived_RR == 0)) {
				s = format("", 11) + "Interaction rejected";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T11,'Interaction rejected',T36,':');
			} else if (NP > NPold) {// "have not cleared the stack with rus rou"
				for (int IP = NPold; IP <= NPold + 1; IP++) {
					KE = E[IP - 1] - Math.abs(IQ[IP - 1]) * RM;
					if (IQ[IP - 1] != 0) {
						s = format("", 11) + "compton electron created";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'compton electron
										// created',T36,':');
					} else {
						s = format("", 11) + "compton scattered photon";
						ll = s.length();
						ll = 36 - ll;
						s = s + format(":", ll);
						CNTOUT(IP, s);// (T11,'compton scattered
										// photon',T36,':');
					}
				}
			}
			if (i_survived_RR > 0) {// "whether the stack has been cleared or not"
									// OUTPUT i_survived_rr,prob_rr;
									// (T10,'Russian Roulette eliminated ',I2,
									// ' particle(s) with probability ',F8.5)
				seqStr = format("", 10) + "Russian Roulette eliminated "
						+ format(i_survived_RR, 2)
						+ " particle(s) with probability "
						+ format(prob_RR, 8, true);// +"\n";
				if (iprint > 0)
					eq.printSequence(seqStr);

				s = format("", 10) + "Now on top of stack";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T10,'Now on top of stack',T36,':');
			}
		}

		else if (IARG == 19) {
			s = " Photoelectric about to occur";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Photoelectric about to occur',T36,':');
		} else if (IARG == 20) {
			if ((NPold == NP) && ((IQ[NP - 1] == 0) && (i_survived_RR == 0))) {
				String s0 = format("", 11) + "Photon energy below N-shell"
						+ "\n";
				s = format("", 11) + "Photon discarded";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				s = s0 + s;

				CNTOUT(NP, s);// (T11,'Photon energy below N-shell',/,
				// T11,'Photon discarded',T36,':');
			} else if ((IQ[NPold - 1] == -1) && (i_survived_RR == 0)) {
				KE = E[NPold - 1] - RM;
				s = format("", 10) + "Resulting photoelectron";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NPold, s);// (T10,'Resulting photoelectron',T36,':');
			} else if (i_survived_RR > 0) {// "done some russian roulette"
				if ((NP == NPold - 1) || (IQ[NPold - 1] != -1)) {
					if (i_survived_RR > 1) {// "eliminated more than the photoelectron"
											// OUTPUT i_survived_rr-1,prob_rr;
											// (T10,'Russian Roulette eliminated
											// ',I4,
											// ' particle(s) with probability
											// ',F8.5,' plus');
						seqStr = format("", 10)
								+ "Russian Roulette eliminated "
								+ format(i_survived_RR - 1, 2)
								+ " particle(s) with probability "
								+ format(prob_RR, 8, true);// +"\n";
						if (iprint > 0)
							eq.printSequence(seqStr);

					}
					// OUTPUT prob_rr;
					// (T10,'Russian Roulette eliminated resulting
					// photoelectron',
					// ' with probability ',F8.5);
					seqStr = format("", 10)
							+ "Russian Roulette eliminated resulting photoelectron"
							+ " with probability " + format(prob_RR, 8, true);// +"\n";
					if (iprint > 0)
						eq.printSequence(seqStr);

				} else {// "NPold could hold the photoelectron"
					KE = E[NPold - 1] - RM;
					s = format("", 10) + "Resulting photoelectron?";
					ll = s.length();
					ll = 36 - ll;
					s = s + format(":", ll);
					CNTOUT(NPold, s);// (T10,'Resulting
										// photoelectron?',T36,':');

					// OUTPUT i_survived_rr,prob_rr;
					// (T10,'Russian Roulette eliminated ',I4,
					// ' particle(s) with probability ',F8.5);
					seqStr = format("", 10) + "Russian Roulette eliminated "
							+ format(i_survived_RR, 2)
							+ " particle(s) with probability "
							+ format(prob_RR, 8, true);// +"\n";
					if (iprint > 0)
						eq.printSequence(seqStr);

				}
				s = format("", 10) + "Now on top of stack";
				ll = s.length();
				ll = 36 - ll;
				s = s + format(":", ll);
				CNTOUT(NP, s);// (T10,'Now on top of stack',T36,':');
			}
		}

		else if (IARG == 24) {
			s = " Rayleigh scattering occured";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(NP, s);// (' Rayleigh scattering occured',T36,':');
		}

		else if (IARG == 25) {
			s = format("", 10) + "Fluorescent X-ray created";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);

			CNTOUT(NP, s);// (T10,'Fluorescent X-ray created',T36,':');
		}

		else if (IARG == 26) {
			s = format("", 10) + "Coster-Kronig e- created";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);

			CNTOUT(NP, s);// (T10,'Coster-Kronig e- created',T36,':');
		}

		else if (IARG == 27) {
			s = format("", 10) + "Auger electron created";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);

			CNTOUT(NP, s);// (T10,'Auger electron created',T36,':');
		}

		if ((IARG == 0) && (IWATCH == 2)) {
			// OUTPUT USTEP,TUSTEP,VSTEP,TVSTEP,EDEP;
			// (T5,'USTEP,TUSTEP,VSTEP,TVSTEP,EDEP',T36,': ',5(1PE13.4));
			s = format("", 5) + "USTEP,TUSTEP,VSTEP,TVSTEP,EDEP";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);

			seqStr = s + format(USTEP, 13) + format(TUSTEP, 13)
					+ format(VSTEP, 13) + format(TVSTEP, 13) + format(EDEP, 13);// +"\n";
			if (iprint > 0)
				eq.printSequence(seqStr);

			ICOUNT = ICOUNT + 1;
		}

		if ((NP == 1) || (IARG == 0))
			return;
		if (IARG <= 3) {
			N = NP - 1;
			KE = E[N - 1] - Math.abs(IQ[N - 1]) * RM;
			s = format("", 10) + "Now on top of stack";
			ll = s.length();
			ll = 36 - ll;
			s = s + format(":", ll);
			CNTOUT(N, s);// (T10,'Now on top of stack',T36,':');
		}

		return;
	}

	// "*******************************************************************************
	// "
	// "
	// " *****************
	// " * *
	// " * SIGMA.MORTRAN *
	// " * *
	// " *****************
	// "
	// "
	// " SIGMA IS A STATISTICAL ANALYSIS ROUTINE DESIGNED TO BE USED BY EGS
	// " USER PROGRAMS TO GIVE THE TOTALS OR AVERAGES AND THEIR UNCERTAINTIES
	// " OF THE DATA CALCULATED BY THE MONTE CARLO CODE.
	// " THE UNCERTAINTIES ARE RETURNED AS PERCENTS.
	// "
	// " VARIABLES
	// " =========
	// "
	// " DATA(NDATA,ISTAT) THE TWO DIMENSIONAL ARRAY OF DATA TO BE
	// " ANALYZED. ISTAT IS THE NUMBER OF STATISTICAL
	// " BATCHES AND NDATA IS THE NUMBER OF ERRORS TO
	// " BE CALCULATED. AFTER THE END OF THE CALCULATION,
	// " DATA(N,1) CONTAINS THE TOTAL OR AVERAGE AND
	// " DATA(N,2) CONTAINS THE ERROR. NDATA SHOULD
	// " BE < OR = $MAXDATA AND ISTAT SHOULD BE < OR =
	// " $STAT WHCH MUST BE DEFINED IN THE MAIN ROUTINE.
	// " Note $STAT must be 2 or greater, even if istat=1
	// "
	// " MODE = 0 => ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS
	// " IGNORED. (eg. STOPPING POWER RATIO)
	// " = 1 => ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS NOT
	// " IGNORED. (e.g. DOSE)
	// " = 2 => ANALYSIS ON TOTAL VALUES (eg. TOTAL EDEP)
	// "
	// " IERR = 0 => NORMAL COMPLETION.
	// " = 1 => WARNING: MODE OUT OF RANGE, DEFAULTED TO 0
	// " = 10 => ERROR: ONLY ONE BATCH INPUT, QUICK CALCULATION
	// " DONE. ERROR=99.9%
	// " = 11 => ERROR: NO NON-ZERO DATA FOUND IN A GIVEN SET,
	// " ERROR=99.9%
	// " = -1 => FATAL ERROR: NDATA OR ISTAT OUT OF RANGE, NO
	// " CALCULATION DONE.
	// "
	// "
	// " VERSION 1 A.F.B. 83/7/22
	// " Version 2 IK Jan 6 6000 implemented implicit none
	// "
	// "*******************************************************************************

	/**
	 * SIGMA IS A STATISTICAL ANALYSIS ROUTINE DESIGNED TO BE USED BY EGS USER PROGRAMS TO GIVE THE TOTALS OR 
	 * AVERAGES AND THEIR UNCERTAINTIES OF THE DATA CALCULATED BY THE MONTE CARLO CODE. THE UNCERTAINTIES ARE RETURNED AS PERCENTS. <br>
	 * DATA(NDATA,ISTAT) THE TWO DIMENSIONAL ARRAY OF DATA TO BE ANALYZED. ISTAT IS THE NUMBER OF STATISTICAL BATCHES AND NDATA IS 
	 * THE NUMBER OF ERRORS TO BE CALCULATED. AFTER THE END OF THE CALCULATION, DATA(N,1) CONTAINS THE TOTAL OR AVERAGE AND 
	 * DATA(N,2) CONTAINS THE ERROR. NDATA SHOULD BE LESS OR EQUAL THAN MAXDATA AND ISTAT SHOULD BE LESS OR EQUAL THAN STAT WHCH MUST BE DEFINED IN THE MAIN ROUTINE. <br>
	 * MODE = 0 MEANS ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS IGNORED. (eg. STOPPING POWER RATIO); MODE =1 ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS NOT 
	 * IGNORED. (e.g. DOSE); MODE = 2 ANALYSIS ON TOTAL VALUES (eg. TOTAL EDEP).<br>
	 * 
	 * 
	 * @param NDATA NDATA
	 * @param ISTAT ISTAT
	 * @param MODE MODE
	 */
	public static void SIGMA(int NDATA, int ISTAT, int MODE)// ,int IERR)
	{

		// ;Copyright NRC;

		// $INTEGER NDATA,ISTAT,MODE,IERR;

		// REPLACE {;COMIN/ERROR/;} WITH {
		// ;COMMON/ERROR/DATA($MXDATA,$STAT);
		// $REAL data;
		// }
		// double[][] DATA=new DATA[$MXDATA][$STAT];
		// ;COMIN/ERROR/;

		int N = 0;
		int NON0 = 0;
		int I = 0;
		double STAT = 0.0;
		double SDENOM = 0.0;
		// double emax = 0.0;
		double AVG = 0.0;
		double ERROR = 0.0;
		double DATUM = 0.0;
		double ARGMNT = 0.0;
		// "It is a good idea to use double precision"
		// "in cases with very low stat. uncertainties"

		double EMAX = 99.9;

		IERRsgm = 0; // "ASSUME NORMAL COMPLETION"
		boolean transferb = false;

		// "TEST INPUTS AND SET ERROR CODES AND RETURN IF NEEDED."

		if ((MODE < 0) || (MODE > 2)) {
			MODE = 2;
			IERRsgm = 1;
		}

		if (((NDATA <= 0) || (NDATA > $MXDATA))
				|| ((ISTAT <= 0) || (ISTAT > $STAT))) {
			IERRsgm = -1;
			return;// "FATAL INPUT ERROR, RETURN IMMEDIATELY"
		}
		if (ISTAT == 1) {
			IERRsgm = 10;// "ONLY ONE STATISTICAL BATCH, QUICK CALCULATION"
			// DO N=1,NDATA[DATA(N,2)=EMAX;]
			for (N = 1; N <= NDATA; N++) {
				// DATA(N,2)=EMAX;
				DATA[N - 1][1] = EMAX;
			}
			return;
		}

		// "MOST ANOMALIES HAVE BEEN HANDLED. NOW DO THE ANALYSIS"

		if (MODE != 0) {
			STAT = ISTAT;
			SDENOM = STAT * (STAT - 1.);
		}
		for (N = 1; N <= NDATA; N++) {
			transferb = false;
			NON0 = 0;// "NON-ZERO COUNTER"
			AVG = 0.0;
			ERROR = 0.0;
			for (I = 1; I <= ISTAT; I++) {
				DATUM = DATA[N - 1][I - 1];// DATA(N,I);
				if (DATUM != 0.0) {
					NON0 = NON0 + 1;
					AVG = AVG + DATUM;
					ERROR = ERROR + DATUM * DATUM;
				}
			}
			if (NON0 == 0) {
				IERRsgm = 11;
				ERROR = EMAX;
				// GOTO :TRANSFER:;//"NO NON-ZERO DATA "
				transferb = true;

			} else if ((NON0 == 1) && (MODE == 0)) {
				ERROR = EMAX;
				// GOTO:TRANSFER:;//"ONLY ONE DATUM"
				transferb = true;
			} else {
				if (MODE == 0) {
					STAT = NON0;
					SDENOM = STAT * (STAT - 1.);
				}
			}

			if (!transferb) {
				AVG = AVG / STAT;
				ARGMNT = ERROR - STAT * AVG * AVG;
				// "FLAG -VE SQUARE ROOTS THAT CAN ONLY OCCUR DUE TO ROUND-OFF ERRORS"
				if (ARGMNT < 0.0) {
					// OUTPUT ARGMNT,ERROR,STAT,AVG,SDENOM;
					// (' ***** - SQ RT IN SIGMA.
					// ARGMNT,ERROR,STAT,AVG,SDENOM='/' ',5E12.4);
					ARGMNT = 0.0;
				}
				ERROR = Math.sqrt(ARGMNT / SDENOM);

				if (AVG == 0.) {
					ERROR = EMAX;
				} else {
					ERROR = 100. * ERROR / Math.abs(AVG);
				}

				if (MODE == 2)
					AVG = AVG * STAT;
			}

			// :TRANSFER:;
			// DATA(N,1)=AVG;DATA(N,2)=MIN(EMAX,ERROR);
			DATA[N - 1][0] = AVG;
			DATA[N - 1][1] = Math.min(EMAX, ERROR);
		}// "END OF NDATA LOOP"
		return;
	}// "END OF SIGMA"

	// "*****************************************************************************
	// " The following are routines that implement
	// " the alias sampling technique for sampling from a histogram
	// " distribution coded for use with EGSnrc
	// "
	// " subroutine prepare_alias_sampling
	// " function alias_sample
	// "
	// " I. Kawrakow, January 2000
	// "
	// "*****************************************************************************

	/**
	 * Other useful routines that implement the alias sampling technique for sampling from a histogram.
	 * @param nsbin nsbin, number of bins in the histogram
	 * @param fs_array fs_array, bin probabilities
	 * @param ws_array ws_array, output: alias table ready for sampling
	 * @param ibin_array ibin_array, outpu: for alias table ready for sampling
	 */
	public static void prepare_alias_sampling(int nsbin, double[] fs_array,
			double[] ws_array, int[] ibin_array) {
		// "====================================================================
		// "
		// " inputs: nsbin: number of bins in the histogram
		// " fs_array: bin probabilities
		// "
		// " Note that we don't need the bin limits at this point, they
		// " are needed for the actual sampling (in alias_sample)
		// "
		// " outputs: ws_array, ibin_array: alias table ready for sampling
		// "
		// "====================================================================

		// ;Copyright NRC;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL fs_array(nsbin),ws_array(nsbin);

		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;

		boolean exit1b = false;
		boolean exit2b = false;

		sum = 0;
		for (i = 1; i <= nsbin; i++) {
			if (fs_array[i - 1] < 1.e-30) {
				fs_array[i - 1] = 1.e-30;
			}
			ws_array[i - 1] = -fs_array[i - 1];
			ibin_array[i - 1] = 1;
			sum = sum + fs_array[i - 1];
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			exit1b = false;
			exit2b = false;
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (Math.abs(ws_array[j_h - 1]) > sum) {
						// GOTO :AT_EXIT1:;
						exit1b = true;
						break;
					}
				}
			}
			if (!exit1b)
				j_h = nsbin;
			// :AT_EXIT1:

			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (Math.abs(ws_array[j_l - 1]) < sum) {
						// GOTO :AT_EXIT2:;
						exit2b = true;
						break;
					}
				}
			}
			if (!exit2b)
				j_l = nsbin;
			// :AT_EXIT2:

			aux = sum - Math.abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;

			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1;
			}

		}

		return;
	}

	/**
	 * Samples from an alias table which must have been prepared using prepare_alias_sampling.
	 * @param nsbin nsbin
	 * @param xs_array xs_array
	 * @param ws_array ws_array
	 * @param ibin_array ibin_array
	 * @return the result
	 */
	public static double alias_sample(int nsbin, double[] xs_array,// )
			double[] ws_array, int[] ibin_array) {
		// "===============================================================
		// "
		// " samples from an alias table which must have been prepared
		// " using prepare_alias_table
		// "
		// "===============================================================

		// ;Copyright NRC;
		// implicit none;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL xs_array(0:nsbin),ws_array(nsbin);

		// ;COMIN/RANDOM/;
		double alias_sample = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		double aj = 0.0;
		int j = 0;

		v1 = random01();
		v2 = random01();
		aj = 1.0 + v1 * nsbin;
		Double dbl = new Double(aj);
		j = dbl.intValue();// minim1 maxim nsbin!!
		if (j > nsbin)
			j = nsbin; // " this happens only if $RANDOMSET produces
						// " numbers in (0,1]--------> is not the DEFAULT
						// case!!!
		aj = aj - j;
		if (aj > ws_array[j - 1]) {
			j = ibin_array[j - 1];
		}
		alias_sample = (1.0 - v2) * xs_array[j - 1] + v2 * xs_array[j];// ok, xs
																		// =0
																		// biased
		return alias_sample;
	}

	// ********************************AUX****************************************************
	/**
	 * Initialize data for radiative Compton corrections.
	 */
	public static void RADC_HATCH() {
		int ne = 0;
		int nu = 0;
		int ixtot = 1;
		int ibtot = 1;
		int nx = 0;
		int nbox = 0;
		// int ny1 = 0;
		// int ny2 = 0;
		// int nw = 0;
		// int jh = 0;
		// int jl = 0;
		int lgle = 0;
		int icheck = 0;
		double aux1 = 0.0;
		double gmfp = 0.0;
		double gmfp_old = 0.0;
		double gbr1 = 0.0;
		double gbr1_old = 0.0;
		double gbr2 = 0.0;
		double gbr2_old = 0.0;
		double gle = 0.0;
		double frad = 0.0;
		double acheck = 0.0;
		double acheck1 = 0.0;
		// double sum = 0.0;
		double aux = 0.0;

		int iread = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		int lnr = 0;

		int isigs = 0;
		boolean sigsB = false;
		boolean frejB = false;
		int jfrej = 0;
		boolean sigdB = false;
		boolean multiB = false;
		boolean nxB = false;
		boolean nboxB = false;
		int jx = 0;
		boolean xB = false;
		int jset = 0;
		boolean setB = false;
		boolean maxB = false;
		boolean fdatB = false;
		boolean binsB = false;
		boolean ixmin1B = false;
		boolean ixmax1B = false;
		boolean ixmin2B = false;
		boolean ixmax2B = false;
		boolean ixmin3B = false;
		boolean ixmax3B = false;
		boolean ixmin4B = false;
		boolean ixmax4B = false;

		// *************Check if radiative corrections are turned on
		// *******************
		if (radc_flag == 0) {
			// $egs_info(*,'Radiative Compton corrections not requested');
			seqStr = " Radiative Compton corrections not requested ";// + " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			return;
		}
		// ************ Inform the user that radiative corrections are being
		// used ******
		// $egs_info(*,' ');
		// $egs_info(*,'Radiative Compton corrections requested:');
		// $egs_info('(a,$)',' Reading radiative Compton data ...');
		seqStr = " Radiative Compton corrections requested: ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		seqStr = "    Reading radiative Compton data ...";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		// ************ Construct the name of the data file
		// ****************************
		String filename = datas + file_sep + egsDatas + file_sep + radcompton
				+ defaultext;
		// ************ Read the data
		// **************************************************
		FileInputStream in = null;
		try {
			//FileInputStream in = new FileInputStream(filename);
			in = new FileInputStream(filename);
			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if (lnr == 1)// read radc_emin
						{
							String s = desc.toString();
							radc_emin = stringToDouble(s);// System.out.println(radc_emin);
						}
						if (lnr == 2)// read radc_emax
						{
							String s = desc.toString();
							radc_emax = stringToDouble(s);// System.out.println(radc_emax);
						}
						if (lnr == 3)// read radc_dw
						{
							String s = desc.toString();
							radc_dw = stringToDouble(s);// System.out.println(radc_dw);
						}
						if (lnr == 4)// read ne
						{
							String s = desc.toString();
							ne = stringToInt(s);// System.out.println(ne);
						}
						if (lnr == 5)// read nu
						{
							String s = desc.toString();
							nu = stringToInt(s);// System.out.println(nu);

							if ((ne != $RADC_NE) || (nu != $RADC_NU)) {
								// $egs_fatal(*,'radc_init: inconsistent data');
								STOPPROGRAM = true;

								seqStr = "radc_init: inconsistent data";// +" \n";
								// if(iprint>1)
								eq.printSequence(seqStr);

								return;
							}
							radc_dle = Math.log(radc_emax / radc_emin) / ne;
							radc_dlei = 1.0 / radc_dle;
						}
						if (lnr > 5) {
							// read(radc_unit,*,err=:radc-read-error:)
							// (radc_sigs(i),i=0,ne);
							if (!sigsB) {
								String s = desc.toString();
								radc_sigs[isigs] = stringToDouble(s);// System.out.println(radc_sigs[isigs]);
								isigs++;
								if (isigs > ne) {
									sigsB = true;
									isigs = 0;// reset contor
								}
							}
							// DO i=0,ne [
							// read(radc_unit,*,err=:radc-read-error:)
							// (radc_frej(i,j),j=0,nu);
							// ]
							else if (!frejB) {
								String s = desc.toString();
								radc_frej[isigs][jfrej] = stringToDouble(s);// System.out.println(" i "+isigs+" j "+jfrej+" :  "+radc_frej[isigs][jfrej]);
								jfrej++;
								if (jfrej > nu) {
									jfrej = 0;
									isigs++;
								}
								if (isigs > ne) {
									frejB = true;
									isigs = 0;// reset contor
								}
							}
							// read(radc_unit,*,err=:radc-read-error:)
							// (radc_sigd(i),i=0,ne);
							else if (!sigdB) {
								String s = desc.toString();
								radc_sigd[isigs] = stringToDouble(s);// System.out.println(radc_sigd[isigs]);
								isigs++;
								if (isigs > ne) {
									sigdB = true;
									isigs = 0;// reset contor
									ixtot = 1;
									ibtot = 1;
								}
							} else if (!multiB) {
								// read(radc_unit,*) nx,nbox;
								if (!nxB) {
									String s = desc.toString();
									nx = stringToInt(s);// System.out.println(nx);
									nxB = true;
								} else if (!nboxB) {
									String s = desc.toString();
									nbox = stringToInt(s);// System.out.println(nbox);
									nboxB = true;
									if (ixtot - 1 + nx > $RADC_NX) {
										// $egs_fatal(*,'Incosistent number of
										// box boundaries');
										STOPPROGRAM = true;
										seqStr = "radc_init: Incosistent number of box boundaries";// +" \n";
										// if(iprint>1)
										eq.printSequence(seqStr);

										return;

									}
									if (ibtot - 1 + nbox > $RADC_NBOX) {
										// $egs_fatal(*,'Incosistent number of
										// boxes');
										STOPPROGRAM = true;
										seqStr = "radc_init: Incosistent number of boxes";// +" \n";
										// if(iprint>1)
										eq.printSequence(seqStr);

										return;
									}
									radc_startx[isigs] = ixtot;
									radc_startb[isigs] = ibtot;
									jx = ixtot;
									jset = ibtot;
								}
								// read(radc_unit,*,err=:radc-read-error:)
								// (radc_x(j),j=ixtot,ixtot-1+nx);
								else if (!xB) {
									String s = desc.toString();
									radc_x[jx - 1] = stringToDouble(s);// System.out.println(radc_x[jx-1]);
									jx++;
									if (jx > ixtot - 1 + nx) {
										xB = true;
									}

								} else if (!setB) {
									if (!maxB)// fdatBbinsBixmin1Bixmax1B
									{
										String s = desc.toString();
										radc_Smax[jset - 1] = stringToDouble(s);
										maxB = true;
									} else if (!fdatB) {
										String s = desc.toString();
										radc_fdat[jset - 1] = stringToDouble(s);
										fdatB = true;
									} else if (!binsB) {
										String s = desc.toString();
										radc_bins[jset - 1] = stringToInt(s);
										binsB = true;
									} else if (!ixmin1B) {
										String s = desc.toString();
										radc_ixmin1[jset - 1] = stringToInt(s);
										ixmin1B = true;
									} else if (!ixmax1B) {
										String s = desc.toString();
										radc_ixmax1[jset - 1] = stringToInt(s);
										ixmax1B = true;
									} else if (!ixmin2B) {
										String s = desc.toString();
										radc_ixmin2[jset - 1] = stringToInt(s);
										ixmin2B = true;
									} else if (!ixmax2B) {
										String s = desc.toString();
										radc_ixmax2[jset - 1] = stringToInt(s);
										ixmax2B = true;
									} else if (!ixmin3B) {
										String s = desc.toString();
										radc_ixmin3[jset - 1] = stringToInt(s);
										ixmin3B = true;
									} else if (!ixmax3B) {
										String s = desc.toString();
										radc_ixmax3[jset - 1] = stringToInt(s);
										ixmax3B = true;
									} else if (!ixmin4B) {
										String s = desc.toString();
										radc_ixmin4[jset - 1] = stringToInt(s);
										ixmin4B = true;
									} else if (!ixmax4B) {
										String s = desc.toString();
										radc_ixmax4[jset - 1] = stringToInt(s);
										ixmax4B = true;

										// System.out.println(radc_Smax[jset-1]+" "+radc_fdat[jset-1]+" "+radc_bins[jset-1]+
										// " "+radc_ixmin1[jset-1]+" "+radc_ixmax1[jset-1]+" "+radc_ixmin2[jset-1]+" "+radc_ixmax2[jset-1]+
										// " "+radc_ixmin3[jset-1]+" "+radc_ixmax3[jset-1]+" "+radc_ixmin4[jset-1]+" "+radc_ixmax4[jset-1]);

										// reset flags
										maxB = false;
										fdatB = false;
										binsB = false;
										ixmin1B = false;
										ixmax1B = false;
										ixmin2B = false;
										ixmax2B = false;
										ixmin3B = false;
										ixmax3B = false;
										ixmin4B = false;
										ixmax4B = false;
										//
										jset++;
										if (jset > ibtot - 1 + nbox) {
											setB = true;
											ixtot = ixtot + nx;
											ibtot = ibtot + nbox;
											// restart
											setB = false;
											xB = false;
											nboxB = false;
											nxB = false;
											//
											isigs++;
											if (isigs > ne) {
												multiB = true;
												isigs = 0;// reset contor
											}
										}

									}

								}// if(!setB)
							}// if(!multiB)

						}

					}// have data
					desc = new StringBuffer();
				}// else
			}// main while
			//in.close();
		}// try
		catch (Exception exc) {

		}
		finally {
	        if( null != in ) {
	            try 
	            {
	                in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
		}
		seqStr = " OK !";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		radc_le1 = -Math.log(radc_emin * 0.5110034) * radc_dlei;
		// "Change the cross section and branching ratios"
		seqStr = "    Modifying cross sections and branching ratios ...";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		for (int medium = 1; medium <= NMED; medium++) {
			// "write(17,*) '************* Doing medium ',medium;
			for (int i = 1; i <= MGE[medium - 1]; i++)// mge(medium) [
			{
				gle = (i - GE0[medium - 1]) / GE1[medium - 1];
				lgle = i;
				// $EVALUATE gmfp USING gmfp(gle);
				gmfp = GMFP1[lgle - 1][medium - 1] * gle
						+ GMFP0[lgle - 1][medium - 1];
				// $EVALUATE gbr1 USING gbr1(gle);
				gbr1 = GBR11[lgle - 1][medium - 1] * gle
						+ GBR10[lgle - 1][medium - 1];
				// $EVALUATE gbr2 USING gbr2(gle);
				gbr2 = GBR21[lgle - 1][medium - 1] * gle
						+ GBR20[lgle - 1][medium - 1];

				if (Math.exp(gle) > radc_emin * 0.5110034) {
					acheck = radc_dlei * gle + radc_le1;
					Double dbl = new Double(acheck);
					icheck = dbl.intValue();
					acheck = acheck - icheck;
					acheck1 = 1.0 - acheck;
					frad = (radc_sigs[icheck] + radc_sigd[icheck]) * acheck1
							+ (radc_sigs[icheck + 1] + radc_sigd[icheck + 1])
							* acheck;
				} else {
					frad = 1.0;
				}
				// "write(17,*) i,gle,exp(gle),gmfp,gbr1,gbr2,frad;
				aux = 1.0 / (1.0 + (gbr2 - gbr1) * (frad - 1.0));
				aux1 = (gbr2 - gbr1) * frad;
				gmfp = gmfp * aux;
				gbr1 = gbr1 * aux;
				gbr2 = gbr1 + aux1 * aux;
				// "write(17,*) 'New: ',gmfp,gbr1,gbr2;
				if (i > 1) {
					// "write(17,*) 'Old gmfp:
					// ',gmfp1(i-1,medium),gmfp0(i-1,medium);
					GMFP1[i - 2][medium - 1] = (gmfp - gmfp_old)
							* GE1[medium - 1];
					GMFP0[i - 2][medium - 1] = gmfp - GMFP1[i - 2][medium - 1]
							* gle;
					// "write(17,*) 'New gmfp:
					// ',gmfp1(i-1,medium),gmfp0(i-1,medium);
					// "write(17,*) 'Old gbr1:
					// ',gbr11(i-1,medium),gbr10(i-1,medium);
					GBR11[i - 2][medium - 1] = (gbr1 - gbr1_old)
							* GE1[medium - 1];
					GBR10[i - 2][medium - 1] = gbr1 - GBR11[i - 2][medium - 1]
							* gle;
					// "write(17,*) 'New gbr1:
					// ',gbr11(i-1,medium),gbr10(i-1,medium);
					// "write(17,*) 'Old gbr2:
					// ',gbr21(i-1,medium),gbr20(i-1,medium);
					GBR21[i - 2][medium - 1] = (gbr2 - gbr2_old)
							* GE1[medium - 1];
					GBR20[i - 2][medium - 1] = gbr2 - GBR21[i - 2][medium - 1]
							* gle;
					// "write(17,*) 'New gbr2:
					// ',gbr21(i-1,medium),gbr20(i-1,medium);
				}
				gmfp_old = gmfp;
				gbr1_old = gbr1;
				gbr2_old = gbr2;
			}
			GMFP1[MGE[medium - 1] - 1][medium - 1] = GMFP1[MGE[medium - 1] - 2][medium - 1];
			GMFP0[MGE[medium - 1] - 1][medium - 1] = GMFP0[MGE[medium - 1] - 2][medium - 1];
			GBR11[MGE[medium - 1] - 1][medium - 1] = GBR11[MGE[medium - 1] - 2][medium - 1];
			GBR10[MGE[medium - 1] - 1][medium - 1] = GBR10[MGE[medium - 1] - 2][medium - 1];
			GBR21[MGE[medium - 1] - 1][medium - 1] = GBR21[MGE[medium - 1] - 2][medium - 1];
			GBR20[MGE[medium - 1] - 1][medium - 1] = GBR20[MGE[medium - 1] - 2][medium - 1];
		}

		seqStr = " OK !";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

	}

	/**
	 * Initialize data for triplet production process.
	 */
	public static void init_triplet() {
		// implicit none;
		// $declare_max_medium;
		// $COMIN-INIT-TRIPLET;

		double[] energies = new double[$N_TRIPLET_DATA];
		double[][] sig_pair = new double[$N_ELEMENT][$N_TRIPLET_DATA];
		double[][] sig_triplet = new double[$N_ELEMENT][$N_TRIPLET_DATA];
		double[] f_triplet = new double[$N_TRIPLET_DATA];
		double[] sigp = new double[$N_TRIPLET_DATA];
		double[] sigt = new double[$N_TRIPLET_DATA];
		double[] as = new double[$N_TRIPLET_DATA];
		double[] bs = new double[$N_TRIPLET_DATA];
		double[] cs = new double[$N_TRIPLET_DATA];
		double[] ds = new double[$N_TRIPLET_DATA];
		double logE = 0.0;
		double f_new = 0.0;
		double f_old = 0.0;
		// double spline = 0.0;

		int ntrip = 0;
		int ifirst = 0;
		int imed = 0;
		int iz1 = 0;
		int izi = 0;

		char lineSep = '\n';
		String mark = "#";
		int lnr = 0;
		int iread = 0;
		boolean haveData = false;
		StringBuffer desc = new StringBuffer();
		int i = 0;
		boolean enB = false;
		boolean sigpB = false;
		boolean sigtB = false;
		int iel = 0;

		if (itriplet == 0)
			return;
		String filename = datas + file_sep + egsDatas + file_sep + tripletfile
				+ defaultext;

		seqStr = " init_triplet: reading triplet data ... ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);
		FileInputStream in = null;
		try {
			//FileInputStream in = new FileInputStream(filename);
			in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						// read(triplet_unit,*) ntrip;
						if (lnr == 1)// read radc_emin
						{
							String s = desc.toString();
							ntrip = stringToInt(s);// System.out.println(ntrip);
							if (ntrip > $N_TRIPLET_DATA) {
								// $egs_fatal(*,'Max. number of data points per
								// element is ',$N_TRIPLET_DATA);
								STOPPROGRAM = true;

								seqStr = "Max. number of data points per element is "
										+ $N_TRIPLET_DATA;// +" \n";
								// if(iprint>1)
								eq.printSequence(seqStr);

								return;
							}

						}

						if (lnr > 1) {
							String ss = desc.toString();
							// read(triplet_unit,*,err=:error_triplet_data:)
							// (energies(i),i=1,ntrip);
							if (!enB) {
								String s = desc.toString();
								energies[i] = stringToDouble(s);// System.out.println(energies[i]);
								i++;
								if (i == ntrip) {
									enB = true;
									i = 0;// reset contor
								}
							} else if (ss.compareTo(mark) == 0) {
								// System.out.println("y");
								boolean b = true;
								@SuppressWarnings("unused")
								int itst = 0;
								while (b) {
									itst++;
									int indx = in.read();
									if ((char) indx == lineSep) {
										// System.out.println(itst);
										itst = 0;
										break;
									}
								}
							} else if (!sigpB) {
								String s = desc.toString();
								sig_pair[iel][i] = stringToDouble(s);// System.out.println(sig_pair[iel][i]);
								i++;
								if (i == ntrip) {
									sigpB = true;
									i = 0;// reset contor
								}

							} else if (!sigtB) {
								String s = desc.toString();
								sig_triplet[iel][i] = stringToDouble(s);// System.out.println(sig_triplet[iel][i]);
								i++;
								if (i == ntrip) {
									sigtB = true;
									i = 0;// reset contor
									// -----------
									iel++;
									sigpB = false;
									sigtB = false;
									if (iel == $N_ELEMENT)
										break;
								}

							}

						}

					}// have data
					desc = new StringBuffer();
				}// else
			}// main while
			//in.close();
		}// try
		catch (Exception exc) {

		}
		finally {
	        if( null != in ) {
	            try 
	            {
	                in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
	    }
		seqStr = " OK ";// +" \n";
		if (iprint > 1)
			eq.printSequence(seqStr);

		ifirst = 0;
		for (i = 1; i <= ntrip; i++) {
			if ((ifirst == 0) && (energies[i - 1] > 4.01 * RM))
				ifirst = i;
			energies[i - 1] = Math.log(energies[i - 1]);
		}
		log_4rm = Math.log(4 * RM);
		energies[ifirst - 2] = log_4rm;// energies[ifirst-1] = log_4rm;
		dl_triplet = (energies[ntrip - 1] - log_4rm) / $MAX_TRIPLET;
		dli_triplet = 1.0 / dl_triplet;
		bli_triplet = 1.0 - log_4rm / dl_triplet;

		for (imed = 1; imed <= NMED; imed++) {

			seqStr = "   Preparing triplet fraction data for medium " + imed;// +
																				// " \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

			Double iz1d = new Double(ZELEM[imed - 1][0] + 0.1);
			iz1 = iz1d.intValue();// ZELEM[imed][1] + 0.1;
			for (i = 1; i <= ntrip; i++)// DO i=1,ntrip [
			{
				sigp[i - 1] = PZ[imed - 1][0] * sig_pair[iz1 - 1][i - 1];
				sigt[i - 1] = PZ[imed - 1][0] * sig_triplet[iz1 - 1][i - 1];
				for (iel = 2; iel <= NNE[imed - 1]; iel++) // DO iel=2,nne(imed)
															// [
				{
					Double izid = new Double(ZELEM[imed - 1][iel - 1] + 0.1);
					izi = izid.intValue();// zelem(imed,iel) + 0.1;
					sigp[i - 1] = sigp[i - 1] + PZ[imed - 1][iel - 1]
							* sig_pair[izi - 1][i - 1];
					sigt[i - 1] = sigt[i - 1] + PZ[imed - 1][iel - 1]
							* sig_triplet[izi - 1][i - 1];
				}
			}

			for (i = ifirst; i <= ntrip; i++) {
				// f_triplet(i-ifirst+2) = sigt[i-1]/(sigp(i) + sigt(i));

				f_triplet[i - ifirst + 1] = sigt[i - 1]
						/ (sigp[i - 1] + sigt[i - 1]);
			}
			f_triplet[0] = 0.0;// f_triplet(1) = 0;
			// #####for now we have ntrip-ifirst+2 elements: restrict the
			// energies array:
			double[] ener = new double[ntrip - ifirst + 2];
			for (i = 0; i < ntrip - ifirst + 2; i++) {
				ener[i] = energies[i + ifirst - 2];
			}
			// #######################################################################
			// set_spline(energies[ifirst-1],f_triplet,as,bs,cs,ds,ntrip-ifirst+2);
			// set_spline(energies,f_triplet,as,bs,cs,ds,ntrip-ifirst+2);
			set_spline(ener, f_triplet, as, bs, cs, ds, ntrip - ifirst + 2);

			logE = log_4rm;
			f_old = 0.0;
			for (i = 1; i <= $MAX_TRIPLET - 1; i++) {
				logE = logE + dl_triplet;
				// f_new =
				// spline(logE,energies(ifirst-1),as,bs,cs,ds,ntrip-ifirst+2);

				// f_new = spline(logE,energies,as,bs,cs,ds,ntrip-ifirst+2);
				f_new = spline(logE, ener, as, bs, cs, ds, ntrip - ifirst + 2);
				// f_new = spline(logE,ntrip-ifirst+2);

				a_triplet[i - 1][imed - 1] = (f_new - f_old) * dli_triplet;
				b_triplet[i - 1][imed - 1] = f_new - a_triplet[i - 1][imed - 1]
						* logE;
				f_old = f_new;
			}

			seqStr = " OK ";// +" \n";
			if (iprint > 1)
				eq.printSequence(seqStr);

		}

		return;

	}

	/**
	 * Reading byte file utility. Convert a 4 byte object from little endian to big endian byte order or vice versa.
	 * @param c c
	 */
	// some reading bytes file utillities
	public static void egs_swap_4(byte[] c) {
		// "============================================================================"
		// " Convert a 4 byte object from little endian to big endian byte order        "
		// " or vice versa                                                              "
		// character c(4),tmp;
		byte tmp;
		tmp = c[3];
		c[3] = c[0];
		c[0] = tmp;
		tmp = c[2];
		c[2] = c[1];
		c[1] = tmp;
	}

	/**
	 * Reading byte file utility. Convert a 2 byte object from little endian to big endian byte order or vice versa.
	 * @param c c
	 */
	public static void egs_swap_2(byte[] c) {
		// "============================================================================"
		// " Convert a 2 byte object from little endian to big endian byte order        "
		// " or vice versa                                                              "
		// character c(4),tmp;
		byte tmp;
		tmp = c[1];
		c[1] = c[0];
		c[0] = tmp;
	}

	// readLong
	// (((long)(a & 0xff) << 56) | ((long)(b & 0xff) << 48) | ((long)(c & 0xff)
	// << 40) | ((long)(d & 0xff) << 32) | ((long)(e & 0xff) << 24) | ((long)(f
	// & 0xff) << 16) | ((long)(g & 0xff) << 8) | ((long)(h & 0xff)))

	/**
	 * Convert a 2 byte array to short.
	 * @param c c
	 * @return the result
	 */
	public static short byte2toShort(byte[] c) {
		int result = ((c[0] & 0xff) << 8) | (c[1] & 0xff);// or:
		// short result=(short)((a << 8) * | (b & 0xff))
		return (short) result;
	}

	/**
	 * Convert a 4 byte array to int.
	 * @param c c
	 * @return the result
	 */
	public static int byte4toInt(byte[] c) {
		// (((a & 0xff) << 24) | ((b & 0xff) << 16) | ((c & 0xff) << 8) | (d &
		// 0xff))
		int result = (((c[0] & 0xff) << 24) | ((c[1] & 0xff) << 16)
				| ((c[2] & 0xff) << 8) | (c[3] & 0xff));
		return result;
	}

	/**
	 * Convert a 4 byte array to float.
	 * @param c c
	 * @return the result
	 */
	public static float byte4toFloat(byte[] c) {
		float result = 0;
		int num = ((c[0] << 24) & 0xFF000000) | ((c[1] << 16) & 0xFF0000)
				| ((c[2] << 8) & 0xFF00) | (c[3] & 0xFF);
		result = Float.intBitsToFloat(num);

		return result;
	}

	// ====================================

	/**
	 * Reads in spin rejection data for multiple elastic scattering and initializes interpolation arrays for the screening parameter, 
	 * elastic cross section, first and second MS moments. Called by mscati.
	 */
	public static void init_spin() {

		// "======================================================================="
		// "                                                                       "
		// " Reads in spin rejection data for multiple elastic scattering and      "
		// " initializes interpolation arrays for the screening parameter,         "
		// " elastic cross section, first and second MS moments                    "
		// "                                                                       "
		// " I. Kawrakow, NRC                                                      "
		// "======================================================================="

		// ; Copyright NRC;
		// implicit none;
		// $declare_max_medium;
		// COMIN/Spin-Data,ELECIN,MEDIA,BREMPR,THRESH,EGS-IO/;
		// "BREMPR is needed for the elemental composition"

		double[][] eta_array = new double[2][$MAXE_SPI1 + 1];
		double[][] c_array = new double[2][$MAXE_SPI1 + 1];
		double[][] g_array = new double[2][$MAXE_SPI1 + 1];
		double[] earray = new double[$MAXE_SPI1 + 1];
		// double[] tmp_array = new double[$MAXE_SPI1 + 1];
		double sum_Z2 = 0.0;
		double sum_Z = 0.0;
		double sum_A = 0.0;
		double sum_pz = 0.0;
		double Z = 0.0;
		double tmp = 0.0;
		double Z23 = 0.0;
		double g_m = 0.0;
		double g_r = 0.0;
		double sig = 0.0;
		double dedx = 0.0;
		double dum1 = 0.0;
		double dum2 = 0.0;
		double dum3 = 0.0;
		double aux_o = 0.0;
		double tau = 0.0;
		double tauc = 0.0;
		double beta2 = 0.0;
		double eta = 0.0;
		double gamma = 0.0;
		double fmax = 0.0;
		double eil = 0.0;
		double e = 0.0;
		double si1e = 0.0;
		double si2e = 0.0;
		double si1p = 0.0;
		double si2p = 0.0;
		double aae = 0.0;
		double etap = 0.0;
		double[] elarray = new double[$MAXE_SPI1 + 1];
		double[] farray = new double[$MAXE_SPI1 + 1];
		double[] af = new double[$MAXE_SPI1 + 1];
		double[] bf = new double[$MAXE_SPI1 + 1];
		double[] cf = new double[$MAXE_SPI1 + 1];
		double[] df = new double[$MAXE_SPI1 + 1];
		// double spline = 0.0;
		double fine = 137.03604;
		double TF_constant = 0.88534138;
		// String endzs = "";
		// int Zi = 0;
		int n_ener = 0;
		// int n_q = 0;
		int n_point = 0;
		int je = 0;
		int ndata = 0;
		int leil = 0;

		double dloge = 0.0;
		double eloge = 0.0;

		float[] fmax_array = new float[$MAXQ_SPIN + 1];// (0:$MAXQ_SPIN);
		short[] i2_array = new short[512];
		short ii2 = 0;
		int ii4 = 0;

		int medium = 0;
		int iq = 0;
		int i = 0;
		int j = 0;
		int k = 0;
		int i_ele = 0;
		// int iii = 0;
		int iZ = 0;
		// int iiZ = 0;
		int neke = 0;

		// ,length,ii4,//irec;
		// character spin_file*256;
		// character*6 string;
		// integer*4 lnblnk1; "used to be lnblnk but changed for compilers that
		// "do not have lnblnk, we supply lnblnk1 in this file

		// $INTEGER spin_unit, rec_length, want_spin_unit;
		// integer egs_get_unit;
		// character data_version*32,endianess*4;
		boolean swap = false;

		// $declare_write_buffer;

		// " stupid GNU compiler complains that the arguments to egs_swap_ are of "
		// " one kind here, but of some other kind there => need character arrays and "
		// " equivalence statements"
		// real*4 tmp_4;
		// character c_2(2), c_4(4);
		// equivalence (ii2,c_2), (tmp_4,c_4);

		byte[] b_4 = new byte[4];
		byte[] b_2 = new byte[2];
		StringBuffer desc = new StringBuffer();
		String little_endian = "1234";

		// " First construct the path to the spin dbase directory "
		String file = datas + file_sep + egsDatas + file_sep + spinfile;
		try {
			FileInputStream data_in = new FileInputStream(file);
			try {
				// read(spin_unit,rec=1,err=:spin-read-error:)
				// data_version,endianess,
				// espin_min,espin_max,b2spin_min,b2spin_max;

				// swap?
				// The fortran interpreter read data in closely relationship
				// with current CPU endianess
				// In Java, all "read" methods is realted to byte by byte
				// reading, so DO NOT CARE
				// ABOUT current CPU endianess BUT machine generated data
				// endianess!!!!!
				// So, if endianess of machine generated data is little_endian
				// ==>MUST swap!!
				data_in.skip(32);// //version=32
				int iread = 0;
				for (int iend = 1; iend <= 4; iend++)// 4 bytes
				{
					iread = data_in.read();// read next byte of data
					desc.append((char) iread);
				}
				String endianess = desc.toString();
				if (endianess.compareTo(little_endian) == 0)
					swap = true;
				else
					swap = false;

				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				espin_min = byte4toFloat(b_4);// System.out.println(espin_min+"  "+swap);
				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				espin_max = byte4toFloat(b_4);// System.out.println(espin_max+"  "+swap);
				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				b2spin_min = byte4toFloat(b_4);// System.out.println(b2spin_min+"  "+swap);
				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				b2spin_max = byte4toFloat(b_4);// System.out.println(b2spin_max+"  "+swap);

				seqStr = "Reading spin data base from  " + file;
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "The endianess of machine generated data is "
						+ endianess + ",   swap:" + swap;
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "Ranges: " + espin_min + ", " + espin_max + ", "
						+ b2spin_min + ", " + b2spin_max;
				if (iprint > 1)
					eq.printSequence(seqStr);

				n_ener = $MAXE_SPIN;
				// n_q = $MAXQ_SPIN;
				n_point = $MAXU_SPIN;
				dloge = Math.log(espin_max / espin_min) / n_ener;
				eloge = Math.log(espin_min);
				earray[0] = espin_min;// earray(0) = espin_min;
				for (i = 1; i <= n_ener; i++) {
					eloge = eloge + dloge;
					earray[i] = Math.exp(eloge);
				}
				dbeta2 = (b2spin_max - b2spin_min) / n_ener;
				beta2 = b2spin_min;
				earray[n_ener + 1] = espin_max;
				for (i = n_ener + 2; i <= 2 * n_ener + 1; i++) {
					beta2 = beta2 + dbeta2;
					if (beta2 < 0.999) {
						earray[i] = 511.0034 * (1.0 / Math.sqrt(1.0 - beta2) - 1.0);
					} else {
						earray[i] = 50585.1;
					}
				}
				// " Convert to MeV and set interpolation interavals"
				// *
				// * IK: Moved here Feb 4 2004.
				// * Bug found and reported by Dr Helmut Schlattl.
				// * /
				espin_min = espin_min / 1000.0;
				espin_max = espin_max / 1000.0;
				dlener = Math.log(espin_max / espin_min) / $MAXE_SPIN;
				dleneri = 1.0 / dlener;
				espml = Math.log(espin_min);
				dbeta2 = (b2spin_max - b2spin_min) / $MAXE_SPIN;
				dbeta2i = 1.0 / dbeta2;
				dqq1 = 0.5 / $MAXQ_SPIN;
				dqq1i = 1.0 / dqq1;

				for (medium = 1; medium <= NMED; medium++) {
					seqStr = "  medium " + medium + " ..................... ";
					if (iprint > 1)
						eq.printSequence(seqStr);

					for (iq = 0; iq <= 1; iq++) {
						for (i = 0; i <= $MAXE_SPI1; i++) {
							// /eta_array(iq,i),c_array(iq,i),g_array(iq,i)/=0;
							eta_array[iq][i] = 0.0;
							c_array[iq][i] = 0.0;
							g_array[iq][i] = 0.0;
							for (j = 0; j <= $MAXQ_SPIN; j++) {
								for (k = 0; k <= $MAXU_SPIN; k++) {
									spin_rej[medium - 1][iq][i][j][k] = 0.0;
								}
							}
						}
					}
					sum_Z2 = 0.0;
					sum_A = 0.0;
					sum_pz = 0.0;
					sum_Z = 0.0;

					for (i_ele = 1; i_ele <= NNE[medium - 1]; i_ele++) {
						Z = ZELEM[medium - 1][i_ele - 1];
						Double dbl = new Double(Z + 0.5);
						iZ = dbl.intValue();// int(Z+0.5);

						tmp = PZ[medium - 1][i_ele - 1] * Z * (Z + 1);
						// "For now, we take into account the contribution of atomic"
						// "electrons to elastic scattering by replacing Z**2 with  "
						// "Z*(Z+1). The part of the scattering power that is taken "
						// "into account by discrete Moller/Bhabha events is        "
						// "substracted below => bc is energy dependent. We will    "
						// "worry about better approaches in the future (a realistic"
						// "inelastic scattering model is needed first)             "
						sum_Z2 = sum_Z2 + tmp;
						sum_Z = sum_Z + PZ[medium - 1][i_ele - 1] * Z;
						sum_A = sum_A + PZ[medium - 1][i_ele - 1]
								* WA[medium - 1][i_ele - 1];
						sum_pz = sum_pz + PZ[medium - 1][i_ele - 1];
						Z23 = Math.pow(Z, 0.6666667);// Z23 = Z**0.6666667;

						// ##############################
						data_in = new FileInputStream(file);// reset
						// based on comparission with old data, bytes must be
						// skipped!!
						int nskip = (iZ - 1) * 70656 + 1104;// 70656=======see
															// next
						data_in.skip(nskip);
						// ###############################
						for (iq = 0; iq <= 1; iq++) {
							for (i = 0; i <= $MAXE_SPI1; i++) {
								// irec = 1 + (iz-1)*4*(n_ener+1) +
								// 2*iq*(n_ener+1) + i+1;
								// $FOOL-INTEL-OPTIMIZER(25) '**** energy
								// ',i,earray(i),irec;
								// read(spin_unit,rec=irec,err=:spin-read-error:)
								// dum1,dum2,dum3,aux_o,fmax_array,i2_array;
								// dum1,2,3,aux forms 4*4 bytes
								// fmax is 4*($MAXQ_SPIN+1) and i2 is 2*512==>
								// nskip_base is:16+4*16+1024=1104. this must be
								// multiplied by:
								// 2 (from iq loop) and by ($MAXE_SPI1+1) from i
								// loop!!=>
								// nskip_base=2*(2*$MAXE_SPIN+1+1)*1104;//$MAXE_SPIN=15=>
								// nskip_base=2*32*1104=70656!!!!!
								// read 4 normal floats
								data_in.read(b_4);
								if (swap)
									egs_swap_4(b_4);
								dum1 = byte4toFloat(b_4);
								// System.out.println(dum1+" dum1 -----iq=    "+iq);
								data_in.read(b_4);
								if (swap)
									egs_swap_4(b_4);
								dum2 = byte4toFloat(b_4);
								// System.out.println(dum2+" dum2 -----iq=    "+iq);
								data_in.read(b_4);
								if (swap)
									egs_swap_4(b_4);
								dum3 = byte4toFloat(b_4);
								// System.out.println(dum3+" dum3 -----iq=    "+iq);
								data_in.read(b_4);
								if (swap)
									egs_swap_4(b_4);
								aux_o = byte4toFloat(b_4);
								// System.out.println(aux_o+" aux_o -----iq=    "+iq);

								// fmax
								for (int ij = 1; ij <= $MAXQ_SPIN + 1; ij++)// 16;ij++)
								{
									data_in.read(b_4);
									if (swap)
										egs_swap_4(b_4);
									fmax_array[ij - 1] = byte4toFloat(b_4);
									// System.out.println(fmax_array[ij-1]+" fmax------iq=    "+iq);
								}
								// i2
								for (int ij = 1; ij <= 512; ij++) {
									data_in.read(b_2);
									if (swap)
										egs_swap_2(b_2);
									i2_array[ij - 1] = byte2toShort(b_2);// 1
																			// biased!!!
									// int ifi2=ifi1;
									// if(ifi2<0)
									// {
									// ifi2=ifi2+ 65536;
									// System.out.println(i2_array[ij-1]+" i2_change------iq=    "+iq);
									// }
									// else
									// {
									// System.out.println(ifi2+" i2------iq=    "+iq);
									// }
								}
								eta_array[iq][i] = eta_array[iq][i] + tmp
										* Math.log(Z23 * aux_o);
								tau = earray[i] / 511.0034; // "energy in the file is in keV"
								beta2 = tau * (tau + 2.0)
										/ ((tau + 1.0) * (tau + 1.0));
								eta = Z23
										/ ((fine * TF_constant) * (fine * TF_constant))
										* aux_o / 4.0 / tau / (tau + 2.0);
								c_array[iq][i] = c_array[iq][i]
										+ tmp
										* (Math.log(1.0 + 1.0 / eta) - 1.0 / (1.0 + eta))
										* dum1 * dum3;
								g_array[iq][i] = g_array[iq][i] + tmp * dum2;
								for (j = 0; j <= $MAXQ_SPIN; j++) {
									// tmp_4 = fmax_array[j];
									for (k = 0; k <= $MAXU_SPIN; k++) {
										// ii2 = i2_array((n_point+1)*j + k+1);
										ii2 = i2_array[(n_point + 1) * j + k];// 1
																				// biased
										ii4 = ii2;
										if (ii4 < 0)
											ii4 = ii4 + 65536;
										dum1 = ii4;
										dum1 = dum1 * fmax_array[j] / 65535.0;// tmp_4/65535.0;
										spin_rej[medium - 1][iq][i][j][k] = spin_rej[medium - 1][iq][i][j][k]
												+ tmp * dum1;
									}
								}

							}// for( i=0;i<=$MAXE_SPI1;i++)
						}// for( iq=0;iq<=1;iq++)
					}// for( i_ele=1;i_ele<=NNE[medium-1];i_ele++)

					// " spin_rej will be used as a rejection function in MS sampling, "
					// " so scale maximum to unity"
					for (iq = 0; iq <= 1; iq++) {
						for (i = 0; i <= $MAXE_SPI1; i++) {
							for (j = 0; j <= $MAXQ_SPIN; j++) {
								fmax = 0;
								for (k = 0; k <= $MAXU_SPIN; k++) {
									if (spin_rej[medium - 1][iq][i][j][k] > fmax) {
										fmax = spin_rej[medium - 1][iq][i][j][k];
									}
								}
								for (k = 0; k <= $MAXU_SPIN; k++) {
									spin_rej[medium - 1][iq][i][j][k] = spin_rej[medium - 1][iq][i][j][k]
											/ fmax;
								}
							}
						}
					}
					// " Process eta_array, c_array and g_array to their final form "
					for (i = 0; i <= $MAXE_SPI1; i++) {
						tau = earray[i] / 511.0034;
						beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
						for (iq = 0; iq <= 1; iq++) {
							aux_o = Math.exp(eta_array[iq][i] / sum_Z2)
									/ ((fine * TF_constant) * (fine * TF_constant));
							eta_array[iq][i] = 0.26112447 * aux_o
									* BLCC[medium - 1] / XCC[medium - 1];
							eta = aux_o / 4.0 / tau / (tau + 2.0);
							gamma = 3.0
									* (1.0 + eta)
									* (Math.log(1.0 + 1.0 / eta)
											* (1.0 + 2.0 * eta) - 2.0)
									/ (Math.log(1.0 + 1.0 / eta) * (1.0 + eta) - 1.0);
							g_array[iq][i] = g_array[iq][i] / sum_Z2 / gamma;
							c_array[iq][i] = c_array[iq][i]
									/ sum_Z2
									/ (Math.log(1.0 + 1.0 / eta) - 1.0 / (1.0 + eta));
						}
					}
					// " Prepare interpolation table for the screening parameter "
					eil = (1 - EKE0[medium - 1]) / EKE1[medium - 1];
					e = Math.exp(eil);
					if (e <= espin_min) {
						si1e = eta_array[0][0];
						si1p = eta_array[1][0];
					} else {
						if (e <= espin_max) {
							aae = (eil - espml) * dleneri;
							Double jed = new Double(aae);
							je = jed.intValue();
							aae = aae - je;
						} else {
							tau = e / 0.5110034;
							beta2 = tau * (tau + 2.0)
									/ ((tau + 1.0) * (tau + 1.0));
							aae = (beta2 - b2spin_min) * dbeta2i;
							Double jed = new Double(aae);
							je = jed.intValue();
							aae = aae - je;
							je = je + $MAXE_SPIN + 1;
						}
						si1e = (1.0 - aae) * eta_array[0][je] + aae
								* eta_array[0][je + 1];
						si1p = (1.0 - aae) * eta_array[1][je] + aae
								* eta_array[1][je + 1];
					}
					neke = MEKE[medium - 1];
					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						e = Math.exp(eil);
						if (e <= espin_min) {
							si2e = eta_array[0][0];
							si2p = eta_array[1][0];
						} else {
							if (e <= espin_max) {
								aae = (eil - espml) * dleneri;
								Double jed1 = new Double(aae);
								je = jed1.intValue();
								aae = aae - je;
							} else {
								tau = e / 0.5110034;
								beta2 = tau * (tau + 2.0)
										/ ((tau + 1.0) * (tau + 1.0));
								aae = (beta2 - b2spin_min) * dbeta2i;
								Double jed1 = new Double(aae);
								je = jed1.intValue();
								aae = aae - je;
								je = je + $MAXE_SPIN + 1;
							}
							si2e = (1.0 - aae) * eta_array[0][je] + aae
									* eta_array[0][je + 1];
							si2p = (1.0 - aae) * eta_array[1][je] + aae
									* eta_array[1][je + 1];
						}

						etae_ms1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						etae_ms0[i - 1][medium - 1] = si2e
								- etae_ms1[i - 1][medium - 1] * eil;
						etap_ms1[i - 1][medium - 1] = (si2p - si1p)
								* EKE1[medium - 1];
						etap_ms0[i - 1][medium - 1] = si2p
								- etap_ms1[i - 1][medium - 1] * eil;

						si1e = si2e;
						si1p = si2p;
					}
					etae_ms1[neke - 1][medium - 1] = etae_ms1[neke - 2][medium - 1];
					etae_ms0[neke - 1][medium - 1] = etae_ms0[neke - 2][medium - 1];
					etap_ms1[neke - 1][medium - 1] = etap_ms1[neke - 2][medium - 1];
					etap_ms0[neke - 1][medium - 1] = etap_ms0[neke - 2][medium - 1];
					// "Prepare correction to the first MS moment due to spin effects"
					// "first electrons"

					for (i = 0; i <= $MAXE_SPIN; i++) {
						elarray[i] = Math.log(earray[i] / 1000.0);
						farray[i] = c_array[0][i];
					}
					for (i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
						elarray[i] = Math.log(earray[i + 1] / 1000.0);
						farray[i] = c_array[0][i + 1];
					}
					ndata = $MAXE_SPI1 + 1;
					if (UE[medium - 1] > 1.e5) {
						elarray[ndata - 1] = Math.log(UE[medium - 1]);
					} else {
						elarray[ndata - 1] = Math.log(1e5);
					}
					farray[ndata - 1] = 1.0;
					set_spline(elarray, farray, af, bf, cf, df, ndata);
					eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
					si1e = spline(eil, elarray, af, bf, cf, df, ndata);

					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						si2e = spline(eil, elarray, af, bf, cf, df, ndata);
						q1ce_ms1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						q1ce_ms0[i - 1][medium - 1] = si2e
								- q1ce_ms1[i - 1][medium - 1] * eil;
						si1e = si2e;
					}
					q1ce_ms1[neke - 1][medium - 1] = q1ce_ms1[neke - 2][medium - 1];
					q1ce_ms0[neke - 1][medium - 1] = q1ce_ms0[neke - 2][medium - 1];
					// "now positrons"

					for (i = 0; i <= $MAXE_SPIN; i++) {
						farray[i] = c_array[1][i];
					}
					for (i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
						farray[i] = c_array[1][i + 1];
					}
					set_spline(elarray, farray, af, bf, cf, df, ndata);
					eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
					si1e = spline(eil, elarray, af, bf, cf, df, ndata);

					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						si2e = spline(eil, elarray, af, bf, cf, df, ndata);
						q1cp_ms1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						q1cp_ms0[i - 1][medium - 1] = si2e
								- q1cp_ms1[i - 1][medium - 1] * eil;
						si1e = si2e;
					}
					q1cp_ms1[neke - 1][medium - 1] = q1cp_ms1[neke - 2][medium - 1];
					q1cp_ms0[neke - 1][medium - 1] = q1cp_ms0[neke - 2][medium - 1];

					// "prepare interpolation table for the second MS moment correction"
					// "e-"
					for (i = 0; i <= $MAXE_SPIN; i++) {
						farray[i] = g_array[0][i];
					}
					for (i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
						farray[i] = g_array[0][i + 1];
					}
					set_spline(elarray, farray, af, bf, cf, df, ndata);
					eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
					si1e = spline(eil, elarray, af, bf, cf, df, ndata);

					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						si2e = spline(eil, elarray, af, bf, cf, df, ndata);
						q2ce_ms1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						q2ce_ms0[i - 1][medium - 1] = si2e
								- q2ce_ms1[i - 1][medium - 1] * eil;
						si1e = si2e;
					}
					q2ce_ms1[neke - 1][medium - 1] = q2ce_ms1[neke - 2][medium - 1];
					q2ce_ms0[neke - 1][medium - 1] = q2ce_ms0[neke - 2][medium - 1];
					// "e+"
					for (i = 0; i <= $MAXE_SPIN; i++) {
						farray[i] = g_array[1][i];
					}
					for (i = $MAXE_SPIN + 1; i <= $MAXE_SPI1 - 1; i++) {
						farray[i] = g_array[1][i + 1];
					}
					set_spline(elarray, farray, af, bf, cf, df, ndata);
					eil = (1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
					si1e = spline(eil, elarray, af, bf, cf, df, ndata);

					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						si2e = spline(eil, elarray, af, bf, cf, df, ndata);
						q2cp_ms1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						q2cp_ms0[i - 1][medium - 1] = si2e
								- q2cp_ms1[i - 1][medium - 1] * eil;
						si1e = si2e;
					}
					q2cp_ms1[neke - 1][medium - 1] = q2cp_ms1[neke - 2][medium - 1];
					q2cp_ms0[neke - 1][medium - 1] = q2cp_ms0[neke - 2][medium - 1];

					// "Now substract scattering power that is already taken into account in"
					// "discrete Moller/Bhabha events"
					tauc = TE[medium - 1] / 0.5110034;
					si1e = 1.0;
					for (i = 1; i <= neke - 1; i++) {
						eil = (i + 1.0 - EKE0[medium - 1]) / EKE1[medium - 1];
						e = Math.exp(eil);
						leil = i + 1;
						tau = e / 0.5110034;
						if (tau > 2.0 * tauc) {
							// $EVALUATE sig USING esig(eil);
							sig = ESIG1[leil - 1][medium - 1] * eil
									+ ESIG0[leil - 1][medium - 1];
							// $EVALUATE dedx USING ededx(eil);
							dedx = EDEDX1[leil - 1][medium - 1] * eil
									+ EDEDX0[leil - 1][medium - 1];
							sig = sig / dedx;
							if (sig > 1.e-6) {// "To be sure that this is not a CSDA calc."
												// $EVALUATE etap USING
												// etae_ms(eil);
								etap = etae_ms1[leil - 1][medium - 1] * eil
										+ etae_ms0[leil - 1][medium - 1];
								eta = 0.25 * etap * XCC[medium - 1]
										/ BLCC[medium - 1] / tau / (tau + 2.0);
								g_r = (1.0 + 2.0 * eta)
										* Math.log(1.0 + 1.0 / eta) - 2.0;
								g_m = Math.log(0.5 * tau / tauc)
										+ (1.0 + ((tau + 2.0) / (tau + 1.0))
												* ((tau + 2.0) / (tau + 1.0)))
										* Math.log(2.0 * (tau - tauc + 2.0)
												/ (tau + 4.0))
										- 0.25
										* (tau + 2.0)
										* (tau + 2.0 + 2.0 * (2.0 * tau + 1.0)
												/ ((tau + 1.0) * (tau + 1.0)))
										* Math.log((tau + 4.0) * (tau - tauc)
												/ tau / (tau - tauc + 2.0))
										+ 0.5
										* (tau - 2.0 * tauc)
										* (tau + 2.0)
										* (1.0 / (tau - tauc) - 1.0 / ((tau + 1.0) * (tau + 1.0)));
								if (g_m < g_r) {
									g_m = g_m / g_r;
								} else {
									g_m = 1.0;
								}
								si2e = 1.0 - g_m * sum_Z / sum_Z2;
							} else {
								si2e = 1.0;
							}
						} else {
							si2e = 1.0;
						}
						blcce1[i - 1][medium - 1] = (si2e - si1e)
								* EKE1[medium - 1];
						blcce0[i - 1][medium - 1] = si2e
								- blcce1[i - 1][medium - 1] * eil;
						si1e = si2e;
					}
					blcce1[neke - 1][medium - 1] = blcce1[neke - 2][medium - 1];
					blcce0[neke - 1][medium - 1] = blcce0[neke - 2][medium - 1];

					// "We will not bother to do the same for positrons at this time"

					seqStr = "done!";
					if (iprint > 1)
						eq.printSequence(seqStr);

				}// medium
					// close(spin_unit);
			} catch (EOFException eof) {
				// System.out.println ("End of File");
				// break;
			}
			data_in.close();
		} catch (IOException ioexc) {
			// System.out.println ( "IO Exception =: " + ioexc );
		}

	}

	/**
	 * Initialize data for pair production. Called by HATCH.
	 */
	public static void init_nrc_pair() {

		// implicit none;
		// $declare_max_medium;
		// ;COMIN/MEDIA,BREMPR,ELECIN,NRC-PAIR-DATA,THRESH,USEFUL,EGS-IO/;
		// $declare_write_buffer;

		// character nrcp_file*256, endianess*4;
		// integer egs_get_unit;
		// $INTEGER nrcp_unit, want_nrcp_unit, rec_length;
		int i = 0; // lnblnk1;
		double tmp = 0.0;
		double ddx = 0.0;
		double xx = 0.0;
		double Z = 0.0;
		float emin = 0.0f;
		float emax = 0.0f;
		int ne = 0;
		int nb = 0;
		int ix = 0;
		int ie = 0;
		// int irec,
		int i_ele = 0;
		int nbb = 0;
		int iz = 0;
		// int $cdum_size = 4 * ($NRC_PAIR_NXX - 4) - 1;
		// character endian, cdum($cdum_size);
		int endian = 0;
		String endianess = "";
		boolean swap = false;
		// real*4 tmp_4,
		float[] tarray = new float[$NRC_PAIR_NXX];
		// int itmp_4=0;
		// character c_4(4), ic_4(4);
		// equivalence (tmp_4,c_4), (itmp_4, ic_4);

		byte[] b_4 = new byte[4];
		// byte[] b_2 = new byte[2];

		// " First construct the path to the base directory "
		String file = datas + file_sep + egsDatas + file_sep + pair_nrc_file;

		FileInputStream data_in = null;
		try {
			//FileInputStream data_in = new FileInputStream(file);
			data_in = new FileInputStream(file);
			
			try {
				// swap?
				// The fortran interpreter read data in closely relationship
				// with current CPU endianess
				// In Java, all "read" methods is realted to byte by byte
				// reading, so DO NOT CARE
				// ABOUT current CPU endianess BUT machine generated data
				// endianess!!!!!
				// So, if endianess of machine generated data is little_endian
				// ==>MUST swap!!
				data_in.skip(16);// read endian
				endian = data_in.read();
				if (endian == 0)// little_endian=>mustSWAP
				{
					endianess = "1234";
					swap = true;
				} else {
					endianess = "4321";
					swap = false;
				}
				// read(nrcp_unit,rec=1,err=:nrcp-read-error:) emin, emax, ne,
				// nb, endian, cdum;
				data_in.close();
				data_in = new FileInputStream(file);// reset
				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				emin = byte4toFloat(b_4);
				// System.out.println("emin= "+emin+"  "+swap);

				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				emax = byte4toFloat(b_4);
				// System.out.println("emax= "+emax+"  "+swap);

				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				ne = byte4toInt(b_4);
				// System.out.println("ne= "+ne+"  "+swap);

				data_in.read(b_4);
				if (swap)
					egs_swap_4(b_4);
				nb = byte4toInt(b_4);
				// System.out.println("nb= "+nb+"  "+swap);

				// IF( ichar(endian) = 0 ) [ endianess = '1234'; ] ELSE [
				// endianess = '4321'; ]
				seqStr = "Reading NRC pair data base from " + file;
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "The endianess of the machine generated data is "
						+ endianess + ",   swap:" + swap;
				if (iprint > 1)
					eq.printSequence(seqStr);
				seqStr = "Energy range of the data: " + emin + ",   " + emax;
				if (iprint > 1)
					eq.printSequence(seqStr);

				if (nb != $NRC_PAIR_NXX) {
					// $egs_fatal(*,'Inconsistent x-grid size');
					STOPPROGRAM = true;
					seqStr = "Inconsistent x-grid size";
					eq.printSequence(seqStr);

					return;

				}
				if (ne != $NRC_PAIR_NEE) {
					STOPPROGRAM = true;
					seqStr = "Inconsistent energy grid size";
					eq.printSequence(seqStr);

					return;
					// $egs_fatal(*,'Inconsistent energy grid size');
				}

				nrcp_emin = emin;
				nrcp_emax = emax;
				nrcp_dle = Math.log((emax - 2.0) / (emin - 2.0)) / (ne - 1.0);
				nrcp_dlei = 1.0 / nrcp_dle;

				nbb = nb / 2;
				ddx = Math.sqrt(0.5) / nbb;
				for (ix = 0; ix <= nbb; ix++) {
					xx = ddx * ix;
					nrcp_xdata[ix] = xx * xx;
				}// nrcp_xdata[ix+1] = xx*xx; }
				for (ix = nbb - 1; i >= 0; i--) {
					xx = ddx * ix;
					nrcp_xdata[nb - ix - 1] = 1.0 - xx * xx;
				}// nrcp_xdata(nb-ix) = 1 - xx*xx; }

				for (int medium = 1; medium <= NMED; medium++) {
					seqStr = "  medium " + medium + " ..................... ";
					if (iprint > 1)
						eq.printSequence(seqStr);

					for (ie = 1; ie <= $NRC_PAIR_NEE; ie++) {
						for (ix = 1; ix <= $NRC_PAIR_NXX; ix++) {
							nrcp_fdata[ix - 1][ie - 1][medium - 1] = 0.0;
						}
					}
					for (i_ele = 1; i_ele <= NNE[medium - 1]; i_ele++) {
						Z = ZELEM[medium - 1][i_ele - 1];
						Double izd = new Double(Z + 0.5);
						iz = izd.intValue();// int(Z+0.5);
						tmp = PZ[medium - 1][i_ele - 1] * Z * Z;
						// irec = (iz-1)*ne + 2;
						// ##############################
						data_in = new FileInputStream(file);// reset
						// based on comparission of total bytes vs the next
						// total,some bytes must be skipped!!
						int nskip = (iz - 1) * 21840 + 260;// 21840=======see
															// next
						data_in.skip(nskip);
						// ###############################
						// 21840=$NRC_PAIR_NEE*$NRC_PAIR_NXX*4!!!
						for (ie = 1; ie <= $NRC_PAIR_NEE; ie++) {
							// read(nrcp_unit,rec=irec,err=:nrcp-read-error:)
							// tarray;
							for (int it = 1; it <= $NRC_PAIR_NXX; it++) {
								data_in.read(b_4);
								if (swap)
									egs_swap_4(b_4);
								tarray[it - 1] = byte4toFloat(b_4);
								// System.out.println(tarray[it-1]+"  "+swap);
							}

							for (ix = 1; ix <= $NRC_PAIR_NXX; ix++) {
								// tmp_4 = tarray(ix);
								// IF( swap ) [ call egs_swap_4(c_4); ]
								// nrcp_fdata(ix,ie,medium)=nrcp_fdata(ix,ie,medium)+tmp*tmp_4;
								nrcp_fdata[ix - 1][ie - 1][medium - 1] = nrcp_fdata[ix - 1][ie - 1][medium - 1]
										+ tmp * tarray[ix - 1];
							}
							// irec = irec + 1;
						}
					}
					for (ie = 1; ie <= $NRC_PAIR_NEE; ie++) {
						// prepare_alias_table(nb-1,nrcp_xdata,nrcp_fdata(1,ie,medium),
						// nrcp_wdata(1,ie,medium),nrcp_idata(1,ie,medium));

						double[] patx = new double[nb];
						double[] patf = new double[nb];
						double[] patw = new double[nb - 1];
						int[] pati = new int[nb - 1];
						for (int ll = 0; ll <= nb - 1; ll++) {
							patx[ll] = nrcp_xdata[ll];
							patf[ll] = nrcp_fdata[ll][ie - 1][medium - 1];
						}
						prepare_alias_table(nb - 1, patx, patf, patw, pati);

						for (int ll = 0; ll < nb - 1; ll++) {
							nrcp_wdata[ll][ie - 1][medium - 1] = patw[ll];
							nrcp_idata[ll][ie - 1][medium - 1] = pati[ll];
						}

					}

					seqStr = " done!";
					if (iprint > 1)
						eq.printSequence(seqStr);

				}

			} catch (EOFException eof) {
				// System.out.println ("End of File");
				// break;
			}

			//data_in.close();
		} catch (IOException ioexc) {
			// System.out.println ( "IO Exception =: " + ioexc );
		}finally {
	        if( null != data_in ) {
	            try 
	            {
	                data_in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
	    }

	}

	// if user change some var-> it could be usefull to call reset in order to
	// perform
	// a new hatch, simulation etc.Before hatch call the egs_set_defaults must
	// be also called
	/**
	 * Reset all global variable for re-use.
	 */
	public static void reset() {
		ECUT = new double[$MXREG];
		PCUT = new double[$MXREG];
		E = new double[$MXSTACK];
		X = new double[$MXSTACK];
		Y = new double[$MXSTACK];
		Z = new double[$MXSTACK];
		U = new double[$MXSTACK];
		V = new double[$MXSTACK];
		W = new double[$MXSTACK];
		DNEAR = new double[$MXSTACK];
		WT = new double[$MXSTACK];
		IQ = new int[$MXSTACK];
		IR = new int[$MXSTACK];
		LATCH = new int[$MXSTACK];
		NP = 0;
		NPold = 0;
		RHOR = new double[$MXREG];
		MED = new int[$MXREG];
		IRAYLR = new int[$MXREG];
		MEDIA = new String[$MXMED];
		RLC = new double[$MXMED];
		RLDU = new double[$MXMED];
		RHO = new double[$MXMED];
		MSGE = new int[$MXMED];
		MGE = new int[$MXMED];
		MSEKE = new int[$MXMED];
		MEKE = new int[$MXMED];
		MLEKE = new int[$MXMED];
		MCMFP = new int[$MXMED];
		MRANGE = new int[$MXMED];
		IRAYLM = new int[$MXMED];
		iz_array = new int[$MXTOTSH];
		be_array = new double[$MXTOTSH];
		Jo_array = new double[$MXTOTSH];
		erfJo_array = new double[$MXTOTSH];
		ne_array = new int[$MXTOTSH];
		shn_array = new int[$MXTOTSH];
		shell_array = new int[$MXMDSH][$MXMED];
		eno_array = new double[$MXMDSH][$MXMED];
		n_shell = new int[$MXMED];
		ibcmp = new int[$MXREG];
		binding_energies = new double[$MXSHELL][$MXELEMENT];
		interaction_prob = new double[$MXSHELL][$MXELEMENT];
		relaxation_prob = new double[$MXTRANS][$MXELEMENT];
		edge_energies = new double[$MXEDGE][$MXELEMENT];
		edge_number = new int[$MXELEMENT];
		edge_a = new double[$MXEDGE][$MXELEMENT];
		edge_b = new double[$MXEDGE][$MXELEMENT];
		edge_c = new double[$MXEDGE][$MXELEMENT];
		edge_d = new double[$MXEDGE][$MXELEMENT];
		iedgfl = new int[$MXREG];
		iphter = new int[$MXREG];
		$MAX_EII_BINS = $N_EII_BINS * $MAX_EII_SHELLS;
		eii_xsection_a = new double[$MAX_EII_BINS];
		eii_xsection_b = new double[$MAX_EII_BINS];
		eii_z = new int[$MAX_EII_SHELLS];
		eii_sh = new int[$MAX_EII_SHELLS];
		eii_a = new double[$MAX_EII_SHELLS];
		eii_b = new double[$MAX_EII_SHELLS];
		eii_nshells = new int[$MXELEMENT];
		eii_nsh = new int[$MXMED];
		eii_cons = new double[$MXMED];
		eii_first = new int[$MXMED][$MXEL];
		eii_no = new int[$MXMED][$MXEL];
		u_relax = 0.0;
		ish_relax = 0;
		iZ_relax = 0;
		SMAXIR = new double[$MXREG];
		e_max_rr = new double[$MXREG];
		i_do_rr = new int[$MXREG];
		AP = new double[$MXMED];
		AE = new double[$MXMED];
		UP = new double[$MXMED];
		UE = new double[$MXMED];
		TE = new double[$MXMED];
		THMOLL = new double[$MXMED];
		PZERO = 0.0;
		MEDIUM = 0;
		MEDOLD = 0;
		DL1 = new double[8][$MXMED];
		DL2 = new double[8][$MXMED];
		DL3 = new double[8][$MXMED];
		DL4 = new double[8][$MXMED];
		DL5 = new double[8][$MXMED];
		DL6 = new double[8][$MXMED];
		ALPHI = new double[2][$MXMED];
		BPAR = new double[2][$MXMED];
		DELPOS = new double[2][$MXMED];
		WA = new double[$MXMED][$MXEL];
		PZ = new double[$MXMED][$MXEL];
		ZELEM = new double[$MXMED][$MXEL];
		RHOZ = new double[$MXMED][$MXEL];
		PWR2I = new double[$MXPWR2I];
		DELCM = new double[$MXMED];
		ZBRANG = new double[$MXMED];
		LZBRANG = new double[$MXMED];
		NNE = new int[$MXMED];
		ASYM = new String[$MXMED][$MXEL];
		nb_fdata = new double[$MXBRXS + 1][$MXBRES][$MXMED];
		nb_xdata = new double[$MXBRXS + 1][$MXBRES][$MXMED];
		nb_wdata = new double[$MXBRXS][$MXBRES][$MXMED];
		nb_idata = new int[$MXBRXS][$MXBRES][$MXMED];
		nb_emin = new double[$MXMED];
		nb_emax = new double[$MXMED];
		nb_lemin = new double[$MXMED];
		nb_lemax = new double[$MXMED];
		nb_dle = new double[$MXMED];
		nb_dlei = new double[$MXMED];
		log_ap = new double[$MXMED];
		EDEP = 0.0;
		TSTEP = 0.0;
		TUSTEP = 0.0;
		USTEP = 0.0;
		VSTEP = 0.0;
		TVSTEP = 0.0;
		RHOF = 0.0;
		EOLD = 0.0;
		ENEW = 0.0;
		EKE = 0.0;
		ELKE = 0.0;
		GLE = 0.0;
		E_RANGE = 0.0;
		iausfl = new int[$MXAUS];
		IDISC = 0;
		IROLD = 0;
		IRNEW = 0;
		ums_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
		fms_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
		wms_array = new double[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
		ims_array = new int[$MAXL_MS + 1][$MAXQ_MS + 1][$MAXU_MS + 1];
		llammin = 0.0;
		llammax = 0.0;
		dllamb = 0.0;
		dllambi = 0.0;
		dqms = 0.0;
		dqmsi = 0.0;
		spin_rej = new double[$MXMED][2][$MAXE_SPI1 + 1][$MAXQ_SPIN + 1][$MAXU_SPIN + 1];
		espin_min = 0.0;
		espin_max = 0.0;
		espml = 0.0;
		b2spin_min = 0.0;
		b2spin_max = 0.0;
		dbeta2 = 0.0;
		dbeta2i = 0.0;
		dlener = 0.0;
		dleneri = 0.0;
		dqq1 = 0.0;
		dqq1i = 0.0;
		fool_intel_optimizer = false;
		is_ch_step = false;
		SINC0 = 0.0;
		SINC1 = 0.0;
		SIN0 = new double[$MXSINC];
		SIN1 = new double[$MXSINC];
		THETA = 0.0;
		SINTHE = 0.0;
		COSTHE = 0.0;
		SINPHI = 0.0;
		COSPHI = 0.0;
		esig_e = new double[$MXMED];
		psig_e = new double[$MXMED];
		esige_max = 0.0;
		psige_max = 0.0;
		range_ep = new double[2][$MXEKE][$MXMED];
		E_array = new double[$MXEKE][$MXMED];
		etae_ms0 = new double[$MXEKE][$MXMED];
		etae_ms1 = new double[$MXEKE][$MXMED];
		etap_ms0 = new double[$MXEKE][$MXMED];
		etap_ms1 = new double[$MXEKE][$MXMED];
		q1ce_ms0 = new double[$MXEKE][$MXMED];
		q1ce_ms1 = new double[$MXEKE][$MXMED];
		q1cp_ms0 = new double[$MXEKE][$MXMED];
		q1cp_ms1 = new double[$MXEKE][$MXMED];
		q2ce_ms0 = new double[$MXEKE][$MXMED];
		q2ce_ms1 = new double[$MXEKE][$MXMED];
		q2cp_ms0 = new double[$MXEKE][$MXMED];
		q2cp_ms1 = new double[$MXEKE][$MXMED];
		blcce0 = new double[$MXEKE][$MXMED];
		blcce1 = new double[$MXEKE][$MXMED];
		expeke1 = new double[$MXMED];
		EKE0 = new double[$MXMED];
		EKE1 = new double[$MXMED];
		XR0 = new double[$MXMED];
		TEFF0 = new double[$MXMED];
		BLCC = new double[$MXMED];
		XCC = new double[$MXMED];
		ESIG0 = new double[$MXEKE][$MXMED];
		ESIG1 = new double[$MXEKE][$MXMED];
		PSIG0 = new double[$MXEKE][$MXMED];
		PSIG1 = new double[$MXEKE][$MXMED];
		EDEDX0 = new double[$MXEKE][$MXMED];
		EDEDX1 = new double[$MXEKE][$MXMED];
		PDEDX0 = new double[$MXEKE][$MXMED];
		PDEDX1 = new double[$MXEKE][$MXMED];
		EBR10 = new double[$MXEKE][$MXMED];
		EBR11 = new double[$MXEKE][$MXMED];
		PBR10 = new double[$MXEKE][$MXMED];
		PBR11 = new double[$MXEKE][$MXMED];
		PBR20 = new double[$MXEKE][$MXMED];
		PBR21 = new double[$MXEKE][$MXMED];
		TMXS0 = new double[$MXEKE][$MXMED];
		TMXS1 = new double[$MXEKE][$MXMED];
		IUNRST = new int[$MXMED];
		EPSTFL = new int[$MXMED];
		IAPRIM = new int[$MXMED];
		EBINDA = new double[$MXMED];
		GE0 = new double[$MXMED];
		GE1 = new double[$MXMED];
		GMFP0 = new double[$MXGE][$MXMED];
		GMFP1 = new double[$MXGE][$MXMED];
		GBR10 = new double[$MXGE][$MXMED];
		GBR11 = new double[$MXGE][$MXMED];
		GBR20 = new double[$MXGE][$MXMED];
		GBR21 = new double[$MXGE][$MXMED];
		RCO0 = new double[$MXMED];
		RCO1 = new double[$MXMED];
		RSCT0 = new double[$MXRAYFF][$MXMED];
		RSCT1 = new double[$MXRAYFF][$MXMED];
		COHE0 = new double[$MXGE][$MXMED];
		COHE1 = new double[$MXGE][$MXMED];
		DPMFP = 0.0;
		MPGEM = new int[$MXSGE][$MXMED];
		NGR = new int[$MXMED];
		IRAYL = 0;
		NEKE = 0;
		I1ST = 1;
		got_data = false;
		tperp = 0.0;
		irl = 0;
		STOPPROGRAM = false;
		startSimulationTime = 0;
		rng_array = new double[24];
		rng_array1 = new double[$NRANMAR];
		seeds = new int[24];
		seedin = 0;
		luxury_level = 0;
		state = new int[25];
		carry = 0;
		i24 = 0;
		j24 = 0;
		next = new int[24];
		nskipRnd = 0;
		status = 0;
		jseed = 0;
		icarry = 0;
		kRnd = 0;
		jRnd = 0;
		not_initialized = true;
		uni = 0;
		twom24 = 0.0;
		twop24 = 0.0;
		urndm = new int[97];
		crndm = 0;
		cdrndm = 0;
		cmrndm = 0;
		i4opt = 0;
		ixx = 0;
		jxx = 0;
		fool_optimizer = 0;

		DATA = new double[$MXDATA][$STAT];
		IERRsgm = 0;
		ICOUNT = 0;
		JHSTRY = 1;
		KE = 0.0;
		nf = NumberFormat.getInstance(Locale.US);
		nff = new DecimalFormat(pattern, dfs);
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
		radc_sigs = new double[$RADC_NE + 1];
		radc_sigd = new double[$RADC_NE + 1];
		radc_frej = new double[$RADC_NE + 1][$RADC_NU + 1];
		radc_x = new double[$RADC_NX];
		radc_fdat = new double[$RADC_NBOX];
		radc_Smax = new double[$RADC_NBOX];
		radc_emin = 0.0;
		radc_emax = 0.0;
		radc_lemin = 0.0;
		radc_dw = 0.0;
		radc_dle = 0.0;
		radc_dlei = 0.0;
		radc_le1 = 0.0;
		radc_bins = new int[$RADC_NBOX];
		radc_ixmin1 = new int[$RADC_NBOX];
		radc_ixmax1 = new int[$RADC_NBOX];
		radc_ixmin2 = new int[$RADC_NBOX];
		radc_ixmax2 = new int[$RADC_NBOX];
		radc_ixmin3 = new int[$RADC_NBOX];
		radc_ixmax3 = new int[$RADC_NBOX];
		radc_ixmin4 = new int[$RADC_NBOX];
		radc_ixmax4 = new int[$RADC_NBOX];
		radc_startx = new int[$RADC_NE + 1];
		radc_startb = new int[$RADC_NE + 1];
		a_triplet = new double[$MAX_TRIPLET][$MXMED];
		b_triplet = new double[$MAX_TRIPLET][$MXMED];
		dl_triplet = 0.0;
		dli_triplet = 0.0;
		bli_triplet = 0.0;
		log_4rm = 0.0;
		nrcp_idata = new int[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
		nrcp_xdata = new double[$NRC_PAIR_NXX];
		nrcp_fdata = new double[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
		nrcp_wdata = new double[$NRC_PAIR_NXX][$NRC_PAIR_NEE][$MXMED];
		nrcp_emin = 0.0;
		nrcp_emax = 0.0;
		nrcp_dle = 0.0;
		nrcp_dlei = 0.0;

		photon_xsections = "";
		eiifile = "eii_ik";
		sig_ismonotone = new boolean[2][$MXMED];
	}

	// ====================
	public static String photon_xsections = "";
	public static final String photon_xsections_epdl = "epdl";
	public static final String photon_xsections_xcom = "xcom";
	// "If photon_xsections is not empty, photon cross sections will be"
	// "re-initialized using data files  "
	// "  'photon_xsection'_photo.data   "
	// "  'photon_xsection'_pair.data    "
	// "  'photon_xsection'_triplet.data "
	// "  'photon_xsection'_rayleigh.data"
	// "that must be placed in $HEN_HOUSE/data"
	// REPLACE {$MXINPUT} WITH {2000};
	public static int $MXINPUT = 2000;

	// "============================================================================="
	/**
	 * (Re)-initialize the photon cross section data using data files from the series prefix. Called by HATCH. 
	 * @param prefix prefix
	 * @param out out must be 1 (for now).
	 */
	public static void egs_init_user_photon(String prefix, int out) {
		// "============================================================================="
		// implicit none;
		// $declare_max_medium;//REPLACE {$declare_max_medium;} WITH {;};
		// character*(*) prefix;
		// $INTEGER out;

		int ndat = 0;
		int[] sorted = new int[$MXEL];
		int i = 0;
		int j = 0;
		int k = 0;
		int iz = 0;
		// int iz_old = 0;
		int nge = 0;
		double[] z_sorted = new double[$MXEL];
		double[] pz_sorted = new double[$MXEL];
		double[] sig_photo = new double[$MXGE];
		double[] sig_pair = new double[$MXGE];
		double[] sig_triplet = new double[$MXGE];
		double[] sig_rayleigh = new double[$MXGE];
		double sigma = 0.0;
		double cohe = 0.0;
		double gmfp = 0.0;
		double gbr1 = 0.0;
		double gbr2 = 0.0;
		double sig_KN = 0.0;
		double gle = 0.0;
		double e = 0.0;
		double sig_p = 0.0;
		double cohe_old = 0.0;
		double gmfp_old = 0.0;
		double gbr1_old = 0.0;
		double gbr2_old = 0.0;
		double[] etmp = new double[$MXINPUT];
		double[] ftmp = new double[$MXINPUT];
		double sumZ = 0.0;
		double sumA = 0.0;
		// double con1 = 0.0;
		double con2 = 0.0;
		double egs_KN_sgma0 = 0.0;

		seqStr = "(Re)-initializing photon cross section data" + " \n"
				+ "Using data files from the series " + prefix;
		if (iprint > 1)
			eq.printSequence(seqStr);
		// System.out.println(seqStr);

		String data_dir = datas + file_sep + egsDatas + file_sep;// defaultext=.data;
		String photo_file = data_dir + prefix + "_photo.data";
		String pair_file = data_dir + prefix + "_pair.data";
		String triplet_file = data_dir + prefix + "_triplet.data";
		String rayleigh_file = data_dir + prefix + "_rayleigh.data";

		String s = "";
		int ll = 0;
		if (out == 1) {
			seqStr = "Photon cross sections initialized from " + prefix
					+ " data files";
			if (iprint > 1)
				eq.printSequence(seqStr);
			seqStr = "============================================================================";
			if (iprint > 1)
				eq.printSequence(seqStr);
			seqStr = "Grid energies and cross sections are output";
			if (iprint > 1)
				eq.printSequence(seqStr);
			seqStr = "Grid energies and cross sections are output";
			if (iprint > 1)
				eq.printSequence(seqStr);

			// write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a)')
			// 'Energy',' GMFP(cm) ',' Pair ','Compton',' GMFP(cm) ';

			s = format("", 5) + "Energy";
			ll = s.length();
			ll = 19 - ll;
			s = s + format(" ", ll);

			s = s + " GMFP(cm) ";
			ll = s.length();
			ll = 34 - ll;
			s = s + format(" ", ll);

			s = s + " Pair ";
			ll = s.length();
			ll = 49 - ll;
			s = s + format(" ", ll);

			s = s + "Compton";
			ll = s.length();
			ll = 64 - ll;
			s = s + format(" ", ll);

			s = s + " GMFP(cm) ";

			seqStr = s;
			if (iprint > 1)
				eq.printSequence(seqStr);
			// System.out.println(seqStr);
			// write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a/)')
			// '(MeV)','no Rayleigh','(fraction)','(fraction)','with Rayleigh';

			s = format("", 5) + "(MeV)";
			ll = s.length();
			ll = 19 - ll;
			s = s + format(" ", ll);

			s = s + "no Rayleigh";
			ll = s.length();
			ll = 34 - ll;
			s = s + format(" ", ll);

			s = s + "(fraction)";
			ll = s.length();
			ll = 49 - ll;
			s = s + format(" ", ll);

			s = s + "(fraction)";
			ll = s.length();
			ll = 64 - ll;
			s = s + format(" ", ll);

			s = s + "with Rayleigh";

			seqStr = s;
			if (iprint > 1)
				eq.printSequence(seqStr);
			// System.out.println(seqStr);

		}
		// Replace binding energies with the edges in the photo-absorption file
		int iread = 0;// int lnr =0;//data number
		// int lnrr =0;//line number
		// int indx=0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";char comma=',';
		// char lineSep = '\n';//System.getProperty("line.separator").charAt(0);

		boolean ndatB = true;
		k = 1;
		boolean neB = true;
		iz = 1;
		try {
			FileInputStream in = new FileInputStream(photo_file);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						// lnr++;

						if (ndatB)// new dat line!!
						{
							s = desc.toString();
							ndat = stringToInt(s);// System.out.println(ndat);

							ndatB = false;
						} else {
							s = desc.toString();

							if (neB) {
								etmp[k - 1] = stringToDouble(s);
								neB = false;
							} else {
								ftmp[k - 1] = stringToDouble(s);
								// System.out.println(etmp[k-1]+ "   "+
								// ftmp[k-1]);
								k = k + 1;
								neB = true;
							}

							if (k == ndat + 1) {
								// ===========
								k = 0;
								for (j = ndat; j >= 2; j--) {
									if (etmp[j - 1] - etmp[j - 2] < 1.0e-5) {
										k = k + 1;
										binding_energies[k - 1][iz - 1] = Math
												.exp(etmp[j - 1]);
										// System.out.println("k "+k+" iz "+iz+
										// "  be "+
										// binding_energies[k-1][iz-1]);
										if (k >= 4)
											break;
									}
								}
								// ============
								iz = iz + 1;
								ndatB = true;
								k = 1;
							}
						}

					}// have data
						// if ((char)iread == lineSep)
						// lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
		}// try
		catch (Exception exc) {

		}

		for (int medium = 1; medium <= NMED; medium++) {
			MGE[medium - 1] = $MXGE;
			nge = $MXGE;
			GE1[medium - 1] = nge - 1.0;
			GE1[medium - 1] = GE1[medium - 1]
					/ Math.log(UP[medium - 1] / AP[medium - 1]);
			GE0[medium - 1] = 1.0 - GE1[medium - 1] * Math.log(AP[medium - 1]);

			seqStr = "Working on medium " + medium + " which is "
					+ MEDIA[medium - 1] + " ... ";
			if (iprint > 1)
				eq.printSequence(seqStr);

			sumZ = 0.0;
			sumA = 0.0;
			for (i = 1; i <= NNE[medium - 1]; i++) {
				z_sorted[i - 1] = ZELEM[medium - 1][i - 1];// (medium,i);
				sumZ = sumZ + PZ[medium - 1][i - 1] * ZELEM[medium - 1][i - 1];// zelem(medium,i);
				sumA = sumA + PZ[medium - 1][i - 1] * WA[medium - 1][i - 1];
			}

			// con1 = sumZ * RHO[medium - 1] / (sumA * 1.6605655);
			con2 = RHO[medium - 1] / (sumA * 1.6605655);

			egs_heap_sort(NNE[medium - 1], z_sorted, sorted);

			for (i = 1; i <= NNE[medium - 1]; i++) {
				pz_sorted[i - 1] = PZ[medium - 1][sorted[i - 1] - 1];
			}

			egsi_get_data(0, photo_file, nge, NNE[medium - 1], z_sorted,
					pz_sorted, GE1[medium - 1], GE0[medium - 1], sig_photo);

			egsi_get_data(0, rayleigh_file, nge, NNE[medium - 1], z_sorted,
					pz_sorted, GE1[medium - 1], GE0[medium - 1], sig_rayleigh);
			egsi_get_data(1, pair_file, nge, NNE[medium - 1], z_sorted,
					pz_sorted, GE1[medium - 1], GE0[medium - 1], sig_pair);
			egsi_get_data(2, triplet_file, nge, NNE[medium - 1], z_sorted,
					pz_sorted, GE1[medium - 1], GE0[medium - 1], sig_triplet);

			for (i = 1; i <= nge; i++) {
				gle = (i - GE0[medium - 1]) / GE1[medium - 1];
				e = Math.exp(gle);
				egs_KN_sgma0 = egs_KN_sigma0(e);
				sig_KN = sumZ * egs_KN_sgma0;
				sig_p = sig_pair[i - 1] + sig_triplet[i - 1];
				sigma = sig_KN + sig_p + sig_photo[i - 1];
				gmfp = 1.0 / (sigma * con2);
				gbr1 = sig_p / sigma;
				gbr2 = gbr1 + sig_KN / sigma;
				cohe = sigma / (sig_rayleigh[i - 1] + sigma);
				if (out == 1) {
					// write(ounit,'(5(1pe15.6))')
					// e,gmfp,gbr1,gbr2-gbr1,gmfp*cohe;
					seqStr = format(e, 15) + format(gmfp, 15)
							+ format(gbr1, 15) + format(gbr2 - gbr1, 15)
							+ format(gmfp * cohe, 15);
					if (iprint > 1)
						eq.printSequence(seqStr);
				}
				if (i > 1) {
					GMFP1[i - 2][medium - 1] = (gmfp - gmfp_old)
							* GE1[medium - 1];
					GMFP0[i - 2][medium - 1] = gmfp - GMFP1[i - 2][medium - 1]
							* gle;
					GBR11[i - 2][medium - 1] = (gbr1 - gbr1_old)
							* GE1[medium - 1];
					GBR10[i - 2][medium - 1] = gbr1 - GBR11[i - 2][medium - 1]
							* gle;
					GBR21[i - 2][medium - 1] = (gbr2 - gbr2_old)
							* GE1[medium - 1];
					GBR20[i - 2][medium - 1] = gbr2 - GBR21[i - 2][medium - 1]
							* gle;
					COHE1[i - 2][medium - 1] = (cohe - cohe_old)
							* GE1[medium - 1];
					COHE0[i - 2][medium - 1] = cohe - COHE1[i - 2][medium - 1]
							* gle;
				}
				gmfp_old = gmfp;
				gbr1_old = gbr1;
				gbr2_old = gbr2;
				cohe_old = cohe;
				GMFP1[nge - 1][medium - 1] = GMFP1[nge - 2][medium - 1];
				GMFP0[nge - 1][medium - 1] = gmfp - GMFP1[nge - 1][medium - 1]
						* gle;
				GBR11[nge - 1][medium - 1] = GBR11[nge - 2][medium - 1];
				GBR10[nge - 1][medium - 1] = gbr1 - GBR11[nge - 1][medium - 1]
						* gle;
				GBR21[nge - 1][medium - 1] = GBR21[nge - 2][medium - 1];
				GBR20[nge - 1][medium - 1] = gbr2 - GBR21[nge - 1][medium - 1]
						* gle;
				COHE1[nge - 1][medium - 1] = COHE1[nge - 2][medium - 1];
				COHE0[nge - 1][medium - 1] = cohe - COHE1[nge - 1][medium - 1]
						* gle;

			}

			seqStr = "OK";
			if (iprint > 1)
				eq.printSequence(seqStr);
		}

		return; // end;
	}

	// "============================================================================="
	/**
	 * Sort the real array rarray of dimension n in ascending order and at the same time put into the integer array 
	 * jarray the original position of the elements, e.g. if rarray was on input (5,14,8,2), it will be after completion 
	 * of heap_sort (2,5,8,14) and jarray will be (4,1,3,2). Heap_sort uses the heap sort algorithm, the implementation is 
	 * based on hpsort from Numerical Recipies with a couple of modifications.
	 * @param n n
	 * @param rarray rarray
	 * @param jarray jarray
	 */
	public static void egs_heap_sort(int n, double[] rarray, int[] jarray) {
		// "************************************************************************
		// " egs_heap_sort will sort the real array rarray of dimension n in
		// " ascending order and at the same time put into the integer array
		// " jarray the original position of the elements, e.g.
		// " if rarray was on input (5,14,8,2), it will be after completion
		// " of heap_sort (2,5,8,14) and jarray will be (4,1,3,2).
		// " heap_sort uses the heap sort algorithm, the implementation is
		// " based on hpsort from Numerical Recipies with a couple of
		// " modifications.
		// "
		// " Iwan Kawrakow, NRC, July 2001
		// "*************************************************************************

		// implicit none;

		// $INTEGER n,jarray(*);
		// $REAL rarray(*);

		int i = 0;
		int ir = 0;
		int j = 0;
		int l = 0;
		int ira = 0;
		double rra = 0.0;

		for (i = 1; i <= n; i++) {
			jarray[i - 1] = i;
		}// { jarray(i)=i; }
		if (n < 2)
			return;
		l = n / 2 + 1;
		ir = n;

		while (true) {
			if (l > 1) {
				l = l - 1;
				rra = rarray[l - 1];// rarray(l);
				ira = l;
			} else {
				rra = rarray[ir - 1];// rarray(ir);
				ira = jarray[ir - 1];// jarray(ir);
				rarray[ir - 1] = rarray[0];// rarray(ir)=rarray(1);
				jarray[ir - 1] = jarray[0];// jarray(ir)=jarray(1);
				ir = ir - 1;
				if (ir == 1) {
					rarray[0] = rra;// rarray(1)=rra;
					jarray[0] = ira; // jarray(1)=ira;
					return;
				}
			}
			i = l;
			j = l + l;
			while (true) {
				if (j > ir)
					break;
				if (j < ir) {
					// IF (rarray(j) < rarray(j+1) ) j=j+1;
					if (rarray[j - 1] < rarray[j])
						j = j + 1;
				}
				if (rra < rarray[j - 1])// if (rra < rarray(j))
				{
					// rarray(i)=rarray(j); jarray(i)=jarray(j);
					rarray[i - 1] = rarray[j - 1];
					jarray[i - 1] = jarray[j - 1];
					i = j;
					j = j + j;
				} else {
					j = ir + 1;
				}
			}
			// rarray(i)=rra; jarray(i)=ira;
			rarray[i - 1] = rra;
			jarray[i - 1] = ira;
		}
		// return; //end;
	}

	/**
	 * Internally used. Called by egs_init_user_photon.
	 * @param flag flag
	 * @param filename filename
	 * @param n n
	 * @param ne ne
	 * @param zsorted zsorted
	 * @param pz_sorted pz_sorted
	 * @param ge1 ge1
	 * @param ge0 ge0
	 * @param data data
	 */
	public static void egsi_get_data(int flag, String filename, int n, int ne,
			double[] zsorted, double[] pz_sorted, double ge1, double ge0,
			double[] data) {
		// "=========================================================================="
		// implicit none;
		// $INTEGER flag,iunit,n,ne;
		// $REAL ge1,ge0,zsorted(*),pz_sorted(*),data(*);
		double[] etmp = new double[$MXINPUT];
		double[] ftmp = new double[$MXINPUT];
		double gle = 0.0;
		double sig = 0.0;
		double p = 0.0;
		double e = 0.0;
		int i = 0;
		// int j = 0;
		int k = 0;
		int kk = 0;
		int iz = 0;
		int ndat = 0;
		int iiz = 0;

		double rm = 0.5110034;
		double eth = 0.0;

		// rewind(iunit);
		int iz_old = 0;

		for (k = 1; k <= n; k++) {
			data[k - 1] = 0.0;
		}

		i = 1;
		Double iizd = new Double(zsorted[i - 1] + 0.5);// iiz =
														// int(zsorted(i)+0.5);
		iiz = iizd.intValue();

		// @@@@@@@@@@@@@@@@@@@@@@@READ===============================
		int iread = 0;// int lnr =0;//data number
		String s = "";
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";char comma=',';
		// char lineSep = '\n';//System.getProperty("line.separator").charAt(0);

		boolean ndatB = true;
		k = 1;
		boolean neB = true;
		iz = iz_old + 1;
		
		FileInputStream in = null;//new FileInputStream(filename);
		try {
			//FileInputStream in = new FileInputStream(filename);
			in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset

						if (ndatB)// new dat line!!
						{
							s = desc.toString();
							ndat = stringToInt(s);
							// System.out.println(ndat);
							if (ndat > $MXINPUT) {
								// $egs_fatal(*,'Too many input data points.
								// Max. is ',$MXINPUT);
								STOPPROGRAM = true;
								seqStr = "Too many input data points. Max. is "
										+ $MXINPUT;
								eq.printSequence(seqStr);
								return;
							}

							ndatB = false;
						} else {
							s = desc.toString();

							if (neB) {
								if (flag == 0)
									etmp[k - 1] = stringToDouble(s);
								else
									etmp[k] = stringToDouble(s);

								neB = false;
							} else {
								if (flag == 0) {
									ftmp[k - 1] = stringToDouble(s);
									// System.out.println(etmp[k-1]+ "   "+
									// ftmp[k-1]);
								} else {
									ftmp[k] = stringToDouble(s);
									// System.out.println(etmp[k]+ "   "+
									// ftmp[k]);
								}
								k = k + 1;
								neB = true;
							}

							if (k == ndat + 1) {
								// ===========
								if (flag != 0) {
									if (flag == 1) {
										eth = 2. * rm;
									} else {
										eth = 4. * rm;
									}
									ndat = ndat + 1;
									for (k = 2; k <= ndat; k++) {
										ftmp[k - 1] = ftmp[k - 1]
												- 3.
												* Math.log(1.0 - eth
														/ Math.exp(etmp[k - 1]));
									}
									ftmp[0] = ftmp[1];
									etmp[0] = Math.log(eth);
								}

								// ============
								iz = iz + 1;
								ndatB = true;
								k = 1;

								if (iz > iiz) // we have it
								{
									iz_old = iiz;
									// ==========rest
									for (k = 1; k <= n; k++) {
										gle = (k - ge0) / ge1;
										e = Math.exp(gle);
										if (gle < etmp[0]
												|| gle >= etmp[ndat - 1]) {
											if (flag == 0) {
												// $egs_fatal(*,'Energy
												// ',exp(gle),
												// ' is outside the available
												// data range of ',
												// exp(etmp(1)),exp(etmp(ndat)));
												STOPPROGRAM = true;
												seqStr = "Energy "
														+ Math.exp(gle)
														+ " is outside the available data range of "
														+ Math.exp(etmp[0])
														+ " ,  "
														+ Math.exp(etmp[ndat - 1]);
												eq.printSequence(seqStr);
												return;
											} else {
												if (gle < etmp[0]) {
													sig = 0.0;
												} else {
													sig = Math
															.exp(ftmp[ndat - 1]);
												}
											}
										} else {
											for (kk = 1; kk <= ndat - 1; kk++) {
												if (gle >= etmp[kk - 1]
														&& gle < etmp[kk])
													break;
											}
											p = (gle - etmp[kk - 1])
													/ (etmp[kk] - etmp[kk - 1]);
											sig = Math.exp(p * ftmp[kk]
													+ (1.0 - p) * ftmp[kk - 1]);
										}
										if (flag != 0 && e > eth)
											sig = sig * (1. - eth / e)
													* (1. - eth / e)
													* (1. - eth / e);
										data[k - 1] = data[k - 1]
												+ pz_sorted[i - 1] * sig;
									}
									// ===============
									k = 1;// reset
									iz = iz_old + 1;// reset

									// i stuff
									i = i + 1;
									iizd = new Double(zsorted[i - 1] + 0.5);// iiz
																			// =
																			// int(zsorted(i)+0.5);
									iiz = iizd.intValue();

									if (i > ne)
										return;// break;
								}// if(iz>iiz)
							}
						}

					}// have data
						// if ((char)iread == lineSep)
						// lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			//in.close();
		}// try
		catch (Exception exc) {

		}finally {
	        if( null != in ) {
	            try 
	            {
	                in.close();
	            } catch(IOException ex) {
	                // log or fail if you like
	            }
	        }
	    }
		// @@@@@@@@@@@@@@@@@@@@@@@@@READ==========================
		return;

	}

	/**
	 * Internally used. Called by egs_init_user_photon.
	 * @param e e
	 * @return the result
	 */
	public static double egs_KN_sigma0(double e) {
		// "=========================================================================="
		// implicit none;
		// $REAL e;
		// $REAL rm,con,ko,c1,c2,c3,eps1,eps2;
		// data rm/0.5110034/,con/0.1274783851/;

		double rm = 0.0;
		double con = 0.0;
		double ko = 0.0;
		double c1 = 0.0;
		double c2 = 0.0;
		double c3 = 0.0;
		double eps1 = 0.0;
		double eps2 = 0.0;

		rm = 0.5110034;
		con = 0.1274783851;
		double egs_KN_sigma0 = 0.0;
		ko = e / rm;
		if (ko < 0.01) {
			egs_KN_sigma0 = 8. * con / 3.
					* (1 - ko * (2 - ko * (5.2 - 13.3 * ko))) / rm;
			return egs_KN_sigma0;
		}
		c1 = 1. / (ko * ko);
		c2 = 1. - 2.0 * (1.0 + ko) * c1;
		c3 = (1.0 + 2.0 * ko) * c1;
		eps2 = 1.0;
		eps1 = 1. / (1.0 + 2.0 * ko);
		egs_KN_sigma0 = (c1 * (1. / eps1 - 1. / eps2) + c2
				* Math.log(eps2 / eps1) + eps2 * (c3 + 0.5 * eps2) - eps1
				* (c3 + 0.5 * eps1))
				/ e * con;
		return egs_KN_sigma0; // end;
	}

	/**
	 * Not used, egs_KN_sigma0 is better.
	 * @param e e
	 * @return the result
	 */
	public static double egs_KN_sigma1(double e) {
		// "=========================================================================="
		// implicit none;
		// $REAL e;
		double rm = 0.0;
		double con = 0.0;
		double ko = 0.0;
		double c1 = 0.0;
		double c2 = 0.0;
		double c3 = 0.0;
		double eps1 = 0.0;
		double eps2 = 0.0;

		rm = 0.5110034;
		con = 0.1274783851;

		ko = e / rm;
		c1 = 1. / (ko * ko);
		c2 = 1. - 2.0 * (1.0 + ko) * c1;
		c3 = (1.0 + 2.0 * ko) * c1;
		eps2 = 1.;
		eps1 = 1. / (1.0 + 2.0 * ko);
		double egs_KN_sigma1 = c1 * (1. / eps1 - 1. / eps2);
		egs_KN_sigma1 = egs_KN_sigma1 + Math.log(eps2 / eps1) * (c2 - c1) - c2
				* (eps2 - eps1);
		egs_KN_sigma1 = egs_KN_sigma1 + c3 * (eps2 - eps1)
				* (1.0 - 0.5 * (eps1 + eps2));
		egs_KN_sigma1 = egs_KN_sigma1
				+ (eps2 - eps1)
				* (0.5 * (eps1 + eps2) - (eps1 * eps1 + eps2 * eps2 + eps1
						* eps2) / 3.0);
		egs_KN_sigma1 = egs_KN_sigma1 * con;
		return egs_KN_sigma1; // end;
	}
}
