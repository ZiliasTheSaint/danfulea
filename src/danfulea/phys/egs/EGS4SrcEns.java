package danfulea.phys.egs;

import java.io.FileInputStream;

/**
 * Utility class for handling particle sources and particle energy (mono-energetic or spectrum). 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 07 NOV. 2005
 */
public class EGS4SrcEns {
	
		public static EgsQuestion eq;
		public static double $ONE_EPS=0.9999;// "USED TO KEEP THE BEAM INSIDE THE TARGET      "
		////"ENFLAG   FLAG THAT SOURCE ENERGY IS READ HERE (ISOURC=21,22):
		public static int ENFLAG=0;
		public static int iMONOENERGETIC=0;
		public static int iSPECTRUM=1;
		public static int monoindex=0;//0=MONOENERGETIC;1=SPECTRUM
		public static int MONOEN=0;
		public static double ein=0.0;
		public static double ikemev=0.0;//default//INCIDENT KINETIC ENERGY(MEV)
		private static double ikemev_default=1.25;
		private static double ikemev_min=0.001;//"databases for EGSnrc stop at 1 keV"
		private static double ikemev_max=200000.0;//"not sure what real upper limit is"
		public static String enerFilename="";
		private static final String datas="Data";
		private static final String egsDatas="egs";
		private static final String dataspec="spectra";
		private static String file_sep=System.getProperty("file.separator");
		private static String defaultext=".spectrum";
		public static int $NENSRC=300;//     "MAX # OF POINTS IN ENERGY DISTRIBUTION  "
		public static int NENSRC=0;
		public static double[] ENSRCD=new double[$NENSRC+1];//(0:$NENSRC)
		public static double[] SRCPDF=new double[$NENSRC];//($NENSRC)

		public static boolean is_mono=false;
		public static double[] srcpdf_at=new double[$NENSRC];//($NENSRC),
		public static int[] srcbin_at=new int[$NENSRC];//($NENSRC)
		public static double enmin=0.0;
		public static double enmax=0.0;
		public static double Ek_max=0.0;
		public static double Emono=0.0;
		public static int mode=0;
		public static int iNONE=0;
		public static int iINCLUDE=1;
		public static int IOUTSP=0;//0=NONE,1=INCLUDE

		public static double sume1=0.0;public static double sume=0.0;
		//used for read egsinputs file---> We want to read it from GUIs!!
		//@@@@@@@@@@@@@@@SRC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		//RADIUS OF THE BEAM AT THE FRONT OF THE TARGET or OUTER RADIUS of SOURCE for ISOURC = 3:
		public static double RBEAM=0.0;//RADIUS OF THE BEAM AT THE FRONT OF THE TARGET
		public static double RBEAM2=0.0;// RBEAM**2
		public static double xin=0.0;//    X-COORDINATE SELECTED
		public static double yin=0.0;//    Y-COORDINATE SELECTED
		public static double zin=0.0;//    Z-COORDINATE SELECTED
		public static double uin=0.0;//    DIRECTION COSINE WITH RESPECT TO THE X-AXIS
		public static double vin=0.0;//    DIRECTION COSINE WITH RESPECT TO THE Y-AXIS
		public static double win=0.0;//    DIRECTION COSINE WITH RESPECT TO THE Z-AXIS
		public static double wtin=0.0;
		public static double WEIGHT=0.0;// WEIGHT OF THE PARTICLE
		public static int NRCFLG = 10;// => ENTRY IS FROM A FRONT FACE
		public static double R2=0.0;//     XIN**2+yin**2
		public static int IXIN=0;//   RADIUS NUMBER IN WHICH THE SELECTED POINT RESIDES
		public static int irin=0;//REGION NUMBER OF THE ENTRANCE ELEMENT IN THE GEOMETRY
		public static int iqin=0;

		//MINIMUM BEAM RADIUS FOR ISOURC=14 or INNER RADIUS of SOURCE,ISOURCE = 3:
		public static double RMINBM=0.0;// IS MINIMUM BEAM RADIUS FOR ISOURC=14
		//DISTANCE OF THE SOURCE FROM THE MIDDLE OF THE TARGET:
		public static double DISTRH=0.0;// PERPENCICULAR DISTANCE OF THE POINT SOURCE OFF THE
		//"                   SYMMETRY (Z) AXIS OF THE GEOMETRY
		//DISTANCE OF THE POINT SOURCE FROM THE FRONT OF THE TARGET:
		public static double DISTZ=0.0;//  PERPENCICULAR DISTANCE OF THE POINT SOURCE FROM THE
		//"                   FIRST PLANE (ZPLANE(1) IN THE GEOMETRY)
		public static double DISTZ2=0.0;// DISTZ**2
		//FACIAL AREA OF THE TARGET:
		public static double AFACE=0.0;//  AREA TO WHICH THE BEAM IS COLLIMATED ON THE FACE
		public static double RMINSQ=0.0;////   RMINBM**2
		public static int ISOURC=0;//number of the source
		//ACCUMULATES PARTICLE WEIGHTS (WEIGHT^2) FROM THE SOURCE,
		//"          EVENTUALLY GIVES SOLID ANGLE (UNCERTAINTY ON SOLID ANGLE) OF
		//"          THE TARGET
		public static double SCOMEG=0.0;
		public static double SCOMEG2=0.0;
		//PROBABILITY THAT INCIDENT BEAM STRIKES THE FLAT FACE (ISOURC=12):
		public static double PROBFC=0.0;
		//DISTANCE OF THE POINT SOURCE FROM THE BACK OF THE TARGET:
		public static double DISTB=0.0;//  PERPENCICULAR DISTANCE OF THE POINT SOURCE FROM THE
		//"                   LAST PLANE (ZPLANE(NPLANE) IN THE GEOMETRY)
		public static double DISTB2=0.0;// DISTB**2
		//PROBABILITY THAT INCIDENT BEAM STRIKES THE FLAT FACE (ISOURC=12):
		public static double PROBBK=0.0;
		//HALF-HEIGHT OF THE BEAM AT THE CENTER OF THE TARGET:
		public static double XBEAM=0.0;//  HALF-WIDTH WITHIN WHICH TO CHOOSE THE POINT
		//HALF-WIDTH OF THE BEAM AT THE CENTER OF THE TARGET:
		public static double ZBEAM=0.0;//  HALF-DEPTH WITHIN WHICH TO CHOOSE THE POINT
		//Z-AXIS OFFSET OF THE CENTER OF THE TARGET:
		public static double ZCOFST=0.0;// Z OFFSET OF THE CENTER OF THE GEOMETRY
		public static double RCYL1=0.0;//RADIUS OF THE TARGET
		public static int IZ1=0;
		//Z-AXIS OFFSET OF THE SOURCE:
		public static double ZSOFST=0.0;// Z OFFSET OF THE SOURCE
		public static double RCYL2=0.0;//"RCYL2    = RCYL1**2
		//PROBABILITY THAT INCIDENT BEAM STRIKES THE CURVED SIDE (ISOURC=12):
		public static double PROBSD=0.0;
		//RECTANGULAR AREA OF THE TARGET SEEN BY THE BEAM:
		public static double ASIDE=0.0;
		public static double DSTRH2=0.0;//= DISTRH**2
		//INVERSE OF THE RADIAL CUMULATIVE PROBABILITY FUNCTION (ISOURC=20):
		public static int $MXRDIST=1000;//   "# OF POINTS IN RADIAL DISTRIBUTION FIT       "
		public static double[][] RCDFIN=new double[$MXRDIST][2];
		//RADIAL MESH POINTS FOR THE RADIAL PROBABILITY FUNCTION (ISOURC=20):
		public static double[] RDISTF=new double[100];
		//RADIAL PROBABILITY FUNCTION (ISOURC=20):
		public static double[] RPDF=new double[100];
		//RADIAL CUMULATIVE PROBABILITY FUNCTION (ISOURC=20)
		public static double[] RCDF=new double[100];
		public static double AINFLU=0.0;//   INCIDENT FLUENCE
		public static double UINC=0.0;//     INCIDENT X-AXIS DIRECTION COSINE (ISOURC=0)
		public static double VINC=0.0;//     INCIDENT Y-AXIS DIRECTION COSINE (ISOURC=0)
		public static double WINC=0.0;//     INCIDENT Z-AXIS DIRECTION COSINE (ISOURC=0)
		public static int NRDIST=0;//   # RADIAL BINS IN DISTRIBUTION HISTOGRAM (ISOURC=20)
		//SOURCE ENERGY INPUT FROM FULL PHASE-SPACE SOURCE FILE (ISOURC=21,22):
		public static double EINSRC=0.0;
		public static int NSHIST=0;//   NUMBER OF PARTICLES IN SOURCE FILE (ISOURC=21,22)
		public static int NSLEFT=0;//   NUMBER OF PARTICLES REMAINING IN SOURCE FILE (ISOURC=21,22)
		public static double EKSRCM=0.0;//   MAXIMUM ENERGY OF SOURCE (ISOURC=21,22)
		public static String SPCNAM="";//   FILE NAME OF SOURCE DATA FILE (ISOURC=21,22)
		public static double ZSMAX=0.0;//    MAXIMUM of SOURCE Z extent for ISOURC = 3
		public static double ZSMIN=0.0;//    MINIMUM of SOURCE Z extent for ISOURC = 3
		public static int IMODE    =0;// FOR DEFAULT BEAM FILE FORMAT WITH 7 VARIABLES PER RECORD"
		//"         =2 FOR PH-SP FILES CREATED BY BEAM WITH 8 VARIABLES/RECORD, +ZLAST"
		public static int NCASE_PHSP=0;//  TOTAL NUMBER OF PARTICLES STORED IN THE PH-SP FILE"
		//   Stores particle no. to be read from phase space source (ISOURC=21,22)"
		public static int NPHSPN=0;
		//   Stores no. of primary histories read from phsp source (ISOURC=21,22)"
		public static int NHSTRY=0;
		//   NO. of times a particle has crossed the phase space plane (21,22)"
		public static int NPASS =0;
		//   Total number of particles read from phase space source (ISOURC=21,22)"
		public static int NNREAD=0;
		// For ISOURC=21,22, set to 1 if we cannot read no. primary histories"
		public static int DOSE_STAT=0;
		//  number of particles incident from primary source (ISOURC=21,22)"
		public static double NINCSRC=0.0;
		//   no of times to recycle each source particle (ISOURC=21,22)"
		public static int NRCYCL=0;
		// no. of parallel jobs into which a simulation is split(ISOURC=21,22)"
		public static int IPARALLEL=0;
		//   set to a different integer value in the range 1<=PARNUM<=IPARALLEL
		//"         for the IPARALLEL parallel jobs (ISOURC=21,22)
		public static int PARNUM=0;
		//  keeps track of how many times a particle has been recycled
		//"         (ISOURC=21,22)
		public static int CYCLNUM=0;
		//   keeps track of how many times a phsp source has been restarted
		//"         (ISOURC=21,22)
		public static int OUTCNT=0;
		public static double XINOLD=0.0;//   holds XIN for recycled particles (ISOURC=21,22)
		public static double YINOLD=0.0;//   holds yin for recycled particles (ISOURC=21,22)
		public static double ZINOLD=0.0;//   holds zin for recycled particles (ISOURC=21,22)
		public static double UINOLD=0.0;//   holds uin for recycled particles (ISOURC=21,22)
		public static double VINOLD=0.0;//   holds vin for recycled particles (ISOURC=21,22)
		public static double WINOLD=0.0;//   holds win for recycled particles (ISOURC=21,22)
		public static double WEIGHTOLD=0.0;//   holds WEIGHT for recycled particles (ISOURC=21,22)
		public static int IRINOLD=0;//  holds IRIN for recycled particles (ISOURC=21,22)
		public static int NRCFLGOLD=0;//   holds NRCFLG for recycled particles (ISOURC=21,22)
		//   no. of histories left in a parallel run excluding current run
		//"         used for parallel runs with ISOURC=21,22 only
		public static int N_LEFT=0;
		//  no. of previous chunk of parallel run.  Used for
		//"         parallel runs with ISOURC=21,22 only
		public static int N_RUN_CHUNK_OLD=0;
		//  min. value for NPHSPN in chunk of phsp source.  Used for
		//"         parallel runs with ISOURC=21,22 only
		public static int NPHSPN_MIN=0;
		//  max. value for NPHSPN in chunk of phsp source.  Used for
		//"         parallel runs with ISOURC=21,22 only
		public static int NPHSPN_MAX=0;
		// no. of particles in phsp source per chunk of parallel run
		//"         For parallel runs with ISOURC=21,22 only
		public static int P_PER_PHSP_CHUNK=0;
		public static int IPHATIN=0;//  set to 1 if this is a fat photon (from DBS)--ISOURC=23 only
		public static double SVTMP1=0.0;
		public static double SVTMP2=0.0;
		public static double SVTMP3=0.0;
		public static double SVTMP4=0.0;
		public static int NSRCRG=0;
		public static int ISRPT=0;
		public static int ISRPAR=0;

		public static final int parallel_beam_incident_from_the_front=0;
		public static final int point_source_on_axis_incident_from_the_front=1;
		public static final int broad_parallel_beam_incident_from_the_front=2;
		public static final int uniform_isotropically_radiating_disk=3;
		public static final int axial_fluence_for_various_beam_radii=4;

		public static final int parallel_beam_incident_from_the_side=10;
		public static final int point_source_incident_from_the_side=11;
		public static final int point_source_incident_from_off_axis=12;
		public static final int broad_parallel_beam_incident_from_any_angle=13;
		public static final int point_source_on_axis_from_front_but_restricted=14;//FOR COLLIMATOR STUDIES
		public static final int point_source_incident_from_off_axis_with_weights=15;//for ion chamber
		public static final int circular_or_rectangular_source_off_axis=16;

		public static final int parallel_beam_incident_from_the_front_with_radial_distribution=20;
		public static final int full_phase_space_on_front=21;
		public static final int full_phase_space_from_any_angle=22;
		public static final int BEAM_treatment_head_simulation_from_any_angle=23;

		public static int ipart=0;
		public static final int ipart_photon=0;
		public static final int ipart_positron=1;
		public static final int ipart_all=2;
		public static final int ipart_charged=3;
		public static final int ipart_electron=4;

		public static int nsource_option=0;
		public static double[] source_option=new double[9];
		public static final double source_option_min=-999999.;
		public static final double source_option_max=999999.;
		public static final double source_option_default=1.;

		public static final int MODEIN_LOCAL=0;
		public static final int MODEIN_EXTERNAL=1;
		public static final int NRDIST_MIN=1;
		public static final int NRDIST_MAX=100;
		public static final int NRDIST_DEFAULT=1;
		public static final double RDISTF_MIN=0.0;
		public static final double RDISTF_MAX=999999.0;
		public static final double RDISTF_DEFAULT=1.0;
		public static final double RPDF_MIN=0.0;
		public static final double RPDF_MAX=999999.0;
		public static final double RPDF_DEFAULT=1.0;
		public static int RDIST_IOUTSP=0;
		public static final int RDIST_IOUTSP_NONE=0;
		public static final int RDIST_IOUTSP_INCLUDE=1;

		public static int IFPB=0;
		//####################INSIDE
		public static double TEMP1=0.0;
		public static double TEMP2=0.0;
		public static double TEMP3=0.0;
		public static double TEMP4=0.0;
		public static double TEMP5=0.0;
		public static double TEMP6=0.0;
		public static double TEMP7=0.0;
		public static double TEMP8=0.0;
		public static double TEMP9=0.0;
		public static double FNORM=0.0;
		public static int IBNSOK=0;
		public static int IRDIST=0;
		public static double GRIDSZ=0.0;
		public static double AK=0.0;
		public static int MODEIN=0;
		public static double dist_phsp=0.0;
		public static double theta_phsp=0.0;
		public static double cost_phsp=0.0;
		public static double sint_phsp=0.0;
		public static double chamber_c=0.0;
		public static int count_phsp=0;
		public static double xoffset=0.0;
		public static double yoffset=0.0;
		public static String the_beam_code="";
		public static String the_pegs_file="";
		public static String the_input_file="";
		public static int last_nhstry=0;
		public static double min_weight_23=0.0;
		public static double max_weight_23=0.0;
		public static int secret_option_23=0;
		//------------
		public static double cost=0.0;
		public static double sint=0.0;
		public static double d=0.0;
		public static double d2=0.0;
		public static double yo=0.0;
		public static double zo=0.0;
		public static double dz=0.0;
		public static double R=0.0;
		//public static double R2=0.0;
		public static double xmin=0.0;
		public static double delx=0.0;
		public static double ymin=0.0;
		public static double dely=0.0;
		public static double zc=0.0;
		public static double area=0.0;
		public static double w0=0.0;
		public static double count=0.0;
		public static double sumw=0.0;
		public static double sumw2=0.0;
		public static boolean just_fb=false;
		public static boolean just_side=false;
		public static double pi=0.0;
		public static double ro2=0.0;

		public static double d0=0.0;
		public static double d02=0.0;
		public static double dr=0.0;
		public static double theta=0.0;
		public static double xo=0.0;
		public static double a_side=0.0;
		public static double a_fb=0.0;
		public static double atot=0.0;
		public static double wfb=0.0;
		public static double rs=0.0;
		public static double delxs=0.0;
		public static double delys=0.0;
		public static double rs2=0.0;
		public static boolean point_source=false;
		public static boolean do_fb=false;
		public static boolean do_side=false;
		public static boolean is_circle=false;
		public static boolean do_both=false;
		public static double COTANG=0.0;
		public static double FACTOR=0.0;
		//#########################
		///------------------------SRCSPH
		public static double ABEAM=0.0;//    AREA OF THE BEAM ON THE SURFACE
		public static double DISTR=0.0;//   DISTANCE OF THE SOURCE FROM THE MIDDLE OF THE TARGET
		public static double DISTR2=0.0;//   = DISTR**2
		//   = 0 => PARALLEL SOURCE (ISOURC=0,10)
		//"         = 1 => POINT    SOURCE (ISOURC=1,11)
		public static int ISRCTY=0;
		public static int IDSTON=0;//   = 0 RADIAL DISTRIBUTION OFF, = 1 ON (ISOURC=10,11)
		public static double RSHADW=0.0;
		public static double YINP=0.0;
		public static double RHO2=0.0;
		public static double RRHO=0.0;
		public static double ZINP=0.0;
		public static double TNTST2=0.0;
		public static double D=0.0;

		//$$$$$$$$$$MAIN
		public static int $STAT=10;        //"# BINS FOR UNCERTAINTY ANALYSIS"
		public static double[] OMEGIS=new double[$STAT];
		public static int IS=0;
		//-----------------------------------
		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		//from main
		public static int NCASET=0;

		/**
		 * Reset all global variables for re-use.
		 */
		public static void reset()
		{
			ENFLAG=0;iMONOENERGETIC=0;iSPECTRUM=1;monoindex=0;ikemev=0.0;MONOEN=0;
			ikemev_default=1.25;ikemev_min=0.001;ikemev_max=200000.0;
			enerFilename="";//$NENSRC=300
			;NENSRC=0;ENSRCD=new double[$NENSRC+1];//(0:$NENSRC)
			SRCPDF=new double[$NENSRC];is_mono=false;srcpdf_at=new double[$NENSRC];//($NENSRC),
			srcbin_at=new int[$NENSRC];//($NENSRC)
			enmin=0.0;enmax=0.0;Ek_max=0.0;Emono=0.0;mode=0;iNONE=0;iINCLUDE=1;IOUTSP=0;
			sume1=0.0;sume=0.0;
			//=========================SRC===============================================
			RBEAM=0.0;RBEAM2=0.0;xin=0.0;yin=0.0;zin=0.0;uin=0.0;vin=0.0;win=0.0;wtin=0.0;
			WEIGHT=0.0;NRCFLG = 10;R2=0.0;IXIN=0;irin=0;RMINBM=0.0;DISTRH=0.0;ein=0.0;
			DISTZ=0.0;DISTZ2=0.0;AFACE=0.0;RMINSQ=0.0;ISOURC=0;SCOMEG=0.0;SCOMEG2=0.0;
			PROBFC=0.0;DISTB=0.0;DISTB2=0.0;PROBBK=0.0;XBEAM=0.0;ZBEAM=0.0;ZCOFST=0.0;
			RCYL1=0.0;IZ1=0;ZSOFST=0.0;RCYL2=0.0;PROBSD=0.0;ASIDE=0.0; DSTRH2=0.0;
			//$MXRDIST=1000;
			RCDFIN=new double[$MXRDIST][2];RDISTF=new double[100];
			RPDF=new double[100];RCDF=new double[100];AINFLU=0.0;UINC=0.0;VINC=0.0;
			WINC=0.0;NRDIST=0;EINSRC=0.0;NSHIST=0;NSLEFT=0;EKSRCM=0.0;SPCNAM="";
			ZSMAX=0.0;ZSMIN=0.0;IMODE    =0;NCASE_PHSP=0;NPHSPN=0;NHSTRY=0;
			NPASS =0;NNREAD=0;DOSE_STAT=0;NINCSRC=0.0;NRCYCL=0;IPARALLEL=0;PARNUM=0;
			CYCLNUM=0;OUTCNT=0;XINOLD=0.0;YINOLD=0.0;ZINOLD=0.0;UINOLD=0.0;VINOLD=0.0;
			WINOLD=0.0;WEIGHTOLD=0.0;IRINOLD=0;NRCFLGOLD=0;N_LEFT=0;N_RUN_CHUNK_OLD=0;
			NPHSPN_MIN=0;NPHSPN_MAX=0;P_PER_PHSP_CHUNK=0;IPHATIN=0;SVTMP1=0.0;SVTMP2=0.0;
			SVTMP3=0.0;SVTMP4=0.0;NSRCRG=0;ISRPT=0;ISRPAR=0;
			ABEAM=0.0;DISTR=0.0;DISTR2=0.0;ISRCTY=0;IDSTON=0;
			IFPB=0;

			TEMP1=0.0;TEMP2=0.0;TEMP3=0.0;TEMP4=0.0;TEMP5=0.0;TEMP6=0.0;TEMP7=0.0;TEMP8=0.0;
			TEMP9=0.0;FNORM=0.0;IBNSOK=0;MODEIN=0;dist_phsp=0.0;theta_phsp=0.0;cost_phsp=0.0;
			sint_phsp=0.0;chamber_c=0.0;count_phsp=0;xoffset=0.0;yoffset=0.0;the_beam_code="";
			the_pegs_file="";the_input_file="";last_nhstry=0;min_weight_23=0.0;max_weight_23=0.0;
			secret_option_23=0;

			ipart=0;ISOURC=0;source_option=new double[9];nsource_option=0;iqin=0;RDIST_IOUTSP=0;
			cost=0.0;sint=0.0;d=0.0;d2=0.0;yo=0.0;zo=0.0;dz=0.0;R=0.0;xmin=0.0;
			delx=0.0;ymin=0.0;dely=0.0;zc=0.0;area=0.0;w0=0.0;count=0.0;sumw=0.0;
			sumw2=0.0;just_fb=false;just_side=false;pi=0.0;ro2=0.0;
			d0=0.0;d02=0.0;dr=0.0;theta=0.0;xo=0.0;a_side=0.0;a_fb=0.0;atot=0.0;
			wfb=0.0;rs=0.0;delxs=0.0;delys=0.0;rs2=0.0;point_source=false;do_fb=false;
			do_side=false;is_circle=false;do_both=false;

			NCASET=0;COTANG=0.0;FACTOR=0.0;GRIDSZ=0.0;AK=0.0;IRDIST=0;
			RSHADW=0.0;YINP=0.0;RHO2=0.0;RRHO=0.0;ZINP=0.0;TNTST2=0.0;D=0.0;

			OMEGIS=new double[$STAT];IS=0;
		}

		/**
		 * ENSRC uses the alias sampling technique to sample the source energy from a histogrammed input energy spectrum. 
		 * All data needed are stored in the static arrays ensrcd,srcpdf,srcpdf_at and srcbin_at. 
		 * Call this routine from anywhere in your user code (but before the first SHOWER call) to get input data.
		 */
		public static void ENSRC()
		{

			//ENFLAG must be known here-> via SRC!!!! or $NENSRC=300 can be modified && eq!!!!!
			if(ENFLAG==1)
			{//"FULL PHASE-SPACE INFORMATION READ PREVIOUSLY"
			    MONOEN = 2;
			    return;
			}

			//before call ENSRC EGS4SrcEns.monoindex=# must be set!!!
			MONOEN = monoindex;
			if(MONOEN == 0)
			{
				//if here before call ENSRC, EGS4SrcEns.ikemev=# must be set!!!
				if (ikemev<ikemev_min || ikemev>ikemev_max)
				{
					//use default
					ikemev=ikemev_default;
				}
				is_mono=true;ein=ikemev;Emono=ein;
			}
			else//energy spectrum
			{
				//if here before call ENSRC, EGS4SrcEns.enerFilename=# must be set!!!IOUTSP!!
				is_mono=false;
				readSpectrum(enerFilename);

		      	EGS4.seqStr="    HAVE READ"+EGS4.format(NENSRC,5)+" INPUT ENERGY BINS FROM FILE";
			  	if(EGS4.iprint>1)
			  		eq.printSequence(EGS4.seqStr);

			    if(mode==0)
			    {
			      	EGS4.seqStr="      Counts/bin assumed";
				  	if(EGS4.iprint>1)
				  		eq.printSequence(EGS4.seqStr);
				}
				else if(mode == 1)
				{
			      	EGS4.seqStr="      Counts/MeV assumed";
				  	if(EGS4.iprint>1)
				  		eq.printSequence(EGS4.seqStr);

				    for(int IB=1;IB<=NENSRC;IB++)
				    {
						SRCPDF[IB-1] = SRCPDF[IB-1]*(ENSRCD[IB]-ENSRCD[IB-1]);
					}
				}// "end mode = 1 block"
				else
				{
			      	EGS4.seqStr="*****MODE not 0 or 1 in spectrum file? **";
				  	if(EGS4.iprint>1)
				  		eq.printSequence(EGS4.seqStr);
			    }

			    ein=ENSRCD[NENSRC];//"SET TO MAX ENERGY FOR SOME CHECKS"

		      	EGS4.seqStr="    ENERGY RANGES FROM"+EGS4.format(enmin,10,true)+
		      	" MeV TO"+EGS4.format(ein,12,true)+" MeV";
			  	if(EGS4.iprint>1)
			  		eq.printSequence(EGS4.seqStr);

				enmax = ein;

				if(IOUTSP!=0 && IOUTSP!=1)
				{
					IOUTSP=0;//default value
				}
			}
		}

		/**
		 * Initializatin routine (if spectrum, prepare alias sampling). To be called after NSRC and before the first SHOWER call.
		 */
		public static void ENSRC1()
		{
			if( is_mono ) return;

			//" Rewritten by IK: to guarantee that the input spectrum is exactly sampled,"
			//"                  The alias sampling technique is employed
			//"                  needs prepare_alias_sampling which is in nrcaux.mortran

			//" Check that enmin < ensrcd(1) "
			if( enmin >= ENSRCD[1] )
			{
			    EGS4.STOPPROGRAM=true;
				EGS4.seqStr=" Bad spectrum: minimum energy is > top of first bin! ";
				eq.printSequence(EGS4.seqStr);
				return;
			    //stop;
			}

			//prepare_alias_sampling(nensrc,srcpdf,srcpdf_at,srcbin_at);
			EGS4.prepare_alias_sampling(NENSRC,SRCPDF,srcpdf_at,srcbin_at);

		}

		/**
		 * Produces a summary of the input data.
		 */
		public static void ENSRCO()
		{
			if( is_mono )
			{
		      	EGS4.seqStr="Mono-energy: "+EGS4.format(Emono,10,true)+" MeV";
			  	if(EGS4.iprint>1)
			  		eq.printSequence(EGS4.seqStr);

			  return;
			}
			if(ENFLAG == 1)
			{// "this is a phase space input, print nothing here"
				return;
			}
			//WRITE(IOUT,105) FILNAM(:lnblnk1(FILNAM)),SPEC_TITLE(:lnblnk1(SPEC_TITLE));
			//105 FORMAT(18x,' Spectrum file and title:'/18x,A/18x,A);
	      	EGS4.seqStr=" Spectrum file:  "+enerFilename;
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

			sume1 = 0.0; sume = 0.0;
			for(int ib=1;ib<=NENSRC;ib++)
			{
			  sume = sume + SRCPDF[ib-1];
			  sume1 = sume1 + SRCPDF[ib-1]*(ENSRCD[ib]+ENSRCD[ib-1])/2.0;
			}
			//WRITE(IOUT,'(20x,a,f10.4)') ' Average spectrum energy is ',sume1/sume;
			double rap=sume1/sume;
	      	EGS4.seqStr=" Average spectrum energy is "+EGS4.format(rap,10,true);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

			if(IOUTSP==1)
			{
			  	if(mode==0)
			  	{
			      	EGS4.seqStr=" Counts/bin assumed";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			   	else
			   	{
			      	EGS4.seqStr=" Counts/MeV assumed";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			   	//WRITE(IOUT,110)NENSRC,enmin;
			   	for(int ib=1;ib<=NENSRC;ib++)
			   	{
			      	EGS4.seqStr=EGS4.format(ib,5)+" "+EGS4.format(SRCPDF[ib-1],10,true)+" "+
			      	EGS4.format(srcpdf_at[ib-1],10,true)+" "+EGS4.format(srcbin_at[ib-1],10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
	//System.out.println(EGS4.seqStr);
			   	  //WRITE(IOUT,120)IB,ensrcd(ib),srcpdf(ib),srcpdf_at(ib),srcbin_at(ib);
				}
			}
			return;

		}

		/**
		 * Called from main for each history to sample the spectrum (MONOEN=1) or return the source energy if mono-energetic beam.
		 * @return the source energy
		 */
		public static double ENSRCH()
		{
			double result =0.0;

			if( is_mono )
			{
				result = Emono; return result;
			}

			//result = EGS4.alias_sample(nensrc,ensrcd,srcpdf_at,srcbin_at);
			result = EGS4.alias_sample(NENSRC,ENSRCD,srcpdf_at,srcbin_at);
			return result;
		}

		/**
		 * Puts maximum spectrum energy (source energy for mono-energetic beams) into EK_MAX variable an also return it.
		 * @return the result
		 */
		public static double ENSRC_EMAX()//(Ek_max);
		{
			double result =0.0;

			if( is_mono ) { Ek_max = Emono; }
			else { Ek_max = ENSRCD[NENSRC]; }

			result=Ek_max;
			return result;
		}

		/**
		 * Given the energy file, read the spectrum.
		 * @param enerfile enerfile
		 */
		private static void readSpectrum(String enerfile)
		{
			String filename=datas+file_sep+egsDatas+file_sep+dataspec+file_sep+enerfile+defaultext;
			int iread =0;int lnr =0;//data number
			int lnrr =0;//line number
			int indx=1;
			StringBuffer desc=new StringBuffer();
			boolean haveData=false;
			//String equals="=";
			char comma=',';
			char lineSep = '\n';//System.getProperty("line.separator").charAt(0);

			boolean enB=false;boolean srB=false;

		try
		{
			FileInputStream in = new FileInputStream(filename);

	       	while ((iread = in.read()) != -1)
			{
				if (!Character.isWhitespace((char)iread)&&((char)iread!=comma))
				{
					desc.append((char)iread);
					haveData=true;
				}
				else
				{
					if (haveData)// we have data
					{
						//System.out.println("ok");
						haveData=false;//reset
						if(lnrr!=0)//for skip first line
						{
							lnr++;
							//READ(9,*) nensrc,ensrcd(0),mode;
							if(lnr==1)
							{
								String s=desc.toString();
								NENSRC=EGS4.stringToInt(s);//System.out.println(NENSRC);
							}
							if(lnr==2)
							{
								String s=desc.toString();
								ENSRCD[0]=EGS4.stringToDouble(s);//System.out.println(ENSRCD[0]);
							}
							if(lnr==3)
							{
								String s=desc.toString();
								mode=EGS4.stringToInt(s);//System.out.println(mode);
							    enmin = ENSRCD[0];
							    if(NENSRC > $NENSRC)
							    {
							       //OUTPUT NENSRC,$NENSRC;
							       //(//' ********** Asked for too many energy bins******'/
							       //' NENSRC =',I4, ' reduced to max allowed =',I4/1x,30('*')//);
					      	EGS4.seqStr=" ********** Asked for too many energy bins******";
						  	if(EGS4.iprint>1)
						  		eq.printSequence(EGS4.seqStr);
					      	EGS4.seqStr=" NENSRC ="+NENSRC+" reduced to max allowed ="+$NENSRC;
						  	if(EGS4.iprint>1)
						  		eq.printSequence(EGS4.seqStr);


							       NENSRC = $NENSRC;
								}

							}
							else if(lnr>3)
							{
					//READ(9,*)(ENSRCD(IB),SRCPDF(IB),IB=1,NENSRC);
								String s=desc.toString();
								if(!enB)
								{
									ENSRCD[indx]=EGS4.stringToDouble(s);
									enB=true;
								}
								else if(!srB)
								{
									SRCPDF[indx-1]=EGS4.stringToDouble(s);
	 //System.out.println(ENSRCD[indx]+"    "+SRCPDF[indx-1]);

	 								indx++;
	 								if(indx==NENSRC+1)
	 									break;
									srB=true;enB=false;srB=false;
								}
							}

						}

					}//have data
					if ((char)iread == lineSep)
						lnrr++;
					desc=new StringBuffer();
				}//else
			}//main while
			in.close();
		}//try
		catch (Exception exc)
		{

		}

		}
	
	//"RANDOMLY SELECTS A POINT WITHIN A CIRCLE IN THE NRCC RZ GEOMETRY
	//"    INPUTS  RBEAM  RADIUS WITHIN WHICH TO CHOOSE THE POINT
	//"            RBEAM2 RBEAM**2
	//"    OUTPUTS XIN    X-COORDINATE SELECTED
	//"            yin    Y-COORDINATE SELECTED
	//"            R2     XIN**2+yin**2
	//"            IXIN   RADIUS NUMBER IN WHICH THE SELECTED POINT RESIDES
	//"*** THIS MACRO IS DESIGNED TO WORK WITH THE STANDARD NRCC RZ GEOMETRY ***
	//"*** ROUTINE ONLY                                                      ***
		/**
		 * RANDOMLY SELECTS A POINT WITHIN A CIRCLE IN THE RZ GEOMETRY. <p>
		 * INPUTS (GLOBAL VARIABLES):  RBEAM  RADIUS WITHIN WHICH TO CHOOSE THE POINT; RBEAM2 RBEAM**2<p>
		 * OUTPUTS(GLOBAL VARIABLES):  xin    X-COORDINATE SELECTED; yin    Y-COORDINATE SELECTED; R2     XIN**2+yin**2; 
		 * IXIN   RADIUS NUMBER IN WHICH THE SELECTED POINT RESIDES
		 */
		public static void CHOOSE_POINT_IN_CIRCLE()
		{
		  while(true)
		  {
		    xin=EGS4.random01();xin=(2.0*xin-1.0)*RBEAM;
		    yin=EGS4.random01();yin=(2.0*yin-1.0)*RBEAM;
		    R2=xin*xin+yin*yin;

		    if(R2<=RBEAM2)
		    	break;
		  }// UNTIL R2.LE.RBEAM2;
		  for(int IX=1;IX<=EGS4Geom.NR;IX++)
		  {
			  IXIN=IX;
			  if(R2<=EGS4Geom.CYRAD2[IX-1])break;
		  }
		}
	//"SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING A
	//"FRONT FACE IN THE NRCC RZ GEOMETRY
	//"    INPUTS  RBEAM  RADIUS WITHIN WHICH TO CHOOSE THE POINT
	//"            RBEAM2 RBEAM**2
	//"            RMINBM IS MINIMUM BEAM RADIUS FOR ISOURC=14
	//"            DISTRH PERPENCICULAR DISTANCE OF THE POINT SOURCE OFF THE
	//"                   SYMMETRY (Z) AXIS OF THE GEOMETRY
	//"            DISTZ  PERPENCICULAR DISTANCE OF THE POINT SOURCE FROM THE
	//"                   FIRST PLANE (ZPLANE(1) IN THE GEOMETRY)
	//"            DISTZ2 DISTZ**2
	//"            AFACE  AREA TO WHICH THE BEAM IS COLLIMATED ON THE FACE
	//"    OUTPUTS XIN    X-COORDINATE SELECTED
	//"            yin    Y-COORDINATE SELECTED
	//"            zin    Z-COORDINATE SELECTED
	//"            uin    DIRECTION COSINE WITH RESPECT TO THE X-AXIS
	//"            vin    DIRECTION COSINE WITH RESPECT TO THE Y-AXIS
	//"            win    DIRECTION COSINE WITH RESPECT TO THE Z-AXIS
	//"            IRIN   REGION NUMBER OF THE ENTRANCE ELEMENT IN THE GEOMETRY
	//"            WEIGHT WEIGHT OF THE PARTICLE
	//"            NRCFLG = 10 => ENTRY IS FROM A FRONT FACE
	//"*** THIS MACRO IS DESIGNED TO WORK WITH THE STANDARD NRCC RZ GEOMETRY ***
	//"*** ROUTINE ONLY                                                      ***
		/**
		 * SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING A FRONT FACE IN THE RZ GEOMETRY.
		 */
		public static void ENTRY_FRONT_FACE()
		{
			  double D=0.0;

			  CHOOSE_POINT_IN_CIRCLE();
			  irin=2+(IXIN-1)*EGS4Geom.NZ;
			  zin=EGS4Geom.ZPLANE[0];//ZPLANE(1);
			  if(ISOURC==14 && R2<=RMINSQ)
			  {
			   	WEIGHT=0.0;//"I.E. GET HOWFAR TO TERMINATE HISTORY"
			  }
			  else
			  {
			   	if(DISTRH==0.0)
			   	{
					D=Math.sqrt(R2+DISTZ*DISTZ);//point source on symmetry z axis!!
				}
			   	else
			   	{
					//point source at DISTRH distance from symmetry axis,
					//DISTRH is formaly put on Y direction!!!
					D=Math.sqrt(R2+DISTRH*(DISTRH-2.0*yin)+DISTZ*DISTZ);//OK!
				}
			   	uin=xin/D;vin=(yin-DISTRH)/D;win=DISTZ/D;
			   	NRCFLG=10;
			   	WEIGHT=AFACE*DISTZ/(D*D*D)/PROBFC;
			   	SCOMEG=SCOMEG+WEIGHT;
			   	SCOMEG2=SCOMEG2+WEIGHT*WEIGHT;
		      }
		}
	//"SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING A
	//"BACK FACE IN THE NRCC RZ GEOMETRY
	//"    INPUTS  RBEAM  RADIUS WITHIN WHICH TO CHOOSE THE POINT
	//"            RBEAM2 RBEAM**2
	//"            DISTRH PERPENCICULAR DISTANCE OF THE POINT SOURCE OFF THE
	//"                   SYMMETRY (Z) AXIS OF THE GEOMETRY
	//"            DISTB  PERPENCICULAR DISTANCE OF THE POINT SOURCE FROM THE
	//"                   LAST PLANE (ZPLANE(NPLANE) IN THE GEOMETRY)
	//"            DISTB2 DISTB**2
	//"            AFACE  AREA TO WHICH THE BEAM IS COLLIMATED ON THE FACE
	//"    OUTPUTS XIN    X-COORDINATE SELECTED
	//"            yin    Y-COORDINATE SELECTED
	//"            zin    Z-COORDINATE SELECTED
	//"            uin    DIRECTION COSINE WITH RESPECT TO THE X-AXIS
	//"            vin    DIRECTION COSINE WITH RESPECT TO THE Y-AXIS
	//"            win    DIRECTION COSINE WITH RESPECT TO THE Z-AXIS
	//"            IRIN   REGION NUMBER OF THE ENTRANCE ELEMENT IN THE GEOMETRY
	//"            WEIGHT WEIGHT OF THE PARTICLE
	//"            NRCFLG = 30 => ENTRY IS FROM THE BACK FACE
	//"*** THIS MACRO IS DESIGNED TO WORK WITH THE STANDARD NRCC RZ GEOMETRY ***
	//"*** ROUTINE ONLY  ***
		/**
		 * SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING A BACK FACE IN THE RZ GEOMETRY.
		 */
		public static void ENTRY_BACK_FACE()
		{
			double D=0.0;

			  CHOOSE_POINT_IN_CIRCLE();
			  irin=1+IXIN*EGS4Geom.NZ;
			  zin=EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1];//ZPLANE(NPLANE);
			  if(DISTRH==0.0)
			  {
				  D=Math.sqrt(R2+DISTB2);
			  }
			  else
			  {
				  D=Math.sqrt(R2+DISTRH*(DISTRH-2.0*yin)+DISTB2);
			  }
			  uin=xin/D;vin=(yin-DISTRH)/D;win=-DISTB/D;
			  NRCFLG=30;
			  WEIGHT=AFACE*DISTB/(D*D*D)/PROBBK;
			  SCOMEG=SCOMEG+WEIGHT;
			  SCOMEG2=SCOMEG2+WEIGHT*WEIGHT;
		}
	//"RANDOMLY SELECTS A POINT WITHIN A RECTANGLE IN THE NRCC RZ GEOMETRY
	//"    INPUTS  XBEAM  HALF-WIDTH WITHIN WHICH TO CHOOSE THE POINT
	//"            ZBEAM  HALF-DEPTH WITHIN WHICH TO CHOOSE THE POINT
	//"            ZCOFST Z OFFSET OF THE CENTER OF THE GEOMETRY
	//"    OUTPUTS XIN    X-COORDINATE SELECTED
	//"            yin    Y-COORDINATE SELECTED
	//"            zin    Z-COORDINATE SELECTED
	//"            IRIN   REGION NUMBER IN WHICH THE SELECTED POINT RESIDES
	//"*** THIS MACRO IS DESIGNED TO WORK WITH THE STANDARD NRCC RZ GEOMETRY ***
	//"*** ROUTINE ONLY                                                      ***
		/**
		 * RANDOMLY SELECTS A POINT WITHIN A RECTANGLE IN THE RZ GEOMETRY.<p>
		 * INPUTS(GLOBAL VARIABLES):  XBEAM  HALF-WIDTH WITHIN WHICH TO CHOOSE THE POINT; ZBEAM  HALF-DEPTH WITHIN WHICH TO CHOOSE THE POINT; 
		 * ZCOFST Z OFFSET OF THE CENTER OF THE GEOMETRY;<p>
		 * OUTPUTS(GLOBAL VARIABLES) xin    X-COORDINATE SELECTED; yin    Y-COORDINATE SELECTED; zin    Z-COORDINATE SELECTED; 
		 * IRIN   REGION NUMBER IN WHICH THE SELECTED POINT RESIDES.
		 */
		public static void CHOOSE_POINT_IN_RECTANGLE()
		{
			  xin=EGS4.random01();xin=(2.0*xin-1.0)*XBEAM;
			  zin=EGS4.random01();zin=(2.0*zin-1.0)*ZBEAM+ZCOFST;
			  yin=Math.sqrt((RCYL1-xin)*(RCYL1+xin));

			  for(int IZ=2;IZ<=EGS4Geom.NPLANE;IZ++)
			  {
				  IZ1=IZ;
				  if(zin<=EGS4Geom.ZPLANE[IZ-1])
				  	break;
			  }
			  irin=(EGS4Geom.NR-1)*EGS4Geom.NZ+IZ1;
		}
	//"SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING THE
	//"SIDE IN THE NRCC RZ GEOMETRY
	//"    INPUTS  XBEAM  HALF-WIDTH WITHIN WHICH TO CHOOSE THE POINT
	//"            ZBEAM  HALF-DEPTH WITHIN WHICH TO CHOOSE THE POINT
	//"            DISTRH Y OFFSET OF POINT SOURCE
	//"            ZCOFST Z OFFSET OF THE CENTER OF THE GEOMETRY
	//"            ZSOFST Z OFFSET OF THE SOURCE
	//"    OUTPUTS XIN    X-COORDINATE SELECTED
	//"            yin    Y-COORDINATE SELECTED
	//"            zin    Z-COORDINATE SELECTED
	//"            uin    DIRECTION COSINE WITH RESPECT TO THE X-AXIS
	//"            vin    DIRECTION COSINE WITH RESPECT TO THE Y-AXIS
	//"            win    DIRECTION COSINE WITH RESPECT TO THE Z-AXIS
	//"            IRIN   REGION NUMBER OF THE ENTRANCE ELEMENT IN THE GEOMETRY
	//"            WEIGHT WEIGHT OF THE PARTICLE
	//"            NRCFLG = 20 => ENTRY IS FROM A FRONT FACE
	//"*** THIS MACRO IS DESIGNED TO WORK WITH THE STANDARD NRCC RZ GEOMETRY ***
	//"*** ROUTINE ONLY                                                      ***
		/**
		 * SPECIFIES ALL THE INPUT PHASE SPACE FOR A POINT SOURCE IRRADIATING THE SIDE IN THE RZ GEOMETRY.
		 */
		public static void ENTRY_SIDE()
		{
			double D=0.0;
			  CHOOSE_POINT_IN_RECTANGLE();
			  D=Math.sqrt(RCYL2+DSTRH2-2.0*DISTRH*yin+(zin-ZSOFST)*(zin-ZSOFST));
			  uin=xin/D;vin=(yin-DISTRH)/D;win=(zin-ZSOFST)/D;
			  WEIGHT=ASIDE*(DISTRH-RCYL2/yin)/(D*D*D)/PROBSD;
			  SCOMEG=SCOMEG+WEIGHT;
			  SCOMEG2=SCOMEG2+WEIGHT*WEIGHT;
			  NRCFLG=20;
		}
	//"*******************************************************************************
	//"                            SOURCE INPUT (Rev 1.14)
	//"                            **********************
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
	//"  (for SOURCE 21,22,23)  all    (2)  include all of the particles
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
	//;
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
	//;
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
	//"-------------------------------------------------------------------------------
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
	//;
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
	//"  SOURCE OPTIONS                (M4)  IMODE, NRCYCL, IPARALLEL, PARNUM
	//"
	//"               IMODE              0=> 7 variables/record: X,Y,U,V,E,WT,LATCH
	//"                                  2=> 8 variables/record: the above + ZLAST
	//"              NRCYCL     Number of times to recycle each particle in a phase
	//"                         space source.  Each particle in the phase space
	//"                         file is used a total of NRCYCL+1 times before
	//"                         going on to the next particle.
	//"                         If NCASE > no. of particles in the phase space file,
	//"                         then use of NRCYCL is essential for accurate
	//"                         statistics.
	//"                         If NRCYCL is set > 0, then the user-input value is
	//"                         used.
	//"                         If NRCYCL is set <=0 then NRCYCL is automatically
	//"                         calculated to use the entire phase space file with no
	//"                         restarts, unless INCIDENT PARTICLE= positron.
	//"                         The calculated NRCYCL does not take into
	//"                         account particles that are rejected because they
	//"                         miss the geometry.  If the automatically-calculated
	//"                         value of still results in restarts, then use the
	//"                         following guideline for accurate statistics:
	//"                         1. If there is only one restart and only a small
	//"                            fraction of the source is re-used on the second
	//"                            pass, the effect on statistics is unlikely to
	//"                            be significant.
	//"                         2. If there is one restart and a significant portion
	//"                            of the source is re-used on the second pass, set
	//"                            NRCYCL=2*NRCYCL_prev+1, where NRCYCL_prev is the
	//"                            previous automatically-calculated value
	//"                         3. If the source is restarted more than once, try
	//"                            setting NRCYCL=(no. of restarts+1)*NRCYCL_prev+1
	//"            IPARALLEL   set >1 if you are distributing the job among
	//"                        IPARALLEL machines.  IPARALLEL is used with PARNUM
	//"                        (see below) to partition a phase space source into
	//"                        IPARALLEL equal parts.
	//"            PARNUM      For each of the IPARALLEL parallel jobs, PARNUM
	//"                        should have a different integer value in the range
	//"                        1<=PARNUM<=IPARALLEL.  The partition of the phase
	//"                        space source that is used for a particular job is
	//"                        then given by:
	//"                           (PARNUM-1)*(NCASE_PHSP/IPARALLEL)<NPHSPN<=
	//"                                            (PARNUM)*(NCASE_PHSP/IPARALLEL)
	//"                        where NCASE_PHSP is the total number of particles in
	//"                        the phsp source and NPHSPN is the particle no. chosen.
	//"
	//"  FILSPC                        (C)   filename (with ext) contains
	//"                                      phase space information
	//"                                      (maximum of 80 characters)
	//"                                      (assigned to unit 42)
	//"
	//;
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
	//"  SOURCE OPTIONS                (M9)  IMODE, DIST, ANGLE, ZOFFSET, NRCYCL,
	//"                                      IPARALLEL, PARNUM, XOFFSET, YOFFSET
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
	//"
	//"              NRCYCL     Number of times to recycle each particle in a phase
	//"                         space source.  Each particle in the phase space
	//"                         file is used a total of NRCYCL+1 times before
	//"                         going on to the next particle.
	//"                         If NCASE > no. of particles in the phase space file,
	//"                         then use of NRCYCL is essential for accurate
	//"                         statistics.
	//"                         If NRCYCL is set > 0, then the user-input value is
	//"                         used.
	//"                         If NRCYCL is set <=0 then NRCYCL is automatically
	//"                         calculated to use the entire phase space file with no
	//"                         restarts, unless INCIDENT PARTICLE= positron.
	//"                         The calculated NRCYCL does not take into
	//"                         account particles that are rejected because they
	//"                         miss the geometry.  If the automatically-calculated
	//"                         value of still results in restarts, then use the
	//"                         following guideline for accurate statistics:
	//"                         1. If there is only one restart and only a small
	//"                            fraction of the source is re-used on the second
	//"                            pass, the effect on statistics is unlikely to
	//"                            be significant.
	//"                         2. If there is one restart and a significant portion
	//"                            of the source is re-used on the second pass, set
	//"                            NRCYCL=2*NRCYCL_prev+1, where NRCYCL_prev is the
	//"                            previous automatically-calculated value
	//"                         3. If the source is restarted more than once, try
	//"                            setting NRCYCL=(no. of restarts+1)*NRCYCL_prev+1
	//"            IPARALLEL   set >1 if you are distributing the job among
	//"                        IPARALLEL machines.  IPARALLEL is used with PARNUM
	//"                        (see below) to partition a phase space source into
	//"                        IPARALLEL equal parts.
	//"            PARNUM      For each of the IPARALLEL parallel jobs, PARNUM
	//"                        should have a different integer value in the range
	//"                        1<=PARNUM<=IPARALLEL.  The partition of the phase
	//"                        space source that is used for a particular job is
	//"                        then given by:
	//"                           (PARNUM-1)*(NCASE_PHSP/IPARALLEL)<NPHSPN<=
	//"                                            (PARNUM)*(NCASE_PHSP/IPARALLEL)
	//"                        where NCASE_PHSP is the total number of particles in
	//"                        the phsp source and NPHSPN is the particle no. chosen.
	//"    XOFFSET,YOFFSET     X and Y offset of phase space plane (cm).
	//"                        Offset will be applied before rotating the source.
	//"
	//"  FILSPC                        (C)   filename (with ext) contains
	//"                                      phase space information
	//"                                      (maximum of 80 characters)
	///"                                      (assigned to unit 42)
	//"
	//"-------------------------------------------------------------------------------
	//"
	//"                     >>>>>>>> SOURCE  23 <<<<<<<<
	//"
	//"    BEAM TREATMENT HEAD SIMULATION AS SOURCE INCIDENT FROM AND ANGLE,  ""toc:
	//"    INSIDE OR OUTSIDE PHANTOM                                          ""toc:
	//"
	//"    PARTICLES ARE READ DIRECTLY FROM A BEAM SIMULATION COMPILED AS A
	//"    SHARED LIBRARY.  Particles are read at the scoring plane in
	//"    the BEAM simulation (although no phase space file is scored) and are
	//"    tranlated/rotated by the inputs DIST, ANGLE, XOFFSET, YOFFSET, ZOFFSET,
	//"    described below.  Then it is checked
	//"    whether they are already inside the geometry. If yes, the region index
	//"    is determined and the shower intiated. If not, it is checked whether
	//"    the particle trajectory will intersect the geometry (assuming that the
	//"    geometry is surrounded by vacuum). If not, the particle is rejected and
	//"    the next one taken from the BEAM simulation (more histories are run in
	//"    the BEAM simulation if required).  If yes, the particle
	//"    is placed on the entry point and the shower is initiated.
	//"
	//"  BEAM CODE                     (C)  The name of the accelerator code being
	//"                                     used as a source including the BEAM_
	//"                                     prefix (ie BEAM_accelname).  This code
	//"                                     must have been compiled as a shared
	//"                                     library (see the BEAM manual for more
	//"                                     details) and exist as
	//"                                     libBEAM_accelname.so (for Linux/Unix) or
	//"                                     libBEAM_accelname.dll (for Windows) in
	//"                                     directory $EGS_HOME/bin/config.
	//"
	//"  INPUT FILE                    (C)  The name of a working input file
	//"                                     (no .egsinp extension) for
	//"                                     the BEAM code BEAM_accelname.  This
	//"                                     input file must specify output of a
	//"                                     phase space file at one scoring plane.
	//"                                     Particles that would have been scored
	//"                                     in the phase space file are extracted
	//"                                     and used as the incident particles in
	//"                                     the DOSXYZ simulation instead.  The
	//"                                     input file must exist in the directory
	//"                                     $EGS_HOME/BEAM_accelname.
	//"
	//"  PEGS FILE                     (C)  The name of the pegs4 data set (no
	//"                                     .pegs4dat extension) used
	//"                                     by BEAM_accelname with the input file
	//"                                     specified by INPUT FILE.  The pegs4
	//"                                     data set must exist in either
	//"                                     $HEN_HOUSE/pegs4/data or in
	//"                                     $EGS_HOME/pegs4/data.
	//"
	//"  WEIGHT WINDOW                 (M2)  MIN_WEIGHT_23, MAX_WEIGHT_23
	//"
	//"               MIN_WEIGHT_23         Min. weight of particles to use from
	//"                                     the BEAM simulation (defaults to -1E30)
	//"               MAX_WEIGHT_23         Max. weight of particles to use from
	///"                                     the BEAM simulation (defaults to 1E30)
	//"
	//"  SOURCE OPTIONS                (M5)  DIST, ANGLE, ZOFFSET, XOFFSET, YOFFSET
	//"
	//"               DIST               Perpendicular distance of the phase-space
	//"                                  plane to the point of rotation in cm.
	//"               ANGLE              Angle of rotation in degrees. The rotation
	//"                                  is performed around an axis that is parallel
	//"                                  to the x-axis and passes through the point
	//"                                  (x,y,z)=(0,0,ZOFFSET).
	//"               ZOFFSET            Point of rotation. If |ZOFFSET| > 1e4,
	//"                                  the centre of the geometry is taken as
	//"                                  the point of rotation (but note that
	///"                                  the maximum value allowed by the input
	//"                                  routine is 1e6, so that |ZOFFSET| must
	//"                                  be between 1e4 and 1e6 to use the centre
	//"                                  of the geometry automatically).
	//"               XOFFSET,YOFFSET    X and Y offset of scoring plane in BEAM
	//"                                  simulation (cm).  Offsets are applied before
	//"                                  rotating the source.
	//"
	//"    Examples:
	//"       - to have BEAM simulation incident on the upper z-face of the geometry,
	//"         use DIST=0, ANGLE=0, ZOFFSET=zplane(1)
	//"         This is the same as source 21
	//"       - to have BEAM simulation incident on the lower z-face of the geometry,
	//"         use DIST=0, ANGLE=180, ZOFFSET=zplane(n)
	//"       - to have BEAM simulation incident from, say, 60 degrees with
	//"         a distance to the centre of the geometry of 30 cm, use
	//"         DIST=30, ANGLE=60, ZOFFSET=9999.
	//"
	//"*******************************************************************************

		/**
		 * Handles radiation source in RZ geometry. 
		 */
		public static void SRCRZ()
		{
			ENFLAG = 0;

			if(ipart<0 || ipart>4){ipart=0;}//default
			//NO 21,22,23 since there are related to additional files!! ->No example and we will include these in the future (maybe)!
			if(ISOURC<0 || (ISOURC>=5 && ISOURC<=9)
			  ||(ISOURC>=17 && ISOURC<=19) || ISOURC>20)//ISOURC>23)
			{
				ISOURC=0;//default
			}

			if (ISOURC!=20)
			{
			   if(ISOURC==22)
			   {
				     nsource_option=9;
			   }
			   else if( ISOURC==23 )
			   {
				     nsource_option=5;
			   }
			   else
			   {
				     nsource_option=4;
			   }

			    for(int i=1;i<=nsource_option;i++)
			    {
				   if(source_option[i-1]<source_option_min || source_option[i-1]>source_option_max)
			   		source_option[i-1]=source_option_default;
			    }

				TEMP1=source_option[0];TEMP2=source_option[1];TEMP3=source_option[2];
				TEMP4=source_option[3];
				if(ISOURC==22)
				{
					TEMP5=source_option[4];TEMP6=source_option[5];TEMP7=source_option[6];
					TEMP8=source_option[7];TEMP9=source_option[8];
				}
				if( ISOURC==23 ) TEMP5=source_option[4];
		    }

			iqin=ipart;
			if (iqin==4) {iqin=-1;}
			//EGS4.seqStr=" IQIN: "+iqin;
			//if(EGS4.iprint>1)
			//	eq.printSequence(EGS4.seqStr);

			//"--------------------------------------------------------------------------"
			//"|                           optional inputs                              |"
			//"--------------------------------------------------------------------------"
			if (ISOURC==20)
			{
				if(MODEIN<0 || MODEIN>0){MODEIN=0;}//default no external file!!!GUI

				if (MODEIN==1)//not allowed
				{
				      //VALUES_SOUGHT(IVAL)='RDIST FILENAME';
			    }
			    else
			    {
					if(NRDIST<NRDIST_MIN || NRDIST>NRDIST_MAX)
						NRDIST=NRDIST_DEFAULT;
					for(int i=1;i<=NRDIST;i++)
					{
						if(RDISTF[i-1]<RDISTF_MIN || RDISTF[i-1]>RDISTF_MAX)
							RDISTF[i-1]=RDISTF_DEFAULT;
						if(RPDF[i-1]<RPDF_MIN || RPDF[i-1]>RPDF_MAX)
							RPDF[i-1]=RPDF_DEFAULT;
					}
				}

				if(RDIST_IOUTSP<0 || RDIST_IOUTSP>1){RDIST_IOUTSP=0;}//default

			}

			//"Check IQIN is OK and set 0 if not acceptable value"
			if((iqin < -1) || (iqin > 1 && (ISOURC < 21 || ISOURC > 23))
			    || (iqin > 3 && (ISOURC < 21 | ISOURC > 23))) iqin = 0;

			//"these three lines to effect source correlation"
			//"SVTMP1=TEMP1;SVTMP2=TEMP2;SVTMP2=TEMP2;SVTMP2=TEMP2;"
			//"after jans"
			SVTMP1=TEMP1;SVTMP2=TEMP2;SVTMP3=TEMP3;SVTMP4=TEMP4;
		}

		/**
		 * Initialize and print information for some pre-defined sources.
		 */
		public static void SRCINI()
		{
			if (EGS4.iGeom==EGS4.iCavity)//RZ GEOM
			{

			NHSTRY=0; //"initialize no. of primary histories"
			          //"will have to change this if we start storing this in .egsdat"
			          //"files"
			last_nhstry = 0;

			if(ISOURC == 0)
			{// "FRONTAL PARALLEL BEAM"
			    RBEAM=TEMP1;
			    FNORM=TEMP2*TEMP2+TEMP3*TEMP3+TEMP4*TEMP4;
			    if(FNORM==0.0){ UINC=0.0;VINC=0.0;WINC=1.0; }
			    else
			    {
			        FNORM=Math.sqrt(FNORM);
			        UINC=TEMP2/FNORM;VINC=TEMP3/FNORM;WINC=TEMP4/FNORM;
				}
			    TEMP5=EGS4Geom.RCYL[EGS4Geom.NR];//0 biased
			    if((RBEAM<=0.0)||(RBEAM>TEMP5))RBEAM=TEMP5;

				EGS4.seqStr=" Electric charge of the source:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Parallel beam incident from the front(+ve Z-axis)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Radius of beam at the front face of the target:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" X-axis direction cosine:"+
				EGS4.format(UINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Y-axis direction cosine:"+
				EGS4.format(VINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Z-axis direction cosine:"+
				EGS4.format(WINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			}

			else if(ISOURC==1)
			{// "FRONTAL POINT SOURCE"
			    DISTZ=TEMP1;RBEAM=TEMP2;TEMP3=EGS4Geom.RCYL[EGS4Geom.NR];//0 biased
			    if(DISTZ<=0.0)DISTZ=100.;
			    if((RBEAM<=0.0)||(RBEAM>TEMP3))RBEAM=TEMP3;

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" POINT SOURCE ON AXIS INCIDENT FROM THE FRONT(+VE Z-AXIS)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" SOURCE DISTANCE TO THE FRONT FACE OF THE TARGET:"+
				EGS4.format(DISTZ,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE FRONT FACE OF THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			}

			else if(ISOURC==2)
			{// "BROAD FRONTAL PARALLEL BEAM"
			    RBEAM=0.0;
			    EGS4Geom.NR=1;EGS4Geom.nreg=EGS4Geom.NPLANE;
			    EGS4Geom.RCYL[1]=1000.0;EGS4Geom.CYRAD2[0]=1.0E6;//CYRAD2(1)=1.0E6;

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" BROAD PARALLEL BEAM INCIDENT FROM THE FRONT (+VE Z-AXIS)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE FRONT FACE OF THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			    if(EGS4.NMED!=1)
			    {
			        //"EACH SLAB IN THE PHANTOM MUST BE HOMOGENEOUS"
					EGS4.seqStr=" ALL REGIONS IN A SLAB MUST BE SAME MATERIAL!!";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			}

			else if(ISOURC == 3)
			{// "ISOTROPICALLY RADIATING DISK"
			    RBEAM=TEMP2;
			    if(RBEAM > EGS4Geom.RCYL[EGS4Geom.NR])
			    {
					EGS4.seqStr=" Source radius too big at "+EGS4.format(RBEAM,10,true)+
					" ,Reduced to "+EGS4.format(EGS4Geom.RCYL[EGS4Geom.NR],10,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

			       RBEAM = EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
				}
			    RBEAM2 = RBEAM*RBEAM;
			    RMINBM = TEMP1; if(RMINBM > RBEAM){RMINBM=RBEAM;}
			    RMINSQ = RMINBM*RMINBM;
			    ZSMIN=TEMP3;ZSMAX=TEMP4;
			    if(ZSMIN < EGS4Geom.ZPLANE[0])//ZPLANE(1))
			    { ZSMIN= EGS4Geom.ZPLANE[0];}//ZPLANE(1);}
			    if(ZSMAX > EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])
			    { ZSMAX = EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1];}

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" UNIFORM ISOTROPICALLY RADIATING RING: INNER,OUTER Radius"+
				EGS4.format(RMINBM,10,true)+","+EGS4.format(RBEAM,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=EGS4.format("",50)+"FRONT,Back DEPTH:"+
				EGS4.format(ZSMIN,10,true)+","+EGS4.format(ZSMAX,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}

			else if(ISOURC == 4)
			{// "POINT SOURCE WITH VARIOUS BEAM RADII"
			    RBEAM=TEMP1;
			    if(EGS4Geom.RCYL[EGS4Geom.NR]<650.)
			    {
			        //"MUST HAVE A BROAD PHANTOM changed limit to 650 cm April 1990 DR"
			        if(EGS4Geom.NR<EGS4Geom.$MAXRADII)
			        {
						EGS4Geom.NR=EGS4Geom.NR+1;
						EGS4Geom.nreg=EGS4Geom.NZ*EGS4Geom.NR+1;
					}
			        else
			        {
						EGS4.seqStr=" WARNING: LAST RADIAL BIN INCREASED TO 1000cm";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
					}
			        EGS4Geom.RCYL[EGS4Geom.NR]=1000.0;
			        EGS4Geom.CYRAD2[EGS4Geom.NR-1]=1.0E6;
				}
			    if((RBEAM<=0.0)||(RBEAM>EGS4Geom.RCYL[EGS4Geom.NR]))
			    	RBEAM=EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" CENTRAL AXIS FLUENCE VS BEAM RADIUS";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADII INPUT ABOVE WILL BE CONSIDERED AS BEAM RADII";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF CENTRAL AXIS ZONE: "+EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			    if(EGS4.NMED!=1)
			    {
			        //"EACH SLAB IN THE PHANTOM MUST BE HOMOGENEOUS"
					EGS4.seqStr=" ALL REGIONS IN A SLAB MUST BE SAME MATERIAL!!";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			}

			else if(ISOURC == 10)
			{// "SIDE PARALLEL BEAM"
			    XBEAM=TEMP1;ZBEAM=TEMP2;
			    TEMP3=EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
			    if((XBEAM<=0.0)||(XBEAM>TEMP3))XBEAM=TEMP3;
			    //TEMP3=0.5*(ZPLANE(NPLANE)-ZPLANE(1));
			    TEMP3=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
			    if((ZBEAM<=0.0)||(ZBEAM>TEMP3))ZBEAM=TEMP3;

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" PARALLEL BEAM INCIDENT FROM THE SIDE";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" HALF-WIDTH(X-AXIS) OF THE BEAM AT THE TARGET MIDPOINT"+
				EGS4.format(XBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" HALF-HEIGHT(Z-AXIS/) OF THE BEAM AT THE TARGET MIDPOINT"+
				EGS4.format(ZBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}

			else if(ISOURC == 11)
			{// "Side point source centered at target middle"
			    DISTRH = TEMP1;
			    if(DISTRH < 0.0) DISTRH=-DISTRH;
			    if(DISTRH == 0.0) DISTRH=100.0;
			    XBEAM = TEMP2;   ZBEAM = TEMP3;
			    TEMP4 = EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
			    if(DISTRH <= TEMP4)
			    {
					EGS4.STOPPROGRAM=true;
					EGS4.seqStr="ERROR: SOURCE FROM WITHIN TARGET, INPUT IGNORED!!";
					eq.printSequence(EGS4.seqStr);
					return;
				}
			    //"TEMP4 is the maximum 1/2 width of the beam by simple geometry"
			    TEMP4=TEMP4*Math.sqrt((DISTRH+TEMP4)*(DISTRH-TEMP4))/DISTRH; //"SHADOWING"
			    if((XBEAM <= 0.0) || (XBEAM > TEMP4)) XBEAM=TEMP4;
			    //TEMP4=0.5*(ZPLANE(NPLANE)-ZPLANE(1));
			    TEMP4=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
			    if((ZBEAM <= 0.0) || (ZBEAM > TEMP4)) ZBEAM=TEMP4;
			    //IF(ICORRL = 0)[

				EGS4.seqStr=" Electric charge of the source:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Point source on axis incident from the side(Y-axis)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Source distance to target midpoint:"+
				EGS4.format(DISTRH,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Half-width(X-axis) of the beam at the target midpoint"+
				EGS4.format(XBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Half-height(Z-axis/) of the beam at the target midpoint"+
				EGS4.format(ZBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			   //]
			}

			else if(ISOURC == 12)
			{// "point source off axis"
			    DISTRH=TEMP1;       //"DISTRH is distance of source from center of geomety"
			    if(DISTRH < 0.0)DISTRH=-DISTRH;
			    DISTZ=TEMP2;        //"DISTZ is distance of source from front of geometry"
			                        //"DISTZ > 0     point located in front of front face
			                        //"0 > DISTZ > -(ZPLANE(NPLANE)-ZPLANE(1))
			                        //"              point located between front, rear faces
			                        //"DISTZ < -(ZPLANE(NPLANE)-ZPLANE(1))
			                        //"              point located rear of rear plane
			    RBEAM=EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
			    if(DISTRH > RBEAM)
			    {
					//"the source is further away than the outer radius"
			        //"XBEAM    half-height of the beam at the center of the target"
			        //"         following is simple geometry-its really 1/2 width if you"
			        //"         think of Z as the height direction"
			        XBEAM=RBEAM*Math.sqrt((DISTRH-RBEAM)*(DISTRH+RBEAM))/DISTRH;
			        //"ZBEAM is half-width (height) of the beam at the center of the target"
			        //ZBEAM=0.5*(ZPLANE(NPLANE)-ZPLANE(1));
			        ZBEAM=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
				}
			    else
			    {
					// "source either inside geometry (error caught below)"
			        // "or it just hits the top or bottom of the geometry, not the side"
			        XBEAM=0.0; ZBEAM=0.0;
				}
			    if((DISTRH <= RBEAM) &&
			    	//DISTZ >= ZPLANE(1)-ZPLANE(NPLANE) &&
			    	DISTZ >= EGS4Geom.ZPLANE[0]-EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1] &&
			    	DISTZ<=0.0)
			    {
			        //"first condition => inside cylinder outer boundary"
			        //"2nd condition => it is above lower plane"
			        //"3rd condition => it is below front face"
					EGS4.STOPPROGRAM=true;
					EGS4.seqStr=" Source from within target, input ignored!!";
					eq.printSequence(EGS4.seqStr);
					return;

			        //OUTPUT;(' Source from within target, input ignored');
			        //ERROR_FLAG=1;
				}

				EGS4.seqStr=" Electric charge of the source:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Point source off axis";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Distance of source off the central axis"+
				EGS4.format(DISTRH,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Perpendicular distance of source from front plane"+
				EGS4.format(DISTZ,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Radius of the beam on the front face of the target:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Half-width(X-axis) of the beam between the tangent points"+
				EGS4.format(XBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Half-height(z-axis) of the beam between the tangent points"+
				EGS4.format(ZBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			}

			else if(ISOURC == 13)
			{// "Frontal parallel beam at any angle wrt the target"
			    FNORM=TEMP1*TEMP1+TEMP2*TEMP2+TEMP3*TEMP3;
			    if(FNORM == 0.0){ UINC=0.0;VINC=0.0;WINC=1.0; }
			    else
			    {
			        //"AT THIS POINT THE INCIDENT BEAM IS ROTATED SO THAT THE"
			        //"DIRECTION IS IN THE Y-Z PLANE AND POINTED IN THE      "
			        //"NEGATIVE Y DIRECTION. THE NEGATIVE DIRECTION IS CHOSEN"
			        //"BECAUSE THE ENTRANCE POINT IS CHOSEN TO BE ON THE     "
			        //"POSITIVE Y-SIDE OF THE CYLINDER. (SEE MACRO $CHOOSE-  "
			        //"POINT-IN-RECTANGLE.) THIS ROTATION IS PERMITTED       "
			        //"OF THE AZIMUTHAL SYMMETRY OF THE TARGET.THE REST OF   "
			        //"CODING IS LEFT GENERAL WHERE POSSIBLE TO ALLOW THE    "
			        //"SYMMETRY TO BE RELAXED.                               "
			        FNORM=Math.sqrt(FNORM);
			        WINC=TEMP3/FNORM;
			        VINC=-Math.sqrt((1.0-WINC)*(1.0+WINC));
			        UINC=0.0;
				}
			    RBEAM=EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
			    if(Math.abs(WINC)!=1.0)
			    {
			        XBEAM=RBEAM;
			        //ZBEAM=0.5*(ZPLANE(NPLANE)-ZPLANE(1));
			        ZBEAM=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
				}
			    else{ XBEAM=0.0;ZBEAM=0.0; }

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" PARALLEL BEAM AT ANY ANGLE WITH RESPECT TO THE TARGET";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" X-AXIS DIRECTION COSINE:"+
				EGS4.format(UINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Y-AXIS DIRECTION COSINE:"+
				EGS4.format(VINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Z-AXIS DIRECTION COSINE:"+
				EGS4.format(WINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE FRONT FACE OF THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" HALF-WIDTH(X-AXIS) OF THE BEAM BETWEEN THE TANGENT POINTS"+
				EGS4.format(XBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" HALF-HEIGHT(Z-AXIS) OF THE BEAM BETWEEN THE TANGENT POINTS"+
				EGS4.format(ZBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}

			else if(ISOURC == 14)
			{
			    //"FRONTAL POINT SOURCE with all particles below some radius skipped"
			    DISTZ=TEMP1;RBEAM=TEMP2;RMINBM=TEMP3;
			    if(DISTZ<=0.0)DISTZ=100.;
			    if((RBEAM<=0.0)||(RBEAM>EGS4Geom.RCYL[EGS4Geom.NR]))
			    	RBEAM=EGS4Geom.RCYL[EGS4Geom.NR];//RCYL(NR);
			    if(RMINBM==0.0)
			    {
					EGS4.seqStr=" WHY USE ISOURC=14 WITH RMINB=0.0?";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					//OUTPUT;(/' ***WHY USE ISOURC=14 WITH RMINB=0.0?***');
				}

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" POINT SOURCE ON AXIS INCIDENT FROM THE FRONT(+VE Z-AXIS)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" SOURCE DISTANCE TO THE FRONT FACE OF THE TARGET:"+
				EGS4.format(DISTZ,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE FRONT FACE OF THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" INNER RADIUS OF BEAM ON FRONT FACE OF THE TARGET:"+
				EGS4.format(RMINBM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}

			else if(ISOURC == 15)
			{
			    //"Point source from an arbitrary distance and direction."
			    //"This uses an alternative coding to source 12 to avoid "
			    //"strongly varying weights in some cases "
			    src15_ini();//(temp1,temp2,temp3,temp4,ERROR_FLAG);
			}

			else if(ISOURC == 16)
			{
			    //"A disk or a rectangle irradiating the geometry "
			    //"from an arbitrary distance or angle            "
			    src16_ini();//(temp1,temp2,temp3,temp4,ERROR_FLAG);
			}

			else if(ISOURC == 20)
			{//"Radial distribution"
			    //OUTPUT;(' LOCAL INPUT(0) OR EXTERNAL FILE(1): ',$);
			    //"IF(MODEIN.NE.1) MODEIN=0;"  "DEFAULT"
				EGS4.seqStr=" RADIAL DISTRIBUTION:";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			    if(MODEIN==0)
			    {
			        //"INPUT FROM KEYBOARD OR .INP FILE"
			        //OUTPUT;(' NUMBER OF RADIAL BINS: ',$);
			        if((NRDIST<1)||(NRDIST>100))
			        {
						EGS4.seqStr="NUMBER RADIAL BINS OUT OF RANGE (<1 OR >100)=> RESET TO 100";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);

			            NRDIST=100;
					}
					EGS4.seqStr="NRDIST: "+NRDIST;
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
			        //OUTPUT NRDIST;(' INPUT',I4,' SETS OF RDISTF,RPDF IN 2F20.0 FORMAT');
			        //OUTPUT;('   RDISTF INCREASING IN SIZE, RPDF NON-NEGATIVE');
			        for(int IB=1;IB<=NRDIST;IB++)
			        {
			        	//"What is this?  DR"
					}
				}
			    //ELSE[
			    //    "EXTERNAL FILE INPUT"
			    //    OUTPUT;
			    //    (' INPUT NAME OF FILE WITH SPECTRUM (up to 80 CHAR WITH .EXT): '/);
			    //    READ(5,600,END=:EOF_SRCRZ:)FILNAM;
			    //    OUTPUT FILNAM;(/' READ INPUT RADIAL DISTRIBUTION FROM: '/80A1);
			    //    OPEN(UNIT=9,file=filnam,STATUS='OLD');
			    //    READ(9,*)NRDIST;
			    //    IF((NRDIST.LT.1).OR.(NRDIST.GT.100))[
			    //        OUTPUT;
			    //        (' *** NUMBER RADIAL BINS OUT OF RANGE (<1 OR >100),',
			    //        ' RESET TO 100 ***');
			    //        NRDIST=100;
			    //    ]
			    //    READ(9,*)(RDISTF(IB),RPDF(IB),IB=1,NRDIST);
			    //    CLOSE(UNIT=9);
			    //    OUTPUT NRDIST;('    HAVE READ',I5,' INPUT RADIAL BINS FROM FILE');
			    //]

			    //"DO A CHECK ON THE RADIAL DISTRIBUTION"
			    int ICOUNT=0;
			    double RLAST=0.0;
			    int IERROR=0;
			    R_DIST_INPUT:
			    while(true)
			    {
			        ICOUNT=ICOUNT+1;
			        if(ICOUNT>NRDIST){break R_DIST_INPUT;}
			        if(RDISTF[ICOUNT-1]<=RLAST)
			        {
			            IERROR=1;

						EGS4.seqStr="RDISTF>=LAST ONE. NOT ALLOWED=>TERMINATING RADIAL DISTRIBUTION INPUT";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
					}
			        else if(RDISTF[ICOUNT-1]>EGS4Geom.RCYL[EGS4Geom.NR])
			        {
			            IERROR=1;

						EGS4.seqStr="RDISTF>"+EGS4.format(EGS4Geom.RCYL[EGS4Geom.NR],14,true)+
						" ,GEOMETRY SIZE.=>TERMINATING RADIAL DISTRIBUTION INPUT";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
					}
			        else if(RPDF[ICOUNT-1]<0.0)
			        {
			            IERROR=1;

						EGS4.seqStr=" PDF < 0 NOT ALLOWED,TERMINATING RADIAL DISTRIBUTION INPUT";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
					}

			        if(IERROR==1)
			        {
			            ICOUNT=ICOUNT-1;
			            if(ICOUNT==0)
			            {
							EGS4.STOPPROGRAM=true;
							EGS4.seqStr=" NO RADIAL DISTRIBUTION DEFINED=>STOPPING EXECUTION!";
							eq.printSequence(EGS4.seqStr);
							return;
						}

						EGS4.seqStr=" RADIAL DITRIBUTION INPUT APPEARS TO BE INCOMPLETE, NRDIST RESET TO "+
						EGS4.format(ICOUNT,12);
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);

			            NRDIST=ICOUNT;
			            break R_DIST_INPUT;
					}

			        RLAST=RDISTF[ICOUNT-1];
				}

				EGS4.seqStr="    RADIAL DISTRIBUTION RANGES FROM 0 TO"+
				EGS4.format(RDISTF[NRDIST-1],12,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			    //OUTPUT;(' PRINT DISTRIBUTION DATA IN OUTPUT SUMMARY, YES(1) OF NO(0): ',$);
			    if(RDIST_IOUTSP!=1) RDIST_IOUTSP=0;//RDIST_IOUTSP
			    //WRITE(6,'('' '')')
			}
	/*
	//NOT USED HERE
			ELSEIF( ISOURC = 23 ) [ "A full treatment head simulation source using BEAM"
			    dist_phsp = temp1; theta_phsp = temp2;
			    cost_phsp = cos(theta_phsp*0.017453292222);
			    sint_phsp = sin(theta_phsp*0.017453292222);
			    chamber_c = temp3; xoffset = temp4; yoffset = temp5;
			    ENFLAG=1;
			    iqinc = iqin;
			    OUTPUT iqinc,min_weight_23,max_weight_23,
			           $cstring(the_beam_code),$cstring(the_pegs_file),
			           $cstring(the_input_file),
			           dist_phsp,theta_phsp,chamber_c,xoffset,yoffset;
			    (/' Full treatment head simulation from an angle'/,
			      '    Particle charge                               : ',i3,/
			      '    Weight window                                 : ',2g15.6,/
			      '    BEAMnrc user code                             : ',a,/
			      '    PEGS data file                                : ',a,/
			      '    Input file                                    : ',a,/
			      '    Rotation point to phsp-plane midpoint distance: ',f10.4,/
			      '    Incident angle (degrees)                      : ',f10.4,/
			      '    Rotation point (will be set to chamber center'/,
			      '      if greater than 1e4 or less than -1e4):     : ',g14.4,
			      ' cm.'/
			      '    X offset of phsp plane (before rotation)      : ',f10.4,' cm'/
			      '    Y offset of phsp plane (before rotation)      : ',f10.4,' cm');
			     write(6,'(//a)') 'About to call init_beamsource';
			     call init_beamsource(i_parallel,$CONFIGURATION_NAME,
			                         hen_house,egs_home,the_beam_code,
			                         the_pegs_file,the_input_file);
			     call maxenergy_beamsource(EKSRCM);
			     EIN=EKSRCM;
			     OUTPUT EKSRCM;
			     ('    Max. kinetic energy of simulation             : ',f10.4,' MeV'/);
			]
	//NOT USED HERE
			ELSEIF((ISOURC = 21) | (ISOURC = 22))["full phase space for each particle"
			    IQINC=IQIN;"IQINC is used in SRCHST"
			    IF(ISOURC = 22) [
			        dist_phsp=TEMP2; theta_phsp = temp3;
			        cost_phsp = cos(theta_phsp*0.017453292222); "convert to radians"
			        sint_phsp = sin(theta_phsp*0.017453292222);
			        chamber_c = temp4;
			        count_phsp = 0;
			        xoffset = TEMP8;
			        yoffset = TEMP9;
			    ]
			    IMODE=TEMP1;"MODE OF THE PH-SP FILE"
			    IF(ISOURC=21)[
			       NRCYCL=TEMP2;
			       IPARALLEL=TEMP3;
			       PARNUM=TEMP4;
			    ]
			    ELSEIF(ISOURC=22)[
			       NRCYCL=TEMP5;
			       IPARALLEL=TEMP6;
			       PARNUM=TEMP7;
			    ]
			    IF(IMODE.NE.2)IMODE=0;"DEFAULT TO MODE0"
			    OUTPUT IMODE;
			    (/' MODE',I2,' Phase-space file to be read from unit 42...'/
			      ' Input name of file with phase space data (1 - 80 CHARS,',
			      'with .EXT): ');
			    OUTPUT FILSPC;(/' Reading phase space information from: '/10x,80A1);
			    $OPEN_PHSP_FOR_READ(IMODE,42,SPCNAM,LINE1,NCASE_PHSP,
			                        TEMP1,TEMP2,TEMP3,NINCSRC);

			       "NCASE_PHSP: the total # of particles in the ph-sp file"
			       "the 3rd variable is the number of photons"
			       "the 4th variable is the maximum kinetic energy of all particles"
			       "the 5th variable is the minimum kinetic energy of electrons"
			       "the 6th variable is the # of incident partices from original source"

			     OUTPUT NCASE_PHSP,TEMP1,TEMP2,TEMP3,NINCSRC;
			     (' Total number of particles in file      :',I10/
			      ' Total number of photons                :',I10/
			      ' (the rest are electrons and positrons)'//
			      ' Maximum kinetic energy of the particles:',F10.3,' MeV'/
			      ' Minimum kinetic energy of the electrons:',F10.3,' MeV'/
			      ' # particles incident when phase space created :',F15.0/);
			      EKSRCM=TEMP2;EIN=TEMP2;"EIN WILL BE USED BY PRESTA LATER"

			     IF (NCASE_PHSP < 0 | TEMP1 < 0 | TEMP2 <= 0 | TEMP3 <0 |
			         NINCSRC <=0.0) ["This means something wrong "
			            OUTPUT;( 3(1x,79('*')/)/' Something is wrong in the above'
			            /'Likely the file is the wrong format (need to swap bytes)'/
			            3(1x,79('*')/));
			            STOP;
			     ]

			     "read record 2 of phsp source, if NHSTRYM gets incremented, then we"
			     "can read NHSTRY from this source."
			     NHSTRYM=0;
			     $READ_PHSP_RECORD(IMODE,42,2:
			           NHSTRYM,ZLASTM,LATCHM,EINM,WEIGHTM,XINM,YINM,UINM,VINM);
			     IF(NHSTRYM~=1)[
			       DOSE_STAT=1;
			       OUTPUT;(//' ***WARNING***'/
			       ' Cannot read no. of primary (non-phsp) histories from ph-sp source.'/
			       ' Dose and fluence will be analyzed assuming each particle read from'/
			       ' the ph-sp file is an independent history.  May result in an'/
			       ' underestimate of uncertainties.'//);
			     ]

			     IF( isourc = 22 ) [
			         OUTPUT dist_phsp,theta_phsp,chamber_c,xoffset,yoffset;
			         (/' Phase space file from an angle'/,
			           '    Rotation point to phsp-plane midpoint distance: ',f10.4,/
			           '    Incident angle (degrees)                      : ',f10.4,/
			           '    Rotation point (will be set to chamber center'/,
			           '      if greater than 1e4 or less than -1e4):     : ',g14.4,
			           ' cm.'/
			           '    X offset of phsp plane (before rotation)      : ',f10.4,' cm'/
			           '    Y offset of phsp plane (before rotation)      : ',f10.4,' cm'/);
			     ]

			    IF(IPARALLEL>1 & n_parallel>0)[
			       OUTPUT IPARALLEL, n_parallel;
			       (/' You have set IPARALLEL in the input file to ',I4,/
			         ' But you are also running a C compiled code with n_parallel=',I4,/
			         ' IPARALLEL will be reset to 1, and control of the parallel run '/
			         ' will be from the code.'/);
			       IPARALLEL=1;
			    ]

			    IF(IPARALLEL<=0) IPARALLEL=1;
			    IF(PARNUM<0) PARNUM=0;
			    IF(IPARALLEL>1)[
			       IF(PARNUM>=1 & PARNUM<=IPARALLEL)[
			         OUTPUT IPARALLEL,INT((PARNUM-1)*NCASE_PHSP/IPARALLEL)+1,
			                INT(PARNUM*NCASE_PHSP/IPARALLEL);
			                  (/' This is one of ',I4,' parallel jobs.'/
			                    ' It will use from particle ',I12,' to particle ',I12,/
			                    ' from the phase space source in the simulation.'/);
			       ]
			       ELSE[
			         OUTPUT IPARALLEL;
			     (/' IPARALLEL input indicates that this is one of ',I4,' parallel jobs.'/
			       ' But PARNUM is out of range (<1 or >IPARALLEL).  Therefore, phsp '/
			       ' source will not be partitioned.'/);
			         PARNUM=0;
			         IPARALLEL=1;
			       ]
			    ]
			    IF(NRCYCL<=0)[
			     IF(IQIN=1)["cannot estimate NRCYCL"
			       OUTPUT;(/' NRCYCL cannot be calculated automatically because'/
			                ' INCIDENT PARTICLE= positrons '/);
			       NRCYCL=0;
			     ]
			     ELSE[
			       OUTPUT;(/' NRCYCL will be calculated automatically'/);
			       IF(IQIN=-1 | IQIN=3)[
			         TEMPDIV=NCASE_PHSP-TEMP1;
			       ]
			       ELSEIF(IQIN=0)[
			         TEMPDIV=TEMP1;
			       ]
			       ELSEIF(IQIN=2)[
			         TEMPDIV=NCASE_PHSP;
			       ]
			       IF(NINT(dble(IPARALLEL*NCASE)/dble(TEMPDIV))<=1)[
			         "IPARALLEL*NCASE slightly > TEMPDIV or < TEMPDIV"
			         NRCYCL=0;
			       ]
			       ELSEIF(MOD(IPARALLEL*NCASE,TEMPDIV)=0)[
			         "IPARALLEL*NCASE is an exact multiple of TEMPDIV"
			         NRCYCL=(IPARALLEL*NCASE)/TEMPDIV-1;
			       ]
			       ELSE[
			         NRCYCL=(IPARALLEL*NCASE)/TEMPDIV;
			       ]
			     ]
			    ]
			    OUTPUT NRCYCL;
			  (/' Particles will be recycled ',I4,' times before moving on to next one.'/);

			    CYCLNUM=0; "counts number of times a particle has been recycled"
			    ENFLAG=1; "flag to tell ensrc that no questions needed re energy"
			    N_RUN_CHUNK_OLD=0; "initialize for parallel runs"
			]
			ELSE[ "defaults to parallel beam from the front"
			    ISOURC=0;RBEAM=RCYL(NR);
			    OUTPUT IQIN,RBEAM;
			    (/
			    ' Electric charge of the source:',T60,I12/
			    ' Parallel beam incident from the front(+ve Z-axis)'/
			    ' Radius of beam at the front face of the target:',T60,F10.4,' cm'/);
			]

			RETURN;  "normal return"

			:EOF_SRCRZ:;  "bad input"
			ERROR_FLAG=1;RETURN;

			"all data has now been fed into the routine"

			"************************************************************************
			*/
		    }
		    else if(EGS4.iGeom==EGS4.iCavitySPH)//2=SPH Geom!!!
		    {
			//" The arguments to this call aren't used here.  They were probably used
			//" for a source type not yet defined for sph, so I'll leave it in -- JT
			if((ISOURC!=0)&&(ISOURC!=1)&&(ISOURC!=4)&&(ISOURC!=10)&&(ISOURC!=11))
			{
			    //"DEFAULT SOURCE TYPE- PARALLEL BEAM FROM THE FRONT"
			    ISOURC=0;TEMP1=0.0;TEMP2=0.0;TEMP3=0.0;TEMP4=1.0;
			}

			if(ISOURC==0)
			{
			    //"PARALLEL BEAM AT ANY ANGLE"
			    RBEAM=TEMP1;
			    if((RBEAM<=0.0)||(RBEAM>EGS4Geom.RSPH[EGS4Geom.NR]))
			    	RBEAM=EGS4Geom.RSPH[EGS4Geom.NR];//RSPH(NR);
			    FNORM=TEMP2*TEMP2+TEMP3*TEMP3+TEMP4*TEMP4;
			    if(FNORM==0.0)
			    {
			        UINC=0.0;VINC=0.0;WINC=1.0;
				}
			    else
			    {
			        //"TO SIMPLIFY, THE BEAM IS ROTATED SO THAT THERE IS NO         "
			        //"U-COMPONENT. THIS IS ALLOWED BY THE SYMMETRY OF THE GEOMETRY."
			        FNORM=Math.sqrt(FNORM);
			        WINC=TEMP4/FNORM;
			        VINC=Math.sqrt((1.0-WINC)*(1.0+WINC));
			        UINC=0.0;
				}

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" PARALLEL BEAM AT ANY ANGLE";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" X-AXIS DIRECTION COSINE:"+
				EGS4.format(UINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Y-AXIS DIRECTION COSINE:"+
				EGS4.format(VINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Z-AXIS DIRECTION COSINE:"+
				EGS4.format(WINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}

			else if(ISOURC==1)
			{
			    //"POINT SOURCE AT ANY ANGLE"
			    DISTR=Math.abs(TEMP1);
			    if(DISTR<=EGS4Geom.RSPH[EGS4Geom.NR])DISTR=100.;          //"DEFAULT"
			    if(DISTR<=EGS4Geom.RSPH[EGS4Geom.NR])DISTR=100.+
			    EGS4Geom.RSPH[EGS4Geom.NR]; //"DEFAULT FOR BIG TARGET"
			    RBEAM=TEMP2;
			    //"RSHADW=MAXIMUM BEAM RADIUS ALLOWED BY SHADOWING"
			    RSHADW=(EGS4Geom.RSPH[EGS4Geom.NR]/DISTR)*
			    Math.sqrt((DISTR+EGS4Geom.RSPH[EGS4Geom.NR])*
			    (DISTR-EGS4Geom.RSPH[EGS4Geom.NR]));
			    if((RBEAM<=0.0)||(RBEAM>RSHADW))RBEAM=RSHADW;
			    FNORM=TEMP3*TEMP3+TEMP4*TEMP4+TEMP5*TEMP5;
			    if(FNORM==0.0)
			    {
			        UINC=0.0;VINC=0.0;WINC=1.0;
				}
			    else
			    {
			        //"TO SIMPLIFY, THE BEAM IS ROTATED SO THAT THERE IS NO         "
			        //"U-COMPONENT. THIS IS ALLOWED BY THE SYMMETRY OF THE GEOMETRY."
			        FNORM=Math.sqrt(FNORM);
			        WINC=TEMP5/FNORM;
			        VINC=Math.sqrt((1.0-WINC)*(1.0+WINC));
			        UINC=0.0;
				}

				EGS4.seqStr=" ELECTRIC CHARGE OF THE SOURCE:"+EGS4.format(iqin,12);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" POINT SOURCE BEAM AT ANY ANGLE";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" SOURCE DISTANCE TO THE MIDDLE OF THE TARGET:"+
				EGS4.format(DISTR,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" RADIUS OF THE BEAM ON THE TARGET:"+
				EGS4.format(RBEAM,10,true)+" cm";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" X-AXIS DIRECTION COSINE:"+
				EGS4.format(UINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Y-AXIS DIRECTION COSINE:"+
				EGS4.format(VINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Z-AXIS DIRECTION COSINE:"+
				EGS4.format(WINC,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			}

			else if(ISOURC==4)
			{//"point spread function calculation"

			    //"DO NOTHING, I.E.:"
			    //"XIN = YIN = ZIN = 0.0"
			    //"UIN = VIN = 0.0; WIN = 1.0"

				EGS4.seqStr="PARTICLE LOCATED AT ORIGIN,MOVING PARALLEL TO Z-AXIS";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="USED FOR POINT SPREAD FUNCTION CALCULATIONS";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

			}
	/*
	//NO
			ELSEIF(ISOURC.EQ.10)[
			    "PARALLEL BEAM AT ANY ANGLE WITH RADIAL DISTRIBUTION"
			    FNORM=TEMP1**2+TEMP2**2+TEMP3**2;
			    IF(FNORM.EQ.0.0)[
			        UINC=0.0;VINC=0.0;WINC=1.0;
			        ]
			    ELSE[
			        "TO SIMPLIFY, THE BEAM IS ROTATED SO THAT THERE IS NO         "
			        "U-COMPONENT. THIS IS ALLOWED BY THE SYMMETRY OF THE GEOMETRY."
			        FNORM=SQRT(FNORM);
			        WINC=TEMP3/FNORM;
			        VINC=SQRT((1.0-WINC)*(1.0+WINC));
			        UINC=0.0;
			        ]
			    OUTPUT IQIN,UINC,VINC,WINC;
			    (' ELECTRIC CHARGE OF THE SOURCE:',T60,I12/
			    ' PARALLEL BEAM AT ANY ANGLE WITH RADIAL DISTRIBUTION'/
			    ' X-AXIS DIRECTION COSINE:',T60,F10.4/
			    ' Y-AXIS DIRECTION COSINE:',T60,F10.4/
			    ' Z-AXIS DIRECTION COSINE:',T60,F10.4/);
			    ]

			ELSEIF(ISOURC.EQ.11)[
			    "POINT SOURCE AT ANY ANGLE WITH RADIAL DISTRIBUTION"
			    DISTR=ABS(TEMP1);
			    IF(DISTR.LE.RSPH(NR))DISTR=100.;          "DEFAULT"
			    IF(DISTR.LE.RSPH(NR))DISTR=100.+RSPH(NR); "DEFAULT FOR BIG TARGET"
			    FNORM=TEMP2**2+TEMP3**2+TEMP4**2;
			    IF(FNORM.EQ.0.0)[
			        UINC=0.0;VINC=0.0;WINC=1.0;
			        ]
			    ELSE[
			        "TO SIMPLIFY, THE BEAM IS ROTATED SO THAT THERE IS NO         "
			        "U-COMPONENT. THIS IS ALLOWED BY THE SYMMETRY OF THE GEOMETRY."
			        FNORM=SQRT(FNORM);
			        WINC=TEMP4/FNORM;
			        VINC=SQRT((1.0-WINC)*(1.0+WINC));
			        UINC=0.0;
			        ]
			    OUTPUT IQIN,DISTR,UINC,VINC,WINC;
			    (' ELECTRIC CHARGE OF THE SOURCE:',T60,I12/
			    ' POINT SOURCE AT ANY ANGLE WITH RADIAL DISTRIBUTION'/
			    ' SOURCE DISTANCE TO THE MIDDLE OF THE TARGET:',T60,F10.4,' cm'/
			    ' X-AXIS DIRECTION COSINE:',T60,F10.4/
			    ' Y-AXIS DIRECTION COSINE:',T60,F10.4/
			    ' Z-AXIS DIRECTION COSINE:',T60,F10.4/);
			    ]

			OUTPUT IOUTSP;
			(' PRINT DISTRIBUTION DATA IN OUTPUT SUMMARY, YES(1) OR NO(0): ',I2);

			IF((ISOURC.EQ.10).OR.(ISOURC.EQ.11))[ "INCLUDE RADIAL DISTRIBUTION"

			    IF(MODEIN.NE.1) MODEIN=0; "DEFAULT"

			    IF(MODEIN.EQ.0)[ "HAS READ RADIAL DIST LOCALLY (FROM .INP)
			        OUTPUT NRDIST;(' NUMBER OF RADIAL BINS: ',I4);
			        IF((NRDIST.LT.1).OR.(NRDIST.GT.100))[
			            OUTPUT;
			            (' *** NUMBER RADIAL BINS OUT OF RANGE (<1 OR >100),',
			            ' RESET TO 100 ***');
			            NRDIST=100;
			            ]
			        OUTPUT NRDIST;('   INPUT',I4,' SETS OF RDISTF,RPDF IN 2F20.0 FORMAT');
			        OUTPUT;('   RDISTF INCREASING IN SIZE, RPDF NON-NEGATIVE');
			        DO I=1,NRDIST [
			           RDISTF(I)=VALUE(NUM_RDISTF,I);
			           RPDF(I)=VALUE(NUM_RPDF,I);
			           ]
			        ]
			    ELSE[ "READ INPUT FROM A SEPARATE SOURCE FILE"
			        OUTPUT SPCFIL;(' FILE WITH DISTRIBUTION DATA READ AS:',A80);
			        OPEN(UNIT=9,STATUS='OLD',file=SPCFIL);
			        READ(9,*)NRDIST;
			        IF((NRDIST.LT.1).OR.(NRDIST.GT.100))[
			            OUTPUT;
			            (' *** NUMBER RADIAL BINS OUT OF RANGE (<1 OR >100),',
			            ' RESET TO 100 ***');
			            NRDIST=100;
			            ]
			        READ(9,*)(RDISTF(IB),RPDF(IB),IB=1,NRDIST);
			        CLOSE(UNIT=9);
			        OUTPUT NRDIST;('    HAVE READ',I5,' INPUT RADIAL BINS FROM FILE');
			        ]

			    "DO A CHECK ON THE RADIAL DISTRIBUTION"
			    ICOUNT=0;
			    RLAST=0;
			    IERROR=0;
			    :R-DIST-INPUT:LOOP[
			        ICOUNT=ICOUNT+1;
			        IF(ICOUNT.GT.NRDIST)[EXIT:R-DIST-INPUT:;]
			        IF(RDISTF(ICOUNT).LE.RLAST)[
			            IERROR=1;
			            OUTPUT;
			            (' *** RDISTF<=LAST ONE. NOT ALLOWED,',
			            ' TERMINATING RADIAL DISTRIBUTION INPUT.');
			            ]
			        ELSEIF(RDISTF(ICOUNT).GT.RSPH(NR))[
			            IERROR=1;
			            OUTPUT RSPH(NR);
			            (' *** RDISTF>',G14.7,', GEOMETRY SIZE.',
			            ' TERMINATING RADIAL DISTRIBUTION INPUT');
			            ]
			        ELSEIF(RPDF(ICOUNT).LT.0.0)[
			            IERROR=1;
			            OUTPUT;(' PDF < 0 NOT ALLOWED,',
			            ' TERMINATING RADIAL DISTRIBUTION INPUT');
			            ]
			        IF(IERROR.EQ.1)[
			            ICOUNT=ICOUNT-1;
			            IF(ICOUNT.EQ.0)[
			                OUTPUT;(' *** NO RADIAL DISTRIBUTION DEFINED,',
			                ' STOPPING EXECUTION ***');
			                STOP;
			                ]
			            OUTPUT ICOUNT;
			            (' RADIAL DITRIBUTION INPUT APPEARS TO BE INCOMPLETE,',
			            ' NRDIST RESET TO ',I12);
			            NRDIST=ICOUNT;
			            EXIT:R-DIST-INPUT:;
			            ]
			        RLAST=RDISTF(ICOUNT);
			        ]

			    OUTPUT RDISTF(NRDIST);
			    ('    RADIAL DISTRIBUTION RANGES FROM 0 TO',F12.3,' cm');

			    ]

			*/
			}
		}

		/**
		 * For a given int, this method covert it to a double. Useful to avoid integer division.
		 * @param n n
		 * @return the result
		 */
		public static double dble(int n)
		{
			Integer in = new Integer(n);
			return in.doubleValue();
		}

		/**
		 * Calculates one time only constants (such as WEIGHT or incident fluence AINFLU) that may vary with source type. 
		 */
		public static void SRCOTO()
		{
			if (EGS4.iGeom==EGS4.iCavity)//1=RZ Geom
			{
			//"Calculation of one time only constants that may vary with source type"
			IFPB=1;//"default flag to not being frontal parallel beam"
			pi=Math.PI;//just in case!!!

			if((ISOURC==0)||(ISOURC==2)||(ISOURC==4))
			{
			    //"frontal parallel beam source"
			    IFPB=0;//"set flag- passed in common source"
			    if((ISOURC==0)&&(WINC!=1.)) IFPB=1; //"not for angled incidence"
			    RBEAM=RBEAM*$ONE_EPS; //"scale down to insure a strike on the target"
			    RBEAM2=RBEAM*RBEAM;
			    WEIGHT=1.0; //"incident weight"
			    //"incident fluence"
			    AFACE=pi*RBEAM2;
			    if(ISOURC == 0){AINFLU=dble(NCASET)/AFACE;}
			    else{AINFLU=dble(NCASET);}
			}

			else if(ISOURC == 1)
			{// "frontal point source"
			    RBEAM=RBEAM*$ONE_EPS; //"scale down to insure a strike on the target"
			    RBEAM2=RBEAM*RBEAM;
			    DISTZ2=DISTZ*DISTZ;
			    DISTRH=0.0; //"on-axis"
			    ZCOFST=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			    ZSOFST=EGS4Geom.ZPLANE[0]-DISTZ;
			    AFACE=pi*RBEAM*RBEAM;
			    PROBFC=1.0;
			    //"incident fluence times solid angle"
			    //"  the solid angle factor is cancelled by the weighting routine"
			    //"  see bielajew,rogers and nahum, phys med biol 1985 for details"
			    AINFLU=dble(NCASET)/DISTZ2;
			}

			else if(ISOURC == 3)
			{// "isotropically radiating disk"
			    WEIGHT=1.0;
			    AINFLU=dble(NCASET);
				//"   Following are used to sample Z easily
			    ZSOFST=0.5*(ZSMAX - ZSMIN)+ ZSMIN;
			    ZBEAM=0.5*(ZSMAX - ZSMIN);
			}

			else if(ISOURC == 10)
			{// "SIDE PARALLEL BEAM"
			    XBEAM=XBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    ZBEAM=ZBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    ZCOFST=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			    RCYL1=EGS4Geom.RCYL[EGS4Geom.NR];
			    RCYL2=RCYL1*RCYL1;
			    ASIDE=4.0*XBEAM*ZBEAM;
			    WEIGHT=1.; //"INCIDENT WEIGHT"
			    AINFLU=dble(NCASET)/ASIDE; //"INCIDENT FLUENCE"
			}

			else if(ISOURC == 11)
			{// "SIDE POINT SOURCE CENTERED AT TARGET MIDDLE"
			    XBEAM=XBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    ZBEAM=ZBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    ZCOFST=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			    ZSOFST=ZCOFST; //"BEAM CENTERED ON TARGET MID-POINT"
			    RCYL1=EGS4Geom.RCYL[EGS4Geom.NR];
			    RCYL2=RCYL1*RCYL1;
			    DSTRH2=DISTRH*DISTRH;
			    ASIDE=4.0*XBEAM*ZBEAM;
			    PROBSD=1.0;
			    //"INCIDENT FLUENCE AT THE CENTER OF THE GEOMETRY TIMES SOLID ANGLE"
			    AINFLU = dble(NCASET)/DSTRH2;
			}

			else if(ISOURC == 12)
			{// "point source off axis"
			    RBEAM = RBEAM*$ONE_EPS; //"scale down to insure a strike on the target"
			    XBEAM = XBEAM*$ONE_EPS; //"scale down to insure a strike on the target"
			    ZBEAM = ZBEAM*$ONE_EPS; //"scale down to insure a strike on the target"
			    //"Z-axis offset of center of geom"
			    ZCOFST = 0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			    ZSOFST = EGS4Geom.ZPLANE[0]-DISTZ;  //"Z-axis offset of the source"
			    RBEAM2 = RBEAM*RBEAM;
			    RCYL1 = EGS4Geom.RCYL[EGS4Geom.NR];
			    RCYL2 = RCYL1*RCYL1;
			    DSTRH2 = DISTRH*DISTRH;
			    DISTZ2 = DISTZ*DISTZ;
			    DISTB = DISTZ-(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
			    //"DISTB is distance of point source from back of geomery"
			    DISTB2 = DISTB*DISTB;
			    AFACE = pi*RBEAM2;
			    ASIDE = 4.0*XBEAM*ZBEAM;
			    //"Incident fluence referred to midpoint"
			    AINFLU = NCASET/(DSTRH2+(ZCOFST-ZSOFST)*(ZCOFST-ZSOFST));
			    if(DISTRH <= EGS4Geom.RCYL[EGS4Geom.NR])//RCYL(NR))
			    {// "either front or back face exclusively"
			        PROBSD = 0.0;
			        if(DISTZ >= 0.0)
			        {//"front face only"
			        	PROBFC = 1.0;PROBBK = 0.0;
					}
			        else if(DISTZ <= (EGS4Geom.ZPLANE[0]-EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]))
			        {//"back face only"
			            PROBFC = 0.0;PROBBK = 1.0;
					}
				}
			    else if((DISTZ >= (EGS4Geom.ZPLANE[0]-EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]) &&
			    		(DISTZ <= 0.0)))
			    {
			        //"hits side exclusively"
			        PROBFC = 0.0;PROBBK = 0.0;PROBSD = 1.0;
				}
			    else
			    {// "intermediate case - can hit two surfaces"
			        COTANG = ZSOFST/DISTRH;
			        if(COTANG < 0.0)
			        {//"incident from front and side"
			            FACTOR = AFACE*Math.abs(COTANG)/ASIDE;
			            PROBSD = 1.0/(1.0+FACTOR);
			            PROBFC = FACTOR/(1.0+FACTOR);
			            PROBBK = 0.0;
					}
			        else
			        {//"incident from back and side"
			            FACTOR = AFACE*COTANG/ASIDE;
			            PROBSD = 1.0/(1.0+FACTOR);
			            PROBBK = FACTOR/(1.0+FACTOR);
			            PROBFC = 0.0;
					}
			    }//"end two surface case"
			}

			else if(ISOURC == 15) { src15_oto();}//(iout); }
			else if(ISOURC == 16) { src16_oto();}//(iout); ]

			else if(ISOURC == 13)
			{// "BROAD PARALLEL BEAM AT ANY ANGLE wrt THE TARGET"
			    RBEAM=RBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    XBEAM=XBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    ZBEAM=ZBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    RBEAM2=RBEAM*RBEAM;
			    RCYL1=EGS4Geom.RCYL[EGS4Geom.NR];
			    RCYL2=RCYL1*RCYL1;
			    ZCOFST=0.5*(EGS4Geom.ZPLANE[0]+EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]);
			    AFACE=pi*RBEAM2;
			    ASIDE=4.0*XBEAM*ZBEAM;
			    WEIGHT=1.0; //"INCIDENT WEIGHT"
			    AINFLU=NCASET/(Math.abs(WINC)*AFACE+Math.sqrt(UINC*UINC+VINC*VINC)*ASIDE);
			    if(WINC==1.0)
			    {
			        //"FRONT FACE EXCLUSIVELY"
			        PROBFC=1.0;PROBSD=0.0;PROBBK=0.0;
				}
			    else if(WINC==0.0)
			    {
			        //"HITS SIDE EXCLUSIVELY"
			        PROBFC=0.0;PROBSD=1.0;PROBBK=0.0;
				}
			    else if(WINC==-1.0)
			    {
			        //"BACK FACE EXCLUSIVELY"
			        PROBFC=0.0;PROBSD=0.0;PROBBK=1.0;
				}
			    else
			    {
			        //"INTERMEDIATE CASE - CAN HIT TWO SURFACES"
			        if(WINC<0.0)
			        {//"INCIDENT FROM BACK AND SIDE"
			            FACTOR=AFACE*Math.abs(WINC)/(ASIDE*Math.sqrt(UINC*UINC+VINC*VINC));
			            PROBFC=0.0;
			            PROBSD=1.0/(1.0+FACTOR);
			            PROBBK=FACTOR/(1.0+FACTOR);
					}
			        else
			        {//"INCIDENT FROM FRONT AND SIDE"
			            FACTOR=AFACE*WINC/(ASIDE*Math.sqrt(UINC*UINC+VINC*VINC));
			            PROBFC=FACTOR/(1.0+FACTOR);
			            PROBSD=1.0/(1.0+FACTOR);
			            PROBBK=0.0;
					}
				}
			}

			else if(ISOURC == 14)
			{// "FRONTAL POINT SOURCE restricted radius"
			    RBEAM=RBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    RBEAM2=RBEAM*RBEAM;
			    RMINSQ=RMINBM*RMINBM;
			    DISTZ2=DISTZ*DISTZ;
			    DISTRH=0.0; //"ON-AXIS"
			    ZCOFST=0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			    ZSOFST=EGS4Geom.ZPLANE[0]-DISTZ;
			    AFACE=pi*RBEAM*RBEAM;
			    PROBFC=1.0;
			    //"INCIDENT FLUENCE TIMES SOLID ANGLE"
			    //"  THE SOLID ANGLE FACTOR IS CANCELLED BY THE WEIGHTING ROUTINE"
			    //"  SEE BIELAJEW,ROGERS AND NAHUM, PHYS MED BIOL 1985 FOR DETAILS"
			    AINFLU=dble(NCASET)/DISTZ2;
			}

			else if(ISOURC == 20)
			{
			    //"CALCULATE THE CPDF FROM THE PDF AND NORMALIZE IT"
			    //"NOTE THAT THE RADIAL WEIGHTING IS INCLUDED AT THIS STAGE"
			    RCDF[0]=0.5*RDISTF[0]*RDISTF[0]*RPDF[0];
			    for(int IB=2;IB<=NRDIST;IB++)
			    {
			        RCDF[IB-1]=RCDF[IB-2]+
			        0.5*(RDISTF[IB-1]*RDISTF[IB-1]-RDISTF[IB-2]*RDISTF[IB-2])*RPDF[IB-1];
				}
			    FNORM=1./RCDF[NRDIST-1];
			    IBNSOK=0;
			    GRIDSZ=1./dble($MXRDIST);
			    for(int IB=1;IB<=NRDIST;IB++)
			    {
			        RCDF[IB-1]=FNORM*RCDF[IB-1];
			        if(IB==1){if(RCDF[0]<GRIDSZ)IBNSOK=1;}
			        else if((RCDF[IB-1]-RCDF[IB-2])<GRIDSZ){IBNSOK=1;}
			    }//"END OF LOOP ON IB"
			    if(IBNSOK!=0)
			    {
					EGS4.seqStr="WARNING: SOME OF NORMALIZED BIN PROBABILITIES ,SO SMALL BINS MAY BE MISSED";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}

			    //"CALCULATE RCDFIN - AN ARRAY WHICH ALLOWS THE RAPID SAMPLING FOR THE"
			    //"RADIAL BY PRECOMPUTING THE RESULTS FOR A FINE GRID"
				boolean FOUND_R_BIN=false;
			    for(int K=1;K<=$MXRDIST;K++)
			    {
					FOUND_R_BIN=false;
			        AK=dble(K)*GRIDSZ;
			        for(int I=1;I<=NRDIST;I++)
			        {
						FOUND_R_BIN=false;
						IRDIST=I;
						if(AK<=RCDF[I-1])
						{
							//GOTO :FOUND-R-BIN:;
							FOUND_R_BIN=true;
							break;
						}
					}
			        //"WE SHOULD NEVER FALL THROUGH TO HERE"
			        if(!FOUND_R_BIN)
			        {
						EGS4.seqStr="FELL THROUGH RADIAL SAMPLING ROUTINE!!";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
					}
			        //:FOUND-R-BIN:
			        if(IRDIST!=1)
			        {
						//RCDFIN(K,1)=RDISTF(IRDIST-1)**2;
						RCDFIN[K-1][0]=RDISTF[IRDIST-2]*RDISTF[IRDIST-2];
					}
			        else{RCDFIN[K-1][0]=0.0;}
			        RCDFIN[K-1][1]=RDISTF[IRDIST-1]*RDISTF[IRDIST-1]-RCDFIN[K-1][0];
			    }//"END LOOP OVER K"
			    //"OUTPUT IF IWATCH IS ACTIVE"
			    //if(IWATCH.NE.0)[
			    //    OUTPUT;(/' RDISTF,   RPDF,   RCDF:'//);
			    //    OUTPUT (RDISTF(IB),RPDF(IB),RCDF(IB),IB=1,NRDIST);(3E17.7);
			    //    OUTPUT;(//);
			    //]
				EGS4.seqStr=" RDISTF,   RPDF,   RCDF:";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				for (int IB=1;IB<=NRDIST;IB++)
				{
					EGS4.seqStr=EGS4.format(RDISTF[IB-1],17,false)+
					EGS4.format(RPDF[IB-1],17,false)+EGS4.format(RCDF[IB-1],17,false);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}

			    WEIGHT=1.0; //"INCIDENT WEIGHT"
			    //"calculate incident fluence"
			    AINFLU=dble(NCASET)/(pi*RDISTF[NRDIST-1]*RDISTF[NRDIST-1]);
			}
	/*
	//NOT ALLOWED HERE!!!!!
			ELSEIF( ISOURC = 23 ) [
			    AINFLU=dble(NCASET); WEIGHT=1.0;
			    IF( chamber_c > 1e4 | chamber_c < -1e4 ) [
			        chamber_c = 0.5*(zplane(1) + zplane(nplane));
			    ]
			]
	//NOT ALLOWED HERE!!!!!
			ELSEIF((ISOURC = 21) | (ISOURC = 22))["full phase space"
			     AINFLU=dble(NCASET); "number of incident particles"
			     WEIGHT=1.0; "INCIDENT WEIGHT"
			     IF( chamber_c > 1e4 | chamber_c < -1e4 ) [
			         chamber_c = 0.5*(zplane(1) + zplane(nplane));
			     ]
			     "Note that the weight for each particle is picked up from the phase"
			     "space file on the SRCHST call for these 2 sources"

			     "Now, calculate the initial value of NPHSPN, the no.-1 of the first"
			     "particle to use"
			     IF(IPARALLEL>1 & PARNUM > 0)["one of IPARALLEL runs"
			        IF(IHSTRY<NCASE_PHSP/IPARALLEL)[
			           NPHSPN=INT((PARNUM-1)*NCASE_PHSP/IPARALLEL)+IHSTRY;
			        ]
			        ELSE[
			           tmp_mod = NCASE_PHSP;
			           NPHSPN=tmp_mod*(PARNUM-1)/IPARALLEL+
			                   MOD(IHSTRY,tmp_mod/IPARALLEL);
			        ]
			     ]
			     ELSE["not a parallel run"
			         IF(IHSTRY<NCASE_PHSP)[NPHSPN = IHSTRY; ]
			         ELSE[ tmp_mod = NCASE_PHSP; NPHSPN=MOD(IHSTRY,tmp_mod); ]
			     ]
			]


			RETURN;

			"************************************************************************

			*/
			}
			else if(EGS4.iGeom==EGS4.iCavitySPH)//2=SPH Geom
			{

			//"CALCULATION OF ONE TIME ONLY CONSTANTS THAT MAY VARY WITH SOURCE TYPE"

			if(ISOURC==0)
			{
			    //"PARALLEL BEAM AT ANY ANGLE"
			    IDSTON=0; //"NO RADIAL DISTRIBUTION"
			    ISRCTY=0; //"PARALLEL BEAM"
			    RBEAM=RBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    RBEAM2=RBEAM*RBEAM;
			    ABEAM=Math.PI*RBEAM2; //"BEAM'S AREA ON THE FRONT OF THE TARGET"
			    WEIGHT=1.0; //"INCIDENT WEIGHT"
			    AINFLU=dble(NCASET)/ABEAM; //"INCIDENT FLUENCE"
			}

			else if(ISOURC==1)
			{
			    //"POINT SOURCE AT ANY ANGLE"
			    IDSTON=0; //"NO RADIAL DISTRIBUTION"
			    ISRCTY=1; //"POINT SOURCE"
			    RBEAM=RBEAM*$ONE_EPS; //"SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			    RBEAM2=RBEAM*RBEAM;
			    ABEAM=Math.PI*RBEAM2; //"BEAM'S AREA ON THE FRONT OF THE TARGET"
			    DISTR2=DISTR*DISTR;
			    //"INCIDENT FLUENCE TIMES SOLID ANGLE"
			    //"  THE SOLID ANGLE FACTOR IS CANCELLED BY THE WEIGHTING ROUTINE"
			    //"  SEE BIELAJEW,ROGERS AND NAHUM, PHYS MED BIOL 1985 FOR DETAILS"
			    AINFLU=dble(NCASET)/DISTR2;
			}

			else if(ISOURC==4)
			{//"point spread function calculation"
			    WEIGHT = 1.0;//"incident weight"
			    AINFLU=dble(NCASET); //"INCIDENT FLUENCE"
			                          //"i.e., we give dose/incident particle"
			}
	/*
			ELSEIF((ISOURC.EQ.10).OR.(ISOURC.EQ.11))[ "RADIAL DISTRIBUTION"

			    IDSTON=1; "RADIAL DISTRIBUTION FLAG ON"

			    "CALCULATE THE CPDF FROM THE PDF AND NORMALIZE IT"
			    "NOTE THAT THE RADIAL WEIGHTING IS INCLUDED AT THIS STAGE"
			    RCDF(1)=0.5*RDISTF(1)**2*RPDF(1);
			    DO IB=2,NRDIST[
			        RCDF(IB)=RCDF(IB-1)+0.5*(RDISTF(IB)**2-RDISTF(IB-1)**2)*RPDF(IB);
			        ]
			    FNORM=1./RCDF(NRDIST);
			    IBNSOK=0;
			    GRIDSZ=1./FLOAT($MXRDIST);
			    DO IB=1,NRDIST[
			        RCDF(IB)=FNORM*RCDF(IB);
			        IF(IB.EQ.1)[IF(RCDF(1).LT.GRIDSZ )IBNSOK=1;]
			        ELSEIF((RCDF(IB)-RCDF(IB-1)).LT.GRIDSZ)[IBNSOK=1.0;]
			        ]"END OF LOOP ON IB"
			    IF(IBNSOK.NE.0)[OUTPUT;(//'0*******WARNING******'/
			    T20,'SOME OF NORMALIZED BIN PROBABILITIES SO SMALL BINS MAY BE MISSED'/);]
			    "CALCULATE RCDFIN - AN ARRAY WHICH ALLOWS THE RAPID SAMPLING FOR THE"
			    "RADIAL BY PRECOMPUTING THE RESULTS FOR A FINE GRID"
			    DO K=1,$MXRDIST[
			        AK=FLOAT(K)*GRIDSZ;
			        DO I=1,NRDIST[
			            IRDIST=I;IF(AK.LE.RCDF(I))[GOTO :FOUND-R-BIN:;]
			            ]
			        "WE SHOULD NEVER FALL THROUGH TO HERE"
			        OUTPUT;(' *** FELL THROUGH RADIAL SAMPLING ROUTINE *** ');
			        :FOUND-R-BIN:
			        IF(IRDIST.NE.1)[RCDFIN(K,1)=RDISTF(IRDIST-1)**2;]
			        ELSE[RCDFIN(K,1)=0.0;]
			        RCDFIN(K,2)=RDISTF(IRDIST)**2-RCDFIN(K,1);
			        ]"END LOOP OVER K"
			    "OUTPUT IF IWATCH IS ACTIVE"
			    IF(IWATCH.NE.0)[
			        OUTPUT;(' RDISTF,   RPDF,   RCDF:'//);
			        OUTPUT (RDISTF(IB),RPDF(IB),RCDF(IB),IB=1,NRDIST);(3E17.7);
			        OUTPUT;(//);
			        ]

			    IF(ISOURC.EQ.10)[ "PARALLEL BEAM"
			        ISRCTY=0; "PARALLEL BEAM"
			        WEIGHT=1.0; "INCIDENT WEIGHT"
			        ABEAM=PI*RDISTF(NRDIST)**2; "BEAM'S AREA ON THE FRONT OF THE TARGET"
			        AINFLU=FLOAT(NCASET)/ABEAM;
			        ]
			    ELSE[ "POINT SOURCE"
			        ISRCTY=1; "POINT SOURCE"
			        DISTR2=DISTR**2;
			        "INCIDENT FLUENCE TIMES SOLID ANGLE"
			        "  THE SOLID ANGLE FACTOR IS CANCELLED BY THE WEIGHTING ROUTINE"
			        "  SEE BIELAJEW,ROGERS AND NAHUM, PHYS MED BIOL 1985 FOR DETAILS"
			        AINFLU=FLOAT(NCASET)/DISTR2;
			        "RSHADW=MAXIMUM BEAM RADIUS ALLOWED BY SHADOWING"
			        RSHADW=(RSPH(NR)/DISTR)*SQRT((DISTR+RSPH(NR))*(DISTR-RSPH(NR)));
			        RBEAM=RSHADW; "NEEDED TO REJECT LATER IF R SELECTED IS TOO BIG"
			        RBEAM=RBEAM*$ONE-EPS; "SCALE DOWN TO INSURE A STRIKE ON THE TARGET"
			        ABEAM=PI*RBEAM2; "BEAM'S AREA ON THE FRONT OF THE TARGET"
			        ]
			    ]

			RETURN;
			*/
			}
		}

		/**
		 * Calculate the history dependent constants that may change depending on the source type.
		 */
		public static void SRCHST()
		{
			if (EGS4.iGeom==EGS4.iCavity)//1=RZ geom
			{
			int IZ=0;int IX=0;
			/*
	//NOT ALLOWED HERE!!!!!
			IF(n_parallel>0 & (ISOURC=21 | ISOURC=22))[
			"set up chunk of phase space file to sample"
			    "if this is the first chunk in this run, calculate no. of particles in a"
			    "phsp chunk"
			    IF(N_RUN_CHUNK_OLD=0) P_PER_PHSP_CHUNK=NCASE_PHSP/(n_parallel*$N_CHUNKS);
			    N_RUN_CHUNK=(NCASE-N_LEFT)*n_parallel*$N_CHUNKS/NCASE;
			    IF(N_RUN_CHUNK ~= N_RUN_CHUNK_OLD)["have moved on to a new run chunk"
			       N_RUN_CHUNK_OLD=N_RUN_CHUNK;
			       NPHSPN_MIN=(N_RUN_CHUNK-1)*P_PER_PHSP_CHUNK+1;
			       IF(N_LEFT=0)["this is the last run just use up the rest of the"
			                     "phsp source"
			           NPHSPN_MAX=NCASE_PHSP;
			       ]
			       ELSE["calculate the max value of INPHSP"
			           NPHSPN_MAX=NPHSPN_MIN+P_PER_PHSP_CHUNK-1;
			       ]
			       NPHSPN=NPHSPN_MIN-1; "srchst later adds 1 to NPHSPN"
			       CYCLNUM=0; "reset NRCYCL counter"
			        write(6,'(/a/,a,i12,a,i12/,a//)')
			  '      This simulation uses a phase space source.',
			  '      This run will use from particle',NPHSPN_MIN,' to particle ',
			         NPHSPN_MAX,
			  '      in the source file.';
			    ]
			]
	*/

			//"Calculate the history dependent constants that may change depending on"
			//"the source type"
			if((ISOURC == 0) || (ISOURC == 2) ||(ISOURC == 4))
			{//"frontal parallel beam source"
			    if(RBEAM == 0.0)
			    { //"pencil source"
			        xin=0.0;yin=0.0;//XIN,YIN/=0.0;
			        irin=2;
				}
			    else
			    { //"choose a point randomly in a circle"
			        CHOOSE_POINT_IN_CIRCLE();
			        irin=2+(IXIN-1)*EGS4Geom.NZ;
				}
			    zin=EGS4Geom.ZPLANE[0]; //"INCIDENT Z POSITION"
			    if(IFPB == 0)
			    {
					uin=0.0;vin=0.0;win=1.0;
				}
				else
				{
					uin=UINC;vin=VINC;win=WINC;
				}
			    NRCFLG=10; //"geometrical flag"
			    WEIGHT=1.0;// "this version needs to set the weight in case correlation"
			               // "changed it"
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 1 || ISOURC == 14)
			{//"frontal point source on axis"
			    ENTRY_FRONT_FACE();
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 3)
			{// "isotropically radiating disk"
			 //   "choose z-coordinate for the generation of a particle within region NSRCRG"
			    zin=EGS4.random01();
			    zin=ZSOFST+(2.0*zin-1.0)*ZBEAM;
			    for(IZ=1;IZ<=EGS4Geom.NPLANE-1;IZ++)
			    {
					//IF(ZIN <= ZPLANE(IZ+1) & ZIN >= ZPLANE(IZ))
					if(zin <= EGS4Geom.ZPLANE[IZ] && zin >= EGS4Geom.ZPLANE[IZ-1])
						break;
				}
			    //"choose a point randomly in a circle that has the radius of the source"
			    //"choose a point randomly in a ring" "this is not the fastest way"
			    while(true)
			    {
			        xin=EGS4.random01();xin=(2.0*xin-1.0)*RBEAM;
			        yin=EGS4.random01();yin=(2.0*yin-1.0)*RBEAM;
			        R2=xin*xin+yin*yin;

			        if(R2 <= RBEAM2 && R2 >= RMINSQ)
			        	break;
			    }// UNTIL (R2 <= RBEAM2 & R2 >= RMINSQ);

			    for(IX=1;IX<=EGS4Geom.NR;IX++)
			    {
				   if(R2 <= EGS4Geom.CYRAD2[IX-1])
				   		break;
			    }

			    irin = IZ + EGS4Geom.NZ*(IX-1)+1;

			    //"now determine initial direction cosines"
			    EGS4.COSTHE=EGS4.random01();   EGS4.COSTHE=2.*EGS4.COSTHE-1;
			    EGS4.SINTHE=Math.sqrt(1.0-EGS4.COSTHE*EGS4.COSTHE);
			    double PHI=EGS4.random01();PHI=EGS4.TWOPI*PHI;
			    uin=EGS4.SINTHE*Math.cos(PHI);
			    vin=EGS4.SINTHE*Math.sin(PHI);
			    win=EGS4.COSTHE;
			    NRCFLG=50;  //"geometrical flag meaning particle generated within region"
			    WEIGHT=1.0;
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 10)
			{// "side parallel beam"
			    //"choose a point randomly in the rectangle"
			    CHOOSE_POINT_IN_RECTANGLE();
			    uin=0.0;vin=-1.0;win=0.0; //"incident angles"
			    NRCFLG=20; //"geometrical flag"
			    WEIGHT=1.0;// "this version needs to set the weight in case correlation"
			               // "changed it"
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 11)
			{// "SIDE POINT SOURCE CENTERED AT TARGET MIDDLE"
			 //   "SIDE POINT SOURCE CENTERED AT TARGET MIDDLE"
			    ENTRY_SIDE();
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 12)
			{// "POINT SOURCE OFF AXIS"
			 //    "POINT SOURCE OFF AXIS"
			 //    "DECIDE WHICH SURFACE THE INCIDENT PARTICLE WILL HIT"
			    if(PROBFC==1.0){ENTRY_FRONT_FACE();}
			    else if(PROBBK==1.0){ENTRY_BACK_FACE();}
			    else if(PROBSD==1.0){ENTRY_SIDE();}
			    else
			    {
			        //"MAY HIT EITHER FACE OR SIDE"
			        double WHICH=EGS4.random01();
			        if(WHICH<=PROBSD){ENTRY_SIDE();}
			        else if(WHICH<=(PROBSD+PROBFC)){ENTRY_FRONT_FACE();}
			        else{ENTRY_BACK_FACE();}
				}
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 15)
			{
			    src15_hst();//(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 16)
			{
			    src16_hst();//(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
			    NHSTRY=NHSTRY+1;
			}

			else if(ISOURC == 13)
			{// "BROAD PARALLEL BEAM AT ANY ANGLE wrt THE TARGET"
			 //   "BROAD PARALLEL BEAM AT ANY ANGLE WITH RESPECT TO THE TARGET"
			    uin=UINC;vin=VINC;win=WINC;
			    //"DECIDE WHICH SURFACE THE INCIDENT PARTICLE WILL HIT"
			    if(PROBFC==1.0)
			    {
			        CHOOSE_POINT_IN_CIRCLE();
			        irin=2+(IXIN-1)*EGS4Geom.NZ;
			        zin=EGS4Geom.ZPLANE[0]; //"INCIDENT Z POSITION"
			        NRCFLG=10; //"GEOMETRICAL FLAG"
				}
			    else if(PROBBK==1.0)
			    {
			        CHOOSE_POINT_IN_CIRCLE();
			        irin=1+IXIN*EGS4Geom.NZ;
			        zin=EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]; //"INCIDENT Z POSITION"
			        NRCFLG=30; //"GEOMETRICAL FLAG"
				}
			    else if(PROBSD==1.0)
			    {
			        CHOOSE_POINT_IN_RECTANGLE();
			        NRCFLG=20; //"GEOMETRICAL FLAG"
				}
			    else
			    {
			        //"MAY HIT EITHER FACE OR SIDE"
			        double WHICH=EGS4.random01();
			        if(WHICH<=PROBSD)
			        {
			            CHOOSE_POINT_IN_RECTANGLE();
			            NRCFLG=20; //"GEOMETRICAL FLAG"
					}
			        else if(WHICH<=(PROBSD+PROBFC))
			        {
			            CHOOSE_POINT_IN_CIRCLE();
			            irin=2+(IXIN-1)*EGS4Geom.NZ;
			            zin=EGS4Geom.ZPLANE[0];// "INCIDENT Z POSITION"
			            NRCFLG=10; //"GEOMETRICAL FLAG"
					}
			        else
			        {
			            CHOOSE_POINT_IN_CIRCLE();//$CHOOSE-POINT-IN-CIRCLE;
			            irin=1+IXIN*EGS4Geom.NZ;
			            zin=EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1];// "INCIDENT Z POSITION"
			            NRCFLG=30;// "GEOMETRICAL FLAG"
					}
				}
			    WEIGHT=1.0; //"THIS VERSION NEEDS TO SET THE WEIGHT IN CASE CORRELATION"
			                //"CHANGED IT"
			    NHSTRY=NHSTRY+1;
			}
			//"ISOURC=14 DONE WITH ISOURC=1"

			else if(ISOURC == 20)
			{//"RADIAL DISTRIBUTION, FRONT PARALLEL BEAM"
			 //   "THIS ENTRY DOES THE ACTUAL SAMPLING OF THE INCIDENT RADIAL DISTRIBUTION"
			 //   "THIS WILL RETURN A HISTOGRAM OF VALUES"
			    double RNNO1=EGS4.random01();double RNNO2=EGS4.random01();
			 //   "IN NEXT STATEMENT MIN AVOIDS K=$MXRDIST+1, ALMOST NEVER IMPLEMENTED"
			    Double dbl1=new Double(dble($MXRDIST)*RNNO1+1.);
			    int dbl1i=dbl1.intValue();
			    Double dbl2=new Double(dble($MXRDIST));
			    int dbl2i=dbl2.intValue();
			    int K=Math.min(dbl1i,dbl2i);
			 //   "THIS VERSION ONLY FOR CYLINDRICAL SYMMETRY"
			    //xin=Math.sqrt(RCDFIN(K,1)+RNNO2*RCDFIN(K,2));
			    xin=Math.sqrt(RCDFIN[K-1][0]+RNNO2*RCDFIN[K-1][1]);
			    yin=0.0;
			    zin=EGS4Geom.ZPLANE[0]; //"INCIDENT POSITION"
			 //   "NOW DO A SEQUENTIAL SEARCH FOR THE ENTRANCE REGION"
			    for(IX=1;IX<=EGS4Geom.NR;IX++)
			    {
					IXIN=IX;
					if(xin<=EGS4Geom.RCYL[IX-1])
						break;//EXIT;
				}
			    irin=2+(IXIN-1)*EGS4Geom.NZ;
			    win=1.0;uin=0.0;vin=0.0; //"INCIDENT ANGLES"
			    NRCFLG=10; //"GEOMETRICAL FLAG"
			    WEIGHT=1.0;// "THIS VERSION NEEDS TO SET THE WEIGHT IN CASE CORRELATION"
			               // "CHANGED IT"
			    NHSTRY=NHSTRY+1;
			}
	/*
	//NOT ALLOWED HERE!!!!!
			ELSEIF(ISOURC = 21)["Full phase space for each particle"

			 IF(NRCYCL>0 & CYCLNUM>0 & CYCLNUM<=NRCYCL)["recycle this particle"
			    CYCLNUM=CYCLNUM+1;
			    XIN=XINOLD;
			    YIN=YINOLD;
			    ZIN=ZINOLD;
			    UIN=UINOLD;
			    VIN=VINOLD;
			    WIN=WINOLD;
			    IRIN=IRINOLD;
			    NRCFLG=NRCFLGOLD;
			    WEIGHT=WEIGHTOLD;
			    "the rest of the parameters are global and unchanged"
			    NNREAD=NNREAD+1;
			 ]
			 ELSE["get a new particle from the file"
			    :READ-PARTICLE:;

			    NPHSPN=NPHSPN+1;

			    IF(IPARALLEL>1 & PARNUM>0 & NPHSPN > INT(PARNUM*NCASE_PHSP/IPARALLEL))[
			       NPHSPN=INT((PARNUM-1)*NCASE_PHSP/IPARALLEL)+1;
			       OUTCNT=OUTCNT+1;
			       OUTPUT;
			         (' ***WARNING*** Used all particles from partition in source file.'/
			          '               Restarting from first particle in partition.');
			       WRITE(IOUT,
			        '('' ***WARNING*** Used all particles from partition in source file.''/
			          ''               Restarting from first particle in partition.'')');
			    ]
			    ELSEIF(n_parallel>0 & NPHSPN>NPHSPN_MAX)[
			       "have reached the end of the phsp chunk in this parallel run"
			       NPHSPN=NPHSPN_MIN;
			       OUTCNT=OUTCNT+1;
			       OUTPUT ;
			     (///' ***WARNING*** USED ALL PARTICLES FROM CHUNK IN SOURCE FILE!'/
			       '               RESTARTING FROM FIRST PARTICLE IN THIS CHUNK.'//
			       /1x,79('*')// );
			    ]
			    ELSEIF(NPHSPN>NCASE_PHSP)[
			       NPHSPN=1;
			       OUTCNT=OUTCNT+1;
			       OUTPUT;
			         (' ***WARNING*** Used all particles from source file.'/
			          '               Restarting from first particle.');
			       WRITE(IOUT,
			        '('' ***WARNING*** Used all particles from source file.''/
			          ''               Restarting from first particle.'')');
			    ]

			    IF(OUTCNT>1000)["too many restart warnings"
			       OUTPUT;(' ***ERROR*** >1000 restart warnings.');
			       WRITE(IOUT,
			        '('' ***ERROR*** >1000 restart warnings.'')');
			       STOP;
			    ]

			    $READ_PHSP(IMODE,42,NPHSPN+1:NHSTRY,NPASS,IQIN,WIN,ZLAST,
			               LATCHI,EIN,WEIGHT,XIN,YIN,UIN,VIN);

			    NNREAD=NNREAD+1;

			    IF(DOSE_STAT=1 | NHSTRY=0) NHSTRY=NHSTRY+1;
			    "NHSTRY may also remain zero if this is a partition of a phase space"
			    "file and it does not begin with a new primary history"

			    "We only use the forward particles."
			    IF( NPASS=1 ) GOTO :READ-PARTICLE:; "...signals not 1st time crossing"

			    "check charge of this particle, IQIN, is what the user asked for, IQINC"
			    "IQINC = 0 => photons only"
			    "IQINC = -1 =>electrons only"
			    "IQINC = +1 => positrons only"
			    "IQINC = 2  => all particles"
			    "IQINC = 3 => charged particles"

			    IF(IQIN =-1 & IQINC~=-1 & IQINC~=2 & IQINC ~=3)[ GOTO :READ-PARTICLE:;]
			    IF(IQIN = 1 & IQINC~= 1 & IQINC~=2 & IQINC ~=3)[ GOTO :READ-PARTICLE:;]
			    IF(IQIN = 0 & IQINC~= 0 & IQINC~=2)[ GOTO :READ-PARTICLE:;]
			    "If it gets here, it is acceptable"

			    "Now do a sequential search for the entrance region"
			     R2 = XIN**2 +YIN**2;
			    "Check if particle is outside geometry"
			    IF(R2 > CYRAD2(NR))[ GOTO :READ-PARTICLE:;]
			    DO IX=1,NR[
			       IXIN=IX;
			       IF( R2 <=  CYRAD2(IX) ) EXIT;  "CYRAD2(IX) is RCYL(IX)**2 set in GEOMRZ"
			    ]
			    IRIN=2+(IXIN-1)*NZ;
			    ZIN=ZPLANE(1);
			    NRCFLG = 10; "GEOMETRICAL FLAG, SOURCE INCIDENT ON FRONT FACE"
			    IF(NRCYCL>0)["store non-global parameters for recycling"
			       CYCLNUM=1;
			       XINOLD=XIN;
			       YINOLD=YIN;
			       ZINOLD=ZIN;
			       UINOLD=UIN;
			       VINOLD=VIN;
			       WINOLD=WIN;
			       IRINOLD=IRIN;
			       NRCFLGOLD=NRCFLG;
			       WEIGHTOLD=WEIGHT;
			    ]
			 ]
			]
	//NOT ALLOWED HERE!!!!!
			ELSEIF( ISOURC = 23 ) [ "full BEAM simulation"

			    LOOP [;
			        :retry_sample_beamsource:;
			        call sample_beamsource(ein,xin,yin,zin,uin,vin,win,weight,iqin,latchi,
			                               nhstry,iphatin);
			        IF( iqinc < 2 & iqin ~= iqinc ) goto :retry_sample_beamsource:;
			        IF( iqinc = 3 & iqin = 0 ) goto :retry_sample_beamsource:;
			        IF( weight < min_weight_23 | weight > max_weight_23 ) [
			            "write(6,*) 'Not in weight window: ',weight,iqin,ein,xin,yin;"
			            goto :retry_sample_beamsource:;
			        ]

			        IF( secret_option_23 = 1 ) [
			            IF( xin > 0 ) [ xin = -xin; uin = -uin; ]
			        ]

			        xin_tmp = xin; yin_tmp = yin; zin_tmp = zin;

			        "apply the offset"
			        xin=xin+xoffset; yin=yin+yoffset;

			        "Rotate the position"
			        zin = -dist_phsp*cost_phsp + yin*sint_phsp + chamber_c;
			        yin = dist_phsp*sint_phsp + yin*cost_phsp;

			        "Rotate the direction"
			        vtemp = vin;
			        vin = -win*sint_phsp + vin*cost_phsp;
			        win =  win*cost_phsp + vtemp*sint_phsp;

			        "write(6,*) 'Got ',iqin,ein,xin,yin,zin,uin,vin,win,weight;
			        "Check if particle is within the geometry or will enter it"
			        radp = xin*xin + yin*yin;
			        IF( radp > CYRAD2(nr) | zin < zplane(1) | zin > zplane(nz+1) ) [
			            " outside of geometry. check whether that particle will get in "
			            check = .true.;
			            IF( zin < zplane(1) & win > 0 ) [ tf = (zplane(1)-zin)/win; ]
			            ELSE IF( zin > zplane(nz+1) & win < 0 ) [
			                          tf = (zplane(nz+1) - zin)/win;
			            ]
			            ELSE [ check = .false.; ]
			            IF( check ) [
			                xtmp = xin + uin*tf; ytmp = yin + vin*tf;
			                radp = xtmp*xtmp + ytmp*ytmp;
			                IF( radp <= CYRAD2(nr) ) [
			                    xin = xtmp; yin = ytmp;
			                    IF( zin < zplane(1) ) [ zin = zplane(1); iz = 1; ]
			                    ELSE                  [ zin = zplane(nz+1); iz = nz; ]
			                    DO ix=1,nr [ IF( radp <= CYRAD2(ix) ) EXIT; ]
			                    EXIT;
			                ]
			                radp = xin*xin + yin*yin;
			            ]
			            IF( radp < CYRAD2(nr) ) NEXT;
			            phbb = uin*uin + vin*vin;
			            IF( phbb < 1e-15 ) NEXT;
			            pha = (xin*uin + yin*vin)/phbb;
			            IF( pha > 0 ) NEXT;
			            phb = (radp - CYRAD2(nr))/phbb;
			            phd = pha*pha - phb;
			            IF( phd < 0 ) NEXT;
			            tf = -pha-sqrt(phd);
			            IF( tf < 0 ) NEXT;
			            zin = zin + tf*win;
			            IF( zin < zplane(1) | zin > zplane(nz+1) ) NEXT;
			            ix = nr; xin = xin + uin*tf; yin = yin + vin*tf;
			            DO iz=2,nplane [ IF( zin <= zplane(iz) ) EXIT; ]
			            iz = iz-1;
			            EXIT;
			        ]
			        ELSE [
			            " already inside, find the region index "
			            DO iz=2,nplane [ IF( zin <= zplane(iz) ) EXIT; ]
			            iz = iz-1;
			            DO ix=1,nr [ IF( radp <= CYRAD2(ix) ) EXIT; ]
			            EXIT;
			        ]
			    ]
			    irin = 1 + (ix-1)*nz + iz;
			    "write(6,*) 'Got ',iqin,ein,xin,yin,zin,uin,vin,win,weight,ix,iz;
			    ihstry = ihstry + nhstry - last_nhstry;
			    last_nhstry = nhstry;
			]
	//NOT ALLOWED HERE!!!!!
			ELSEIF(ISOURC = 22)["full phase space for each particle: from side"

			 IF(NRCYCL>0 & CYCLNUM>0 & CYCLNUM<=NRCYCL)["recycle this particle"
			    CYCLNUM=CYCLNUM+1;
			    XIN=XINOLD;
			    YIN=YINOLD;
			    ZIN=ZINOLD;
			    UIN=UINOLD;
			    VIN=VINOLD;
			    WIN=WINOLD;
			    IRIN=IRINOLD;
			    NRCFLG=NRCFLGOLD;
			    WEIGHT=WEIGHTOLD;
			    "the rest of the parameters are global and unchanged"
			    NNREAD=NNREAD+1;
			 ]
			 ELSE["get a new particle from the file"
			    :READ-PARTICLE-SRC22:;

			    NPHSPN=NPHSPN+1;

			    IF(IPARALLEL>1 & PARNUM>0 & NPHSPN > INT(PARNUM*NCASE_PHSP/IPARALLEL))[
			       NPHSPN=INT((PARNUM-1)*NCASE_PHSP/IPARALLEL)+1;
			       OUTCNT=OUTCNT+1;
			       OUTPUT;
			         (' ***WARNING*** Used all particles from partition in source file.'/
			          '               Restarting from first particle in partition.');
			       WRITE(IOUT,
			        '('' ***WARNING*** Used all particles from partition in source file.''/
			          ''               Restarting from first particle in partition.'')');
			    ]
			    ELSEIF(n_parallel>0 & NPHSPN>NPHSPN_MAX)[
			       "have reached the end of the phsp chunk in this parallel run"
			       NPHSPN=NPHSPN_MIN;
			       OUTCNT=OUTCNT+1;
			       OUTPUT ;
			     (///' ***WARNING*** USED ALL PARTICLES FROM CHUNK IN SOURCE FILE!'/
			       '               RESTARTING FROM FIRST PARTICLE IN THIS CHUNK.'//
			       /1x,79('*')// );
			    ]
			    ELSEIF(NPHSPN>NCASE_PHSP)[
			       NPHSPN=1;
			       OUTCNT=OUTCNT+1;
			       OUTPUT;
			         (' ***WARNING*** Used all particles from source file.'/
			          '               Restarting from first particle.');
			       WRITE(IOUT,
			        '('' ***WARNING*** Used all particles from source file.''/
			          ''               Restarting from first particle.'')');
			    ]

			    IF(OUTCNT>1000)["too many restart warnings"
			       OUTPUT;(' ***ERROR*** >1000 restart warnings.');
			       WRITE(IOUT,
			        '('' ***ERROR*** >1000 restart warnings.'')');
			       STOP;
			    ]

			    $READ_PHSP(IMODE,42,NPHSPN+1:NHSTRY,NPASS,IQIN,WIN,ZLAST,
			               LATCHI,EIN,WEIGHT,XIN,YIN,UIN,VIN);

			    NNREAD=NNREAD+1; "counts total number of reads"

			    IF(DOSE_STAT=1 | NHSTRY=0) NHSTRY=NHSTRY+1;

			    "We only use the forward particles."
			    IF( NPASS=1 ) GOTO :READ-PARTICLE-SRC22:;"...signals not 1st time crossing"

			    IF(IQIN=-1 & IQINC~=-1 & IQINC~=2 & IQINC ~=3)[GOTO :READ-PARTICLE-SRC22:;]
			    IF(IQIN= 1 & IQINC~= 1 & IQINC~=2 & IQINC ~=3)[GOTO :READ-PARTICLE-SRC22:;]
			    IF(IQIN= 0 & IQINC~= 0 & IQINC~=2)[GOTO :READ-PARTICLE-SRC22:;]

			    "IF(IQIN=-1 & IQINC~=-1 & IQINC~=2)[ GOTO :READ-PARTICLE-SRC22:;]"
			    "ELSEIF(IQIN=1 & IQINC~=1 & IQINC~=2)[ GOTO :READ-PARTICLE-SRC22:;]"
			    "ELSEIF(IQINC~=0 & IQINC~=2)[GOTO :READ-PARTICLE-SRC22:;]"

			    "apply the offset"
			    xin=xin+xoffset;
			    yin=yin+yoffset;

			    "Rotate the position"
			    zin = -dist_phsp*cost_phsp + yin*sint_phsp + chamber_c;
			    yin = dist_phsp*sint_phsp + yin*cost_phsp;

			    "Rotate the direction"
			    vtemp = vin;
			    vin = -win*sint_phsp + vin*cost_phsp;
			    win =  win*cost_phsp + vtemp*sint_phsp;

			    "Check if particle is within the geometry or will enter it"
			    radp = xin*xin + yin*yin;
			    IF( radp > CYRAD2(nr) | zin < zplane(1) | zin > zplane(nz+1) ) [
			        " outside of geometry. check whether that particle will get in "
			        check = .true.;
			        IF( zin < zplane(1) & win > 0 ) [ tf = (zplane(1)-zin)/win; ]
			        ELSE IF( zin > zplane(nz+1) & win < 0 ) [
			                      tf = (zplane(nz+1) - zin)/win;
			        ]
			        ELSE [ check = .false.; ]
			        IF( check ) [
			            xtmp = xin + uin*tf; ytmp = yin + vin*tf;
			            radp = xtmp*xtmp + ytmp*ytmp;
			            IF( radp <= CYRAD2(nr) ) [
			                xin = xtmp; yin = ytmp;
			                IF( zin < zplane(1) ) [ zin = zplane(1); iz = 1; ]
			                ELSE                  [ zin = zplane(nz+1); iz = nz; ]
			                DO ix=1,nr [ IF( radp <= CYRAD2(ix) ) EXIT; ]
			                goto :FOUND-ENTRY-22:;
			            ]
			            radp = xin*xin + yin*yin;
			        ]
			        IF( radp < CYRAD2(nr) ) goto :READ-PARTICLE-SRC22:;
			        phbb = uin*uin + vin*vin;
			        IF( phbb < 1e-15 ) goto :READ-PARTICLE-SRC22:;
			        pha = (xin*uin + yin*vin)/phbb;
			        IF( pha > 0 ) goto :READ-PARTICLE-SRC22:;
			        phb = (radp - CYRAD2(nr))/phbb;
			        phd = pha*pha - phb;
			        IF( phd < 0 ) goto :READ-PARTICLE-SRC22:;
			        tf = -pha-sqrt(phd);
			        IF( tf < 0 ) goto :READ-PARTICLE-SRC22:;
			        zin = zin + tf*win;
			        IF( zin < zplane(1) | zin > zplane(nz+1) ) goto :READ-PARTICLE-SRC22:;
			        ix = nr; xin = xin + uin*tf; yin = yin + vin*tf;
			        DO iz=2,nplane [ IF( zin <= zplane(iz) ) EXIT; ]
			        iz = iz-1;
			        goto :FOUND-ENTRY-22:;
			    ]
			    " already inside, find the region index "
			    DO iz=2,nplane [ IF( zin <= zplane(iz) ) EXIT; ]
			    iz = iz-1;
			    DO ix=1,nr [ IF( radp <= CYRAD2(ix) ) EXIT; ]

			  :FOUND-ENTRY-22:

			    count_phsp = count_phsp + 1;
			    irin = 1 + (ix-1)*nz + iz;
			    NRCFLG = 20; "GEOMETRICAL FLAG, SOURCE INCIDENT FROM SIDE(?)"
			    IF(NRCYCL>0)["store non-global parameters for recycling"
			       CYCLNUM=1;
			       XINOLD=XIN;
			       YINOLD=YIN;
			       ZINOLD=ZIN;
			       UINOLD=UIN;
			       VINOLD=VIN;
			       WINOLD=WIN;
			       IRINOLD=IRIN;
			       NRCFLGOLD=NRCFLG;
			       WEIGHTOLD=WEIGHT;
			    ]
			 ]
			]
	*/
			return;
		    }
		    else if(EGS4.iGeom==EGS4.iCavitySPH)//2=SPH Geom
		    {
				int IC=0;//int IX=0;

				//"CALCULATE THE HISTORY DEPENDENT CONSTANTS THAT MAY CHANGE DEPENDING ON"
				//"THE SOURCE TYPE"

				if(ISOURC==4)
				{//"not need to sample position and direction"
				   xin    = 0.0;yin = 0.0;zin = 0.0;//"particle at origin"
				   uin    = 0.0;vin = 0.0;win = 1.0;//"moving along z-axis"
				   irin   = 2; //"initial region always inner most sphere"
				   WEIGHT = 1.0;//"incident weight"
				   return;
				}

				if(IDSTON==0)
				{// "NO RADIAL DISTRIBUTION"
				    while(true)
				    {
				        xin=EGS4.random01(); xin =(2.0*xin -1.0)*RBEAM; //"CHOOSE X"
				        YINP=EGS4.random01();YINP=(2.0*YINP-1.0)*RBEAM; //"CHOOSE Y (UNROTATED)"
				        RHO2=xin*xin+YINP*YINP;

				        if(RHO2<=RBEAM2)
				        	break;
				    }// UNTIL RHO2.LE.RBEAM2;
				}
				else
				{// "RADIAL DISTRIBUTION"
				    while(true)
				    {
				        //"SAMPLE THE INCIDENT RADIAL DISTRIBUTION"
				        double RNNO1=EGS4.random01();double RNNO2=EGS4.random01();
				        //"IN NEXT STATEMENT MIN AVOIDS K=$MXRDIST+1, ALMOST NEVER IMPLEMENTED"
				        Double dbl1=new Double(dble($MXRDIST)*RNNO1+1.);
				        int dbl1i=dbl1.intValue();
				        Double dbl2=new Double(dble($MXRDIST));
				        int dbl2i=dbl2.intValue();
				        int K=Math.min(dbl1i,dbl2i);
				        //RRHO=Math.sqrt(RCDFIN(K,1)+RNNO2*RCDFIN(K,2));
				        RRHO=Math.sqrt(RCDFIN[K-1][0]+RNNO2*RCDFIN[K-1][1]);

				        if(RRHO<=RBEAM)
				        	break;
				     }// UNTIL RRHO.LE.RBEAM;

				    xin=EGS4.random01();xin=(2.0*xin-1.0)*RRHO;// "CHOOSE X"
				    RHO2=RRHO*RRHO;
				    YINP=RHO2-xin*xin;
				    if(YINP<=0.0){YINP=0.0;}else{YINP=Math.sqrt(YINP);}// "FIND Y (UNROTATED)"
				}

				ZINP=EGS4Geom.RSPH2[EGS4Geom.NR-1]-RHO2;
				if(ZINP<=0.0){ZINP=0.0;}else{ZINP=-Math.sqrt(ZINP);}// "FIND Z (UNROTATED)"
				ZINP=ZINP*0.999995; //"ALLOW A LITTLE PENETRATION"

				//"ROTATION YINP, ZINP ABOUT X-AXIS TO ACCOUNT FOR THE BEAM INCIDENCE"
				zin=ZINP*WINC-YINP*VINC;
				yin=YINP*WINC+ZINP*VINC;

				//"FIND THE STARTING REGION'S CONE NUMBER"
				//"NOTE THAT ON THE SUREFACE, IRIN=NREG+IC(CONE NUMBER)-NC"
				if(EGS4Geom.NC==1)
				{//"ONLY ONE CONE"
				    irin=EGS4Geom.nreg;//NREG;
				}
				else if(zin>0.0)
				{// "HITS THE FORWARD HEMISPHERE"
				    if(EGS4Geom.NC==2)
				    {// "UNIQUE SOLUTION IN THIS CASE"
				        irin=EGS4Geom.nreg-1;//NREG-1;
					}
				    else
				    {// "MUST SEARCH FOR ENTRANCE REGION"
				        TNTST2=RHO2/(zin*zin);//ZIN**2;
				        IC=1; //"IF SEARCH FAILS, THE PARTICLE STARTS IN THE FIRST CONE"
				        for(int ICTST=1;ICTST<=EGS4Geom.NPLAN1;ICTST++)
				        {
				            if(TNTST2>=EGS4Geom.TANAL2[ICTST])//0 biased
				            {
								IC=IC+1;
								break;
							}
						}
				        IC=Math.min(EGS4Geom.NPLAN1,IC); //"TO KEEP PARTICLE IN THE FORWARD HEMISPHERE"
				        irin=EGS4Geom.nreg+IC-EGS4Geom.NC;
					}
				}
				else if(zin<0.0)
				{// "HITS THE BACKWARD HEMISPHERE"
				    if(EGS4Geom.NC==2)
				    {// "UNIQUE SOLUTION IN THIS CASE"
				        irin=EGS4Geom.nreg;//NREG;
					}
				    else
				    {// "MUST SEARCH FOR ENTRANCE REGION"
				        TNTST2=RHO2/(zin*zin);//ZIN**2;
				        IC=EGS4Geom.NC;// "IF SEARCH FAILS, THE PARTICLE STARTS IN THE LAST CONE"
				        for(int ICTST=EGS4Geom.NPLAN2;ICTST<=EGS4Geom.NC;ICTST++)
				        {
				            if(TNTST2>=EGS4Geom.TANAL2[ICTST-1])
				            {
								IC=IC-1;
								break;
							}
						}
				        IC=Math.max(EGS4Geom.NPLAN2,IC); //"TO KEEP PARTICLE IN THE BACKWARD HEMISPHERE"
				        irin=EGS4Geom.nreg+IC-EGS4Geom.NC;
					}
				}
				else
				{// "ZIN=0.0, RIGHT ON THE EQUATOR - IMPROBABLE CASE"
				    if    (WINC>0.0)
				    {//"GOING FORWARD"
				    	irin=EGS4Geom.nreg+EGS4Geom.NPLAN1-EGS4Geom.NC;
					}
				    else if(WINC<0.0)
				    {// "GOING BACKWARD"
				    	irin=EGS4Geom.nreg+EGS4Geom.NPLAN2-EGS4Geom.NC;
					}
				    else
				    {// "PARALLEL TO EQUATORIAL PLANE"
				     //   "SHARE EQUALLY BETWEEN FRONT AND BACK HEMISPHERES"
				        double SPLIT=EGS4.random01();
				        if(SPLIT<0.5)
				        {
							irin=EGS4Geom.nreg+EGS4Geom.NPLAN1-EGS4Geom.NC;
						}
						else
						{
							irin=EGS4Geom.nreg+EGS4Geom.NPLAN2-EGS4Geom.NC;
						}
					}
				}

				//"FIND INITIAL DIRECTION COSINES"
				if(ISRCTY==0)
				{// "PARALLEL BEAM"
				    uin=UINC;vin=VINC;win=WINC;
				}
				else
				{// "POINT SOURCE"
				    //"FOR D AND WEIGHT USE UNROTATED FORMS - THEY ARE INVARIANT"
				    D=Math.sqrt((DISTR+ZINP)*(DISTR+ZINP)+RHO2);
				    WEIGHT=ABEAM*(DISTR+EGS4Geom.RSPH2[EGS4Geom.NR-1]/ZINP)/(D*D*D);
				    uin=xin/D;vin=(DISTR*VINC+yin)/D;win=(DISTR*WINC+zin)/D;
				    //OMEGIS[IS-1]=OMEGIS(IS)+WEIGHT;
				    OMEGIS[IS-1]=OMEGIS[IS-1]+WEIGHT;
				}

				return;


			}
		}

		/**
		 * Not implemented.
		 */
		public static void SRCEND()
		{
			/*
			if( ISOURC == 23 )
			{
			    finish_beamsource();
			}
			*/
			return;
		}


		//subroutine src15;
		//implicit none;
		//;COMIN/GEOM,IODAT2,RANDOM,SOURCE/;
		//$REAL       temp1,temp2,temp3,temp4; "from input file"
		//$REAL       XIN,YIN,ZIN,UIN,VIN,WIN,WEIGHT;
		//$INTEGER    NRCFLG,IRIN,iout,error_flag;

		//real*8      theta,
		//real*8      fak,x1,y1,x2,y2,x3,y3,,xmax,,,ymax,;
		//logical     ,enter_fb,enter_side;
		//real*8      eta1,eta2,xx,yy,zz,dist2,dist,uli,vli,wli,
		//            xli,yli,zli,rad2,ts,tf,ur2,ux,ux2,rr2,del,error;
		//integer*4   ix,iz;

		//save        cost,sint,d,d2,yo,zo,dz,R,R2,xmin,delx,ymin,dely,
		//            zc,area,w0,count,sumw,sumw2,just_fb,just_side,pi,ro2;

		/**
		 * Internal use by SRCINI.
		 */
		public static void src15_ini()//(temp1,temp2,temp3,temp4,error_flag);
		{
		//"=================================================="
			double fak=0.0;double x1=0.0;double y1=0.0;double x2=0.0;double y2=0.0;
			double x3=0.0;double y3=0.0;double xmax=0.0;double ymax=0.0;
			pi = Math.PI;//4*datan(1d0);
			double theta = TEMP2;
			if( theta > 180 ) { theta = 360 - theta; }
			theta = theta*pi/180;
			cost = Math.cos(theta); sint = Math.sin(theta);
			d = TEMP1; d2 = d*d;
			yo = d*sint; zo = -d*cost;
			//"write(35,*) ' Source position: ',yo,zo;"
			R = EGS4Geom.RCYL[EGS4Geom.NR];
			R2 = R*R;
			dz = 0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
			zc = 0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			if( yo > -R && yo < R && zo > -dz && zo < dz )
			{// "inside geometry"
				EGS4.STOPPROGRAM=true;
				EGS4.seqStr="ERROR: source inside geometry!!";
				eq.printSequence(EGS4.seqStr);
				return;
		    	//error_flag = 1;
			}
			just_fb = false; just_side = false;
			if( yo >= -R && yo <= R ) { just_fb = true; }
			if( zo >= -dz && zo <= dz ) { just_side = true; }

			if( just_fb )
			{
		  		if( zo < -dz )
		  		{
					fak = 1./(1-(-R*sint + dz*cost)/d);
		  		  	x1 = R*fak; y1 = (-R*cost-dz*sint)*fak;
		  			fak = 1./(1-(R*sint + dz*cost)/d);
		  		  	x2 = R*fak; y2 = (R*cost-dz*sint)*fak;
				}
		  		else
		  		{
		    		fak = 1./(1-(-R*sint - dz*cost)/d);
		    		x1 = R*fak; y1 = (-R*cost+dz*sint)*fak;
		    		fak = 1./(1-(R*sint - dz*cost)/d);
		    		x2 = R*fak; y2 = (R*cost+dz*sint)*fak;
				}
		  		ymin = Math.min(y1,y2); ymax = Math.max(y1,y2);
		  		xmax = Math.max(x1,x2); xmin = -xmax;
			}
			else if( just_side )
			{
		  		fak = 1./(1-(R*sint + dz*cost)/d);
		  		x1 = R*fak; y1 = (R*cost-dz*sint)*fak;
		  		fak = 1./(1-(R*sint - dz*cost)/d);
		  		x2 = R*fak; y2 = (R*cost+dz*sint)*fak;
		  		ymin = Math.min(y1,y2); ymax = Math.max(y1,y2);
		  		xmax = Math.max(x1,x2); xmin = -xmax;
			}
			else
			{
		  		if( zo < -dz )
		  		{
		    		fak = 1./(1-(-R*sint + dz*cost)/d);
		    		x1 = R*fak; y1 = (-R*cost-dz*sint)*fak;
				}
		  		else
		  		{
		    		fak = 1./(1-(-R*sint - dz*cost)/d);
		    		x1 = R*fak; y1 = (-R*cost+dz*sint)*fak;
				}
		  		fak = 1./(1-(R*sint + dz*cost)/d);
		  		x2 = R*fak; y2 = (R*cost-dz*sint)*fak;
		  		fak = 1./(1-(R*sint - dz*cost)/d);
		  		x3 = R*fak; y3 = (R*cost+dz*sint)*fak;
		  		xmax = EGS4.max(x1,x2,x3); xmin = -xmax;
		  		ymin = EGS4.min(y1,y2,y3); ymax = EGS4.max(y1,y2,y3);
			}
			delx = xmax - xmin; dely = ymax - ymin; area = delx*dely;
			w0 = area/d2; ro2 = yo*yo;

			EGS4.seqStr=" Electric charge of the source:"+EGS4.format(iqin,12);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr=" Point source off axis";
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr=" Distance of source from centre of chamber: "+
			EGS4.format(d,10,true)+" cm";
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr=" Angle to z-axis (degrees): "+
			EGS4.format(TEMP2,10,true);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr=" Rectangle in the plane perpendicular to the source-chamber "+
			" \n"+" axis seen from the source: "+EGS4.format(xmin,10,true)+
			", "+EGS4.format(xmax,10,true)+", "+EGS4.format(ymin,10,true)+", "+
			EGS4.format(ymax,10,true);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr="Source position: y = "+EGS4.format(yo,8,true)+
			" , z = "+EGS4.format(zo,8,true);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

			return;
	    }

		/**
		 * Internal use by SRCOTO.
		 */
		public static void  src15_oto()//(iout);
		{
			AINFLU = NCASET/d2;
			return;
		}

		/**
		 * Internal use by SRCOUT.
		 */
		public static void src15_describe()
		{
					EGS4.seqStr=" POINT SOURCE OFF AXIS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" distance to chamber centre: "+
					EGS4.format(d,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" angle: "+
					EGS4.format(theta,8,true)+" degrees ";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" source position: y = "+
					EGS4.format(yo,8,true)+" z = "+EGS4.format(zo,8,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					return;
	    }

		/**
		 * Internally used by SRCHST
		 */
		public static void src15_hst()//(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
		{
			double xli=0.0;double yli=0.0;double ur2=0.0;double ux=0.0;double ux2=0.0;double rr2=0.0;
			double del=0.0;double zli=0.0;int iz=0;double rad2=0.0;int ix=0;
			count = count + 1;
			double eta1=EGS4.random01(); double eta2=EGS4.random01();
			double xx = xmin + delx*eta1;
			double yy = ymin + dely*eta2;
			//"write(35,*) ' xx = ',xx,' yy = ',yy;"
			double dist2 = d2 + xx*xx + yy*yy;
			double dist = Math.sqrt(dist2);
			w0 = d/dist;
			double zz = yy*sint; yy = yy*cost;
			double uli = xx/dist;
			double vli = (yy-yo)/dist;
			double wli = (zz-zo)/dist;
			//"write(35,*) ' direction ',uli,vli,wli;"
			double tf = 1.e15; double ts = 1.e15;
			boolean enter_fb = false;boolean enter_side = false;

			if( !just_side )
			{
		  		if( zo < -dz ) { tf = (-dz-zo)/wli; } else { tf = (dz - zo)/wli; }
		  		xli = uli*tf; yli = vli*tf + yo;
		  		if( xli*xli + yli*yli < R2 ) { enter_fb = true; }
		  		else { tf = 1.e15; }
			}
			if( !just_fb )
			{
		  		ur2 = uli*uli + vli*vli; ux = -yo*vli/ur2; ux2 = ux*ux; rr2 = ro2/ur2;
		  		del = ux2 - rr2 + R2/ur2;
		  		if( del >= 0 )
		  		{
		    		ts = ux - Math.sqrt(del);
		    		if( ts > 0 )
		    		{
		      			zli = zo + wli*ts;
		      			if( zli >= -dz && zli <= dz ) { enter_side = true; }
					}
		    		if( !enter_side ) { ts = 1.e15; }
				}
			}
			//"write(35,*) ' distances: ',tf,ts;"

			if( !enter_fb && !enter_side )
			{
		    	WEIGHT = 0.0;
		    	xin = 0.0; yin = 0.0; zin = EGS4Geom.ZPLANE[0];//zplane(1);
		    	uin = 0.0; vin = 0.0; win = 1.0; irin = 2;
		    	return;
			}
			uin = uli; vin = vli; win = wli; WEIGHT = w0*area/dist2;
			sumw = sumw + WEIGHT; sumw2 = sumw2 + WEIGHT*WEIGHT;
			if( tf < ts )
			{
		  		xin = xli; yin = yli;
		  		if( zo < -dz )
		  		{
					zin = EGS4Geom.ZPLANE[0];//zplane(1);
					iz = 1;
				}
		  		else
		  		{
					zin = EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1];//zplane(nplane);
					iz = EGS4Geom.NPLANE-1;
				}
		  		rad2 = xli*xli+yli*yli;
		  		for(ix=1;ix<=EGS4Geom.NR;ix++)
		  		{
					if( rad2 <= EGS4Geom.CYRAD2[ix-1] )//CYRAD2(IX) )
						break;
				}
			}
			else
			{
		  		zin = zli + zc;
		  		xin = ts*uli; yin = yo + vli*ts; ix = EGS4Geom.NR;//nr;
		  		for(iz=1;iz<=EGS4Geom.NZ;iz++)
		  		{
					if( zin >= EGS4Geom.ZPLANE[iz-1] &&
						zin < EGS4Geom.ZPLANE[iz] )//zin < zplane(iz+1) )
						break;
				}
			}

			irin = 1 + (ix-1)*EGS4Geom.NZ + iz;
			//"write(35,*) ' entering at ',xin,yin,zin,irin;"

			return;
	    }

		/**
		 * Computes the solid angle for source whose code is 15.
		 */
		public static void src15_out()//(iout);
		{
			sumw = sumw/count; sumw2 = sumw2/count;
			double error = (sumw2 - sumw*sumw)/(count-1);
			if( error > 0 ) { error = Math.sqrt(error); }

			EGS4.seqStr="  Source 15: solid angle for detector: "+EGS4.format(sumw,14)+
			" +/- "+EGS4.format(sumw2,14);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

			//write(iout,600) sumw,sumw2;
			//600 format(//'  Source 15: solid angle for detector: ',g14.5,' +/- ',g14.5,//);

			return;

		}//end;

		//"======================================================================"
		// subroutine src16;
		//"======================================================================"
		//implicit none;
		//;COMIN/GEOM,IODAT2,RANDOM,SOURCE/;
		//$REAL       temp1,temp2,temp3,temp4; "from input file"
		//$REAL       XIN,YIN,ZIN,UIN,VIN,WIN,WEIGHT;
		//$INTEGER    NRCFLG,IRIN,iout,error_flag;

		//real*8      ,u,aux,u2,
		//            ,dist,dist2,,
		//            rad2,eta,,,zz;
		//logical     ,do_both,
		//integer*4   ix,iz;

		//save        pi,cost,sint,d0,d02,R,R2,dz,dr,zc,theta,
		//            xo,yo,zo,a_side,a_fb,atot,wfb,rs,delxs,delys,rs2,
		//            sumw,sumw2,count,point_source,do_fb,do_side,is_circle;

		/**
		 * Internally used by SRCINI
		 */
		public static void src16_ini()//(temp1,temp2,temp3,temp4,error_flag);
		{
		//"=================================================="
			double u=0.0;//double u2=0.0;
			double aux=0.0;//double dist=0.0;double dist2=0.0;double rad2=0.0;
			//double eta=0.0;double zz=0.0;int ix=0;int iz=0;

			pi =Math.PI;// 4*datan(1d0);
			count = 0;
			theta = TEMP2; theta = theta*pi/180;
			cost = Math.cos(theta); sint = Math.sin(theta);
			d0 = TEMP1; d02 = d0*d0;
			R = EGS4Geom.RCYL[EGS4Geom.NR];
			R2 = R*R;
			dz = 0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]);
			zc = 0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]+EGS4Geom.ZPLANE[0]);
			if( TEMP3 <= 0 && TEMP4 <= 0 )
			{
		    	//"We assume a point source (i.e. the same as source 12 and 15)"
		    	point_source = true; zo = -d0*cost; yo = d0*sint; xo = 0;
		    	do_fb = false; do_side = false; do_both = false;
		    	a_side = 0.0; a_fb = 0.0;
		    	if( Math.abs(yo) > R )
		    	{
		        	u = R/yo; aux = Math.sqrt(1-u*u); dr = R*aux; do_side = true;
		        	a_side = 4*dz*R*yo*(aux-u*Math.atan(aux/u));
				}
		    	if( Math.abs(zo) > dz )
		    	{
		        	a_fb = pi*R2*(Math.abs(zo)-dz); do_fb = true;
				}
		    	atot = a_fb + a_side;  wfb = a_fb/atot;
		    	if( do_side && do_fb ) { do_both = true; }
			}
			else
			{
		    	point_source = false; is_circle = false;
		    	if( TEMP4 <= 0 || TEMP3 <= 0 )
		    	{
		        	//"We assume that the source is a circle with radius temp3"
		        	if( TEMP4 <= 0 ) { rs = TEMP3; } else { rs = TEMP4; }
		        	rs = TEMP3; rs2 = rs*rs; delxs = rs; delys = rs;
		        	is_circle = true;
				}
		    	else
		    	{
		        	//"We assume that the source is a rectangle with half-widths"
		        	//"temp3 and temp4 (x and y in the plane perpendicular to the
		        	//"source center to chamber center axis)"
		        	delxs = TEMP3; delys = TEMP4; rs2 = 1.1*(delxs*delxs+delys*delys);
				}
			}

			EGS4.seqStr=" Electric charge of the source:"+EGS4.format(iqin,12);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

			if( point_source )
			{
				EGS4.seqStr=" Point source off axis (source 16 implementation)";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Distance of source from centre of chamber: "+
				EGS4.format(d0,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Angle to z-axis (degrees): "+
				EGS4.format(TEMP2,10,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr=" Estimated probability to strike front/back face: "+
				EGS4.format(wfb,12,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				double dbl=1-wfb;
				EGS4.seqStr=" Estimated probability to strike side face: "+
				EGS4.format(dbl,12,true);
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
			}
			else
			{
		    	if( is_circle )
		    	{
					EGS4.seqStr=" Disk irradiating the chamber (source 16 implementation)";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Distance of source centre to chamber centre: "+
					EGS4.format(d0,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Angle to z-axis (degrees): "+
					EGS4.format(TEMP2,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Radius of the source: "+
					EGS4.format(rs,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
		    	else
		    	{
					EGS4.seqStr=" Rectangle irradiating the chamber (source 16 implementation)";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Distance of source centre to chamber centre: "+
					EGS4.format(d0,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Angle to z-axis (degrees): "+
					EGS4.format(TEMP2,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Rectangle half-sizes: "+
					EGS4.format(delxs,10,true)+" , "+EGS4.format(delys,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			}

			return;
		}

		/**
		 * Internally used by SRCOTO
		 */
		public static void src16_oto()//(iout);
		{
			//"===================="
			AINFLU = NCASET/d02;
			return;
		}

		/**
		 * Internally used by SRCOUT
		 */
		public static void src16_describe()
		{
		//"========================="
			if( point_source )
			{
					EGS4.seqStr=" POINT SOURCE OFF AXIS (source 16)";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
		    }
			else
			{
				if( is_circle )
				{
					EGS4.seqStr=" Radiating circle off axis (source 16)";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Source radius: "+
					EGS4.format(rs,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
				else
				{
					EGS4.seqStr=" Radiating rectangle off axis (source 16)";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Source half-sizes: "+
					EGS4.format(delxs,10,true)+EGS4.format(delys,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
				}
			}

			double dbl=theta*180/Math.PI;
					EGS4.seqStr=" Distance between source and chamber centres: "+
					EGS4.format(d0,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Angle to z-axis (degrees): "+
					EGS4.format(dbl,10,true)+" degrees";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					return;
	    }

		/**
		 * Internally used by SRCHST
		 */
		public static void src16_hst()//(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
		{
			double u=0.0;double u2=0.0;
			double aux=0.0;double dist=0.0;double dist2=0.0;double rad2=0.0;
			double eta=0.0;double zz=0.0;int ix=0;int iz=0;

			//"=========================================================="
			count = count + 1;
			if( !point_source )
			{
		    	while(true)
		    	{
		        	xo=EGS4.random01(); xo = delxs*(2*xo-1);
		        	yo=EGS4.random01(); yo = delys*(2*yo-1);

		        	if(xo*xo + yo*yo < rs2)
		        		break;
		    	}// UNTIL ( xo*xo + yo*yo < rs2 );
		    	zo = yo*sint - d0*cost; yo = d0*sint + yo*cost;
		    	yo = Math.sqrt(yo*yo+xo*xo);
		    	do_fb = false; do_side = false; do_both = false;
		    	a_side = 0; a_fb = 0;
		    	if( Math.abs(yo) > R )
		    	{
		        	u = R/yo; do_side = true;
		        	if( u < 0.1 )
		        	{
		            	u2 = u*u; aux = 1 - u2/2 - u2*u2/8;
		            	a_side = 4*dz*R*yo*(aux-pi*u/2+u2+u2*u2/6);
					}
		        	else
		        	{
		            	aux = Math.sqrt(1-u*u);
		            	a_side = 4*dz*R*yo*(aux-u*Math.atan(aux/u));
					}
		        	dr = R*aux;
				}
		    	if( Math.abs(zo) > dz )
		    	{
		        	a_fb = pi*R2*(Math.abs(zo)-dz); do_fb = true;
				}
		    	atot = a_fb + a_side;
		    	if( do_side && do_fb ) { wfb = a_fb/atot; do_both = true; }
			}

			if( do_both )
			{
			    eta=EGS4.random01();
			    if( eta < wfb ) { do_fb = true; do_side = false; }
			    else            { do_fb = false; do_side = true; }
			}
			if( do_fb )
			{
			    while(true)
			    {
			        xin=EGS4.random01(); xin = R*(2*xin-1);
			        yin=EGS4.random01(); yin = R*(2*yin-1);
			        rad2 = xin*xin + yin*yin;

			        if(rad2 <= R2) break;
			    }// UNTIL (rad2 <= R2 );
			    for( ix=1;ix<=EGS4Geom.NR;ix++)
			    {
					if( rad2 <= EGS4Geom.CYRAD2[ix-1] )
						break;
				}
			    if ( zo < -dz ) { zin = -dz; iz = 1; }
			    else            { zin =  dz; iz = EGS4Geom.NZ; }
			}
			else
			{
			    zin=EGS4.random01(); zin = dz*(2*zin-1);
			    while(true)
			    {
			        xin=EGS4.random01();eta=EGS4.random01();
			        xin = dr*(2*xin-1);
			        yin = Math.sqrt(R2 - xin*xin);

			        if(eta < (yo-R2/yin)/(yo-R)) break;
			    }// UNTIL ( eta < (yo-R2/yin)/(yo-R) );
			    ix = EGS4Geom.NR;
			    zz = zin + zc;
			    for( iz=1;iz<=EGS4Geom.NZ;iz++)
			    {
					if( zz >= EGS4Geom.ZPLANE[iz-1] & zz < EGS4Geom.ZPLANE[iz] )
						break;
				}
			}
			uin = xin; vin = yin - yo; win = zin - zo;
			dist2 = uin*uin + vin*vin + win*win;
			dist = Math.sqrt(dist2); uin = uin/dist; vin = vin/dist; win = win/dist;
			WEIGHT = atot/dist2/dist;
			sumw = sumw + WEIGHT; sumw2 = sumw2 + WEIGHT*WEIGHT;
			irin = 1 + (ix-1)*EGS4Geom.NZ + iz;
			zin = zin + zc;
			return;
		}


	/*
		$HAVE_LOAD_DSO(#);
	#ifndef HAVE_LOAD_DSO;
	//==========NOT HERE!!!!!!!!!!!!!!!!!
		subroutine init_beamsource(i_parallel,conf_name,
		                         hen_house,egs_home,the_beam_code,
		                         the_pegs_file,the_input_file);
		$INTEGER i_parallel;
		character*(*) conf_name;
		character*(*) hen_house,egs_home,the_beam_code,the_pegs_file,the_input_file;
		write(6,*) 'You need a working C compiler to use source 23!';
		$CALL_EXIT(1);
		end;

		subroutine sample_beamsource(ein,xin,yin,zin,uin,vin,win,weight,iqin,latchi,
		                             nhstry,iphatin);
		return; end;
	//==========NOT HERE!!!!!!!!!!!!!!!!!
		subroutine finish_beamsource;
		return; end;

		#endif;
	//==========NOT HERE!!!!!!!!!!!!!!!!!
	*/

		//@@@@@@@@@@@@@@@@@@SRCSPH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	//"*************************************************************************
	//"
	//"                       ******************
	//"                       *                *
	//"                       * srcsph.mortran *
	//"                       *                *
	//"                       ******************
	//"
	//"
	//"       This subroutine is dedicated exclusively to calculations related
	//"       to the source configuration. It is accessed from four locations:
	//"               1) from the input subroutine. Control is
	//"                  transferred to this routine to read in source specific
	//"                  information. Control is transferred back to the input
	//"                  subroutine following these inputs
	//"               2) From the 'one time only constants' section of the main
	//"                  routine to scale the geometry of the incident beam
	//"                  and calculate the fluence
	//"               3) From the histories loop of the main routine. A call is
	//"                  made prior to each history to set the point of entry
	//"                  into the target, the initial direction cosines, the
	//"                  initial energy and the statistical weight of the
	//"                  incident particle.
	//"               4) From the input summary routine to print source info
	//"
	//"       Fluence is defined at the point on axis at the front face
	//"       of the geometry.
	//****************************************************************************
	//"The following is the $CALL-HOWNEAR macro for PRESTA-II
	//"This routine is intended to be used to calculate the minimum perpendiclar     "
	//"to the nearest bounding surface. This version is specially designed for the   "
	//"NRCC spherical geometry package. A different version is needed for other      "
	//"geometry packages.                                                            "

	//REPLACE {$CALL-HOWNEAR(#);} WITH
	//{
//	    call hownear({P1},x(np),y(np),z(np),ir(np));
	//}
		/**
		 * This routine is intended to be used to calculate the minimum perpendicular to the nearest bounding surface. 
		 * This version is specially designed for the spherical geometry package.
		 */
		public static void HOWNEAR()
		{
			//tperp, "nearest distance to any boundary (output)

			double xx=EGS4.X[EGS4.NP-1];//     "x-position of the particle (input)
			double yy=EGS4.Y[EGS4.NP-1];//     "y-position of the particle (input)
			double zz=EGS4.Z[EGS4.NP-1];//     "z-position of the particle (input)
			int irl=EGS4.IR[EGS4.NP-1];

			int IX=0;int IC=0;double RHOL=0.0;double RHOL2=0.0;double RL=0.0;

	    	//$GET-IX-IC(IRL);
	    	IX=EGS4Geom.GET_IX(irl);
	    	IC=EGS4Geom.GET_IZC(irl);

	    	RHOL2=xx*xx+yy*yy;
	    	RL=Math.sqrt(RHOL2+zz*zz);
	    	EGS4.tperp=EGS4Geom.RSPH[IX]-RL;//0 biased
	    	if(IX!=1)
	    	{
	    	    EGS4.tperp=Math.min(EGS4.tperp,RL-EGS4Geom.RSPH[IX-1]);
			}
	    	if(EGS4Geom.NC!=1)
	    	{
	    	    RHOL=Math.sqrt(RHOL2);
	    	    if(IC!=EGS4Geom.NC)
	    	    {//0 biased
	    	    	EGS4.tperp=Math.min(EGS4.tperp,
	    	    	Math.abs(EGS4Geom.cosalp[IC]*(RHOL-EGS4Geom.TANALP[IC]*zz)));
				}
	    	    if(IC!=1)
	    	    {
	    	        EGS4.tperp=Math.min(EGS4.tperp,
	    	        Math.abs(EGS4Geom.cosalp[IC-1]*(RHOL-EGS4Geom.TANALP[IC-1]*zz)));
	    	        //REAL(Math.abs(EGS4Geom.cosalp[IC-1]*(RHOL-EGS4Geom.TANALP[IC-1]*ZZ))));
				}
			}

			return;
		}

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
	//"  SOURCE NUMBER                 (I)   number of the source
	//"                                      [ISOURC]
	//"
	//;
	//"------------------------------------------------------------------------------
	//"
	//" SOURCE 0:    FOR PARALLEL BEAM FROM ANY ANGLE
	//"              RBEAM,UINC,VINC,WINC
	//"               RBEAM   RADIUS OF THE BEAM AT THE FRONT OF THE TARGET IN CM
	//"                               DEFAULTS TO MAX RADIUS
	//"               UINC    INCIDENT X-AXIS DIRECTION COSINE
	//"               VINC    INCIDENT Y-AXIS DIRECTION COSINE
	//"               WINC    INCIDENT Z-AXIS DIRECTION COSINE
	//"                       NOTE: (UINC,VINC,WINC) GET AUTOMATICALLY NORMALIZED
	//"                             DEFAULTS TO (0.0,0.0,1.0)
	//;
	//"------------------------------------------------------------------------------
	//"
	//" SOURCE 1:   FOR POINT SOURCE INCIDENT FROM ANY ANGLE
	//"             DISTR,RBEAM,UINC,VINC,WINC
	//"               DISTR   DISTANCE OF THE SOURCE FROM THE MIDDLE OF THE TARGET
	//"                       IN CM (DEFAULTS TO 100.)
	//"               RBEAM   RADIUS OF THE BEAM AT THE FRONT OF THE TARGET IN CM
	//"                               DEFAULTS TO MAX RADIUS
	//"               UINC    INCIDENT X-AXIS DIRECTION COSINE
	//"               VINC    INCIDENT Y-AXIS DIRECTION COSINE
	//"               WINC    INCIDENT Z-AXIS DIRECTION COSINE
	//"                       NOTE: (UINC,VINC,WINC) GET AUTOMATICALLY NORMALIZED
	//"                             DEFAULTS TO (0.0,0.0,1.0)
	//;
	//"------------------------------------------------------------------------------
	//"
	//" SOURCE 10:  FOR PARALLEL BEAM FROM ANY ANGLE WITH RADIAL DISTRIBUTION
	//"             UINC,VINC,WINC
	//"               UINC    INCIDENT X-AXIS DIRECTION COSINE
	//"               VINC    INCIDENT Y-AXIS DIRECTION COSINE
	//"               WINC    INCIDENT Z-AXIS DIRECTION COSINE
	//"                       NOTE: (UINC,VINC,WINC) GET AUTOMATICALLY NORMALIZED
	//"                             DEFAULTS TO (0.0,0.0,1.0)
	//;
	//"------------------------------------------------------------------------------
	//"
	//" SOURCE 11:  FOR POINT SOURCE INCIDENT FROM ANY ANGLE WITH RADIAL DISTRIBUTION
	//"             DISTR,RBEAM,UINC,VINC,WINC
	//"               DISTR   DISTANCE OF THE SOURCE FROM THE MIDDLE OF THE TARGET
	//"                       IN CM (DEFAULTS TO 100.)
	//"               UINC    INCIDENT X-AXIS DIRECTION COSINE
	//"               VINC    INCIDENT Y-AXIS DIRECTION COSINE
	//"               WINC    INCIDENT Z-AXIS DIRECTION COSINE
	//"                       NOTE: (UINC,VINC,WINC) GET AUTOMATICALLY NORMALIZED
	//"                             DEFAULTS TO (0.0,0.0,1.0)
	//;
	//"------------------------------------------------------------------------------
	//"
	//" SOURCE 10 OR 11:
	//"       SPECIFY MODE:
	//"               MODE = local (0)     IF RADIAL DISTRIBUTION IS TO BE INPUT
	//"                                    LOCALLY WHETHER THROUGH THE KEYBOARD
	//"                                    (INTERACTIVE USE) OR THROUGH THE
	//"                                    .INP FILE (DEFAULT)
	//"                    = external (1)  IF THE DISTRIBUTION IS TO BE INPUT
	//"                                    VIA AN EXTERNAL FILE
	//;
	//"       IF MODE=local SPECIFY:
	//"
	//"           NUMBER OF RADIAL BINS
	//"           RDISTF(I),RPDF(I)  I=1,NRDIST
	//"                      TOP OF RADIAL BIN AND PROBABILITY OF
	//"                      INITIAL PARTICLE BEING IN THIS BIN.
	//"                      PROBABILITY DOES NOT NEED TO BE NORMALIZED
	//"                      BUT IT SHOULD BE IN UNITS CM**-2
	//"
	//"       IF MODE=external SPECIFY:
	//"
	//"           RDIST FILNAME = FILENAME(WITH EXT) CONTAINING THE ABOVE INFORMATION
	//"
	//"
	//"       SPECIFY:
	//"           DISTRIBUTION DATA = NONE => NO DISTRIBUTION DATA IN OUTPUT SUMMARY
	//"                             = OUTPUT SUMMARY => INCLUDE DISTRIBUTION
	//"                               DATA IN OUTPUT SUMMARY

		/**
		 * Handles radiation source in spherical geometry. 
		 */
		public static void SRCSPH()
		{
			if(ipart<0 || ipart>4){ipart=0;}//default

			if((ISOURC<0 || ISOURC>1)
			  || (ISOURC<10 || ISOURC>11) || (ISOURC!=4))
			{
				ISOURC=0;//default
			}

			//ERR15 = 15;
			//"************************************************************
			//" IF ISOURC>4, QUIT.  WE'RE LEAVING THE CODE IN, BUT SOURCES
			//" 10 AND 11 SEEM TO GO INTO AN INFINITE LOOP AFTER HATCHING
			//" (On entering srchst) -- JT
			//" ***********************************************************

			if (ISOURC>4)
			{
				EGS4.STOPPROGRAM=true;
				EGS4.seqStr="ERROR: SOURCES WITH RADIAL DISTRIBUTIONS ARE NOT AVAILABLE.";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="THEY WERE FOUND TO HANG AFTER THE CALL TO HATCH. THIS HAS NOT AS YET";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="BEEN FIXED SINCE WE HAVE NOT NEEDED THESE SOURCES.";
				eq.printSequence(EGS4.seqStr);
				return;
			}//;

			nsource_option=5;
		    for(int i=1;i<=nsource_option;i++)
		    {
			   if(source_option[i-1]<source_option_min || source_option[i-1]>source_option_max)
		   		source_option[i-1]=source_option_default;
		    }

			iqin=ipart;
			if (iqin==3) {iqin=-1;}
			TEMP1=source_option[0];TEMP2=source_option[1];TEMP3=source_option[2];
			TEMP4=source_option[3];TEMP5=source_option[4];

			EGS4.seqStr="SRCSPH READ:"+EGS4.format(iqin,5)+EGS4.format(ISOURC,5)+
			EGS4.format(TEMP1,12,true)+EGS4.format(TEMP2,12,true)+EGS4.format(TEMP3,12,true)+
			EGS4.format(TEMP4,12,true)+EGS4.format(TEMP5,12,true);
			if(EGS4.iprint>1)
				eq.printSequence(EGS4.seqStr);

	/*
			IF((ISOURC.EQ.10).OR.(ISOURC.EQ.11))
			{// "INCLUDE RADIAL DISTRIBUTION"
			     IVAL = IVAL + 1;
			    NUM_MODEIN = IVAL;
			    VALUES_SOUGHT(IVAL)='MODE';
			    NVALUE(IVAL)=1;
			    TYPE(IVAL)=3;
			    DEFAULT(IVAL)=1;
			    ALLOWED_INPUTS(IVAL,0)='LOCAL';
			    ALLOWED_INPUTS(IVAL,1)='EXTERNAL';
			    $GET_INPUT(NUM_MODEIN);
			    MODEIN=VALUE(NUM_MODEIN,1);

			    IF(MODEIN.EQ.0)[

			        IVAL = IVAL + 1;
			        NUM_RDISTF = IVAL;
			        VALUES_SOUGHT(IVAL)='RDISTF';
			        TYPE(IVAL)=1;
			        VALUE_MIN(IVAL)=-999999;
			        VALUE_MAX(IVAL)=999999;
			        $GET_INPUT(NUM_RDISTF);

			        IVAL = IVAL + 1;
			        NUM_RPDF = IVAL;
			        VALUES_SOUGHT(IVAL)='RPDF';
			        TYPE(IVAL)=1;
			        VALUE_MIN(IVAL)=0;
			        VALUE_MAX(IVAL)=999999;
			        $GET_INPUT(NUM_RPDF);

			        IF(NVALUE(NUM_RDISTF).NE.NVALUE(NUM_RPDF))[
			          WRITE(ERR15,*);
			          WRITE(ERR15,*)'**********************ERROR**********************';
			          WRITE(ERR15,*)'RDISTF HAS',NVALUE(NUM_RDISTF),' VALUES';
			          WRITE(ERR15,*)'RPDF HAS  ',NVALUE(NUM_RPDF),' VALUES';
			          WRITE(ERR15,*)'>>> THEY MUST HAVE THE SAME NUMBER OF VALUES <<<';
			          WRITE(ERR15,*);
			          ERROR_FLAG=1; RETURN;
			        ] ELSE [
			          NRDIST=NVALUE(NUM_RDISTF);
			        ]
			    ]
			    ELSE[ "INPUT FROM SEPARATE SOURCE FILE"
			        IVAL = IVAL + 1;
			        NUM_RDFIL = IVAL;
			        VALUES_SOUGHT(IVAL)='RDIST FILENAME';
			        NVALUE(IVAL)=1;
			        TYPE(IVAL)=2;
			        $GET_INPUT(NUM_RDFIL);
			        READ (CHAR_VALUE(NUM_RDFIL,1),FMT='(A)') SPCFIL;
			    ]
			    NUM_RDIOUTSP = IVAL;
			    VALUES_SOUGHT(IVAL)='DISTRIBUTION DATA';
			    NVALUE(IVAL)=1;
			    TYPE(IVAL)=3;
			    DEFAULT(IVAL)=0;
			    ALLOWED_INPUTS(IVAL,0)='NONE';
			    ALLOWED_INPUTS(IVAL,1)='OUTPUT SUMMARY';
			    $GET_INPUT(NUM_RDIOUTSP);
			    IOUTSP=VALUE(NUM_RDIOUTSP,1);

			    IF(IOUTSP.NE.1) IOUTSP=0;
			]
	*/
			if((iqin<-1)||(iqin>1)) iqin=0;

			SVTMP1=TEMP1;SVTMP2=TEMP2;SVTMP3=TEMP3;SVTMP4=TEMP4;

			return;
		}

		/**
		 * Print out source information.
		 */
		public static void SRCOUT()
		{
			if (EGS4.iGeom==EGS4.iCavity)//RZ GEOM
			{
				EGS4.seqStr="=========================================================================";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="                   SOURCE PARAMETERS";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="=========================================================================";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
	//ISOURC=20;IOUTSP=1;NRDIST=2;
				if(ISOURC == 0)
				{
					EGS4.seqStr=" Parallel beam on front face, radius="+
					EGS4.format(RBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" X-axis direction cosine="+
					EGS4.format(UINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Y-axis direction cosine="+
					EGS4.format(VINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Z-axis direction cosine="+
					EGS4.format(WINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Incident fluence="+
					EGS4.format(AINFLU,13,false)+" /cm**2";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,510) RBEAM,UINC,VINC,WINC,AINFLU;
				}
				else if(ISOURC == 1)
				{
					EGS4.seqStr=" Point source on cyl axis"+
					EGS4.format(DISTZ,8,true)+" cm away, collimated to R ="+
					EGS4.format(RBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Incident fluence(at front face on axis) * Solid angle="+
					EGS4.format(AINFLU,14,false);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,520) DISTZ,RBEAM,AINFLU;
				}
				else if(ISOURC == 2)
				{
					EGS4.seqStr=" BROAD PARALLEL BEAM ON CYLINDRICAL AXIS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,560);
				}
				else if(ISOURC == 3)
				{
					EGS4.seqStr=" ISOTROPICALLY RADIATING, UNIFORM SOURCE LOCATED";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr="      between radii"+
					EGS4.format(RMINBM,10,true)+
					EGS4.format(RBEAM,10,true)+" and DEPTHS"+
					EGS4.format(ZSMIN,10,true)+
					EGS4.format(ZSMAX,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,570)RMINBM,RBEAM,ZSMIN,ZSMAX;
				}
				else if(ISOURC == 4)
				{
					EGS4.seqStr=" CENTRAL AXIS FLUENCE VS BEAM RADIUS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr="      CENTRAL AXIS RADIUS="+
					EGS4.format(RBEAM,9,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" TREAT RADIAL BINS AS BEAM RADIUS, NOT SCORING REGION RADIUS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,590) RBEAM;
				}
				else if(ISOURC == 10)
				{
					EGS4.seqStr=" PARALLEL BEAM INCIDENT ON CYLINDRICAL WALL";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr="   RADIAL DIMENSION="+
					EGS4.format(XBEAM,8,true)+" cm, HALF-LENGTH="+
					EGS4.format(ZBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,530)XBEAM,ZBEAM;
				    if(((XBEAM/$ONE_EPS)<EGS4Geom.RCYL[EGS4Geom.NR])||
				        ((ZBEAM/$ONE_EPS)<(0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]))))
				    {
						EGS4.seqStr="       NOTE BEAM DOES NOT COVER THE DETECTOR!!";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
				    	//WRITE(IOUT,550);
					}
				}
				else if(ISOURC == 11)
				{
					EGS4.seqStr=" POINT SOURCE INCIDENT FROM SIDE"+
					EGS4.format(DISTRH,8,true)+" cm FROM MID-POINT";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr="   RADIAL DIMENSION="+
					EGS4.format(XBEAM,8,true)+" cm, HALF-LENGTH="+
					EGS4.format(ZBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,540)DISTRH,XBEAM,ZBEAM;
				    if((XBEAM/$ONE_EPS<(EGS4Geom.RCYL[EGS4Geom.NR]*Math.sqrt(DISTRH*DISTRH-
				    	EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR])/DISTRH))||
				       (ZBEAM/$ONE_EPS<(0.5*(EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]-EGS4Geom.ZPLANE[0]))))
				    {
						EGS4.seqStr="       NOTE BEAM DOES NOT COVER THE DETECTOR!!";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
						//WRITE(IOUT,550);
					}
				}
				else if(ISOURC == 12)
				{
					EGS4.seqStr=" POINT SOURCE OFF AXIS, RADIAL COORDINATE "+
					EGS4.format(DISTRH,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" PERPENDICULAR DISTANCE OF SOURCE FROM FRONT FACE "+
					EGS4.format(DISTZ,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,580)DISTRH,DISTZ;
				}
				else if(ISOURC == 15){src15_describe();}
				else if(ISOURC == 16){src16_describe(); }
				else if(ISOURC == 13)
				{
					EGS4.seqStr=" BROAD PARALLEL BEAM FROM ANY ANGLE";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" X-axis direction cosine="+
					EGS4.format(UINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Y-axis direction cosine="+
					EGS4.format(VINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Z-axis direction cosine="+
					EGS4.format(WINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Incident fluence="+
					EGS4.format(AINFLU,13,false)+" /cm**2";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,585) UINC,VINC,WINC,AINFLU;
				}
				else if(ISOURC == 14)
				{
					EGS4.seqStr=" Point source on cyl axis"+
					EGS4.format(DISTZ,8,true)+" cm away";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Collimated to radii ="+
					EGS4.format(RMINBM,8,true)+" cm to "+
					EGS4.format(RBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" Incident fluence(at front face on axis) * SOLID ANGLE="+
					EGS4.format(AINFLU,8,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,525) DISTZ,RMINBM,RBEAM,AINFLU;
				}
				else if(ISOURC == 20)
				{
					EGS4.seqStr=" RADIAL DISTRIBUTION, FRONTAL PARELLEL BEAM";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

				    //WRITE(IOUT,595);
				    if(IOUTSP==1)
				    {
						EGS4.seqStr="  RADIAL DISTBN, # OF INCIDENT RADIAL BINS:"+NRDIST;
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);
						EGS4.seqStr="  BIN    KINETIC RADIAL     PROBABILITY  CUMULATIVE PROB";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);

				        //WRITE(IOUT,610)NRDIST;
				        for(int IB=1;IB<=NRDIST;IB++)
				        {
							String s="  "+EGS4.format(IB,3);
							s=s+EGS4.format("",4)+EGS4.format(RDISTF[IB-1],7,true)+
							EGS4.format("",12)+EGS4.format(RPDF[IB-1],10,false)+
							EGS4.format("",4)+EGS4.format(RCDF[IB-1],7,false);
							EGS4.seqStr=s;
							if(EGS4.iprint>1)
								eq.printSequence(EGS4.seqStr);

							//WRITE(IOUT,620)IB,RDISTF(IB),RPDF(IB),RCDF(IB);
							//620   FORMAT(' ',T20,I3,T28,F7.3,T45,1PE10.3,T60,0PF7.4)
						}
					}
				    if(IBNSOK!=0.0)
				    {
						EGS4.seqStr=" WARNING: SOME OF NORMALIZED BIN PROBABILITIES SO SMALL BINS MAY BE MISSED!";
						if(EGS4.iprint>1)
							eq.printSequence(EGS4.seqStr);

						//WRITE(IOUT,630);
					}
				}
				/*  //NOT ALLOWED
				ELSEIF(ISOURC = 21)[WRITE(IOUT,640) FILSPC,ZPLANE(1),AINFLU;]
				ELSEIF(ISOURC = 22)[
				    WRITE(IOUT,677) FILSPC,dist_phsp,theta_phsp,chamber_c,
				                    xoffset,yoffset,AINFLU;
				]
				ELSEIF(ISOURC=23)[
				    WRITE(IOUT,679) iqinc,min_weight_23,max_weight_23,
				           $cstring(the_beam_code),$cstring(the_pegs_file),
				           $cstring(the_input_file),
				           dist_phsp,theta_phsp,chamber_c,xoffset,yoffset,EKSRCM;
				]

				IF (ISOURC = 21 | ISOURC = 22) [
				   WRITE(IOUT,:abcd:) NCASE_PHSP,TEMP1,TEMP2,TEMP3,NINCSRC;
				     :abcd: FORMAT(T15, ' Total number of particles in file      :',I10/
				      T15, ' Total number of photons                :',I10/
				      T15, ' (the rest are electrons and positrons)'//
				      T15, ' Maximum kinetic energy of the particles:',F10.3,' MeV'/
				      T15, ' Minimum kinetic energy of the electrons:',F10.3,' MeV'/
				      T15, ' # particles incident when phase space created :',F12.0/);
				   IF(DOSE_STAT=1)[" cannot read no. of primary histories from this phsp"
				                   " source"
				      WRITE(IOUT,'(//'' ***WARNING***''/
				       '' Cannot read no. of primary (non-phsp) histories from ph-sp source.''/
				       '' Dose and fluence will be analyzed assuming each particle read from''/
				       '' the ph-sp file is an independent history.  May result in an''/
				       '' underestimate of uncertainties.''//)');
				   ]
				   IF(IPARALLEL>1)[
				       IF(PARNUM>=1 & PARNUM<=IPARALLEL)[
				         WRITE(IOUT,'(/'' This is one of '',I4,'' parallel jobs.''/
				                    '' It will use from particle '',I12,'' to particle '',I12,/
				                    '' from the phase space source in the simulation.''/)')
				              IPARALLEL,INT((PARNUM-1)*NCASE_PHSP/IPARALLEL)+1,
				                INT(PARNUM*NCASE_PHSP/IPARALLEL);
				       ]
				   ]
				   IF(NRCYCL>0)[
				       WRITE(IOUT,'(/'' Particles will be recycled '',
				                   I4,'' times before moving on to next one.''/)') NRCYCL;
				   ]
				]*/

			}
		    else if(EGS4.iGeom==EGS4.iCavitySPH)//2=SPH Geom!!!
		    {

				EGS4.seqStr="=========================================================================";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="                   SOURCE PARAMETERS";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr="=========================================================================";
				if(EGS4.iprint>1)
					eq.printSequence(EGS4.seqStr);

				//WRITE(IOUT,500);
				if(ISOURC==0)
				{
					EGS4.seqStr=" PARALLEL BEAM FROM ANY ANGLE, RADIUS ="+
					EGS4.format(RBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" X-axis direction cosine="+
					EGS4.format(UINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Y-axis direction cosine="+
					EGS4.format(VINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Z-axis direction cosine="+
					EGS4.format(WINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Incident fluence="+
					EGS4.format(AINFLU,13,false)+" /cm**2";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,510) RBEAM,UINC,VINC,WINC,AINFLU;
				}
				else if(ISOURC==1)
				{
					EGS4.seqStr=" POINT SOURCE FROM ANY ANGLE"+
					EGS4.format(DISTR,8,true)+" cm AWAY, RADIUS ="+
					EGS4.format(RBEAM,8,true)+" cm";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" X-axis direction cosine="+
					EGS4.format(UINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Y-axis direction cosine="+
					EGS4.format(VINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Z-axis direction cosine="+
					EGS4.format(WINC,10,true);
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=EGS4.format("",10)+" Incident fluence="+
					EGS4.format(AINFLU,13,false)+" /cm**2";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,520) DISTR,RBEAM,UINC,VINC,WINC,AINFLU;
				}
				else if(ISOURC==4)
				{
					EGS4.seqStr=" PARTICLE LOCATED AT ORIGIN MOVING PARALLEL TO Z-AXIS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);
					EGS4.seqStr=" USED FOR POINT SPREAD FUNCTION CALCULATIONS";
					if(EGS4.iprint>1)
						eq.printSequence(EGS4.seqStr);

					//WRITE(IOUT,524) ;
				}

				/* NOT ALLOWED!!
				ELSEIF(ISOURC.EQ.10)[WRITE(IOUT,530) UINC,VINC,WINC,AINFLU;            ]
				ELSEIF(ISOURC.EQ.11)[WRITE(IOUT,540) DISTR,UINC,VINC,WINC,AINFLU;      ]
				//NEXT IS 1 ONLY FOR 10,11!!
				IF(IDSTON.EQ.1)[ "RADIAL DISTRIBUTION"
				    WRITE(IOUT,550);
				    IF(IOUTSP.EQ.1)[
				        WRITE(IOUT,560)NRDIST;
				        DO IB=1,NRDIST[WRITE(IOUT,570)IB,RDISTF(IB),RPDF(IB),RCDF(IB);]
				        ]
				    IF(IBNSOK.NE.0.0) WRITE(IOUT,580);
				    ]
				*/
				return;

			}
		}
}
