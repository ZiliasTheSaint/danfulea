package danfulea.phys.egs;

/**
 * Utility class for EGS4, EGS4 core engine and for custom user applications. The core engine (EGS4Core) and EGS4 call some user defined methods. Here they are implemented based on the task at hand. 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 25 OCT. 2005
 */
public class EGS4Macro {
	public static EgsQuestion eq;
	// AUTO means it is taken directly from EGS4Core or EGS4!!
	// otherwise it is required from tute!!
	// @raycorr----------------------------------------------------------------------------
	public static boolean use_enhance = false;// cav@@@@@@@
	// //CS ENHANCEMENT FACTOR= 1.0 #Photon cross section scaling factors
	public static double cs_enhance = 0.0;// cav@@@@@@@@@@
	public static int LGLE = 0;// cav=AUTO
	public static double COHFAC = 0.0;// cav=AUTO
	public static double GMFP = 0.0;// cav=AUTO
	public static int irl = 0;// cav=AUTO
	public static int[] iefl = new int[EGS4.$MXREG];// dose->($MXREG
	public static double cs_enhance_current = 0.0;
	public static int ienhance = 0;
	// @END
	// raycorr--------------------------------------------------------------------------
	// @--------electron mfp
	// ---------------------------------------------------------------
	public static int ICSDA = 0;// cav
	public static double demfp = 0.0;// cav=AUTO
	// @-------end-electron mfp
	// ---------------------------------------------------------------
	// @---------user range
	// discard--------------------------------------------------------
	public static int irejct = 0;// tut7&cav
	public static double esave = 0.0;// tut 7
	public static double range = 0.0;// tut7=AUTO
	public static double eker1 = 0.0;// cav
	public static double eker0 = 0.0;// cav
	public static double r_cavity_min = 0.0;
	public static double r_cavity_max = 0.0;// cav
	public static double z_cavity_min = 0.0;// cav
	public static double z_cavity_max = 0.0;// cav
	public static double[] rangerr0 = new double[EGS4.$MXRANGE];// cav
	public static double[] rangerr1 = new double[EGS4.$MXRANGE];// cav

	public static double ZL = 0.0;// flur
	public static double RL = 0.0;
	public static double CDIST = 0.0;
	public static double ZMINR = 0.0;
	public static double ZMAXR = 0.0;
	public static double RMINR = 0.0;
	public static double RMAXR = 0.0;
	public static int IGEOM = 0;
	// @--END-user range
	// discard----------------------------------------------------------
	// @--------photon mfp
	// ---------------------------------------------------------------
	public static int n_split = 0;// cav
	public static int iifano = 0;// cav
	public static int IFULL = 0;// cav
	public static int IWATCH = 0;// cav
	public static int NFTIME = 0;// cav
	public static int NFMIN = 0;// cav
	public static int NFMAX = 0;// cav
	public static int IFORCE = 0;// cav
	public static int IQINC = 0;// cav
	public static double GWAIT = 0.0;// cav
	public static double GWTOLD = 0.0;// cav
	public static double EXPMFP = 0.0;// cav
	public static int irange_rej = 0;
	public static int INEED2 = 0;

	public static double CEXPTR = 0.0;// dose
	public static boolean do_fast_step = false;// dose
	public static int NEWNRC = 0;// dose
	public static int[] NP_INC = new int[EGS4.$MXSTACK];// FLUR->($MXSTACK)
	public static int MAFORC = 0;// flur
	public static double[] GWAITf = new double[EGS4.$MXSTACK];// flur->($MXSTACK)

	public static boolean interact_now = false;// Edk

	public static boolean PCUT_DISCARD = false;// cav=AUTO
	public static boolean returnB = false;// cav=AUTO
	public static double peig = 0.0;// cav=AUTO
	public static double eig = 0.0;// cav=AUTO
	public static double GMFPR0 = 0.0;// cav=AUTO
	public static double GBR1 = 0.0;// cav=AUTO
	public static double GBR2 = 0.0;// cav=AUTO
	// @------end -photon
	// mfp------------------------------------------------------------
	// -----------------related to mfp//select mfp parallel beam also
	public static int ismfpfpb = 0;
	public static double[] GWATE = new double[EGS4Geom.$MAXRADII];// ($MAXRADII);
	// "PHOTON INTERACTION FORCING WEIGHTING FACTORS FOR
	// "NORMALLY INCIDENT PARALLEL BEAMS
	public static double DELTAP = 0.0;
	public static double PATHL = 0.0;
	// -----------------------------
	// /ADITIONAL SHOWER=======================================
	public static double DNEARIN = 0.0;
	public static int LATCHIN = 0;
	public static double ei = 0.0;// "total energy"
	public static double Mgle = 0.0;
	public static int LMgle = 0;// "index for GMFP interpolation"
	public static double Mcohfac = 0.0;
	public static int IARG = 0;
	public static int idopp = 0;
	// ---------
	public static double mXXX = 0.0;
	public static double mQ2 = 0.0;
	public static double mX2 = 0.0;
	public static double mCSQTHE = 0.0;
	public static double mREJF = 0.0;
	public static double mRNNORJ = 0.0;
	public static int LmXXX = 0;
	// /=====================================================
	private static int istartRayCorr = 0;
	private static int istartElec = 0;
	private static int istartFrontal = 0;
	private static int istartPhot = 0;

	// setup=the var...this is necessary for changing arrays dimensions!!
	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
		use_enhance = false;
		cs_enhance = 0.0;
		LGLE = 0;
		COHFAC = 0.0;
		GMFP = 0.0;
		irl = 0;
		iefl = new int[EGS4.$MXREG];
		cs_enhance_current = 0.0;
		ienhance = 0;
		ICSDA = 0;
		demfp = 0.0;
		irejct = 0;
		esave = 0.0;
		range = 0.0;
		eker1 = 0.0;
		eker0 = 0.0;
		r_cavity_min = 0.0;
		r_cavity_max = 0.0;
		z_cavity_min = 0.0;
		z_cavity_max = 0.0;
		rangerr0 = new double[EGS4.$MXRANGE];
		rangerr1 = new double[EGS4.$MXRANGE];
		ZL = 0.0;
		RL = 0.0;
		CDIST = 0.0;
		ZMINR = 0.0;
		ZMAXR = 0.0;
		RMINR = 0.0;
		RMAXR = 0.0;
		IGEOM = 0;
		n_split = 0;
		iifano = 0;
		IFULL = 0;
		IWATCH = 0;
		NFTIME = 0;
		NFMIN = 0;
		NFMAX = 0;
		IFORCE = 0;
		IQINC = 0;
		GWAIT = 0.0;
		GWTOLD = 0.0;
		EXPMFP = 0.0;
		irange_rej = 0;
		INEED2 = 0;
		CEXPTR = 0.0;
		do_fast_step = false;
		NEWNRC = 0;
		NP_INC = new int[EGS4.$MXSTACK];
		MAFORC = 0;
		GWAITf = new double[EGS4.$MXSTACK];
		interact_now = false;
		PCUT_DISCARD = false;
		returnB = false;
		peig = 0.0;
		eig = 0.0;
		GMFPR0 = 0.0;
		GBR1 = 0.0;
		GBR2 = 0.0;
		ismfpfpb = 0;
		GWATE = new double[EGS4Geom.$MAXRADII];
		DELTAP = 0.0;
		PATHL = 0.0;
		DNEARIN = 0.0;
		LATCHIN = 0;
		ei = 0.0;
		Mgle = 0.0;
		LMgle = 0;
		Mcohfac = 0.0;
		IARG = 0;
		idopp = 0;
		mXXX = 0.0;
		mQ2 = 0.0;
		mX2 = 0.0;
		mCSQTHE = 0.0;
		mREJF = 0.0;
		mRNNORJ = 0.0;
		LmXXX = 0;
		istartRayCorr = 0;
		istartElec = 0;
		istartFrontal = 0;
		istartPhot = 0;
	}

	/**
	 * Called by core engine, EGS4Core in ELECTR.
	 * @return the result
	 */
	public static boolean USER_RANGE_DISCARD() {
		if (EGS4.iurd == EGS4.iTut7) {
			if ((irejct == 1) && (EGS4.E[EGS4.NP - 1] > esave))// e(np) > esave
																// )) [
			{
				// "As tperp and range already known, check for a simple"
				// "range rejection in the present region               "
				if (EGS4.tperp >= range)// range) [
				{
					EGS4.IDISC = 50 + 49 * EGS4.IQ[EGS4.NP - 1];// iq(np);
					// idisc = 50 + 49*iq(np);
					// "1 for electrons, 99 for positrons"
					// go to :USER-ELECTRON-DISCARD: ;
					return true;
				}
			}
		} else if ((EGS4.iurd == EGS4.iCavity)
				|| (EGS4.iurd == EGS4.iCavitySPH)) {
			if (irejct == 1) {
				// "Note that in EGSnrc the standard range rejection has already"
				// "rejected a particle if it cannot get out of its local region"
				// "and if it is below e_max_rr for that region"

				// "The following implements the original CAVRZnrc range rejection    "
				// "but has a better range calculation (not just range in graphite)"
				// "Note that initialize_range_rejection() must have been called   "
				if (EGS4.iurd == EGS4.iCavity) {
					range_rejection(EGS4.ELKE, EGS4.X[EGS4.NP - 1],
							EGS4.Y[EGS4.NP - 1], EGS4.Z[EGS4.NP - 1],
							EGS4.IR[EGS4.NP - 1], EGS4.IQ[EGS4.NP - 1]);// ,EGS4.IDISC);

					if (EGS4.IDISC != 0) {
						// go to :USER-ELECTRON-DISCARD:
						return true;
					}

				}
				if (EGS4.iurd == EGS4.iCavitySPH) {
					if (EGS4.tperp >= range)// range) [
					{
						EGS4.IDISC = 50 + 49 * EGS4.IQ[EGS4.NP - 1];// iq(np);
						// idisc = 50 + 49*iq(np);
						// "1 for electrons, 99 for positrons"
						// go to :USER-ELECTRON-DISCARD: ;
						return true;
					}

				}

			}

		} else if (EGS4.iurd == EGS4.iFlur) {
			if ((irejct > 0) && (EGS4.MED[irl - 1] > 0)
					&& (EGS4.IQ[EGS4.NP - 1] != 0)) {
				ZL = EGS4.Z[EGS4.NP - 1];
				RL = Math.sqrt(EGS4.X[EGS4.NP - 1] * EGS4.X[EGS4.NP - 1]
						+ EGS4.Y[EGS4.NP - 1] * EGS4.Y[EGS4.NP - 1]);
				IGEOM = EGS4Geom.ntrack[irl - 1];
				CDIST = 0.0; // "Fixed bug for spectra in central region DR/CMA Jan 94"
								// "need to set CDIST when inside region of interest"
				if (IGEOM == 10) {
					CDIST = ZMINR - ZL;
				} else if (IGEOM == 20) {
					CDIST = RL - RMAXR;
				} else if (IGEOM == 30) {
					CDIST = ZL - ZMAXR;
				} else if (IGEOM == 40) {
					CDIST = EGS4.min(ZL - ZMINR, ZMAXR - ZL, RMINR - RL);
				} else if (IGEOM == 50) {
					CDIST = Math.sqrt((RL - RMAXR) * (RL - RMAXR)
							+ (ZMINR - ZL) * (ZMINR - ZL));
				} else if (IGEOM == 60) {
					CDIST = Math.sqrt((RL - RMAXR) * (RL - RMAXR)
							+ (ZL - ZMAXR) * (ZL - ZMAXR));
				} else if (IGEOM == 70) {
					CDIST = Math.sqrt((RL - RMINR) * (RL - RMINR)
							+ (ZMINR - ZL) * (ZMINR - ZL));
				} else if (IGEOM == 80) {
					CDIST = Math.sqrt((RL - RMINR) * (RL - RMINR)
							+ (ZL - ZMAXR) * (ZL - ZMAXR));
				}
				if (CDIST >= range) {
					EGS4.IDISC = 1;
					// replace with go to :USER-ELECTRON-DISCARD: ;
					return true;

				}
			}
		}

		return false;// no user range discard
	}

	/**
	 * Called by core engine, EGS4Core in ELECTR.
	 */
	public static void USER_CONTROLS_TSTEP_RECURSION() {
		if (EGS4.USER_CONTROLS_TSTEP_RECURSION == EGS4.iSpr)// index
		{
			// System.out.println("VINE@@@@@@@");
			if (EGS4.IRNEW != EGS4.IROLD) {
				EGS4.nextB = true;
				return;
			}
		}

		EGS4.nextB = false;
		return;
	}

	/**
	 * Called by EGS4 in HATCH
	 */
	public static void HATCH_USER_INPUT_INIT() {
		if (EGS4.hatchindex == EGS4.iCavitySPH) {
			// -------------------
			for (int J = 0; J < EGS4.$MXREG; J++) {
				if (EGS4.SMAXIR[J] <= 0.0) {
					EGS4.SMAXIR[J] = 1.E10;
				}
			}
			EGS4.LATCHI = 0;
			// "  IF(ESTEPR(J)<=0.0) [ESTEPR(J)=1;] <=== removed this since ESTEPR is not used"
			// ---------------------
			return;
		}

		// the next is the default one
		for (int J = 0; J < EGS4.$MXREG; J++) {
			if (EGS4.SMAXIR[J] <= 0.0) {
				EGS4.SMAXIR[J] = 1.E10;
			}
		}
		return;
	}

	/**
	 * Rayleigh correction/enhancement method. Called by EGS4Core in PHOTON.
	 */
	public static void RAYLEIGH_CORRECTION() {
		if ((EGS4.iraycorr == EGS4.iCavity)
				|| (EGS4.iraycorr == EGS4.iCavitySPH)) {
			if (EGS4.IRAYLR[irl - 1] == 1) {
				// $EVALUATE COHFAC USING COHE(GLE);
				COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1] * EGS4.GLE
						+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

				GMFP = GMFP * COHFAC;
			}
			if (use_enhance) {
				GMFP = GMFP / cs_enhance;
			}
			// " the above is for the cross section enhancement technique, see "
			// " description around definition of common SCORE                "
		} else if (EGS4.iraycorr == EGS4.iDose) {
			if (EGS4.IRAYLR[irl - 1] == 1) {
				// $EVALUATE COHFAC USING COHE(GLE);
				COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1] * EGS4.GLE
						+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

				GMFP = GMFP * COHFAC;
			}
			// "enhancement flag to beef up a local cross-section"
			if (iefl[irl - 1] == 1) {// "these are all zero if cs_enhance is 1.0"
				cs_enhance_current = cs_enhance * EGS4.RHO[EGS4.MEDIUM - 1]
						* GMFP;
				if (cs_enhance_current > 1) {
					GMFP = GMFP / cs_enhance_current;
					ienhance = 1;
				} else {
					cs_enhance_current = 1;
					ienhance = 0;
				}
			}// "end of block for regions in which enhancement done"
			else {
				cs_enhance_current = 1;
				ienhance = 0;
			}

		} else// default one
		{

			if (istartRayCorr == 0)// print once
			{
				EGS4.seqStr = "RAYLEIGH_CORRECTION is required BUT wrong option is used";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "Thus, the default algorithm will be used!!!";
				eq.printSequence(EGS4.seqStr);
			}

			if (EGS4.IRAYLR[irl - 1] == 1) {
				// $EVALUATE COHFAC USING COHE(GLE);
				COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1] * EGS4.GLE
						+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];
				GMFP = GMFP * COHFAC;
			}
			istartRayCorr = 1;// flag
		}
	}

	/**
	 * Called by EGS4Core in ELECTR.
	 */
	public static void SELECT_ELECTRON_MFP() {
		double RNNE1 = 0.0;
		if ((EGS4.isemfp == EGS4.iCavity) || (EGS4.isemfp == EGS4.iCavitySPH)
				|| (EGS4.isemfp == EGS4.iDose)) {
			if (ICSDA == 0) {
				RNNE1 = EGS4.random01();
				if (RNNE1 == 0.0) {
					RNNE1 = 1.E-30;
				}
				// "DEMFP=AMAX1(-ALOG(RNNE1),$EPSEMFP);"
				demfp = Math.max(-Math.log(RNNE1), EGS4Core.$EPSEMFP);
			} else {
				demfp = EGS4.VACDST;
			}
		} else {
			if (istartElec == 0)// print once
			{
				EGS4.seqStr = "SELECT_ELECTRON_MFP is required BUT wrong option is used";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "Thus, the default algorithm will be used!!!";
				eq.printSequence(EGS4.seqStr);
			}

			RNNE1 = EGS4.random01();
			if (RNNE1 == 0.0) {
				RNNE1 = 1.E-30;
			}
			demfp = Math.max(-Math.log(RNNE1), EGS4Core.$EPSEMFP);

			istartElec = 1;
		}

	}

	/**
	 * Called by user application at initialization step (e.g. GammaDetEffA). 
	 */
	public static void SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM() {
		if ((ismfpfpb == EGS4.iCavity) || (ismfpfpb == EGS4.iSpr)) {
			if ((EGS4SrcEns.IFPB == 0) && (IQINC == 0) && (IFORCE == 1)) {
				EGS4.GLE = Math.log(EGS4SrcEns.ein);
				for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
					PATHL = 0.0;
					for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
						// $GET-IRL(IZ,IX);
						// "MACRO THAT GETS THE GEOMETRY NUMBER FROM THE PLANAR
						// AND RADIAL ZONES
						// REPLACE {$GET-IRL(#,#);} WITH
						// {;IRL={P1}+NZ*({P2}-1)+1;}
						irl = EGS4Geom.GET_IRL(IZ, IX);// IZ+NZ*(IX-1)+1;

						EGS4.MEDIUM = EGS4.MED[irl - 1];
						if (EGS4.MEDIUM == 0) {// "vacuum"
							DELTAP = 0.;
						} else {
							// $SET INTERVAL GLE,GE;
							Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1]
									* EGS4.GLE + EGS4.GE0[EGS4.MEDIUM - 1]);
							LGLE = dbl.intValue();

							// $EVALUATE DELTAP USING GMFP(GLE);
							DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
									* EGS4.GLE
									+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

							if (EGS4.IRAYLR[irl - 1] == 1) {
								// $EVALUATE COHFAC USING COHE(GLE);
								COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

							} else {
								COHFAC = 1.0;
							}
						}
						if (DELTAP != 0) {// "only add to pathl if not a vacuum"
											// DELTAP=(ZPLANE(IZ+1)-ZPLANE(IZ))/(COHFAC*DELTAP);
							DELTAP = (EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1])
									/ (COHFAC * DELTAP);
							PATHL = PATHL + DELTAP;
						}
						if (PATHL <= 1.0E-3) {
							// GWATE(IX)=PATHL*(1.-0.5*PATHL);
							GWATE[IX - 1] = PATHL * (1. - 0.5 * PATHL);
						} else {
							// GWATE(IX)=1.-EXP(-PATHL);
							GWATE[IX - 1] = 1. - Math.exp(-PATHL);
						}
					}
				}
			}
		} else if (ismfpfpb == EGS4.iDose) {
			if ((EGS4SrcEns.IFPB == 0) && (IFORCE == 1)) {
				EGS4.GLE = Math.log(EGS4SrcEns.ein);
				for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
					PATHL = 0.0;
					for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
						// $GET-IRL(IZ,IX);
						// "MACRO THAT GETS THE GEOMETRY NUMBER FROM THE PLANAR
						// AND RADIAL ZONES
						// REPLACE {$GET-IRL(#,#);} WITH
						// {;IRL={P1}+NZ*({P2}-1)+1;}
						irl = EGS4Geom.GET_IRL(IZ, IX);// IZ+NZ*(IX-1)+1;

						EGS4.MEDIUM = EGS4.MED[irl - 1];
						if (EGS4.MEDIUM == 0) {// "vacuum"
							DELTAP = 0.;
						} else {
							// $SET INTERVAL GLE,GE;
							Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1]
									* EGS4.GLE + EGS4.GE0[EGS4.MEDIUM - 1]);
							LGLE = dbl.intValue();

							// $EVALUATE DELTAP USING GMFP(GLE);
							DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
									* EGS4.GLE
									+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

							if (EGS4.IRAYLR[irl - 1] == 1) {
								// $EVALUATE COHFAC USING COHE(GLE);
								COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

							} else {
								COHFAC = 1.0;
							}
						}
						if (DELTAP != 0) {// "only add to pathl if not a vacuum"
											// DELTAP=(ZPLANE(IZ+1)-ZPLANE(IZ))/(COHFAC*DELTAP);
							DELTAP = (EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1])
									/ (COHFAC * DELTAP);
							PATHL = PATHL + DELTAP;
						}
						if (PATHL <= 1.0E-3) {
							// GWATE(IX)=PATHL*(1.-0.5*PATHL);
							GWATE[IX - 1] = PATHL * (1. - 0.5 * PATHL);
						} else {
							// GWATE(IX)=1.-EXP(-PATHL);
							GWATE[IX - 1] = 1. - Math.exp(-PATHL);
						}
					}
				}

			}
		} else {
			if (istartFrontal == 0) {
				// do nothing but inform the user!!
				EGS4.seqStr = "SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM is required"
						+ " but incorrect options is used!!";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "ERROR!!=> This call will DO NOTHING!!";
				eq.printSequence(EGS4.seqStr);
			}
			istartFrontal = 1;
		}

	}

	/**
	 * Called by core engine, EGS4Core in PHOTON.
	 */
	public static void SELECT_PHOTON_MFP() {
		double DUMU = 0.0;
		double DUMX = 0.0;
		double DUMY = 0.0;
		double DUMZ = 0.0;
		double PATHL = 0.0;
		double DELTAP = 0.0;
		double EPSLON = 0.0;
		double ARG = 0.0;
		int IRODUM = 0;
		int IRNDUM = 0;
		int IRDUM = 0;
		int MEDDUM = 0;
		int IDUM = 0;
		int MEDTMP = 0;
		double d_eta = 0.0;
		double eta_prime = 0.0;
		double dpmfp_old = 0.0;
		double x_save = 0.0;
		double y_save = 0.0;
		double z_save = 0.0;
		double u_save = 0.0;
		double v_save = 0.0;
		double w_save = 0.0;
		double e_save = 0.0;
		double wt_save = 0.0;
		//double wt_start = 0.0;
		double a_survive = 0.0;
		int ir_save = 0;
		int ip = 0;
		int i_split = 0;
		//int np_start = 0;
		int latch_save = 0;
		int i_survive = 0;
		int i_survive_s = 0;

		double SPMFP = 0.0;
		double TEMP = 0.0;
		double TEMP1 = 0.0;
		int NRCDUM = 0;
		double PATHLT = 0.0;

		// =============================
		int LXXX = 0;
		double X2 = 0.0;
		double Q2 = 0.0;
		double CSQTHE = 0.0;
		double REJF = 0.0;
		double rnno35 = 0.0;
		double xxx = 0.0;
		double rnno37 = 0.0;
		double RNNORJ = 0.0;
		double RNNO36 = 0.0;
		boolean START_MFP_LOOP = false;
		boolean JUST_RAYLEIGH_EVENT = false;
		boolean SAMPLING_LOOPb = false;
		boolean FASTSTEP = false;
		// ============================
		if ((EGS4.ispmfp == EGS4.iCavity) || (EGS4.ispmfp == EGS4.iCavitySPH)) {
			// MACRO USED FOR FORCING INTERACTIONS IN THE GEOMETRY"
			// "USED BY EGS4 FOR VARIANCE REDUCTION"

			// REPLACE {$SELECT-PHOTON-MFP;} WITH {;
			if (n_split > 1) { // "we use photon splitting instead of interaction"
								// "forcing. This is much more efficient as it "
								// "always gives interactions in the chamber "
								// "(provided n_split is large enough) and does"
								// "not lead to varying weights"

				//np_start = EGS4.NP;
				//wt_start = EGS4.WT[EGS4.NP - 1];

				// START-MFP-LOOP:
				while (true) {
					START_MFP_LOOP = false;
					PCUT_DISCARD = false;

					rnno35 = EGS4.random01();
					rnno35 = rnno35 / n_split;
					d_eta = 1. / n_split;
					x_save = EGS4.X[EGS4.NP - 1];
					y_save = EGS4.Y[EGS4.NP - 1];
					z_save = EGS4.Z[EGS4.NP - 1];
					u_save = EGS4.U[EGS4.NP - 1];
					v_save = EGS4.V[EGS4.NP - 1];
					w_save = EGS4.W[EGS4.NP - 1];
					e_save = EGS4.E[EGS4.NP - 1];
					wt_save = EGS4.WT[EGS4.NP - 1] / n_split;
					ir_save = EGS4.IR[EGS4.NP - 1];
					latch_save = EGS4.LATCH[EGS4.NP - 1];
					EGS4.NP = EGS4.NP - 1;
					i_survive = -1;
					if ((iifano == 1) || ((IFULL == 1) && (latch_save != 3))) {
						xxx = EGS4.random01();
						a_survive = xxx * n_split;
						Double ad = new Double(a_survive);
						i_survive = 1 + ad.intValue();// a_survive;
					}
					if ((iifano == 1) || (latch_save == 2)) {
						i_survive_s = -1;
					} else {
						xxx = EGS4.random01();
						a_survive = xxx * n_split;
						Double ad = new Double(a_survive);
						i_survive_s = 1 + ad.intValue();// a_survive;
					}

					dpmfp_old = 0;

					if (EGS4.ispmfp == EGS4.iCavity)
						eta_prime = 1 - rnno35;
					else
						eta_prime = 1 - rnno35 + d_eta;// iCavitySPH
					// ##############################iCavitySPH remove
					if (EGS4.ispmfp == EGS4.iCavity) {
						if ((IWATCH != 0) && (IWATCH != 4)) {// "output a message about splitting"
							EGS4.seqStr = " Splitting photon into "
									+ EGS4.format(n_split, 4)
									+ " photons with weight "
									+ EGS4.format(wt_save, 10);
							if (EGS4.iprint > 0)// WATCH is next->allways show
								eq.printSequence(EGS4.seqStr);
							// OUTPUT n_split,wt_save;
							// (' Splitting photon into ',I4,' photons with
							// weight ',1PE10.3);
						}
					}
					// ##############################end iCavitySPH remove
					FOR: for (i_split = 1; i_split <= n_split; i_split++) {
						JUST_RAYLEIGH_EVENT = false;
						if (EGS4.ispmfp == EGS4.iCavitySPH)
							eta_prime = eta_prime - d_eta;// iCavitySPH
						if (eta_prime <= 0) {
							break FOR;
							// goto :END-MFP-LOOP:;
						}
						EGS4.DPMFP = -Math.log(eta_prime) - dpmfp_old;
						dpmfp_old = dpmfp_old + EGS4.DPMFP;
						EGS4.NP = EGS4.NP + 1;
						if (EGS4.NP > EGS4.$MXSTACK) {
							EGS4.STOPPROGRAM = true;
							EGS4.seqStr = " Stack overflow in $SELECT-PHOTON-MFP ";

							eq.printSequence(EGS4.seqStr);
							returnB = true;
							return;
							// write(6,*) ' Stack overflow in $SELECT-PHOTON-MFP
							// ';
							// stop;
						}
						EGS4.X[EGS4.NP - 1] = x_save;
						EGS4.Y[EGS4.NP - 1] = y_save;
						EGS4.Z[EGS4.NP - 1] = z_save;
						EGS4.U[EGS4.NP - 1] = u_save;
						EGS4.V[EGS4.NP - 1] = v_save;
						EGS4.W[EGS4.NP - 1] = w_save;
						EGS4.WT[EGS4.NP - 1] = wt_save;
						EGS4.E[EGS4.NP - 1] = e_save;
						EGS4.IQ[EGS4.NP - 1] = 0;
						EGS4.IR[EGS4.NP - 1] = ir_save;
						EGS4.LATCH[EGS4.NP - 1] = latch_save;
						if (EGS4.ispmfp == EGS4.iCavity)
							EGS4.DNEAR[EGS4.NP - 1] = 0;// ########iCavitySPH
														// remove
						irl = EGS4.IR[EGS4.NP - 1];
						EGS4.IROLD = irl;
						EGS4.MEDIUM = EGS4.MED[irl - 1];
						// LOOP
						while (true) {
							if (EGS4.MEDIUM != 0) {
								// $SET INTERVAL GLE,GE;
								Double dbl = new Double(
										EGS4.GE1[EGS4.MEDIUM - 1] * EGS4.GLE
												+ EGS4.GE0[EGS4.MEDIUM - 1]);
								LGLE = dbl.intValue();
								// $EVALUATE GMFPR0 USING GMFP(GLE);
								GMFPR0 = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];
								// $SET-RHOF;
								EGS4.RHOF = EGS4.RHOR[irl - 1]
										/ EGS4.RHO[EGS4.MEDIUM - 1];
								GMFP = GMFPR0 / EGS4.RHOF;

								// $RAYLEIGH-CORRECTION;
								if (EGS4.iraycorr == 0) {
									// $RAYLEIGH-CORRECTION;//
									// "A RAYLEIGH SCATTERING TEMPLATE"
									// REPLACE {$RAYLEIGH-CORRECTION;} WITH {
									if (EGS4.IRAYLR[irl - 1] == 1) {
										// $EVALUATE COHFAC USING COHE(GLE);
										COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
												* EGS4.GLE
												+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];
										GMFP = GMFP * COHFAC;
									}
									// }
								} else {
									RAYLEIGH_CORRECTION();
								}

								EGS4.TSTEP = GMFP * EGS4.DPMFP;
							} else {
								EGS4.TSTEP = EGS4.VACDST;
							}
							EGS4.IRNEW = irl;
							EGS4.IDISC = 0;
							EGS4.USTEP = EGS4.TSTEP;
							EGS4.TUSTEP = EGS4.USTEP;
							eq.HOWFAR();// call howfar;
							if (EGS4.ispmfp == EGS4.iCavity) {
								EGS4.VSTEP = EGS4.USTEP;// ########iCavitySPH
														// remove
								EGS4.TVSTEP = EGS4.VSTEP;// ########iCavitySPH
															// remove
								EGS4.EDEP = EGS4.PZERO;// ########iCavitySPH
														// remove
								if (IWATCH > 0)
									EGS4.WATCH(EGS4.$TRANAUSB, IWATCH);
							}
							EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1]
									+ EGS4.U[EGS4.NP - 1] * EGS4.USTEP;// ustep;
							EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1]
									+ EGS4.V[EGS4.NP - 1] * EGS4.USTEP;
							EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1]
									+ EGS4.W[EGS4.NP - 1] * EGS4.USTEP;
							if (EGS4.IDISC > 0) {
								if (EGS4.ispmfp == EGS4.iCavity)
									if (IWATCH > 0)
										EGS4.WATCH(EGS4.$USERDAUS, IWATCH);// ########iCavitySPH
																			// remove
								EGS4.NP = EGS4.NP - 1;
								if (EGS4.NP == 0) {
									EGS4Core.IRCODE = 2;
									returnB = true;
									return;
								}
								// goto :END-MFP-LOOP:;
								break FOR;
							}
							if (EGS4.MEDIUM != 0)
								EGS4.DPMFP = EGS4.DPMFP - EGS4.USTEP / GMFP;
							if (EGS4.IRNEW != EGS4.IROLD) {
								EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
								irl = EGS4.IRNEW;
								EGS4.IROLD = EGS4.IRNEW;
								EGS4.MEDIUM = EGS4.MED[irl - 1];
							}
							if (EGS4.ispmfp == EGS4.iCavity)
								if (IWATCH > 0)
									EGS4.WATCH(EGS4.$TRANAUSA, IWATCH);// ########iCavitySPH
																		// remove

							if ((EGS4.MEDIUM != 0)
									&& (EGS4.DPMFP < EGS4Core.$EPSGMFP))
								break;
						} // until ((EGS4.MEDIUM != 0) && (EGS4.DPMFP <
							// EGS4Core.$EPSGMFP));
						x_save = EGS4.X[EGS4.NP - 1];
						y_save = EGS4.Y[EGS4.NP - 1];
						z_save = EGS4.Z[EGS4.NP - 1];
						ir_save = EGS4.IR[EGS4.NP - 1];
						if (EGS4.IRAYLR[irl - 1] == 1) {
							rnno37 = EGS4.random01();
							if (rnno37 <= (1.0 - COHFAC)) {
								if (i_split != i_survive_s) {
									EGS4.NP = EGS4.NP - 1;
									// goto :JUST-RAYLEIGH-EVENT:;
									JUST_RAYLEIGH_EVENT = true;
								}

								if (!JUST_RAYLEIGH_EVENT) {
									if (IWATCH > 0)
										EGS4.WATCH(EGS4.$RAYLAUSB, IWATCH);
									EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1]
											* n_split;
									EGS4.LATCH[EGS4.NP - 1] = 3;// LAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCH

									RAYLEIGH_SAMPLING_LOOP: while (true)// LOOP
									{
										SAMPLING_LOOPb = false;
										xxx = EGS4.random01();// $RANDOMSET XXX;
										// $SET INTERVAL XXX,RCO;
										Double dbl1 = new Double(
												EGS4.RCO1[EGS4.MEDIUM - 1]
														* xxx
														+ EGS4.RCO0[EGS4.MEDIUM - 1]);
										LXXX = dbl1.intValue();// THIS CAN BE
																// LOCAL!!!!!!

										// $EVALUATE X2 USING RSCT(XXX);
										X2 = EGS4.RSCT1[LXXX - 1][EGS4.MEDIUM - 1]
												* xxx
												+ EGS4.RSCT0[LXXX - 1][EGS4.MEDIUM - 1];// //THIS
																						// CAN
																						// BE
																						// LOCAL!!!!!!

										Q2 = X2 * EGS4.RMSQ
												/ (20.60744 * 20.60744);// //THIS
																		// CAN
																		// BE
																		// LOCAL!!!!!!
										EGS4.COSTHE = 1.
												- Q2
												/ (2. * EGS4.E[EGS4.NP - 1] * EGS4.E[EGS4.NP - 1]);
										if (Math.abs(EGS4.COSTHE) > 1.0) {
											SAMPLING_LOOPb = true;
											// GO TO :RAYLEIGH-SAMPLING-LOOP:;
										}
										if (!SAMPLING_LOOPb) {
											CSQTHE = EGS4.COSTHE * EGS4.COSTHE;// THIS
																				// CAN
																				// BE
																				// LOCAL!!!!!!
											REJF = (1.0 + CSQTHE) / 2.0;// THIS
																		// CAN
																		// BE
																		// LOCAL!!!!!!
											RNNORJ = EGS4.random01();
											// ------------------------------------
											if (RNNORJ <= REJF) {
												break RAYLEIGH_SAMPLING_LOOP;
											}
										}
									} // UNTIL (RNNORJ <= REJF);
									EGS4.SINTHE = Math.sqrt(1.0 - CSQTHE);
									EGS4Core.UPHI(2, 1);// //$SELECT-AZIMUTHAL-ANGLE
														// and OLD-PARTICLE:

									if (EGS4.ispmfp == EGS4.iCavity)
										if (IWATCH > 0)
											EGS4.WATCH(EGS4.$RAYLAUSA, IWATCH);// ########iCavitySPH
																				// remove
									// goto :JUST-RAYLEIGH-EVENT:;
									JUST_RAYLEIGH_EVENT = true;
								}// JUST_RAYLEIGH_EVENT

							}
						}

						if (!JUST_RAYLEIGH_EVENT) {
							RNNO36 = EGS4.random01();// $RANDOMSET RNNO36;
							// $EVALUATE GBR1 USING GBR1(GLE);
							GBR1 = EGS4.GBR11[LGLE - 1][EGS4.MEDIUM - 1]
									* EGS4.GLE
									+ EGS4.GBR10[LGLE - 1][EGS4.MEDIUM - 1];

							if ((RNNO36 <= GBR1)
									&& (EGS4.E[EGS4.NP - 1] > EGS4.RMT2)) {
								if (EGS4.ispmfp == EGS4.iCavity)
									if (IWATCH > 0)
										EGS4.WATCH(EGS4.$PAIRAUSB, IWATCH);// ########iCavitySPH
																			// remove
								EGS4Core.PAIR();
								if (EGS4.ispmfp == EGS4.iCavity)
									if (IWATCH > 0)
										EGS4.WATCH(EGS4.$PAIRAUSA, IWATCH);// ########iCavitySPH
																			// remove
							} else {
								// $EVALUATE GBR2 USING GBR2(GLE);
								GBR2 = EGS4.GBR21[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.GBR20[LGLE - 1][EGS4.MEDIUM - 1];

								if (RNNO36 > GBR2) {
									if (EGS4.ispmfp == EGS4.iCavity)
										if (IWATCH > 0)
											EGS4.WATCH(EGS4.$COMPAUSB, IWATCH);// ########iCavitySPH
																				// remove
									EGS4Core.COMPT();// compt;
									if (EGS4.ispmfp == EGS4.iCavity)
										if (IWATCH > 0)
											EGS4.WATCH(EGS4.$COMPAUSA, IWATCH);// ########iCavitySPH
																				// remove
								} else {
									if (EGS4.ispmfp == EGS4.iCavity)
										if (IWATCH > 0)
											EGS4.WATCH(EGS4.$PHOTOAUSB, IWATCH);// ########iCavitySPH
																				// remove
									EGS4Core.PHOTO();// photo;
									if (EGS4.ispmfp == EGS4.iCavity)
										if (IWATCH > 0)
											EGS4.WATCH(EGS4.$PHOTOAUSA, IWATCH);// ########iCavitySPH
																				// remove
								}
							}
							// "If ifano is on, we traw away all scattered photons"
							// "if the latch of the photons is 2, it is a regenerated"
							// "primary photon => we need to remove resulting scatter"
							// "as well. To save some time, we also discard on the spot"
							// "all electrons that can not get into the cavity, if"
							// "rejection is on"
							ip = EGS4.NPold;
							while (true)// LOOP
							{
								if (EGS4.IQ[ip - 1] == 0) {
									if (i_split != i_survive_s) {
										if (ip < EGS4.NP) {
											// ########iCavitySPH remove
											if (EGS4.ispmfp == EGS4.iCavity) {
												if ((IWATCH != 0)
														&& (IWATCH != 4)) {
													// OUTPUT
													// ip,e(ip),iq(ip),ir(ip),x(ip),y(ip),
													// z(ip),u(ip),v(ip),w(ip),latch(ip),wt(ip);
													// (' Eliminating scattered
													// photon',T36,':',
													// I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
													String s = " Eliminating scattered photon";
													int ll = s.length();
													ll = 36 - ll;
													s = s
															+ EGS4.format(":",
																	ll);

													EGS4.seqStr = s
															+ EGS4.format(ip, 5)
															+ EGS4.format(
																	EGS4.E[ip - 1],
																	9, true)
															+ EGS4.format(
																	EGS4.IQ[ip - 1],
																	4)
															+ EGS4.format(
																	EGS4.IR[ip - 1],
																	4)
															+ EGS4.format(
																	EGS4.X[ip - 1],
																	8, true)
															+ EGS4.format(
																	EGS4.Y[ip - 1],
																	8, true)
															+ EGS4.format(
																	EGS4.Z[ip - 1],
																	8, true)
															+ EGS4.format(
																	EGS4.U[ip - 1],
																	7, true)
															+ EGS4.format(
																	EGS4.V[ip - 1],
																	7, true)
															+ EGS4.format(
																	EGS4.W[ip - 1],
																	7, true)
															+ EGS4.format(
																	EGS4.LATCH[ip - 1],
																	10)
															+ EGS4.format(
																	EGS4.WT[ip - 1],
																	10, false);
													if (EGS4.iprint > 0)
														eq.printSequence(EGS4.seqStr);

												}
											}
											// ########iCavitySPH remove
											EGS4.E[ip - 1] = EGS4.E[EGS4.NP - 1];
											EGS4.IQ[ip - 1] = EGS4.IQ[EGS4.NP - 1];
											EGS4.U[ip - 1] = EGS4.U[EGS4.NP - 1];
											EGS4.V[ip - 1] = EGS4.V[EGS4.NP - 1];
											EGS4.W[ip - 1] = EGS4.W[EGS4.NP - 1];
											EGS4.WT[ip - 1] = EGS4.WT[EGS4.NP - 1];
											EGS4.LATCH[ip - 1] = EGS4.LATCH[EGS4.NP - 1];
										}
										EGS4.NP = EGS4.NP - 1;
									} else {
										EGS4.WT[ip - 1] = EGS4.WT[ip - 1]
												* n_split;
										EGS4.LATCH[ip - 1] = 3;
										ip = ip + 1;
									}
								} else {
									if (irejct == 1) {
										EGS4.EDEP = EGS4.E[ip - 1] - EGS4.PRM;
										EGS4.ELKE = Math.log(EGS4.EDEP);// elke
																		// =
																		// log(edep);
										range_rejection(EGS4.ELKE,
												EGS4.X[ip - 1], EGS4.Y[ip - 1],
												EGS4.Z[ip - 1],
												EGS4.IR[ip - 1],
												EGS4.IQ[ip - 1]);// ,EGS4.IDISC);
										if (EGS4.IDISC == 0) {
											ip = ip + 1;
										} else {
											if (ip < EGS4.NP) {
												EGS4.E[ip - 1] = EGS4.E[EGS4.NP - 1];
												EGS4.IQ[ip - 1] = EGS4.IQ[EGS4.NP - 1];
												EGS4.U[ip - 1] = EGS4.U[EGS4.NP - 1];
												EGS4.V[ip - 1] = EGS4.V[EGS4.NP - 1];
												EGS4.W[ip - 1] = EGS4.W[EGS4.NP - 1];
												EGS4.WT[ip - 1] = EGS4.WT[EGS4.NP - 1];
												EGS4.LATCH[ip - 1] = EGS4.LATCH[EGS4.NP - 1];
											}
											EGS4.NP = EGS4.NP - 1;
											if (EGS4.ispmfp == EGS4.iCavity)
												if (IWATCH > 0)
													EGS4.WATCH(EGS4.$USERDAUS,
															IWATCH);// ########iCavitySPH
																	// remove
										}
									} else {
										ip = ip + 1;
									}
								}
								// -----------------
								if (ip > EGS4.NP)
									break;
								// ------------------
							} // UNTIL (ip > np);
						}// JUST_RAYLEIGH_EVENT

						// :JUST-RAYLEIGH-EVENT:

						// "IF( iifano = 1 | (ifull = 1 & latch_save ~= 3) ) [ "
						if (i_split == i_survive) {
							// "Re-generate the original photon with prob."
							// "1/n_split, so that we get back the original"
							// "weight. This is necessary when ifano is set "
							// "(i.e. dose calculation with attenuation and scatter removed)"
							// "or ifull is 1 (i.e. user wants to get Awall in addition to"
							// "the dose)"
							// "$RANDOMSET xxx;"
							// " IF(xxx*n_split < 1) ["
							EGS4.NP = EGS4.NP + 1;
							EGS4.X[EGS4.NP - 1] = x_save;
							EGS4.Y[EGS4.NP - 1] = y_save;
							EGS4.Z[EGS4.NP - 1] = z_save;
							EGS4.U[EGS4.NP - 1] = u_save;
							EGS4.V[EGS4.NP - 1] = v_save;
							EGS4.W[EGS4.NP - 1] = w_save;
							EGS4.E[EGS4.NP - 1] = e_save;
							EGS4.IR[EGS4.NP - 1] = ir_save;
							if (IFULL == 1) {
								EGS4.LATCH[EGS4.NP - 1] = 2;
							} else {
								EGS4.LATCH[EGS4.NP - 1] = 0;
							}
							EGS4.WT[EGS4.NP - 1] = wt_save * n_split;
							EGS4.IQ[EGS4.NP - 1] = 0;
							// ########iCavitySPH remove
							if (EGS4.ispmfp == EGS4.iCavity) {
								if ((IWATCH != 0) && (IWATCH != 4)) {
									String s = " Regenerating original photon";
									int ll = s.length();
									ll = 36 - ll;
									s = s + EGS4.format(":", ll);

									EGS4.seqStr = s
											+ EGS4.format(EGS4.NP, 5)
											+ EGS4.format(EGS4.E[EGS4.NP - 1],
													9, true)
											+ EGS4.format(EGS4.IQ[EGS4.NP - 1],
													4)
											+ EGS4.format(EGS4.IR[EGS4.NP - 1],
													4)
											+ EGS4.format(EGS4.X[EGS4.NP - 1],
													8, true)
											+ EGS4.format(EGS4.Y[EGS4.NP - 1],
													8, true)
											+ EGS4.format(EGS4.Z[EGS4.NP - 1],
													8, true)
											+ EGS4.format(EGS4.U[EGS4.NP - 1],
													7, true)
											+ EGS4.format(EGS4.V[EGS4.NP - 1],
													7, true)
											+ EGS4.format(EGS4.W[EGS4.NP - 1],
													7, true)
											+ EGS4.format(
													EGS4.LATCH[EGS4.NP - 1], 10)
											+ EGS4.format(EGS4.WT[EGS4.NP - 1],
													10, false);
									if (EGS4.iprint > 0)
										eq.printSequence(EGS4.seqStr);

									// OUTPUT
									// np,e(np),iq(np),ir(np),x(np),y(np),
									// z(np),u(np),v(np),w(np),latch(np),wt(np);
									// (' Regenerating original photon',T36,':',
									// I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
								}
							}
							// ########iCavitySPH remove
							// "]"
						}
						if (EGS4.ispmfp == EGS4.iCavity)
							eta_prime = eta_prime - d_eta;// ########iCavitySPH
															// remove
					}// for( i_split = 1;i_split<=n_split;i_split++)

					// :END-MFP-LOOP:;
					if (EGS4.NP <= 0) {
						EGS4Core.IRCODE = 2;
						returnB = true;
						return;
					}
					if (EGS4.IQ[EGS4.NP - 1] == 0) {
						peig = EGS4.E[EGS4.NP - 1];
						eig = peig;
						irl = EGS4.IR[EGS4.NP - 1];
						EGS4.MEDIUM = EGS4.MED[irl - 1];
						if (eig <= EGS4.PCUT[irl - 1]) {
							// GO TO :PCUT-DISCARD:;
							PCUT_DISCARD = true;
							return;
						}
						EGS4.GLE = Math.log(eig);

						START_MFP_LOOP = true;
						// goto :START-MFP-LOOP:;
					}

					if (!START_MFP_LOOP) {
						returnB = true;
						return;
					}

				}// while for START_MFP_LOOP
			}// if( n_split > 1 )

			NFTIME = NFTIME + 1;
			rnno35 = EGS4.random01();

			// "------------------------------------------------------------"
			// "CHANGE: J.S. Aug 95                                         "
			// "FOLLOWING CHANGE IS JUST TO PREVENT FAKE TRANSPORT IN CASE  "
			// "OF CORRELATED Afl CALCULATIONS.                             "
			// "PREVIOUSLY THERE WAS:                                       "
			// "IF((IFORCE.EQ.0).OR.(NFTIME.GT.NFMAX).OR.(NFTIME.LT.NFMIN))["
			// "NOW THE WT(NP)==0 IS INCLUDED TO PREVENT THE FORCING MACRO  "
			// "FROM INFINITE LOOPING DURING FAKE TRANSPORT FOR CORRELATED  "
			// "Afl CALCULATIONS WHEN SWITCHED TO PARALLEL BEAMS            "
			// "------------------------------------------------------------"
			if (((IFORCE == 0) || (NFTIME > NFMAX))
					|| ((NFTIME < NFMIN) || (EGS4.WT[EGS4.NP - 1] == 0))) {
				// "------------------------------------------------------------"
				// "END OF CHANGE                                               "
				// "------------------------------------------------------------"
				if (rnno35 == 0.0) {
					rnno35 = 1.E-30;
				}
				EGS4.DPMFP = -Math.log(rnno35);
			} else {
				if (EGS4.ispmfp == EGS4.iCavitySPH) {
					if (EGS4SrcEns.ISOURC == 4) {
						EGS4.DPMFP = 0.0;
						FASTSTEP = true;
						// GOTO :SKIPFORCING:;
					}
				}

				// ########iCavitySPH remove
				if (EGS4.ispmfp == EGS4.iCavity) {
					if (((EGS4SrcEns.IFPB == 0) && (IQINC == 0))
							&& ((NFTIME == 1) && (EGS4SrcEns.MONOEN == 0))) {
						// GOTO :FASTSTEP:;
						FASTSTEP = true;
					}
				}
				// "IFPB, IQINC IN COMIN USER-VARIANCE-REDUCTION"
				// ########iCavitySPH remove
				if (!FASTSTEP) {
					DUMU = EGS4.USTEP;
					DUMX = EGS4.X[EGS4.NP - 1];
					DUMY = EGS4.Y[EGS4.NP - 1];
					DUMZ = EGS4.Z[EGS4.NP - 1];
					IRODUM = EGS4.IROLD;
					IRNDUM = EGS4.IRNEW;
					IRDUM = EGS4.IR[EGS4.NP - 1];
					MEDDUM = EGS4.MEDIUM;
					IDUM = EGS4.IDISC;
					PATHL = 0.0;
					MEDTMP = 0;

					while (true)// LOOP
					{
						EGS4.USTEP = EGS4.VACDST;
						EGS4.IROLD = EGS4.IR[EGS4.NP - 1];
						EGS4.MEDIUM = EGS4.MED[EGS4.IROLD - 1];
						if (EGS4.ispmfp == EGS4.iCavitySPH)
							if (EGS4Geom.ntrack[EGS4.IROLD - 1] == 1)
								INEED2 = 1;
						// ########iCavitySPH remove
						if (EGS4.ispmfp == EGS4.iCavity) {
							if (EGS4.MEDIUM == 0) {// "vacuum"
								DELTAP = 0.;
							} else {
								// ########iCavitySPH remove
								if (MEDTMP != EGS4.MEDIUM) {
									MEDTMP = EGS4.MEDIUM;
									// $SET INTERVAL GLE,GE;
									Double dbl = new Double(
											EGS4.GE1[EGS4.MEDIUM - 1]
													* EGS4.GLE
													+ EGS4.GE0[EGS4.MEDIUM - 1]);
									LGLE = dbl.intValue();

									// $EVALUATE DELTAP USING GMFP(GLE);
									DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
											* EGS4.GLE
											+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

								}
								if (EGS4.IRAYLR[EGS4.IROLD - 1] == 1) {
									// $EVALUATE COHFAC USING COHE(GLE);
									COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
											* EGS4.GLE
											+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

								} else {
									COHFAC = 1.0;
								}
								// ########iCavitySPH remove
							}
							// ########iCavitySPH remove
						}
						if (EGS4.ispmfp == EGS4.iCavitySPH) {
							if (MEDTMP != EGS4.MEDIUM) {
								MEDTMP = EGS4.MEDIUM;
								// $SET INTERVAL GLE,GE;
								Double dbl = new Double(
										EGS4.GE1[EGS4.MEDIUM - 1] * EGS4.GLE
												+ EGS4.GE0[EGS4.MEDIUM - 1]);
								LGLE = dbl.intValue();

								// $EVALUATE DELTAP USING GMFP(GLE);
								DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

							}
							if (EGS4.IRAYLR[EGS4.IROLD - 1] == 1) {
								// $EVALUATE COHFAC USING COHE(GLE);
								COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

							} else {
								COHFAC = 1.0;
							}
						}
						eq.HOWFAR();// CALL HOWFAR;
						if (DELTAP != 0)
							PATHL = PATHL + EGS4.USTEP / (DELTAP * COHFAC);
						// "only add to pathl if not a vacuum"
						if (EGS4.IRNEW == 1)
							break;// EXIT;
						EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.U[EGS4.NP - 1];
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.V[EGS4.NP - 1];
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.W[EGS4.NP - 1];
					}
					EGS4.USTEP = DUMU;
					EGS4.X[EGS4.NP - 1] = DUMX;
					EGS4.Y[EGS4.NP - 1] = DUMY;
					EGS4.Z[EGS4.NP - 1] = DUMZ;
					EGS4.IROLD = IRODUM;
					EGS4.IRNEW = IRNDUM;
					EGS4.IR[EGS4.NP - 1] = IRDUM;
					EGS4.MEDIUM = MEDDUM;
					EGS4.IDISC = IDUM;
					if (PATHL <= 1.0E-3) {
						GWAIT = PATHL * (1. - 0.5 * PATHL);
					} else {
						GWAIT = 1. - Math.exp(-PATHL);
					}
					GWTOLD = EGS4.WT[EGS4.NP - 1];
					EGS4.WT[EGS4.NP - 1] = GWTOLD * GWAIT;

				}// if(!FASTSTEP)

				if (EGS4.ispmfp == EGS4.iCavity) {
					// :FASTSTEP:;
					EPSLON = rnno35 * GWAIT;
					if (EPSLON <= 1.0E-3) {
						if (NFTIME == 1)
							EXPMFP = EPSLON * (1. + EPSLON);
						EGS4.DPMFP = EPSLON * (1. + 0.5 * EPSLON);
					} else {
						ARG = 1. / (1. - EPSLON);
						EGS4.DPMFP = Math.log(ARG);
						if (NFTIME == 1)
							EXPMFP = EPSLON * ARG;
					}

				}

				if (EGS4.ispmfp == EGS4.iCavitySPH) {
					if (!FASTSTEP) {
						// :FASTSTEP:;
						EPSLON = rnno35 * GWAIT;
						if (EPSLON <= 1.0E-3) {
							if (NFTIME == 1)
								EXPMFP = EPSLON * (1. + EPSLON);
							EGS4.DPMFP = EPSLON * (1. + 0.5 * EPSLON);
						} else {
							ARG = 1. / (1. - EPSLON);
							EGS4.DPMFP = Math.log(ARG);
							if (NFTIME == 1)
								EXPMFP = EPSLON * ARG;
						}
					}
				}
				// :SKIPFORCING:;=======//:FASTSTEP:;
			}// else
				// }

		}

		else if ((EGS4.ispmfp == EGS4.iDose) || (EGS4.ispmfp == EGS4.iEdk)) {

			NFTIME = NFTIME + 1;

			if (EGS4.ispmfp == EGS4.iEdk) {
				if (interact_now) {
					interact_now = false;
					NFTIME = NFTIME - 1;// "reset to zero since no actual forcing"
					EGS4.DPMFP = 0.0;
					return;
				}

			}

			rnno35 = EGS4.random01();
			if (rnno35 == 0.0)
				rnno35 = 1.E-30;
			if ((IFORCE == 0) || ((NFTIME > NFMAX) || (NFTIME < NFMIN))) {
				SPMFP = -Math.log(rnno35);
				if ((CEXPTR < 1.0)
						&& (NFTIME >= NFMIN)
						&& (NFTIME <= NFMAX)
						&& (((EGS4.W[EGS4.NP - 1] > 0.0) && (CEXPTR > 0.0)) || (CEXPTR < 0.0))) {
					TEMP = CEXPTR * EGS4.W[EGS4.NP - 1];
					TEMP1 = 1.0 - TEMP;
					SPMFP = SPMFP / TEMP1;
					EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1]
							* Math.exp(-SPMFP * TEMP) / TEMP1;
				}
				EGS4.DPMFP = SPMFP;
			} else {
				if (EGS4.ispmfp == EGS4.iDose) {
					if ((do_fast_step) && (NFTIME == 1) && (CEXPTR == 0.0)) {
						FASTSTEP = true;
						// GOTO :FASTSTEP:;
					}
				}

				if (EGS4.ispmfp == EGS4.iEdk) {
					NFTIME = NFTIME + 1;
					FASTSTEP = false;// no faststep!!
				}
				// "IFPB IN COMIN USER-VARIANCE-REDUCTION"
				if (!FASTSTEP) {
					DUMU = EGS4.USTEP;
					DUMX = EGS4.X[EGS4.NP - 1];
					DUMY = EGS4.Y[EGS4.NP - 1];
					DUMZ = EGS4.Z[EGS4.NP - 1];
					IRODUM = EGS4.IROLD;
					IRNDUM = EGS4.IRNEW;
					IRDUM = EGS4.IR[EGS4.NP - 1];
					MEDDUM = EGS4.MEDIUM;
					IDUM = EGS4.IDISC;

					if (EGS4.ispmfp == EGS4.iDose)
						NRCDUM = NEWNRC;

					PATHL = 0.0;
					MEDTMP = 0;
					while (true) {
						EGS4.USTEP = EGS4.VACDST;
						EGS4.IROLD = EGS4.IR[EGS4.NP - 1];
						EGS4.MEDIUM = EGS4.MED[EGS4.IROLD - 1];
						if (EGS4.MEDIUM == 0) {// "vacuum"
							DELTAP = 0.;
						} else {
							if (MEDTMP != EGS4.MEDIUM) {
								MEDTMP = EGS4.MEDIUM;
								// $SET INTERVAL GLE,GE;
								Double dbl = new Double(
										EGS4.GE1[EGS4.MEDIUM - 1] * EGS4.GLE
												+ EGS4.GE0[EGS4.MEDIUM - 1]);
								LGLE = dbl.intValue();

								// $EVALUATE DELTAP USING GMFP(GLE);
								DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

							}
							if (EGS4.IRAYLR[EGS4.IROLD - 1] == 1) {
								// $EVALUATE COHFAC USING COHE(GLE);
								COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];
							} else {
								COHFAC = 1.0;
							}
						}

						if (EGS4.ispmfp == EGS4.iDose)
							EGS4.IRNEW = EGS4.IROLD;

						eq.HOWFAR();// CALL HOWFAR;
						if (DELTAP != 0)
							PATHL = PATHL + EGS4.USTEP / (DELTAP * COHFAC);
						// "only add to pathl if not going through vacuum"
						if (EGS4.IRNEW == 1)
							break;// EXIT;
						EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.U[EGS4.NP - 1];
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.V[EGS4.NP - 1];
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.W[EGS4.NP - 1];
					}
					EGS4.USTEP = DUMU;
					EGS4.X[EGS4.NP - 1] = DUMX;
					EGS4.Y[EGS4.NP - 1] = DUMY;
					EGS4.Z[EGS4.NP - 1] = DUMZ;
					EGS4.IROLD = IRODUM;
					EGS4.IRNEW = IRNDUM;
					EGS4.IR[EGS4.NP - 1] = IRDUM;
					EGS4.MEDIUM = MEDDUM;
					EGS4.IDISC = IDUM;
					NEWNRC = NRCDUM;
					if (CEXPTR == 0) {
						if (PATHL <= 1.0E-3) {
							GWAIT = PATHL * (1. - 0.5 * PATHL);
						} else {
							GWAIT = 1. - Math.exp(-PATHL);
						}
					}
					GWTOLD = EGS4.WT[EGS4.NP - 1];
					EGS4.WT[EGS4.NP - 1] = GWTOLD * GWAIT;
				}
				// :FASTSTEP:;
				if (CEXPTR == 0) {
					EPSLON = rnno35 * GWAIT;
					if (EPSLON <= 1.0E-3) {
						if (NFTIME == 1)
							EXPMFP = EPSLON * (1. + EPSLON);
						EGS4.DPMFP = EPSLON * (1. + 0.5 * EPSLON);
					} else {
						ARG = 1. / (1. - EPSLON);
						EGS4.DPMFP = Math.log(ARG);
						if (NFTIME == 1)
							EXPMFP = EPSLON * ARG;
					}
				} else {
					TEMP = CEXPTR * EGS4.W[EGS4.NP - 1];
					TEMP1 = 1.0 - TEMP;
					PATHLT = PATHL * TEMP1;
					if (Math.abs(PATHLT) <= 1.0E-3) {
						SPMFP = PATHL * rnno35
								* (1.0 - 0.5 * PATHLT * (1.0 - rnno35));
						GWAIT = PATHL * (1.0 - 0.5 * PATHLT)
								* Math.exp(-SPMFP * TEMP);
						EGS4.DPMFP = SPMFP;
					} else {
						GWAIT = 1. - Math.exp(-PATHLT);
						SPMFP = -Math.log(1.0 - rnno35 * GWAIT) / TEMP1;
						GWAIT = GWAIT * Math.exp(-SPMFP * TEMP) / TEMP1;
						EGS4.DPMFP = SPMFP;
					}
				}
			}

		}
		// ;"MACRO USED FOR FORCING INTERACTIONS IN THE GEOMETRY"
		// "USED BY EGSnrc FOR VARIANCE REDUCTION"
		// "NOTE FORCING OPTION HAS NOW BEEN IMPLEMENTED WHICH GIVES BOTH"
		// "THE PRIMARY PHOTON FLUENCE AND THE ELECTRON FLUENCE CORRECTLY FOR A "
		// "PHOTON BEAM!                                                   CMa"
		else if ((EGS4.ispmfp == EGS4.iFlur)) {
			rnno35 = EGS4.random01();
			if (rnno35 == 0.0)
				rnno35 = 1.E-30;
			if (IFORCE == 0) {// "photon interaction forcing option not chosen"
				EGS4.DPMFP = -Math.log(rnno35);
			} else {
				if (NP_INC[EGS4.NP - 1] == 0) {
					NFTIME = NFTIME + 1;
					// "a new photon. Check if it should be forced to interact. "
					if ((NFTIME > NFMAX) || (NFTIME < NFMIN)) {
						// "photon interaction forcing not quired in these cases"
						EGS4.DPMFP = -Math.log(rnno35);
					} else {// "force this photon to interact in the cylinder"
						NP_INC[EGS4.NP - 1] = 1; // "flag, NOW doing photon interaction forcing   "
						EGS4.NP = EGS4.NP + 1;// "create a photon which carries the remaining weight."
						// "Note that this newly created photon will be transported    "
						// "first in this new photon interaction forcing scheme.       "
						EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NP - 2];
						EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NP - 2];
						EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NP - 2];
						EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
						EGS4.E[EGS4.NP - 1] = EGS4.E[EGS4.NP - 2];
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
						EGS4.IQ[EGS4.NP - 1] = EGS4.IQ[EGS4.NP - 2];
						EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
						EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

						// "           IF(IFPB = 0 & NFTIME = 1 & MONOEN = 0 & NP = 2)[  "
						if (MAFORC == 1) {
							// "changed by CMa because all the checks don't work if we   "
							// "have a parallel beams of electrons or positrons.         "
							// "MAFORC=1 means that the particle is from the source. This"
							// "removes the bug for e-, e+ parallel beams!!!!!!          "
							MAFORC = 0; // "do it only once for the source photon          "
							EGS4.WT[EGS4.NP - 1] = 1. - EGS4.WT[EGS4.NP - 2];// "WT(NP-1) calculated in the MAIN.      "
							// "The photon on top carries the remaining weight.          "
							// GOTO
							// :FASTSTEP:;"IFPB IN COMIN USER-VARIANCE-REDUCTION    "
							FASTSTEP = true;
						}
						if (!FASTSTEP) {
							DUMU = EGS4.USTEP;
							// "        DUMX=X(NP);DUMY=Y(NP);DUMZ=Z(NP);"
							// "        commented out by CMa, no longer needed"
							IRODUM = EGS4.IROLD;
							IRNDUM = EGS4.IRNEW;
							// "        IRDUM=IR(NP);"
							MEDDUM = EGS4.MEDIUM;
							IDUM = EGS4.IDISC;
							PATHL = 0.0;
							MEDTMP = 0;
							while (true) {
								EGS4.USTEP = EGS4.VACDST;
								EGS4.IROLD = EGS4.IR[EGS4.NP - 1];
								EGS4.MEDIUM = EGS4.MED[EGS4.IROLD - 1];
								if (EGS4.MEDIUM == 0) {// "vacuum"
									DELTAP = 0.;
								} else {
									if (MEDTMP != EGS4.MEDIUM) {
										MEDTMP = EGS4.MEDIUM;
										// $SET INTERVAL GLE,GE;
										Double dbl = new Double(
												EGS4.GE1[EGS4.MEDIUM - 1]
														* EGS4.GLE
														+ EGS4.GE0[EGS4.MEDIUM - 1]);
										LGLE = dbl.intValue();
										// $EVALUATE DELTAP USING GMFP(GLE);
										DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
												* EGS4.GLE
												+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

									}
									if (EGS4.IRAYLR[EGS4.IROLD - 1] == 1) {
										// $EVALUATE COHFAC USING COHE(GLE);
										COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
												* EGS4.GLE
												+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];
									} else {
										COHFAC = 1.0;
									}
								}
								eq.HOWFAR();// CALL HOWFAR;
								// "only add to pathl if not a vacuum"
								if (DELTAP != 0)
									PATHL = PATHL + EGS4.USTEP
											/ (DELTAP * COHFAC);
								if (EGS4.IRNEW == 1)
									break;
								EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
								EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1]
										+ EGS4.USTEP * EGS4.U[EGS4.NP - 1];
								EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1]
										+ EGS4.USTEP * EGS4.V[EGS4.NP - 1];
								EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1]
										+ EGS4.USTEP * EGS4.W[EGS4.NP - 1];
							}
							EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NP - 2];
							EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NP - 2];
							EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NP - 2];
							EGS4.E[EGS4.NP - 1] = EGS4.E[EGS4.NP - 2];
							EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
							EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
							EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
							EGS4.IQ[EGS4.NP - 1] = EGS4.IQ[EGS4.NP - 2];
							EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
							EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];
							// "recover the position, etc. of this new photon"
							// "           X(NP)=DUMX;Y(NP)=DUMY;Z(NP)=DUMZ;"
							// "           no longer needed"
							EGS4.USTEP = DUMU;
							EGS4.IROLD = IRODUM;
							EGS4.IRNEW = IRNDUM;
							// "           IR(NP)=IRDUM;"
							EGS4.MEDIUM = MEDDUM;
							EGS4.IDISC = IDUM;
							// //"we now calculate the weighting factor for the photon that is"
							// //"forced to interact in the cylinder.                         "
							if (PATHL <= 1.0E-3) {
								GWAITf[EGS4.NP - 2] = PATHL
										* (1. - 0.5 * PATHL);
							} else {
								GWAITf[EGS4.NP - 2] = 1. - Math.exp(-PATHL);
							}
							GWTOLD = EGS4.WT[EGS4.NP - 2];
							EGS4.WT[EGS4.NP - 2] = GWTOLD * GWAITf[EGS4.NP - 2];
							// "we now calculate the weighting factor for the photon that will"
							// "carry the remaining weight.                                   "
							GWAITf[EGS4.NP - 1] = 1. - GWAITf[EGS4.NP - 2];
							EGS4.WT[EGS4.NP - 1] = GWTOLD * GWAITf[EGS4.NP - 1];
						}// if(!FASTSTEP)
							// :FASTSTEP:;
							// "This photon cannot interact inside the cylinder but will go   "
							// "through the cylinder so as to make the photon fluence right.  "
						EGS4.DPMFP = 1.0E30;// "set the point of interaction at infinity.        "
					}// else
				} // "end of NP_INC(NP) = 0 case"

				else {// "NP_INC(NP)=1, now transport the photon that is forced to         "
						// "interact in the cylinder"
					NP_INC[EGS4.NP - 1] = 0; // "re-set the flag"
					EPSLON = rnno35 * GWAITf[EGS4.NP - 1];
					if (EPSLON <= 1.0E-3) {
						EGS4.DPMFP = EPSLON * (1. + 0.5 * EPSLON);
					} else {
						ARG = 1. / (1. - EPSLON);
						EGS4.DPMFP = Math.log(ARG);
					}
				}// "end of NP_INC(NP) = 1 case"
			}// "end of IFORCE = 1 case"
		}
		// " Don't need any photon transport "
		else if ((EGS4.ispmfp == EGS4.iGe)) {
			EGS4.DPMFP = 0.0;
		} else if ((EGS4.ispmfp == EGS4.iSpr)) {
			NFTIME = NFTIME + 1;
			rnno35 = EGS4.random01();
			if ((IFORCE == 0) || (NFTIME > NFMAX) || (NFTIME < NFMIN)) {
				if (rnno35 == 0.0)
					rnno35 = 1.E-30;
				EGS4.DPMFP = -Math.log(rnno35);
			} else {
				if ((EGS4SrcEns.IFPB == 0) && (NFTIME == 1)
						&& (EGS4SrcEns.MONOEN == 0) && (IQINC == 0)) {
					// GOTO :FASTSTEP:;
					FASTSTEP = true;
				}
				// "IFPB IN COMIN USER-VARIANCE-REDUCTION"
				if (!FASTSTEP) {
					DUMU = EGS4.USTEP;
					DUMX = EGS4.X[EGS4.NP - 1];
					DUMY = EGS4.Y[EGS4.NP - 1];
					DUMZ = EGS4.Z[EGS4.NP - 1];
					IRODUM = EGS4.IROLD;
					IRNDUM = EGS4.IRNEW;
					IRDUM = EGS4.IR[EGS4.NP - 1];
					MEDDUM = EGS4.MEDIUM;
					IDUM = EGS4.IDISC;
					PATHL = 0.0;
					MEDTMP = 0;
					while (true) {
						EGS4.USTEP = EGS4.VACDST;
						EGS4.IROLD = EGS4.IR[EGS4.NP - 1];
						EGS4.MEDIUM = EGS4.MED[EGS4.IROLD - 1];
						if (EGS4.MEDIUM == 0) {
							DELTAP = 0.;
						} else {
							if (MEDTMP != EGS4.MEDIUM) {
								MEDTMP = EGS4.MEDIUM;
								// $SET INTERVAL GLE,GE;
								// $EVALUATE DELTAP USING GMFP(GLE);
								Double dbl = new Double(
										EGS4.GE1[EGS4.MEDIUM - 1] * EGS4.GLE
												+ EGS4.GE0[EGS4.MEDIUM - 1]);
								LGLE = dbl.intValue();

								// $EVALUATE DELTAP USING GMFP(GLE);
								DELTAP = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];

							}
							if (EGS4.IRAYLR[EGS4.IROLD - 1] == 1) {
								// $EVALUATE COHFAC USING COHE(GLE);
								COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
										* EGS4.GLE
										+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];

							} else {
								COHFAC = 1.0;
							}

						}
						eq.HOWFAR();
						// "only add to pathl if not a vacuum"
						if (DELTAP != 0)
							PATHL = PATHL + EGS4.USTEP / (DELTAP * COHFAC);
						if (EGS4.IRNEW == 1)
							break;
						EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.U[EGS4.NP - 1];
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.V[EGS4.NP - 1];
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
								* EGS4.W[EGS4.NP - 1];
					}
					EGS4.USTEP = DUMU;
					EGS4.X[EGS4.NP - 1] = DUMX;
					EGS4.Y[EGS4.NP - 1] = DUMY;
					EGS4.Z[EGS4.NP - 1] = DUMZ;
					EGS4.IROLD = IRODUM;
					EGS4.IRNEW = IRNDUM;
					EGS4.IR[EGS4.NP - 1] = IRDUM;
					EGS4.MEDIUM = MEDDUM;
					EGS4.IDISC = IDUM;
					if (PATHL <= 1.0E-3) {
						GWAIT = PATHL * (1. - 0.5 * PATHL);
					} else {
						GWAIT = 1. - Math.exp(-PATHL);
					}
					GWTOLD = EGS4.WT[EGS4.NP - 1];
					EGS4.WT[EGS4.NP - 1] = GWTOLD * GWAIT;
				}// FASTSTEP
					// FASTSTEP:;
				EPSLON = rnno35 * GWAIT;
				if (EPSLON <= 1.0E-3) {
					if (NFTIME == 1)
						EXPMFP = EPSLON * (1. + EPSLON);
					EGS4.DPMFP = EPSLON * (1. + 0.5 * EPSLON);
				} else {
					ARG = 1. / (1. - EPSLON);
					EGS4.DPMFP = Math.log(ARG);
					if (NFTIME == 1)
						EXPMFP = EPSLON * ARG;
				}
			}

		} else// default one
		{
			if (istartPhot == 0) {
				EGS4.seqStr = "SELECT-PHOTON-MFP is required BUT wrong option is used";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "Thus, the default algorithm will be used!!!";
				eq.printSequence(EGS4.seqStr);
			}

			rnno35 = EGS4.random01();
			if (rnno35 == 0.0) {
				rnno35 = 1.E-30;
			}
			EGS4.DPMFP = -Math.log(rnno35);

			istartPhot = 1;
		}

	}

	// #####################################

	// "*******************************************************************************
	/**
	 * Internally used.
	 * @param elke elke
	 * @param x x
	 * @param y y
	 * @param z z
	 * @param ir ir
	 * @param iq iq
	 */

	private static void range_rejection(double elke, double x, double y,
			double z, int ir, int iq)// ,int idisc)
	{
		// "
		// "
		// "*******************************************************************************

		// $IMPLICIT-NONE;

		// $INTEGER ir,iq,idisc;
		// $REAL elke,x,y,z;

		// COMIN/GEOM,USER/;

		double r2 = 0.0;
		double r = 0.0;
		double zl = 0.0;
		double rl = 0.0;
		double tperp = 0.0;
		double range = 0.0;
		int lelke = 0;

		EGS4.IDISC = 0;// idisc = 0; "No range discard
		if (EGS4Geom.ntrack[ir - 1] == 1) {
			return;
		} // "No additional rejection in the cavity

		if (irange_rej == EGS4.iCavity) {
			// " discard sub-threshold electrons without expensive geometry calculations "
			Double lelked = new Double(eker1 * elke + eker0);
			lelke = lelked.intValue();
			if (lelke < 1) {
				EGS4.IDISC = 50 + 49 * iq;
				return;
			}

			r2 = x * x + y * y;
			r = Math.sqrt(r2);
			rl = r - r_cavity_max;
			if (z <= z_cavity_min) {
				zl = z_cavity_min - z;
				if (r <= r_cavity_max) {
					tperp = zl;// EGS4.tperp = zl;
				} else {
					tperp = Math.sqrt(zl * zl + rl * rl);// EGS4.tperp =
															// Math.sqrt( zl*zl
															// + rl*rl );
				}
			} else if (z <= z_cavity_max) {
				if (r <= r_cavity_max) {
					return;
				}
				tperp = rl;// EGS4.tperp = rl;
			} else {
				zl = z - z_cavity_max;
				if (r <= r_cavity_max) {
					tperp = zl;// EGS4.tperp = zl;
				} else {
					tperp = Math.sqrt(zl * zl + rl * rl);// EGS4.tperp =
															// Math.sqrt( zl*zl
															// + rl*rl );
				}
			}

			range = rangerr1[lelke - 1] * elke + rangerr0[lelke - 1];
			if (range <= tperp)// ( range <= EGS4.tperp )
			{
				EGS4.IDISC = 50 + 49 * iq;
			}
		} else if (irange_rej == EGS4.iCavitySPH) {
			r2 = x * x + y * y + z * z;
			if (r2 >= r_cavity_min * r_cavity_min
					&& r2 <= r_cavity_max * r_cavity_max) {
				return;
			}
			r = Math.sqrt(r2);
			if (r < r_cavity_min) {
				rl = r_cavity_min - r;
			} else {
				rl = r - r_cavity_max;
			}
			Double lelked = new Double(eker1 * elke + eker0);
			lelke = lelked.intValue();// eker1*elke + eker0;
			range = rangerr1[lelke - 1] * elke + rangerr0[lelke - 1];
			if (range <= rl) {
				EGS4.IDISC = 50 + 49 * iq;
			}
		}
		// return; end;
	}

	// ###################################
	// "------------------------------------------------------------------------------"
	// "Macro to force the initial photon interaction at the origin. It finds second- "
	// "ary particles' starting regions according to their direction from this first "
	// "interaction using a binary search algorithm.This scheme has been used to avoid"
	// " round-off errors when the particle is at the origin .                        "
	// "------------------------------------------------------------------------------"
	// "To implement this, one has to:                                                "
	// "   - set flag interact_now = .FALSE.                                          "
	// "   - for edk-source call this macro instead of shower                        "
	// "                                                                 EMH Jan,2003"
	// "-----------------------------------------------------------------------------"
	/**
	 * Macro to force the initial photon interaction at the origin. It finds secondary particles' starting regions according to their direction from this first 
	 * interaction using a binary search algorithm.This scheme has been used to avoid round-off errors when the particle is at the origin. 
	 * Replaces shower routine from EGS4Core. 
	 */
	public static void do_photon_shower() {

		// REPLACE {$DO_PHOTON_SHOWER} WITH {
		EGS4.NP = 1;
		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		EGS4.EDEP = 0.0; // "initially no energy deposition"
		DNEARIN = 0.0;
		LATCHIN = 0;
		EGS4.IQ[0] = EGS4SrcEns.iqin;
		EGS4.E[0] = ei;
		EGS4.U[0] = EGS4SrcEns.uin;
		EGS4.V[0] = EGS4SrcEns.vin;
		EGS4.W[0] = EGS4SrcEns.win;
		// $TRANSFER PROPERTIES TO (1) FROM IN;
		EGS4.X[0] = EGS4SrcEns.xin;
		EGS4.Y[0] = EGS4SrcEns.yin;
		EGS4.Z[0] = EGS4SrcEns.zin;
		EGS4.IR[0] = EGS4SrcEns.irin;
		EGS4.WT[0] = EGS4SrcEns.wtin;
		EGS4.DNEAR[0] = DNEARIN;
		EGS4.LATCH[0] = LATCHIN;

		Mgle = Math.log(ei);
		EGS4.MEDIUM = EGS4.MED[EGS4SrcEns.irin - 1];// "uses medium from region irin, set to 2 in INPUTS"
		// $SET INTERVAL Mgle,GE;
		Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1] * Mgle
				+ EGS4.GE0[EGS4.MEDIUM - 1]);
		LMgle = dbl.intValue();

		Mcohfac = 1.0;
		if (EGS4.IRAYLR[EGS4SrcEns.irin - 1] == 1) {
			// $EVALUATE Mcohfac USING COHE(Mgle);
			Mcohfac = EGS4.COHE1[LMgle - 1][EGS4.MEDIUM - 1] * Mgle
					+ EGS4.COHE0[LMgle - 1][EGS4.MEDIUM - 1];
		}

		double RNNO36 = EGS4.random01();
		// "   GBR1=PAIR/(PAIR+COMPTON+PHOTO)=PAIR/GTOTAL"
		// $EVALUATE GBR1 USING GBR1(Mgle);
		GBR1 = EGS4.GBR11[LMgle - 1][EGS4.MEDIUM - 1] * Mgle
				+ EGS4.GBR10[LMgle - 1][EGS4.MEDIUM - 1];

		if ((RNNO36 <= GBR1 * Mcohfac) && (EGS4.E[EGS4.NP - 1] > EGS4.RMT2)) {// "IT WAS A PAIR PRODUCTION"
			EGS4Core.PAIR();
			// "There should be NPold ... NP particles in the stack"
			for (int i = EGS4.NPold; i <= EGS4.NP; i++) {
				EGS4.IR[i - 1] = ibsearchrev(EGS4.W[i - 1], EGS4Geom.NC,
						EGS4Geom.cosalp) + 1;// remember->1 biased!!
			}
			// $AUSCALL($PAIRAUSA);
			IARG = EGS4.$PAIRAUSA;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}

		} else {// "GBR2=(PAIR+COMPTON)/GTOTAL"
				// $EVALUATE GBR2 USING GBR2(Mgle);
			GBR2 = EGS4.GBR21[LMgle - 1][EGS4.MEDIUM - 1] * Mgle
					+ EGS4.GBR20[LMgle - 1][EGS4.MEDIUM - 1];

			if (RNNO36 < GBR2 * Mcohfac) {
				if (idopp > 0) {
					COMPT_NO_DOPP();
				} else {
					EGS4Core.COMPT();
				}
				// "There should be NPold ... NP particles in the stack"
				for (int i = EGS4.NPold; i <= EGS4.NP; i++)// DO i = NPold,NP[
				{
					EGS4.IR[i - 1] = ibsearchrev(EGS4.W[i - 1], EGS4Geom.NC,
							EGS4Geom.cosalp) + 1;
				}
				// $AUSCALL($COMPAUSA);
				IARG = EGS4.$COMPAUSA;
				if (EGS4.iausfl[IARG] != 0) {
					eq.AUSGAB(IARG);
				}

			} else {
				if (RNNO36 < Mcohfac) {
					EGS4Core.PHOTO();
					// "There should be NPold ... NP particles in the stack"
					for (int i = EGS4.NPold; i <= EGS4.NP; i++)// DO i =
																// NPold,NP[
					{
						EGS4.IR[i - 1] = ibsearchrev(EGS4.W[i - 1],
								EGS4Geom.NC, EGS4Geom.cosalp) + 1;
						// IR(i) = ibsearchrev(W(i),nc,cosalp)+1;
					}
					// $AUSCALL($PHOTOAUSA);
					IARG = EGS4.$PHOTOAUSA;
					if (EGS4.iausfl[IARG] != 0) {
						eq.AUSGAB(IARG);
					}

				} else {
					MY_RAYLEIGH_SCATTERING();
				}
			}
		}

		while (EGS4.NP > 0) {
			if (EGS4.IQ[EGS4.NP - 1] == 0) {
				EGS4Core.PHOTON();
			}// (ircode); }
			else {
				EGS4Core.ELECTR();
			}// (ircode); ]
		}

	}

	// C*******************************************************************************
	// C* *
	// C* Function ibsearchrev(a, nsh, b) *
	// C* *
	// C* binary search for an element l of array b such that *
	// C* b[l] => a > b[l+1], *
	// C* *
	// C* Note: array must be decreasingly monotone *
	// C* *
	// C*******************************************************************************
	/**
	 * Binary search for an element l of array b such that b[l] .GE. a .G. b[l+1]. Array must be decreasingly monotone. Used by do_photon_shower, 
	 * My_RAYLEIGH_SCATTERING...
	 * @param a a
	 * @param nsh nsh
	 * @param b b
	 * @return the result
	 */
	public static int ibsearchrev(double a, int nsh, double[] b) {
		// real*4 a, b(*)
		int min = 0;
		int max = 0;
		int help = 0;// nsh
		double x = 0.0;
		// c array b represents array of energy or of angle distribution
		min = 1;// 0;
		max = nsh;// -1;
		x = a;
		do // while ( min.lt.max-1 )
		{
			help = (max + min) / 2; // !bitweise ganzzahlige Division durch 2
			if (b[help - 1] >= x)// if ( b(help).ge.x) then
				min = help;
			else
				max = help;
			// endif
			// enddo
		} while (min < max - 1);
		int ibsearchrev = min;
		return ibsearchrev;
		// returned value must be substracted by 1 for 0 biased
		// arrays!!!!!!!!!!!!!!!!!!!
	}// end

	// "******************************************************************"
	// "                                                          NRCC    "
	/**
	 * Used by do_photon_shower. Subroutine for sampling incoherent (Compton) scattering. If the flag ibcmp(region) is zero, Klein-Nishina is used.  
	 * Otherwise scattering is modelled in the impuls approximation (see R.Ribberfors and K.F.Berggren, Phys.Rev.A26 (1982) 3325). 
	 * As the total cross section from PEGS4 is not modified (and thus calculated using Klein-Nishina), all rejections leed to an 
	 * unscattered photon and a zero energy electron. If a K,L1,L2,L3,M or N vacancy is created, the subsequent atomic relaxation is treated in RELAX. This has as a 
	 * consequence that more than one particle can be created as a result of an incoherent scattering.
	 */
	protected static void COMPT_NO_DOPP() {
		// "                                VERSION 1.00  --  12 JAN 1999     "
		// "******************************************************************"
		// "                                                                  "
		// "   Subroutine for sampling incoherent (Compton) scattering        "
		// "   If the flag ibcmp(region) is zero, Klein-Nishina is used.      "
		// "   Otherwise scattering is modelled in the impuls approximation   "
		// "    (see R.Ribberfors and K.F.Berggren, Phys.Rev.A26 (1982) 3325) "
		// "   As the total cross section from PEGS4 is not modified (and thus"
		// "   calculated using Klein-Nishina), all rejections leed to an     "
		// "   unscattered photon and a zero energy electron.                 "
		// "   If a K,L1,L2,L3,M or N vacancy is created, the subsequent      "
		// "   atomic relaxation is treated in RELAX. This has as a           "
		// "   consequence that more than one particle can be created as a    "
		// "   result of an incoherent scattering. The user should therefore  "
		// "   check their user codes for possible inconsistencies.           "
		// "                                                                  "
		// "   I.Kawrakow, January 1999                                       "
		// "******************************************************************"

		// ; Copyright NRC;
		// $IMPLICIT-NONE;

		// $COMIN-COMPT; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/COMPTON-DATA,DEBUG,STACK,THRESH,UPHIOT,USEFUL,RANDOM/;

		// $DEFINE-LOCAL-VARIABLES-COMPT;
		double PEIG = 0.0;// "precise energy of incident photon"
		double PESG = 0.0;// "precise energy of scattered photon"
		double PESE = 0.0;// "precise total energy of compton electron"
		double ko = 0.0;// "energy of incident photon in units of RM"
		double broi = 0.0;// "1+2*ko"
		double broi2 = 0.0;// "broi*broi"
		double bro = 0.0;// "1/broi"
		double bro1 = 0.0;// "1-bro"
		double alph1 = 0.0;// "probability for the 1/BR part"
		double alph2 = 0.0;// "probability for the BR part"
		double alpha = 0.0;// "alpha1/(alph1+alph2)"
		double rnno15 = 0.0;
		double rnno16 = 0.0;
		double rnno17 = 0.0;
		//double rnno18 = 0.0;
		double rnno19 = 0.0; // "random numbers"
		double br = 0.0;// "scattered photon energy fraction"
		double temp = 0.0;// "aux. variable for polar angle calculation"
		double rejf3 = 0.0;// "rejection function"
		double rejmax = 0.0;// "max. of rejf3 in thge case of uniform sampling"
		double Uj = 0.0;// "binding energy of the selected shell"
		double Jo = 0.0;// "the Compton profile parameter"
		//double br2 = 0.0;// "br*br"
		double fpz = 0.0;
		double fpz1 = 0.0;// "used for limited pz-range rejection"
		double qc = 0.0;// "momentum transfer corresponding to the Compton line energy"
		double qc2 = 0.0;// "qc squared"
		double af = 0.0;// "for calculating F"
		double Fmax = 0.0;// "maximum of F"
		//double frej = 0.0;// "used for F-rejection"
		double eta_incoh = 0.0;
		//double eta = 0.0;// "random numbers"
		double aux = 0.0;
		//double aux1 = 0.0;
		//double aux2 = 0.0;
		double aux3 = 0.0;
		double aux4 = 0.0;// "aux. variables"
		double pzmax = 0.0;// "max. possible z-component of the initial electron momentum"
		//double pz = 0.0;// "initial electron momentum projection"
		//double pz2 = 0.0;// "pz*pz"
		double rnno_RR = 0.0;// "for playing Russian Roulette"

		int irl = 0;// "local region number"
		int i = 0;// "loop variable for shell sampling (and then shell sampled)"
		int j = 0;// "pointer to the shell in the shell data list"
		//int iarg = 0;// "argument for eq.AUSGAB call"
		int ip = 0;// "a loop variable"
		// }

		boolean finishedComptonSampling = false;
		boolean RETRY_PZ = false;
		boolean RESAMPLEb = false;

		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		PEIG = EGS4.E[EGS4.NP - 1];// "PRECISE ENERGY OF INCIDENT GAMMA"
		ko = PEIG / EGS4.RM;// "Gamma energy in units of electron rest energy"
		broi = 1.0 + 2.0 * ko; // "Needed for scattering angle sampling"

		irl = EGS4.IR[EGS4.NP - 1];
		if (EGS4.ibcmp[irl - 1] == 1) {// "User wants to take into account binding effects"
										// "=>first sample the shell and see whether an    "
										// "  interaction is possible                      "
			rnno17 = EGS4.random01();
			// DO i=1,n_shell(medium) [
			for (i = 1; i <= EGS4.n_shell[EGS4.MEDIUM - 1]; i++) {
				rnno17 = rnno17 - EGS4.eno_array[i - 1][EGS4.MEDIUM - 1];// (i,medium);
				if (rnno17 <= 0)
					break;// EXIT;
			}
			// "j is the shell number in the data list"
			j = EGS4.shell_array[i - 1][EGS4.MEDIUM - 1];// (i,medium);
			// "Uj is the binding energy in units of rm"
			Uj = EGS4.be_array[j - 1];// (j);

			// " Binding energy rejection "
			if (ko <= Uj) {
				// goto :INTERACTION-REJECTED:;
				// :INTERACTION-REJECTED:
				// " Create here a zero energy electron if required (check user codes) "
				return;
			}
			// "Jo is the Compton profile parameter"
			Jo = EGS4.Jo_array[j - 1];// (j);

		}

		// " We always sample the scattering angle from Klein-Nishina"
		while (true)// added
		{// added
			RESAMPLEb = false;
			// :RESAMPLE:
			if (ko > 2.0) {// "At high energies the original EGS4 method is most efficient"
				broi2 = broi * broi;
				alph1 = Math.log(broi);
				alph2 = ko * (broi + 1.0) / broi2;
				alpha = alph1 / (alph1 + alph2);
				// LOOP [
				while (true) {
					rnno15 = EGS4.random01();
					rnno16 = EGS4.random01();
					if (rnno15 < alpha) {// "Use 1/br part"
						br = Math.exp(alph1 * rnno16) / broi;
					} else {// "Use the br part."
						br = Math.sqrt(rnno16 + (1.0 - rnno16) / broi2);
					}
					temp = (1.0 - br) / ko / br;
					EGS4.SINTHE = Math.max(0., temp * (2.0 - temp));
					rejf3 = 1.0 - br * EGS4.SINTHE / (1.0 + br * br);
					rnno19 = EGS4.random01();
					// ------------------------------------------
					if (rnno19 <= rejf3) {
						break;
					}
					// ---------------------------------------------
				}// UNTIL rnno19.le.rejf3;
			} else {// "At low energies it is faster to sample br uniformely"
				bro = 1. / broi;
				bro1 = 1.0 - bro;
				rejmax = broi + bro;
				while (true) {
					rnno15 = EGS4.random01();
					rnno16 = EGS4.random01();
					br = bro + bro1 * rnno15;
					temp = (1.0 - br) / ko / br;
					EGS4.SINTHE = Math.max(0., temp * (2.0 - temp));
					rejf3 = (br + 1. / br - EGS4.SINTHE) / rejmax;
					// ------------------------------------------
					if (rnno16 <= rejf3) {
						break;
					}
					// --------------------------------------
				}// UNTIL rnno16.le.rejf3;
			}
			if ((br < 1. / broi) || (br > 1.0))
			// if((br >= 1./broi) && (br <= 1.0))
			{
				// if( (br < 0.99999/broi) || (br > 1.00001 ))
				// {
				// write(6,*) ' sampled br outside of allowed range!
				// ',ko,1./broi,br;
				// }
				// goto :RESAMPLE: ;
				RESAMPLEb = true;
				// break;
			}

			if (!RESAMPLEb)
				break;

		}// main loop added
		EGS4.COSTHE = 1.0 - temp;// costhe = 1 - temp;

		if (EGS4.ibcmp[irl - 1] == 0) {// "User wants to use Klein-Nishina, so we are done"
			Uj = 0.0;
			// goto :FINISHED-COMPTON-SAMPLING:;
			finishedComptonSampling = true;
		}

		if (!finishedComptonSampling) {
			// " Check for rejection due to the limited range of pzmax "
			//br2 = br * br;
			aux = ko * (ko - Uj) * temp;
			pzmax = (aux - Uj) / Math.sqrt(2.0 * aux + Uj * Uj);// pirs 2.2.17
			if (pzmax <= -1.0) {
				// goto :INTERACTION-REJECTED:;
				// :INTERACTION-REJECTED:
				// " Create here a zero energy electron if required (check user codes) "
				return;
			}
			qc2 = 1.0 + br * br - 2.0 * br * EGS4.COSTHE;
			qc = Math.sqrt(qc2);// ~pirs 2.2.4

			if (pzmax > 1.0) {
				pzmax = 1.0;
				af = 0.0;
				Fmax = 1.0;
				fpz = 1.0;
				// goto :RETRY-PZ:;
				RETRY_PZ = true;
			}
			if (!RETRY_PZ) {
				aux3 = 1.0 + 2.0 * Jo * Math.abs(pzmax);
				aux4 = 0.5 * (1 - aux3 * aux3);// =-b
				fpz = 0.5 * Math.exp(aux4);// si
				af = qc * (1.0 + br * (br - EGS4.COSTHE) / qc2);// 2.2.21;22,23

				if (af < 0.0) {
					if (pzmax > 0.0)
						fpz = 1.0 - fpz;
					eta_incoh = EGS4.random01();
					if (eta_incoh > fpz) {
						// goto :INTERACTION-REJECTED:;
						return;
					}
					af = 0.0;
					Fmax = 1.0;
					// goto :RETRY-PZ:;
					RETRY_PZ = true;
				}
			}// if(!RETRY_PZ)

			if (!RETRY_PZ) {
				if (pzmax < -0.15) {
					Fmax = 1.0 - af * 0.15;// pirs 2.2.20
					fpz1 = fpz * Fmax;
				} else if (pzmax < 0.15) {
					Fmax = 1.0 + af * pzmax;
					aux3 = 1.0 / (1 + 0.33267252734 * aux3);
					// "0.33267252734 is p/sqrt(2), p is the parameter from Eq. 7.1.25"
					// "of Abramowitz and Stegun, needed for approximating Erf        "
					aux4 = fpz
							* aux3
							* (0.3480242 + aux3
									* (-0.0958798 + aux3 * 0.7478556))
							+ EGS4.erfJo_array[j - 1];// (j);
					if (pzmax > 0.0) {
						// "fpz1 = 1 - Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4;"
						// "missing factor 1/2 in the above found by Cerneliu Costescu"
						// "0.62665706866 is sqrt(Pi/8)"
						fpz1 = 1.0 - Fmax * fpz - 0.62665706866 * af
								/ EGS4.Jo_array[j - 1] * aux4;
						fpz = 1.0 - fpz;
					} else {
						// "fpz1 = Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4;"
						// "missing factor 1/2 in the above found by Cerneliu Costescu"
						// "0.62665706866 is sqrt(Pi/8)"
						fpz1 = Fmax * fpz - 0.62665706866 * af
								/ EGS4.Jo_array[j - 1] * aux4;
					}
				} else {
					Fmax = 1.0 + af * 0.15;
					fpz1 = 1.0 - Fmax * fpz;
					fpz = 1.0 - fpz;
				}
				eta_incoh = EGS4.random01();
				if (eta_incoh > fpz1)
					return;// goto :INTERACTION-REJECTED:;
			}// if(!RETRY_PZ)
				// "At this point, all rejections are handled, now we need to sample pz "
				// "between -1 and pzmax using the Compton profile of the selected shell"
				// "and F(pz,cos(theta)) as a rejection function                        "

			// :RETRY-PZ:;
		}// :FINISHED-COMPTON-SAMPLING:
			// :FINISHED-COMPTON-SAMPLING:

		br = 1. + ko * (1 - EGS4.COSTHE);
		br = 1. / br;

		PESG = br * PEIG;
		// PESE = PEIG - PESG - Uj + EGS4.PRM;
		PESE = PEIG - PESG + EGS4.PRM;
		EGS4.SINTHE = Math.sqrt(EGS4.SINTHE);
		EGS4Core.UPHI(2, 1);// //$SELECT-AZIMUTHAL-ANGLE and OLD-PARTICLE:
		EGS4.E[EGS4.NP - 1] = PESG;
		aux = 1.0 + br * br - 2.0 * br * EGS4.COSTHE;
		if (aux > 1.e-8) {
			EGS4.COSTHE = (1.0 - br * EGS4.COSTHE) / Math.sqrt(aux);
			EGS4.SINTHE = (1.0 - EGS4.COSTHE) * (1.0 + EGS4.COSTHE);
			if (EGS4.SINTHE > 0.0) {
				EGS4.SINTHE = -Math.sqrt(EGS4.SINTHE);
			} else {
				EGS4.SINTHE = 0.0;
			}
		} else {
			EGS4.COSTHE = 0.0;
			EGS4.SINTHE = -1.0;
		}
		EGS4.NP = EGS4.NP + 1;// inc NP
		// $CHECK-STACK(np,'COMPT');
		// REPLACE {$CHECK-STACK(#,#);} WITH {;
		if (EGS4.NP > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "COMPT"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}
		EGS4Core.UPHI(3, 2);// NEW-PARTICLE
		EGS4.E[EGS4.NP - 1] = PESE;
		EGS4.IQ[EGS4.NP - 1] = -1;

		// / *
		// IF( ibcmp(irl) = 1 ) [
		//
		// " Shell vacancy "
		// IF( Uj > 1e-3 ) [
		// edep = 0;
		// call relax(Uj,shn_array(j),iz_array(j));
		// "relax will put all particles with energies above ecut,pcut on the "
		// "stack, the remaining energy will be scored in edep and deposited  "
		// "localy (via the call to ausgab below)                             "
		// ]
		// ELSE [ edep = Uj; ]
		// $AUSCALL($PHOTXAUS); "generates IARG = 4 call"
		//
		// ]
		// * /

		// " Now play Russian Roulette with resulting electrons if the user asked for it"
		// $PLAY RUSSIAN ROULETTE WITH ELECTRONS FROM NPold+1; "TO NP"
		EGS4.i_survived_RR = 0; // "flag all survive"
		if (EGS4.i_play_RR == 1) {
			if (EGS4.prob_RR <= 0.0) {
				if (EGS4.n_RR_warning < EGS4.$MAX_RR_WARNING) {
					EGS4.n_RR_warning = EGS4.n_RR_warning + 1;
					// OUTPUT prob_RR;
					EGS4.seqStr = "**** Warning, attempt to play Russian Roulette with prob_RR<0! "
							+ EGS4.prob_RR;// +" \n";
					if (EGS4.iprint > 2)
						eq.printSequence(EGS4.seqStr);

				}
			} else {
				ip = EGS4.NPold + 1;
				while (true) {// "handle all particles from p1 to np"
					if (EGS4.IQ[ip - 1] != 0) {// "i.e. charged particles"
						rnno_RR = EGS4.random01();
						if (rnno_RR < EGS4.prob_RR) {// "particle survives"
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1] / EGS4.prob_RR;
							ip = ip + 1; // "increase local pointer"
						} else {// "particle killed"
							EGS4.i_survived_RR = EGS4.i_survived_RR + 1;
							if (ip < EGS4.NP) {// "=>replace it with last particle"
												// "on stack"
								EGS4.E[ip - 1] = EGS4.E[EGS4.NP - 1];
								EGS4.IQ[ip - 1] = EGS4.IQ[EGS4.NP - 1];
								EGS4.WT[ip - 1] = EGS4.WT[EGS4.NP - 1];
								EGS4.U[ip - 1] = EGS4.U[EGS4.NP - 1];
								EGS4.V[ip - 1] = EGS4.V[EGS4.NP - 1];
								EGS4.W[ip - 1] = EGS4.W[EGS4.NP - 1];
							}
							EGS4.NP = EGS4.NP - 1; // "reduce stack by one=> particle gone"
						}// "end of kill particle block"
					} else {// "this is a photon, leave it. Change pointer"
						ip = ip + 1;
					}

					if (ip > EGS4.NP)
						break;
				} // UNTIL (ip > np);
					// "loops until either np is decreased to ip, or ip increased to np"
				if (EGS4.NP == 0) {// " we need at least one particle on the stack "
									// " so that the transport routines can exit properly"
					EGS4.NP = 1;
					EGS4.E[EGS4.NP - 1] = 0.0;
					EGS4.IQ[EGS4.NP - 1] = 0;
					EGS4.WT[EGS4.NP - 1] = 0.0;
				}
			} // "end of russian roulette block"
		} // "end of flag set block"

		return;

		// :INTERACTION-REJECTED:
		// " Create here a zero energy electron if required (check user codes) "
		// return;
		// end;

		// ;Copyright NRC;
	}

	/**
	 * Used by do_photon_shower.
	 */
	protected static void MY_RAYLEIGH_SCATTERING() {
		boolean SAMPLING_LOOPb = false;

		if (EGS4.IRAYLR[EGS4SrcEns.irin - 1] == 1) {
			// $AUSCALL($RAYLAUSB);
			IARG = EGS4.$RAYLAUSB;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}

			EGS4.NPold = EGS4.NP;

			SAMPLING_LOOP: while (true) {
				SAMPLING_LOOPb = false;
				mXXX = EGS4.random01();
				// $SET INTERVAL mXXX,RCO;
				Double dbl1 = new Double(EGS4.RCO1[EGS4.MEDIUM - 1] * mXXX
						+ EGS4.RCO0[EGS4.MEDIUM - 1]);
				LmXXX = dbl1.intValue();

				// $EVALUATE mX2 USING RSCT(mXXX);
				mX2 = EGS4.RSCT1[LmXXX - 1][EGS4.MEDIUM - 1] * mXXX
						+ EGS4.RSCT0[LmXXX - 1][EGS4.MEDIUM - 1];

				mQ2 = mX2 * EGS4.RMSQ / (20.60744 * 20.60744);
				EGS4.COSTHE = 1. - mQ2
						/ (2. * EGS4.E[EGS4.NP - 1] * EGS4.E[EGS4.NP - 1]);
				if (Math.abs(EGS4.COSTHE) > 1.0) {
					// GO TO :SAMPLING-LOOP:;
					SAMPLING_LOOPb = true;
				}
				if (SAMPLING_LOOPb) {
					mCSQTHE = EGS4.COSTHE * EGS4.COSTHE;
					mREJF = (1.0 + mCSQTHE) / 2.0;
					mRNNORJ = EGS4.random01();
					// ------------------------------------
					if (mRNNORJ <= mREJF) {
						break SAMPLING_LOOP;
					}

				}

			}// UNTIL (mRNNORJ <= mREJF);
			EGS4.SINTHE = Math.sqrt(1.0 - mCSQTHE);
			EGS4Core.UPHI(2, 1);
			EGS4.IR[EGS4.NP - 1] = ibsearchrev(EGS4.W[EGS4.NP - 1],
					EGS4Geom.NC, EGS4Geom.cosalp) + 1;
			// $AUSCALL($RAYLAUSA);
			IARG = EGS4.$RAYLAUSA;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}

		}

	}

	/**
	 * Used by some user application (e.g. cavrz). This sets up the needed parameters to do a global range rejection on 
	 * electrons which cannot reach the cavity. For this range rejection, the cavity is any region of the same material as the cavity as 
	 * designated by the user. The procedure followed is to determine the longest range for any material outside the cavity and then use this value in the range rejection. 
	 * Thus, it would be in-efficient if there were any air regions outside the cavity and it is more efficient to define the cavity as being 
	 * bigger than it actually is.
	 */
	public static void initialize_range_rejection() {
		// " This sets up the needed parameters to do a global range rejection
		// on
		// " electrons which cannot reach the cavity. For this range rejection,
		// " the cavity is any region of the same material as the cavity as
		// designated
		// " by the user.
		// " The procedure followed is to determine the longest range for any
		// material
		// " outside the cavity and then use this value in the range rejection.
		// " Thus, it would be in-efficient if there were any air regions
		// outside
		// " the cavity and it is more efficient to define the cavity as being
		// " bigger than it actually is.
		// "
		// " Note that the EGSnrc region by region range rejection is performed
		// " prior to this more approximate global range rejection.

		int irl = 0;
		int medium = 0;
		int cavity_medium = 0;
		int ix = 0;
		int iz = 0;
		int neke = 0;
		int lelke = 0;
		int lelec = 0;
		int qel = 0;
		//int IC = 0;
		double elke = 0.0;
		double elke_min = 0.0;
		double elke_max = 0.0;
		double rhof = 0.0;
		double rhof_min = 0.0;
		double range = 0.0;
		double max_range = 0.0;
		double si = 0.0;
		double sip = 0.0;
		double eke = 0.0;
		double ekei = 0.0;
		double elkei = 0.0;
		double elktmp = 0.0;
		double dedxmid = 0.0;
		double aux = 0.0;
		double fedep = 0.0;
		int lelktmp = 0;
		int j = 0;

		if (irange_rej == EGS4.iCavity || irange_rej == EGS4.iCavitySPH) {

			EGS4.seqStr = "  initialize_range_rejection ....";
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);

			lelec = -1;
			qel = 0;
			cavity_medium = 0;
			// "Determin first the cavity medium and check for possible errors"
			for (irl = 2; irl <= EGS4Geom.nreg; irl++) {
				medium = EGS4.MED[irl - 1];
				if (EGS4Geom.ntrack[irl - 1] == 1) {// "set 1 for all regions user specifies as in cavity"
					if (cavity_medium == 0) {// "only for first region"
						cavity_medium = medium;
						EGS4.seqStr = "  Cavity medium is med # "
								+ cavity_medium;
						if (EGS4.iprint > 1)// at the start===HATCH
							eq.printSequence(EGS4.seqStr);

						// write(6,*) ' Cavity medium is med # ',cavity_medium;
					} else {
						if (medium != cavity_medium) {
							EGS4.seqStr = " Warning: cavity composed of several media: "
									+ medium;
							if (EGS4.iprint > 1)// at the start===HATCH
								eq.printSequence(EGS4.seqStr);

							EGS4.seqStr = " Turning off range rejection! ";
							if (EGS4.iprint > 1)// at the start===HATCH
								eq.printSequence(EGS4.seqStr);

							irejct = 0;
							return;
						}
					}
				}// "end ntrack = 1 block"
			}

			// " Now determine smallest cylinder enclosing the cavity "
			// " (non-cavity regions composed of cavity material are counted as cavity)"
			if (irange_rej == EGS4.iCavity) {
				z_cavity_max = -1.e10;
				z_cavity_min = 1.e10;
			}
			if (irange_rej == EGS4.iCavitySPH) {
				r_cavity_min = 1.e10;
			}
			r_cavity_max = -1.e10;
			for (irl = 2; irl <= EGS4Geom.nreg; irl++)// DO irl = 2,nreg
			{
				medium = EGS4.MED[irl - 1];// med(irl);
				if (EGS4Geom.ntrack[irl - 1] == 1 || medium == cavity_medium) {
					if (irange_rej == EGS4.iCavity) {
						// $GET-IX-IZ(irl);
						// REPLACE {$GET-IX-IZ(#);} WITH {;IX=({P1}-2)/NZ+1;
						// IZ={P1}-1-NZ*(IX-1);}
						ix = EGS4Geom.GET_IX(irl);// ix=(irl-2)/NZ+1;
						iz = EGS4Geom.GET_IZC(irl);// iz=irl-1-NZ*(ix-1);
						if (z_cavity_min > EGS4Geom.ZPLANE[iz - 1]) {
							z_cavity_min = EGS4Geom.ZPLANE[iz - 1];
						}
						if (z_cavity_max < EGS4Geom.ZPLANE[iz]) {
							z_cavity_max = EGS4Geom.ZPLANE[iz];
						}
						// rcyl=0 biased
						if (r_cavity_max < EGS4Geom.RCYL[ix]) {
							r_cavity_max = EGS4Geom.RCYL[ix];
						}
					}
					if (irange_rej == EGS4.iCavitySPH) {
						// $GET-IX-IC(irl);
						// REPLACE {$GET-IX-IC(#);} WITH {;IX=({P1}-2)/NC+1;
						// IC={P1}-1-NC*(IX-1);}
						ix = EGS4Geom.GET_IX(irl);// ix=(irl-2)/NZ+1;
						//IC = EGS4Geom.GET_IZC(irl);// iz=irl-1-NZ*(ix-1);
						if (EGS4Geom.RSPH[ix - 1] < r_cavity_min) {
							r_cavity_min = EGS4Geom.RSPH[ix - 1];
						}
						if (EGS4Geom.RSPH[ix] > r_cavity_max) {
							r_cavity_max = EGS4Geom.RSPH[ix];
						}
					}
				}
			}
			if (irange_rej == EGS4.iCavity)
				EGS4.seqStr = "  The smallest cylinder enclosing the cavity is: ";
			if (irange_rej == EGS4.iCavitySPH)
				EGS4.seqStr = "  The smallest spheres enclosing the cavity are defined by: ";
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);
			if (irange_rej == EGS4.iCavity)
				EGS4.seqStr = z_cavity_min + " < z < " + z_cavity_max;
			if (irange_rej == EGS4.iCavitySPH)
				EGS4.seqStr = r_cavity_min + " < r < " + r_cavity_max;
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "  r < " + r_cavity_max;
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);

			if (irange_rej == EGS4.iCavitySPH) {
				if (EGS4.NMED == 1) {
					EGS4.seqStr = "  There is only one medium in the simulation => "
							+ "no additional range-rejection can be done!";
					if (EGS4.iprint > 1)// at the start===HATCH
						eq.printSequence(EGS4.seqStr);

					return;
				}

			}
			// " Now prepare an array for range calculation in the medium with the "
			// " smallest stopping power (excluding the cavity material)            "
			// " First determine energy window. We assume that the PEGS data has been "
			// " checked to cover the necessary energy range                          "
			elke_min = -1.e10;
			elke_max = 1.e10;
			for (medium = 1; medium <= EGS4.NMED; medium++) {
				if (medium != cavity_medium) {
					neke = EGS4.MEKE[medium - 1];
					elke = (1.0 - EGS4.EKE0[medium - 1])
							/ EGS4.EKE1[medium - 1];
					if (elke > elke_min) {
						elke_min = elke;
					}
					elke = (neke - EGS4.EKE0[medium - 1])
							/ EGS4.EKE1[medium - 1];
					if (elke < elke_max) {
						elke_max = elke;
					}
				}
			}

			EGS4.seqStr = " Range rejection data will be calculated for ";
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);
			double d1 = Math.exp(elke_min);
			double d2 = Math.exp(elke_max);
			EGS4.seqStr = d1 + " < E < " + d2;
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);

			// "Then determine the minimum density scaling factor (rhof)"
			rhof_min = 1.e10;
			for (irl = 2; irl <= EGS4Geom.nreg; irl++)// DO irl = 2,nreg
			{
				medium = EGS4.MED[irl - 1];
				// $SET-RHOF;
				rhof = EGS4.RHOR[irl - 1] / EGS4.RHO[medium - 1];
				if (rhof < rhof_min) {
					rhof_min = rhof;
				}
			}
			EGS4.seqStr = " rhof_min = " + rhof_min;
			if (EGS4.iprint > 1)// at the start===HATCH
				eq.printSequence(EGS4.seqStr);

			rhof = rhof_min;
			eker1 = (EGS4.$MXRANGE - 1.0) / (elke_max - elke_min);
			eker0 = EGS4.$MXRANGE - eker1 * elke_max;
			// "write(6,*) ' eker1 eker0: ',eker1,eker0;
			// "write(6,*);

			// " Calculate maxium range for 1st table energy "
			elke = (1.0 - eker0) / eker1;
			eke = Math.exp(elke);
			max_range = -1.e10;
			for (medium = 1; medium <= EGS4.NMED; medium++)// DO medium=1,nmed
			{
				if (medium != cavity_medium) {
					// $SET INTERVAL elke,eke;
					Double dbl = new Double(EGS4.EKE1[medium - 1] * elke
							+ EGS4.EKE0[medium - 1]);
					lelke = dbl.intValue();

					if (lelke < 1) {
						lelke = 1;
					}
					// $COMPUTE-RANGE;
					ekei = EGS4.E_array[lelke - 1][medium - 1];// (lelke,medium);
					elkei = (lelke - EGS4.EKE0[medium - 1])
							/ EGS4.EKE1[medium - 1];
					// $COMPUTE-DRANGE(eke,ekei,lelke,elke,elkei,range);

					fedep = 1.0 - ekei / eke;// 1 - {P2}/{P1};
					// elktmp =
					// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
					elktmp = 0.5 * (elke + elkei + 0.25 * fedep * fedep
							* (1.0 + fedep * (1.0 + 0.875 * fedep)));
					// " the above evaluates the logarithm of the midpoint energy"
					// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
					lelktmp = lelke;// lelktmp = {P3};
					if (lelec < 0) {
						// $EVALUATE dedxmid USING ededx(elktmp);
						dedxmid = EGS4.EDEDX1[lelktmp - 1][medium - 1] * elktmp
								+ EGS4.EDEDX0[lelktmp - 1][medium - 1];
						aux = EGS4.EDEDX1[lelktmp - 1][medium - 1] / dedxmid; // (lelktmp,medium)/dedxmid;
					} else {
						// $EVALUATE dedxmid USING pdedx(elktmp);
						dedxmid = EGS4.PDEDX1[lelktmp - 1][medium - 1] * elktmp
								+ EGS4.PDEDX0[lelktmp - 1][medium - 1];
						aux = EGS4.PDEDX1[lelktmp - 1][medium - 1] / dedxmid; // (lelktmp,medium)/dedxmid;
					}
					aux = aux * (1.0 + 2.0 * aux) * (fedep / (2.0 - fedep))
							* (fedep / (2.0 - fedep)) / 6.0;
					// {P6} = fedep*{P1}/dedxmid*(1+aux);
					range = fedep * eke / dedxmid * (1.0 + aux);
					// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"
					range = (range + EGS4.range_ep[qel][lelke - 1][medium - 1])
							/ rhof;

					if (range > max_range) {
						max_range = range;
					}
				}
			}
			si = 1.02 * max_range;

			for (j = 1; j <= EGS4.$MXRANGE - 1; j++) {
				elke = (j + 1.0 - eker0) / eker1;
				eke = Math.exp(elke);
				max_range = -1.e10;
				for (medium = 1; medium <= EGS4.NMED; medium++)// DO
																// medium=1,nmed
				{
					if (medium != cavity_medium) {
						// $SET INTERVAL elke,eke;
						Double dbl = new Double(EGS4.EKE1[medium - 1] * elke
								+ EGS4.EKE0[medium - 1]);
						lelke = dbl.intValue();

						// $COMPUTE-RANGE;
						ekei = EGS4.E_array[lelke - 1][medium - 1];// (lelke,medium);
						elkei = (lelke - EGS4.EKE0[medium - 1])
								/ EGS4.EKE1[medium - 1];
						// $COMPUTE-DRANGE(eke,ekei,lelke,elke,elkei,range);

						fedep = 1.0 - ekei / eke;// 1 - {P2}/{P1};
						// elktmp =
						// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
						elktmp = 0.5 * (elke + elkei + 0.25 * fedep * fedep
								* (1.0 + fedep * (1.0 + 0.875 * fedep)));
						// " the above evaluates the logarithm of the midpoint energy"
						// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
						lelktmp = lelke;// lelktmp = {P3};
						if (lelec < 0) {
							// $EVALUATE dedxmid USING ededx(elktmp);
							dedxmid = EGS4.EDEDX1[lelktmp - 1][medium - 1]
									* elktmp
									+ EGS4.EDEDX0[lelktmp - 1][medium - 1];
							aux = EGS4.EDEDX1[lelktmp - 1][medium - 1]
									/ dedxmid; // (lelktmp,medium)/dedxmid;
						} else {
							// $EVALUATE dedxmid USING pdedx(elktmp);
							dedxmid = EGS4.PDEDX1[lelktmp - 1][medium - 1]
									* elktmp
									+ EGS4.PDEDX0[lelktmp - 1][medium - 1];
							aux = EGS4.PDEDX1[lelktmp - 1][medium - 1]
									/ dedxmid; // (lelktmp,medium)/dedxmid;
						}
						aux = aux * (1.0 + 2.0 * aux) * (fedep / (2.0 - fedep))
								* (fedep / (2.0 - fedep)) / 6.0;
						// {P6} = fedep*{P1}/dedxmid*(1+aux);
						range = fedep * eke / dedxmid * (1.0 + aux);
						// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"
						range = (range + EGS4.range_ep[qel][lelke - 1][medium - 1])
								/ rhof;

						if (range > max_range) {
							max_range = range;
						}
					}
				}
				sip = 1.02 * max_range; // "A safety factor"
				rangerr1[j - 1] = (sip - si) * eker1;
				rangerr0[j - 1] = sip - rangerr1[j - 1] * elke;
				// "write(6,'(4g15.6)') eke,sip,rangerr0(j),rangerr1(j);"
				si = sip;
			}
			rangerr1[EGS4.$MXRANGE - 1] = rangerr1[EGS4.$MXRANGE - 2];
			rangerr0[EGS4.$MXRANGE - 1] = rangerr0[EGS4.$MXRANGE - 2];

		}
		// else if (irange_rej==EGS4.iCavitySPH)
		// {

		// }
	}
}
