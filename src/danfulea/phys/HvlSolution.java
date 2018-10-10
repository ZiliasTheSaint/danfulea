package danfulea.phys;

import java.util.ResourceBundle;

import danfulea.math.Convertor;
import danfulea.math.Sort;
import danfulea.math.StatsUtil;
import danfulea.math.numerical.EvalFunc;
import danfulea.math.numerical.Interpolator;
import danfulea.math.numerical.RootFind;

/**
 * This class try to compute theoretical HVL from data acquired using XCOMP5r. XCOMP5r 
 * is old and a much more complex and reliable database, SRS78, is available. 
 * XRay and XRAYspectrum should be used instead. Therefore, this class 
 * is deprecated. 
 * @author Dan Fulea, 2005
 *
 */
@Deprecated
public class HvlSolution {
	private static boolean valid = false;

	public static boolean validate() {
		return valid;
	}

	// return hvl1 prin diferite metode
	// parametrii de intrare: cele 2 siruri de coordonate x,y, valoarea
	// ordinului polinomului
	// de interpolare direct n si a celui pe coordonate inversate nn, si
	// valoarea de zero
	// adica y la x=0->a se studia metoda hvl din fizica.
	// return result de [0] pt hvl1 si de [1] pt hvl2 iar restul e:
	// [0]-polynomial revert,[1]-spline revert
	// [2]-secant, [3]-iter, [4]-birge, [5]-spline direct
	// OBS. PENTRU REZULTATE BUNE SE ALEGE N SI NN PENTRU FUNCTII MONOTONE
	// (DESCRESCATOARE RESPECTIV
	// CRESCATOARE) ASTFEL CA SE ELIMINA CAZUL BIRGE SI SPLINE DIRECT CU MAI
	// MULTE SOLUTII!!!
	// IN PLUS PE CAT POSIBIL TREBUIE CA NUMARUL DE PUNCTE DE INFLEXIUNE SA FIE
	// NUL
	// METODELE DE REZOLVARE DIRECTE FARA APROXIMARE PE REZOLVARI DE ECUATII DE
	// GRAD 2 SAU 3
	// POT DA GRES IN REZOLVARE!!!<->COEFICIENTII SUNT DE ASA NATURA INCAT NU SE
	// POT GASI ZEROURILE
	// Trebe ca numarul de puncte de retea>=3->pentru metodele spline!!!!!!!!
	// IN CEEA CE PRIVESTE EROAREA DE EVALUARE A HVL-URILOR, ACEASTA NU VA FI
	// DATA DE DIFERENTA
	// PE PUNCTELE DE RETEA A FUNCTIEI DE INTERPOLARE PENTRU CA ACEASTA E PENTRU
	// VALIDAREA METODEI,
	// DE ADUCERE LA O FUNCTIE TEORETIZATA MONOTON DESCRESCATOARE CU PUTINE
	// PUNCTE DE INFLEXIUNE(IDEAL 0)
	// SI NICI DE STATISTICA LA ZERO MM AL PENTRU CA SI ACEASTA ESTE PENTRU
	// VALIDARE-MASURATORILE VOR FI
	// REPETATE LA INSUCCES------>EA VA FI DATA DE O INCREDERE ASOCIATA ASOCIATA
	// CA LA HVL TEORETIC(5%)
	// CARE TINE DE PURITATEA PLACILOR DE AL SI DE ACURATETEA DETECTORULUI,
	// CELELALTE EFECTE NEDORITE CA
	// FLUCTUATIILE DE CURENT INTRAND LA VALIDAREA METODEI.
	public static double[][] hvl12(double[] x, double[] y, int n, int nn,
			double zeroValue) {
		double[][] result = new double[6][2];
		double hvl1 = 0.0;
		double qvl = 0.0;
		result[0][0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[0][1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1][0] = -1.0;
		result[1][1] = -1.0;
		result[2][0] = -1.0;
		result[2][1] = -1.0;
		result[3][0] = -1.0;
		result[3][1] = -1.0;
		result[4][0] = -1.0;
		result[4][1] = -1.0;
		result[5][0] = -1.0;
		result[5][1] = -1.0;

		double[][] a = new double[nn + 1][1];// matricea coeficientilor functiei
												// polinomiale de interpolare
		// sortez sirurile dupa y---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(yi, xi);
		// end sort------------------------------------
		// --revert coordonates for HVL1!!!------------------------------

		Interpolator.polynomial(yi, xi, nn, a);
		double[] d = HvlUtil.convertMatrix(a, 0);// vectorul atasat lui a!!!

		for (int i = 0; i < d.length; i++) {
			hvl1 = hvl1 + d[i] * Math.pow(zeroValue / 2, i);
			qvl = qvl + d[i] * Math.pow(zeroValue / 4, i);
		}
		result[0][0] = hvl1;
		result[0][1] = qvl - result[0][0];
		// --end revert --------------------------------------------------
		// -spline revert-------------------------------------------------

		int ni = x.length;
		double[] as = new double[ni];
		double[] bs = new double[ni];
		double[] cs = new double[ni];
		double[] ds = new double[ni];

		int ns = 2;
		double[] xs = new double[ns];
		double[] ys = new double[ns];
		xs[0] = zeroValue / 4;
		xs[1] = zeroValue / 2;

		double dm = (xi[1] - xi[0]) / (yi[1] - yi[0])
				- ((xi[2] - xi[1]) / (yi[2] - yi[1])) + (xi[2] - xi[0])
				/ (yi[2] - yi[0]);
		double dp = -((xi[ni - 2] - xi[ni - 3]) / (yi[ni - 2] - yi[ni - 3]))
				+ (xi[ni - 1] - xi[ni - 2]) / (yi[ni - 1] - yi[ni - 2])
				+ (xi[ni - 1] - xi[ni - 3]) / (yi[ni - 1] - yi[ni - 3]);

		Interpolator.spline(yi, xi, ni, xs, ys, ns, as, bs, cs, ds, dm, dp);

		result[1][0] = ys[1];
		result[1][1] = ys[0] - result[1][0];
		// -end spline revert--------------------------------------------
		// --secant method for HVL1!!!------------------------------------
		double min1 = 0.0;
		double max1 = 0.0;
		double min2 = 0.0;
		double max2 = 0.0;
		// sortez sirurile dupa x---just in case
		Sort.qSort2(xi, yi);
		// end sort------------------------------------
		a = new double[n + 1][1];
		Interpolator.polynomial(xi, yi, n, a);
		d = HvlUtil.convertMatrix(a, 0);
		// --determin valoarea minima si maxima a intervalului in care caut
		// solutia
		// daca poate se va lua +/-2 puncte de retea

		boolean validhvl1 = false;
		boolean validhvl2 = false;

		boolean hl1 = false;// marcheaza cazul particular de punct de retea
		boolean hl2 = false;// marcheaza cazul particular de punct de retea

		double z1 = zeroValue / 2;
		double z2 = zeroValue / 4;

		// cazul pentru ultimul yi particular
		if (yi[yi.length - 1] == z1) {
			hvl1 = x[yi.length - 1];
			result[2][0] = hvl1;
			result[4][0] = hvl1;
			result[5][0] = hvl1;
			hl1 = true;
		}

		for (int i = 0; i <= yi.length - 2; i++)// pana la ultimul element
		{
			if (yi[i] == z1) {
				hvl1 = x[i];
				result[2][0] = hvl1;
				result[4][0] = hvl1;
				result[5][0] = hvl1;
				hl1 = true;
			}

			if (yi[i] > z1 && z1 > yi[i + 1])// yi e monoton descrescator
			{
				if (i - 1 >= 0)
					min1 = xi[i - 1];
				else
					min1 = xi[i];
				if (i + 2 <= yi.length - 1)
					max1 = xi[i + 2];
				else
					max1 = xi[i + 1];

				hvl1 = HvlUtil.hvlSecant(min1, max1, d, 1, zeroValue);// caut>1
				validhvl1 = HvlUtil.validate();
				if (validhvl1)
					result[2][0] = hvl1;

			}

			if (yi[i] == z2) {
				qvl = x[i];
				result[2][1] = qvl - result[2][0];
				result[4][1] = qvl;// urmeaza -hvl1
				result[5][1] = qvl;// urmeaza -hvl1
				hl2 = true;
			}

			if (yi[i] > z2 && z2 > yi[i + 1])// yi e monoton descrescator
			{
				if (i - 1 >= 0)
					min2 = xi[i - 1];
				else
					min2 = xi[i];
				if (i + 2 <= yi.length - 1)
					max2 = xi[i + 2];
				else
					max2 = xi[i + 1];

				qvl = HvlUtil.hvlSecant(min2, max2, d, 2, zeroValue);// caut>2
				validhvl2 = HvlUtil.validate();
				if (validhvl2)
					result[2][1] = qvl - result[2][0];

			}

		}

		if (yi[yi.length - 1] == z2) {
			qvl = x[yi.length - 1];
			result[2][1] = qvl - result[2][0];
			result[4][1] = qvl;// urmeaza -hvl1
			result[5][1] = qvl;// urmeaza -hvl1
			hl2 = true;
		}

		// --end secant----------------------------------------------------
		// --aproximatii succesive
		validhvl1 = false;
		validhvl2 = false;

		if (result[2][0] != -1) {
			hvl1 = HvlUtil.hvlIter(result[2][0], d, 1, zeroValue);
			validhvl1 = HvlUtil.validate();
			if (validhvl1)
				result[3][0] = hvl1;

		}

		if (result[2][1] != -1) {
			qvl = HvlUtil.hvlIter(result[2][1] + result[2][0], d, 2, zeroValue);
			validhvl2 = HvlUtil.validate();
			if (validhvl2)
				result[3][1] = qvl - result[3][0];

		}

		// ----end aproximatii succesive-------------------------------------
		// birge method------------------------------------------------------
		validhvl1 = false;
		validhvl2 = false;

		if (!hl1) {
			double[] u = HvlUtil.updateToHvlBirge(d, 1, zeroValue);// pentru
																	// hvl1
			double[] cc = HvlUtil.invert(u);
			double[] z = RootFind.birge(cc);

			for (int i = 0; i < z.length; i++) {
				if (min1 < z[i] && z[i] < max1) {
					hvl1 = z[i];
					validhvl1 = true;
					break;
				}
			}

		}

		if (validhvl1)
			result[4][0] = hvl1;

		if (!hl2) {
			double[] u2 = HvlUtil.updateToHvlBirge(d, 2, zeroValue);// pentru
																	// hvl1
			double[] cc2 = HvlUtil.invert(u2);
			double[] zz = RootFind.birge(cc2);

			for (int i = 0; i < zz.length; i++) {
				if (min2 < zz[i] && zz[i] < max2) {
					qvl = zz[i];
					validhvl2 = true;
					break;
				}
			}
		} else
			result[4][1] = result[4][1] - result[4][0];

		if (validhvl2)
			result[4][1] = qvl - result[4][0];
		// end birge----------------------------------------------------------
		// direct spline-------------------------------------------------------
		validhvl1 = false;
		validhvl2 = false;

		as = new double[ni];
		bs = new double[ni];
		cs = new double[ni];
		ds = new double[ni];

		dm = (yi[1] - yi[0]) / (xi[1] - xi[0])
				- ((yi[2] - yi[1]) / (xi[2] - xi[1])) + (yi[2] - yi[0])
				/ (xi[2] - xi[0]);
		dp = -((yi[ni - 2] - yi[ni - 3]) / (xi[ni - 2] - xi[ni - 3]))
				+ (yi[ni - 1] - yi[ni - 2]) / (xi[ni - 1] - xi[ni - 2])
				+ (yi[ni - 1] - yi[ni - 3]) / (xi[ni - 1] - xi[ni - 3]);

		Interpolator.spline(xi, yi, ni, null, null, 0, as, bs, cs, ds, dm, dp);

		if (!hl1) {
			int kh = Sort.findHvlPosition(yi, 1, zeroValue);

			if (Sort.validare()) {
				ds[kh] = ds[kh] - zeroValue / 2;// HVL1!!
				double[] sc = new double[3];
				sc[0] = Math.abs(bs[kh]);
				sc[1] = Math.abs(cs[kh]);
				sc[2] = Math.abs(ds[kh]);

				// cea mai mica valoare pentru trecere la rezolvare quadratica
				// daca e cazul
				double sci = Sort.findValue(sc, sc.length);

				if (Math.abs(as[kh]) < sci * 1e-8) {
					double[] z = EvalFunc.quadratic(bs[kh], cs[kh], ds[kh]);
					for (int i = 0; i < 2; i++) {
						if (min1 < z[i] && z[i] < max1) {
							hvl1 = z[i];
							validhvl1 = true;
							break;
						}
					}
				} else {

					double[] z = EvalFunc.cubic(as[kh], bs[kh], cs[kh], ds[kh]);

					for (int i = 0; i < 3; i++) {
						if (min1 < z[i] && z[i] < max1) {
							hvl1 = z[i];
							validhvl1 = true;
							break;
						}
					}
				}
			}

		}

		if (validhvl1)
			result[5][0] = hvl1;

		// HVL2
		if (!hl2) {
			int kh = Sort.findHvlPosition(yi, 2, zeroValue);
			if (Sort.validare()) {
				ds[kh] = ds[kh] - zeroValue / 4;// HVL2!!
				double[] sc = new double[3];
				sc[0] = Math.abs(bs[kh]);
				sc[1] = Math.abs(cs[kh]);
				sc[2] = Math.abs(ds[kh]);
				// cea mai mica valoare pentru trecere la rezolvare quadratica
				// daca e cazul
				double sci = Sort.findValue(sc, sc.length);
				if (Math.abs(as[kh]) < sci * 1e-8) {
					double[] z = EvalFunc.quadratic(bs[kh], cs[kh], ds[kh]);
					for (int i = 0; i < 2; i++) {
						if (min2 < z[i] && z[i] < max2) {
							qvl = z[i];
							validhvl2 = true;
							break;
						}
					}
				} else {
					double[] z = EvalFunc.cubic(as[kh], bs[kh], cs[kh], ds[kh]);
					for (int i = 0; i < 3; i++) {
						if (min2 < z[i] && z[i] < max2) {
							qvl = z[i];
							validhvl2 = true;
							break;
						}
					}
				}
			}

		} else
			result[5][1] = result[5][1] - result[5][0];

		if (validhvl2)
			result[5][1] = qvl - result[5][0];

		// --------------------------------------------------------------------
		return result;
	}

	// Evaluare maxim->MAB la masa KCl vs Effic
	// Se imparte pentru siguranta intervalele k-1,k si k,k+1 unde k e pozitia y
	// de maxim a punctelor de
	// retea in cate 2 parti egale, desi interpolarea polinomiala bazata pe
	// metoda celor mai mici patrate
	// si cu functii spline asigura o netezire (functie monotona) a functiei de
	// interpolare pe intervalele
	// date de punctele de retea!!!Cazul cel mai rau e sa am 2 puncte de extrem
	// pe un astfel de interval
	// dar acest caz e extrem de rar intalnit, aproape imposibil, cele 2 extreme
	// putand fi doar in jurul
	// celor 2 puncte de retea ce formeaza acest interval, de aceea impartirea
	// in 2 parti egale e justificata..
	// [i][0]-valoarea xmaxim, iar [i][1]-valoarea ymaxim
	public static double[][] maxEvaluation(double[] x, double[] y, int n) {
		double[][] result = new double[4][2];

		result[0][0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1][0] = -1.0;
		result[2][0] = -1.0;
		result[3][0] = -1.0;
		result[0][1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1][1] = -1.0;
		result[2][1] = -1.0;
		result[3][1] = -1.0;

		int ni = x.length;
		double[][] a = new double[n + 1][1];// matricea coeficientilor functiei
											// polinomiale de interpolare
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);
		// end sort------------------------------------
		// --secant method for max!!!------------------------------------
		// sortez sirurile dupa x---just in case
		// Sort.qSort2(xi,yi);
		// end sort------------------------------------
		a = new double[n + 1][1];
		Interpolator.polynomial(xi, yi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);
		// --determin valoarea minima si maxima a intervalului in care caut
		// solutia
		// daca poate se va lua +/-2 puncte de retea
		boolean validmax = false;

		int maxpos = Sort.findPosition(yi, 1);// pozitia maximului

		// System.out.println(""+maxpos);
		if (maxpos > 0 && maxpos < yi.length - 1) {
			double amin = xi[maxpos - 1];
			double amax = xi[maxpos + 1];
			double pas1 = (maxpos - amin) / 2.0;// pentru siguranta impart
												// fiecare subinterval la 2!!
			double pas2 = (amax - maxpos) / 2.0;
			double[] ymaxSearch = new double[5];
			double[] xmaxSearch = new double[5];
			double ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			xmaxSearch[0] = xi[maxpos];
			// y[maxpos-1---y[maxpos+1] intervalul pe care calculez
			// extremul(maximul)
			xmaxSearch[1] = HvlUtil.maxSecant(amin, amin + pas1, d);
			validmax = HvlUtil.validate();
			if (!validmax) {
				ymaxSearch[1] = yi[maxpos] - 1.0;// forteaza un numar arbitrar
													// mai mic ca maximul de
													// retea
				// valoarea x-ului in acest caz nu conteaza!!
			} else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[1], i);
				ymaxSearch[1] = ymx;
			}
			// astfel in final acesta va fi considertat ca maxim!!
			xmaxSearch[2] = HvlUtil.maxSecant(amin + pas1, maxpos, d);
			validmax = HvlUtil.validate();
			if (!validmax)
				ymaxSearch[2] = yi[maxpos] - 1.0;
			else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[2], i);
				ymaxSearch[2] = ymx;
			}

			xmaxSearch[3] = HvlUtil.maxSecant(maxpos, maxpos + pas2, d);
			validmax = HvlUtil.validate();
			if (!validmax)
				ymaxSearch[3] = yi[maxpos] - 1.0;
			else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[3], i);
				ymaxSearch[3] = ymx;
			}

			xmaxSearch[4] = HvlUtil.maxSecant(maxpos + pas2, amax, d);
			validmax = HvlUtil.validate();
			if (!validmax)
				ymaxSearch[4] = yi[maxpos] - 1.0;
			else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[4], i);
				ymaxSearch[4] = ymx;
			}
			Sort.qSort2(ymaxSearch, xmaxSearch);

			result[0][0] = xmaxSearch[xmaxSearch.length - 1];
			result[0][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		if (maxpos == 0) {
			double amin = xi[maxpos];
			double amax = xi[maxpos + 1];
			double pas = (amax - amin) / 2.0;
			double[] ymaxSearch = new double[3];
			double[] xmaxSearch = new double[3];
			double ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			xmaxSearch[0] = xi[maxpos];

			xmaxSearch[1] = HvlUtil.maxSecant(amin, amin + pas, d);
			validmax = HvlUtil.validate();
			if (!validmax) {
				ymaxSearch[1] = yi[maxpos] - 1.0;
			} else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[1], i);
				ymaxSearch[1] = ymx;
			}

			xmaxSearch[2] = HvlUtil.maxSecant(amin + pas, amax, d);
			validmax = HvlUtil.validate();
			if (!validmax) {
				ymaxSearch[2] = yi[maxpos] - 1.0;
			} else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[2], i);
				ymaxSearch[2] = ymx;
			}

			Sort.qSort2(ymaxSearch, xmaxSearch);

			result[0][0] = xmaxSearch[xmaxSearch.length - 1];
			result[0][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		if (maxpos == yi.length - 1) {
			double amin = xi[maxpos - 1];
			double amax = xi[maxpos];
			double pas = (amax - amin) / 2.0;
			double[] ymaxSearch = new double[3];
			double[] xmaxSearch = new double[3];
			double ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			xmaxSearch[0] = xi[maxpos];

			xmaxSearch[1] = HvlUtil.maxSecant(amin, amin + pas, d);
			validmax = HvlUtil.validate();
			if (!validmax) {
				ymaxSearch[1] = yi[maxpos] - 1.0;
			} else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[1], i);
				ymaxSearch[1] = ymx;
			}

			xmaxSearch[2] = HvlUtil.maxSecant(amin + pas, amax, d);
			validmax = HvlUtil.validate();
			if (!validmax) {
				ymaxSearch[2] = yi[maxpos] - 1.0;
			} else {
				ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmaxSearch[2], i);
				ymaxSearch[2] = ymx;
			}

			Sort.qSort2(ymaxSearch, xmaxSearch);

			result[0][0] = xmaxSearch[xmaxSearch.length - 1];
			result[0][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		// --end secant----------------------------------------------------
		// --aproximatii succesive
		validmax = false;
		double xmx = 0.0;
		if (result[0][0] != -1) {
			xmx = HvlUtil.maxIter(result[0][0], d);
			validmax = HvlUtil.validate();
			if (validmax) {
				result[1][0] = xmx;
				double ymx = 0.0;
				for (int i = 0; i < d.length; i++)
					ymx = ymx + d[i] * Math.pow(xmx, i);
				result[1][1] = ymx;
			}

		}

		// ----end aproximatii succesive-------------------------------------
		// birge method------------------------------------------------------
		validmax = false;

		double amin = xi[0];
		double amax = xi[xi.length - 1];

		if (maxpos > 0 && maxpos < yi.length - 1) {
			amin = xi[maxpos - 1];
			amax = xi[maxpos + 1];
		} else {
			if (maxpos == 0) {
				amin = xi[maxpos];
				amax = xi[maxpos + 1];
			} else// maxpos=yi.length-1
			{
				amin = xi[maxpos - 1];
				amax = xi[maxpos];
			}
		}

		double[] u = HvlUtil.updateToMaxBirge(d);
		double[] cc = HvlUtil.invert(u);
		double[] z = RootFind.birge(cc);
		double[] zy = new double[z.length];
		double ymx = 0.0;

		for (int i = 0; i < z.length; i++) {
			if (amin < z[i] && z[i] < amax) {
				ymx = 0.0;
				for (int j = 0; j < d.length; j++)
					ymx = ymx + d[j] * Math.pow(z[i], j);
				zy[i] = ymx;
				validmax = true;
			}
		}

		if (validmax) {
			Sort.qSort2(zy, z);
			if (zy[zy.length - 1] > yi[maxpos]) {
				result[2][0] = z[z.length - 1];
				result[2][1] = zy[zy.length - 1];
			} else {
				result[2][0] = xi[maxpos];
				result[2][1] = yi[maxpos];
			}
		}

		// end birge----------------------------------------------------------
		// direct spline-------------------------------------------------------
		validmax = false;

		double[] as = new double[ni];
		double[] bs = new double[ni];
		double[] cs = new double[ni];
		double[] ds = new double[ni];

		double dm = (yi[1] - yi[0]) / (xi[1] - xi[0])
				- ((yi[2] - yi[1]) / (xi[2] - xi[1])) + (yi[2] - yi[0])
				/ (xi[2] - xi[0]);
		double dp = -((yi[ni - 2] - yi[ni - 3]) / (xi[ni - 2] - xi[ni - 3]))
				+ (yi[ni - 1] - yi[ni - 2]) / (xi[ni - 1] - xi[ni - 2])
				+ (yi[ni - 1] - yi[ni - 3]) / (xi[ni - 1] - xi[ni - 3]);

		Interpolator.spline(xi, yi, ni, null, null, 0, as, bs, cs, ds, dm, dp);

		if (maxpos > 0 && maxpos < yi.length - 1) {
			double a1 = as[maxpos - 1];
			double b1 = bs[maxpos - 1];
			double c1 = cs[maxpos - 1];
			double d1 = ds[maxpos - 1];

			double a2 = as[maxpos];
			double b2 = bs[maxpos];
			double c2 = cs[maxpos];
			double d2 = ds[maxpos];

			// coeficientii derivatelor:
			double am1 = 3 * a1;
			double bm1 = 2 * b1;
			double cm1 = c1;

			double am2 = 3 * a2;
			double bm2 = 2 * b2;
			double cm2 = c2;

			double[] ymaxSearch = new double[5];
			double[] xmaxSearch = new double[5];
			ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			ymaxSearch[1] = yi[maxpos] - 1.0;// forteaza un numar arbitrar mai
												// mic ca maximul de retea
			ymaxSearch[2] = yi[maxpos] - 1.0;
			ymaxSearch[3] = yi[maxpos] - 1.0;
			ymaxSearch[4] = yi[maxpos] - 1.0;
			xmaxSearch[0] = xi[maxpos];

			double[] sc = new double[2];
			sc[0] = Math.abs(bm1);
			sc[1] = Math.abs(cm1);
			// cea mai mica valoare pentru trecere la rezolvare quadratica daca
			// e cazul
			double sci = Sort.findValue(sc, sc.length);
			if (Math.abs(am1) < sci * 1e-8)
				am1 = 0.0;

			double[] ec2 = EvalFunc.quadratic(am1, bm1, cm1);
			validmax = !EvalFunc.failB;// Equations.validate();
			if (validmax) {
				if (ec2[0] <= xi[maxpos - 1] || ec2[0] >= xi[maxpos])
					ymaxSearch[1] = yi[maxpos] - 1.0;// forteaza un numar
														// arbitrar mai mic ca
														// maximul de retea
				else {
					xmaxSearch[1] = ec2[0];
					ymaxSearch[1] = a1 * Math.pow(xmaxSearch[1], 3) + b1
							* Math.pow(xmaxSearch[1], 2) + +c1 * xmaxSearch[1]
							+ d1;
				}

				if (ec2[1] <= xi[maxpos - 1] || ec2[1] >= xi[maxpos])
					ymaxSearch[2] = yi[maxpos] - 1.0;// forteaza un numar
														// arbitrar mai mic ca
														// maximul de retea
				else {
					xmaxSearch[2] = ec2[1];
					ymaxSearch[2] = a1 * Math.pow(xmaxSearch[2], 3) + b1
							* Math.pow(xmaxSearch[2], 2) + +c1 * xmaxSearch[2]
							+ d1;
				}
			}

			sc[0] = Math.abs(bm2);
			sc[1] = Math.abs(cm2);
			// cea mai mica valoare pentru trecere la rezolvare quadratica daca
			// e cazul
			sci = Sort.findValue(sc, sc.length);
			if (Math.abs(am2) < sci * 1e-8)
				am2 = 0.0;

			ec2 = EvalFunc.quadratic(am2, bm2, cm2);
			validmax = !EvalFunc.failB;// Equations.validate();
			if (validmax) {
				if (ec2[0] <= xi[maxpos] || ec2[0] >= xi[maxpos + 1])
					ymaxSearch[3] = yi[maxpos] - 1.0;// forteaza un numar
														// arbitrar mai mic ca
														// maximul de retea
				else {
					xmaxSearch[3] = ec2[0];
					ymaxSearch[3] = a2 * Math.pow(xmaxSearch[3], 3) + b2
							* Math.pow(xmaxSearch[3], 2) + +c2 * xmaxSearch[3]
							+ d2;
				}

				if (ec2[1] <= xi[maxpos] || ec2[1] >= xi[maxpos + 1])
					ymaxSearch[4] = yi[maxpos] - 1.0;// forteaza un numar
														// arbitrar mai mic ca
														// maximul de retea
				else {
					xmaxSearch[4] = ec2[1];
					ymaxSearch[4] = a2 * Math.pow(xmaxSearch[4], 3) + b2
							* Math.pow(xmaxSearch[4], 2) + +c2 * xmaxSearch[4]
							+ d2;
				}
			}

			Sort.qSort2(ymaxSearch, xmaxSearch);
			result[3][0] = xmaxSearch[xmaxSearch.length - 1];
			result[3][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		if (maxpos == 0) {
			double a1 = as[maxpos];
			double b1 = bs[maxpos];
			double c1 = cs[maxpos];
			double d1 = ds[maxpos];

			// coeficientii derivatelor:
			double am1 = 3 * a1;
			double bm1 = 2 * b1;
			double cm1 = c1;

			double[] ymaxSearch = new double[3];
			double[] xmaxSearch = new double[3];
			ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			ymaxSearch[1] = yi[maxpos] - 1.0;
			ymaxSearch[2] = yi[maxpos] - 1.0;

			xmaxSearch[0] = xi[maxpos];

			double[] sc = new double[2];
			sc[0] = Math.abs(bm1);
			sc[1] = Math.abs(cm1);
			// cea mai mica valoare pentru trecere la rezolvare quadratica daca
			// e cazul
			double sci = Sort.findValue(sc, sc.length);
			if (Math.abs(am1) < sci * 1e-8)
				am1 = 0.0;

			double[] ec2 = EvalFunc.quadratic(am1, bm1, cm1);
			validmax = !EvalFunc.failB;// Equations.validate();
			if (validmax) {
				if (ec2[0] <= xi[maxpos] || ec2[0] >= xi[maxpos + 1])
					ymaxSearch[1] = yi[maxpos] - 1.0;
				else {
					xmaxSearch[1] = ec2[0];
					ymaxSearch[1] = a1 * Math.pow(xmaxSearch[1], 3) + b1
							* Math.pow(xmaxSearch[1], 2) + +c1 * xmaxSearch[1]
							+ d1;
				}

				if (ec2[1] <= xi[maxpos] || ec2[1] >= xi[maxpos + 1])
					ymaxSearch[2] = yi[maxpos] - 1.0;
				else {
					xmaxSearch[2] = ec2[1];
					ymaxSearch[2] = a1 * Math.pow(xmaxSearch[2], 3) + b1
							* Math.pow(xmaxSearch[2], 2) + +c1 * xmaxSearch[2]
							+ d1;
				}
			}

			Sort.qSort2(ymaxSearch, xmaxSearch);
			result[3][0] = xmaxSearch[xmaxSearch.length - 1];
			result[3][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		if (maxpos == yi.length - 1) {
			double a1 = as[maxpos - 1];
			double b1 = bs[maxpos - 1];
			double c1 = cs[maxpos - 1];
			double d1 = ds[maxpos - 1];

			// coeficientii derivatelor:
			double am1 = 3 * a1;
			double bm1 = 2 * b1;
			double cm1 = c1;

			double[] ymaxSearch = new double[3];
			double[] xmaxSearch = new double[3];
			ymx = 0.0;
			ymaxSearch[0] = yi[maxpos];
			ymaxSearch[1] = yi[maxpos] - 1.0;
			ymaxSearch[2] = yi[maxpos] - 1.0;
			xmaxSearch[0] = xi[maxpos];

			double[] sc = new double[2];
			sc[0] = Math.abs(bm1);
			sc[1] = Math.abs(cm1);
			// cea mai mica valoare pentru trecere la rezolvare quadratica daca
			// e cazul
			double sci = Sort.findValue(sc, sc.length);
			if (Math.abs(am1) < sci * 1e-8)
				am1 = 0.0;

			double[] ec2 = EvalFunc.quadratic(am1, bm1, cm1);
			validmax = !EvalFunc.failB;// Equations.validate();
			if (validmax) {
				if (ec2[0] <= xi[maxpos - 1] || ec2[0] >= xi[maxpos])
					ymaxSearch[1] = yi[maxpos] - 1.0;
				else {
					xmaxSearch[1] = ec2[0];
					ymaxSearch[1] = a1 * Math.pow(xmaxSearch[1], 3) + b1
							* Math.pow(xmaxSearch[1], 2) + +c1 * xmaxSearch[1]
							+ d1;
				}

				if (ec2[1] <= xi[maxpos - 1] || ec2[1] >= xi[maxpos])
					ymaxSearch[2] = yi[maxpos] - 1.0;
				else {
					xmaxSearch[2] = ec2[1];
					ymaxSearch[2] = a1 * Math.pow(xmaxSearch[2], 3) + b1
							* Math.pow(xmaxSearch[2], 2) + +c1 * xmaxSearch[2]
							+ d1;
				}
			}

			Sort.qSort2(ymaxSearch, xmaxSearch);
			result[3][0] = xmaxSearch[xmaxSearch.length - 1];
			result[3][1] = ymaxSearch[ymaxSearch.length - 1];
		}

		// --------------------------------------------------------------------
		return result;
	}

	// metoda de evaluare a HVL1,2 teoretic-> bazat pe date furnizate de metoda
	// Birch_Marshall
	// rezolvata de softul xcomp5r de la Univ de fizica Viena
	// se utilizeaza datele la filtrarile totale:
	// 0.2 mmAl
	// 0.5 mmAl
	// 1 mmAl
	// 1.5 mmAl
	// 2 mmAl
	// 2.5 mmAl
	// 3 mmAl
	// 3.5 mmAl
	// 4 mmAl
	// pentru un unghi al anodului default de >=7 grade respectiv <=27 si pentru
	// kilovoltaje de 30-150 din 10 in 10kV!!
	// s-a folosit metoda de interpolare liniara pe punctele apropiate!!
	private static double[] hvl12teorAt17(double kv, double filtration) {
		valid = false;

		double[] result = new double[2];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");
		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "117kv";
			String s2 = "217kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;
			for (int i = 0; i < n; i++) {
				eq[i] = filtration * mn[i][0] + mn[i][1];

				if (filtration <= f[i + 1]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = filtration * mn2[i][0] + mn2[i][1];

				if (filtration <= f[i + 1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					break;
				}
			}

		}

		return result;
	}

	// la 7 grade anod
	private static double[] hvl12teorAt7(double kv, double filtration) {
		valid = false;

		double[] result = new double[2];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");
		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "17kv";
			String s2 = "27kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;
			for (int i = 0; i < n; i++) {
				eq[i] = filtration * mn[i][0] + mn[i][1];

				if (filtration <= f[i + 1]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = filtration * mn2[i][0] + mn2[i][1];

				if (filtration <= f[i + 1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					break;
				}
			}

		}

		return result;
	}

	//
	// la 27 grade anod
	private static double[] hvl12teorAt27(double kv, double filtration) {
		valid = false;

		double[] result = new double[2];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");
		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "127kv";
			String s2 = "227kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;
			for (int i = 0; i < n; i++) {
				eq[i] = filtration * mn[i][0] + mn[i][1];

				if (filtration <= f[i + 1]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = filtration * mn2[i][0] + mn2[i][1];

				if (filtration <= f[i + 1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					break;
				}
			}

		}

		return result;
	}

	// hvlt propriuzis
	public static double[] hvl12teor(double kv, double filtration,
			double anodAngle) {
		if (XRay.ICALC == 1) {
			double[] rrr = new double[2];
			XRay.computeHVL1("TISS", true);// tissue equivalent
			XRay.computeHVL2("TISS", true);// tissue equivalent
			XRay.computeHVL1("AL", false);// final of interest in mmAl
			XRay.computeHVL2("AL", false);// final of interest in mmAl
			rrr[0] = XRay.HVL1;
			rrr[1] = XRay.HVL2;
			valid = true;
			return rrr;
		}
		// =============================================================
		valid = false;

		double[] result = new double[2];
		double uAnodMin = 7.0;
		double uAnodMed = 17.0;
		double uAnodMax = 27.0;

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces

		if (anodAngle < uAnodMin || anodAngle > uAnodMax) {
			valid = false;
			return result;
		}

		if (anodAngle == uAnodMin) {
			double[] result7 = hvl12teorAt7(kv, filtration);// aci se
															// valideste!!!
			if (validate())
				return result7;
			else
				return result;

		}

		if (anodAngle == uAnodMed) {
			double[] result17 = hvl12teorAt17(kv, filtration);// aci se
																// valideste!!!
			if (validate())
				return result17;
			else
				return result;

		}

		if (anodAngle == uAnodMax) {
			double[] result27 = hvl12teorAt27(kv, filtration);// aci se
																// valideste!!!
			if (validate())
				return result27;

		}

		if (anodAngle < uAnodMed) {
			double[] result7 = hvl12teorAt7(kv, filtration);// aci se
															// valideste!!!
			// System.out.println("7 "+result7[0]+"    "+result7[1]);
			if (validate()) {
				double[] result17 = hvl12teorAt17(kv, filtration);// aci se
																	// valideste!!!
				// System.out.println("17 "+result17[0]+"    "+result17[1]);
				if (validate()) {
					result[0] = linInt(uAnodMed, result17[0], uAnodMin,
							result7[0], anodAngle);
					result[1] = linInt(uAnodMed, result17[1], uAnodMin,
							result7[1], anodAngle);
					return result;
				} else
					return result;

			} else
				return result;

		}

		if (anodAngle < uAnodMax) {
			double[] result17 = hvl12teorAt17(kv, filtration);// aci se
																// valideste!!!
			if (validate()) {
				double[] result27 = hvl12teorAt27(kv, filtration);// aci se
																	// valideste!!!
				if (validate()) {
					result[0] = linInt(uAnodMax, result27[0], uAnodMed,
							result17[0], anodAngle);
					result[1] = linInt(uAnodMax, result27[1], uAnodMed,
							result17[1], anodAngle);
					return result;
				}
			}

		}

		return result;
	}

	//
	//
	// metoda pentru a extrage in functie de order minimul (order=0) sau maximul
	// (order!=0)
	// dintr-un sir de reali care INTERCALEAZA un real dat
	// daca min=max atunci am punct de retea!!!
	// returneaza in plus si pozitia din sir a acestui maxim sau minim
	public static double[] getNearestValueFromSeriesForValue(double[] series,
			double value, int order) {
		valid = false;
		double[] result = new double[2];
		result[0] = 0.0;
		result[1] = 0.0;
		double max = 0.0;
		double min = 0.0;
		double posmin = 0.0;
		double posmax = 0.0;
		String ss = null;
		double f = value;// din lene de modificat in cod
		// first sortez seria intr-un alt sir nou crescator
		double[] s = Sort.newQSort(series);
		// end sort
		for (int i = 0; i < s.length; i++) {
			if (s[i] <= value) {
				if (i != s.length) {
					if (s[i + 1] >= f) {
						min = s[i];
						ss = Convertor.intToString(i);
						posmin = Convertor.stringToDouble(ss);
						if (s[i] == f) {
							max = s[i];
							ss = Convertor.intToString(i);
							posmax = Convertor.stringToDouble(ss);
						} else {
							max = s[i + 1];
							ss = Convertor.intToString(i + 1);
							posmax = Convertor.stringToDouble(ss);
						}

						if (s[i + 1] == f) {
							min = s[i + 1];
							ss = Convertor.intToString(i + 1);
							posmin = Convertor.stringToDouble(ss);
						}
						break;
					}
					valid = true;
				} else// a ajuns la capatul sirului
				{
					min = s[s.length - 1];
					max = min;
					ss = Convertor.intToString(s.length - 1);
					posmin = Convertor.stringToDouble(ss);
					posmax = posmin;
					valid = true;
				}
			} else// cazul in care prima valoare din sir e deja mai mare
			{
				valid = false;
				break;
			}
		}

		if (order == 0) {
			result[0] = min;
			result[1] = posmin;
		} else {
			result[0] = max;
			result[1] = posmax;
		}

		return result;
	}

	// returneaza hvl-ul experimental la orice alt kV!!!
	public static double[] getHvlFromKv(double kv, double filtration,
			double anodAngle, double hvl1, double hvl2, double newkv) {
		if (XRay.ICALC == 1) {
			double[] rrr = new double[2];
			int ikvn = new Double(Math.floor(newkv)).intValue();
			XRay.computeHVL12Kvp(ikvn);
			valid = true;
			rrr[0] = XRay.HVL1_kv;
			rrr[1] = XRay.HVL2_kv;
			return rrr;
		}
		valid = false;

		double result[] = new double[2];
		result[0] = 0.0;
		result[1] = 0.0;

		double uAnodMin = 7.0;
		double uAnodMax = 27.0;

		if (anodAngle < uAnodMin || anodAngle > uAnodMax) {
			valid = false;
			return result;
		}

		double[] hvlt = hvl12teor(kv, filtration, anodAngle);// aici se face
																// validul!!
		if (validate()) {
			double p1 = hvl1 / hvlt[0];
			double p2 = hvl2 / hvlt[1];

			double[] hvltnew = hvl12teor(newkv, filtration, anodAngle);
			if (validate()) {
				result[0] = p1 * (hvltnew[0] - hvlt[0]) + hvl1;
				result[1] = p2 * (hvltnew[1] - hvlt[1]) + hvl2;
			}

		}

		return result;
	}

	// return 0 daca exp<teor real, 1 daca exp<teor real sau 2 daca statistic
	// ele-s egale
	private static int getConfigForHvlComment(double x1, double x2, double ab1,
			double ab2, int incredere) {
		valid = false;
		int i = 0;
		StatsUtil.confidenceLevel = incredere;
		boolean testt = StatsUtil.ttest(x1, x2, ab1, ab2);
		// Stat.tTestCompForMean(x1, x2, ab1, ab2, incredere);
		if (StatsUtil.failB)// !Stat.validate())
			return i;

		if (x1 < x2 && !testt) {
			i = 0;
			valid = true;
		}

		if (x1 > x2 && !testt) {
			i = 1;
			valid = true;
		}

		if (testt) {
			i = 2;
			valid = true;
		}

		return i;
	}

	// factorul de omogenitate si eroarea sa pentru exp si teoretic!
	// erorile sa nu fie date in procente ci in numar!!!
	public static double[][] homFactor(double hvl1e, double hvl2e,
			double hvl1t, double hvl2t, double errhvl1e, double errhvl2e,
			double errhvl1t, double errhvl2t) {
		valid = false;

		double[][] result = new double[2][2];
		result[0][0] = -1.0;// roexp
		result[0][1] = -1.0;// errroexp
		result[1][0] = -1.0;// roteor
		result[1][1] = -1.0;// errroteor
		if (hvl2e != 0.0 && hvl2t != 0.0) {
			double roe = hvl1e / hvl2e;
			double rot = hvl1t / hvl2t;
			double errroe = Math.sqrt(errhvl1e * errhvl1e / (hvl2e * hvl2e)
					+ errhvl2e * errhvl2e * hvl1e * hvl1e
					/ (hvl2e * hvl2e * hvl2e * hvl2e));
			double errrot = Math.sqrt(errhvl1t * errhvl1t / (hvl2t * hvl2t)
					+ errhvl2t * errhvl2t * hvl1t * hvl1t
					/ (hvl2t * hvl2t * hvl2t * hvl2t));
			valid = true;
			result[0][0] = roe;// roexp
			result[0][1] = errroe;// errroexp
			result[1][0] = rot;// roteor
			result[1][1] = errrot;// errroteor
		}
		return result;
	}

	// comentariu al hvl-ului experimental vs teoretic
	public static String hvlComment(double hvl1e, double hvl2e, double hvl1t,
			double hvl2t, double errhvl1e, double errhvl2e, double errhvl1t,
			double errhvl2t, int incredere) {
		String s = null;
		valid = false;

		int hvl1i = getConfigForHvlComment(hvl1e, hvl1t, errhvl1e, errhvl1t,
				incredere);
		if (!validate()) {
			return s;
		}
		int hvl2i = getConfigForHvlComment(hvl2e, hvl2t, errhvl2e, errhvl2t,
				incredere);
		if (!validate()) {
			return s;
		}
		double[][] ro = homFactor(hvl1e, hvl2e, hvl1t, hvl2t, errhvl1e,
				errhvl2e, errhvl1t, errhvl2t);
		if (!validate()) {
			return s;
		}
		int roi = getConfigForHvlComment(ro[0][0], ro[1][0], ro[0][1],
				ro[1][1], incredere);

		String resursa = Convertor.intToString(hvl1i)
				+ Convertor.intToString(hvl2i) + Convertor.intToString(roi);

		// System.out.println(resursa);
		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		s = resources.getString(resursa);

		return s;
	}

	// metoda de evaluare a filtrarii echivalente corespunzatoare devierii
	// HVL1,2 fata de teoretic
	// este util pentru corectia la calculul dozei(PCXMC care isi ia spectrul x
	// teoretic)
	// pot evalua la diferite plane de interes ale pacientului: pentru intreg
	// corpul
	// variabila plane va lua valoarea 6.6 (corpul uman total e echivalent cu
	// 6.6 mmAl)
	// return filtrarea eq vs HVL1, HVL2 si cea general echivalenta+cele 2
	// ponderi!
	private static double[] eqfiltrationAt17(double kv, double hvl1e,
			double hvl2e, double plane) {
		valid = false;

		double[] result = new double[5];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[2] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[3] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[4] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		// ponderi
		double p1 = 0.5;
		double p2 = 0.5;
		if (hvl1e >= plane)
			p1 = 1.0;
		else
			p1 = 1 - (plane - hvl1e) / plane;
		p2 = 1 - p1;
		result[3] = p1;
		result[4] = p2;
		//

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");

		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "f117kv";
			String s2 = "f217kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;

			// hvlt
			double[][] hvlt = new double[f.length][2];
			for (int ii = 0; ii < f.length; ii++) {
				double[] hvltf = hvl12teorAt17(kvp[j], f[ii]);
				hvlt[ii][0] = hvltf[0];
				hvlt[ii][1] = hvltf[1];
			}
			//

			for (int i = 0; i < n; i++) {
				eq[i] = hvl1e * mn[i][0] + mn[i][1];

				if (hvl1e <= hvlt[i + 1][0]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = hvl2e * mn2[i][0] + mn2[i][1];

				if (hvl2e <= hvlt[i + 1][1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					result[2] = (result[0] * p1 + result[1] * p2) / (p1 + p2);
					break;
				}
			}
			//
		}

		return result;
	}

	private static double[] eqfiltrationAt7(double kv, double hvl1e,
			double hvl2e, double plane) {
		valid = false;

		double[] result = new double[5];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[2] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[3] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[4] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		// ponderi
		double p1 = 0.5;
		double p2 = 0.5;
		if (hvl1e >= plane)
			p1 = 1.0;
		else
			p1 = 1 - (plane - hvl1e) / plane;
		p2 = 1 - p1;
		result[3] = p1;
		result[4] = p2;
		//

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");

		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "f17kv";
			String s2 = "f27kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;

			// hvlt
			double[][] hvlt = new double[f.length][2];
			for (int ii = 0; ii < f.length; ii++) {
				double[] hvltf = hvl12teorAt7(kvp[j], f[ii]);
				hvlt[ii][0] = hvltf[0];
				hvlt[ii][1] = hvltf[1];
			}
			//

			for (int i = 0; i < n; i++) {
				eq[i] = hvl1e * mn[i][0] + mn[i][1];

				if (hvl1e <= hvlt[i + 1][0]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = hvl2e * mn2[i][0] + mn2[i][1];

				if (hvl2e <= hvlt[i + 1][1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					result[2] = (result[0] * p1 + result[1] * p2) / (p1 + p2);
					break;
				}
			}
			//
		}

		return result;
	}

	private static double[] eqfiltrationAt27(double kv, double hvl1e,
			double hvl2e, double plane) {
		valid = false;

		double[] result = new double[5];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[2] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[3] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[4] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		// ponderi
		double p1 = 0.5;
		double p2 = 0.5;
		if (hvl1e >= plane)
			p1 = 1.0;
		else
			p1 = 1 - (plane - hvl1e) / plane;
		p2 = 1 - p1;
		result[3] = p1;
		result[4] = p2;
		//

		ResourceBundle resources;
		String baseName = "jdf.phys.resources.HvlSolutionResources";
		resources = ResourceBundle.getBundle(baseName);
		// pentru kvp=30,40,etc la HVL1 si HVL2;
		double[] f = (double[]) resources.getObject("filtration");
		double[] kvp = (double[]) resources.getObject("kv");

		int m = kvp.length;
		if (kv < kvp[0] || kv > kvp[m - 1]) {
			valid = false;
			return result;
		} else
			valid = true;
		int k = 0;
		double hvl1kv = 0.0;
		double hvl2kv = 0.0;
		double[] hl1 = new double[m];
		double[] hl2 = new double[m];
		double[] hlin1 = new double[m];
		double[] hlin2 = new double[m];
		double[] mm1 = new double[m - 1];
		double[] mm2 = new double[m - 1];
		double[] nn1 = new double[m - 1];
		double[] nn2 = new double[m - 1];
		for (int j = 0; j < m; j++) {
			String s = "f127kv";
			String s2 = "f227kv";
			String ss = Convertor.doubleToString(kvp[j]);
			k = 0;
			char[] ch = { ss.charAt(k) };// pentru conditie din while
			String sk = new String(ch);// pentru conditie din while
			while (!sk.equalsIgnoreCase(".")) {
				s = s + ss.charAt(k);
				s2 = s2 + ss.charAt(k);
				k++;
				char[] c = { ss.charAt(k) };
				sk = new String(c);
			}

			double[][] mn = (double[][]) resources.getObject(s);
			double[][] mn2 = (double[][]) resources.getObject(s2);
			int n = mn.length;
			int n2 = mn2.length;
			double[] eq = new double[n];
			double[] eq2 = new double[n2];

			hvl1kv = 0.0;
			hvl2kv = 0.0;

			// hvlt
			double[][] hvlt = new double[f.length][2];
			for (int ii = 0; ii < f.length; ii++) {
				double[] hvltf = hvl12teorAt27(kvp[j], f[ii]);
				hvlt[ii][0] = hvltf[0];
				hvlt[ii][1] = hvltf[1];
			}
			//

			for (int i = 0; i < n; i++) {
				eq[i] = hvl1e * mn[i][0] + mn[i][1];

				if (hvl1e <= hvlt[i + 1][0]) {
					hvl1kv = eq[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl1kv = eq[n - 1];
				}
			}

			hl1[j] = hvl1kv;

			for (int i = 0; i < n2; i++) {
				eq2[i] = hvl2e * mn2[i][0] + mn2[i][1];

				if (hvl2e <= hvlt[i + 1][1]) {
					hvl2kv = eq2[i];
					break;
				} else// a parcurs sirul si nu a gasit
				{
					hvl2kv = eq2[n - 1];
				}
			}

			hl2[j] = hvl2kv;
			// coeficientii de hvlvskv!!
			if (j > 0) {
				mm1[j - 1] = (hl1[j] - hl1[j - 1]) / (kvp[j] - kvp[j - 1]);
				mm2[j - 1] = (hl2[j] - hl2[j - 1]) / (kvp[j] - kvp[j - 1]);
				nn1[j - 1] = (hl1[j - 1] * kvp[j] - hl1[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				nn2[j - 1] = (hl2[j - 1] * kvp[j] - hl2[j] * kvp[j - 1])
						/ (kvp[j] - kvp[j - 1]);
				hlin1[j - 1] = mm1[j - 1] * kv + nn1[j - 1];
				hlin2[j - 1] = mm2[j - 1] * kv + nn2[j - 1];

				if (kv <= kvp[j]) {
					result[0] = hlin1[j - 1];
					result[1] = hlin2[j - 1];
					result[2] = (result[0] * p1 + result[1] * p2) / (p1 + p2);
					break;
				}
			}
			//
		}

		return result;
	}

	public static double[] eqfiltration(double kv, double anodAngle,
			double hvl1e, double hvl2e, double plane) {
		if (XRay.ICALC == 1) {
			boolean b = true;
			double[] rrr = new double[5];
			XRay.computeFiltrationFromHVL1(hvl1e);
			if (!XRay.ISOK) {
				b = XRay.ISOK;
			}
			XRay.computeFiltrationFromHVL2(hvl1e, hvl2e);
			XRay.ISOK = b;// /notify the error!!!
			// XRay.PLANE=20.0;
			XRay.computeFiltration();
			rrr[0] = XRay.eqFiltr_HVL1;
			rrr[1] = XRay.eqFiltr_HVL2;
			rrr[2] = XRay.eqFiltr;
			rrr[3] = XRay.p1_HVL;
			rrr[4] = XRay.p2_HVL;

			valid = true;
			return rrr;
		}

		valid = false;

		double[] result = new double[5];

		result[0] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[1] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[2] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[3] = -1.0;// daca nu se reuseste asta e valoarea de insucces
		result[4] = -1.0;// daca nu se reuseste asta e valoarea de insucces

		double uAnodMin = 7.0;
		double uAnodMed = 17.0;
		double uAnodMax = 27.0;

		if (anodAngle < uAnodMin || anodAngle > uAnodMax) {
			valid = false;
			return result;
		}

		if (anodAngle == uAnodMin) {
			double[] result7 = eqfiltrationAt7(kv, hvl1e, hvl2e, plane);// aci
																		// se
																		// valideste!!!
			if (validate())
				return result7;
			else
				return result;
		}

		if (anodAngle == uAnodMed) {
			double[] result17 = eqfiltrationAt17(kv, hvl1e, hvl2e, plane);// aci
																			// se
																			// valideste!!!
			if (validate())
				return result17;
			else
				return result;
		}

		if (anodAngle == uAnodMax) {
			double[] result27 = eqfiltrationAt27(kv, hvl1e, hvl2e, plane);// aci
																			// se
																			// valideste!!!
			if (validate())
				return result27;
		}

		if (anodAngle < uAnodMed) {
			double[] result7 = eqfiltrationAt7(kv, hvl1e, hvl2e, plane);// aci
																		// se
																		// valideste!!!
			if (validate()) {
				double[] result17 = eqfiltrationAt17(kv, hvl1e, hvl2e, plane);// aci
																				// se
																				// valideste!!!
				if (validate()) {
					result[0] = linInt(uAnodMed, result17[0], uAnodMin,
							result7[0], anodAngle);
					result[1] = linInt(uAnodMed, result17[1], uAnodMin,
							result7[1], anodAngle);
					result[2] = linInt(uAnodMed, result17[2], uAnodMin,
							result7[2], anodAngle);
					result[3] = linInt(uAnodMed, result17[3], uAnodMin,
							result7[3], anodAngle);
					result[4] = linInt(uAnodMed, result17[4], uAnodMin,
							result7[4], anodAngle);
					return result;
				} else
					return result;
			} else
				return result;
		}

		if (anodAngle < uAnodMax) {
			double[] result17 = eqfiltrationAt17(kv, hvl1e, hvl2e, plane);// aci
																			// se
																			// valideste!!!
			if (validate()) {
				double[] result27 = eqfiltrationAt27(kv, hvl1e, hvl2e, plane);// aci
																				// se
																				// valideste!!!
				if (validate()) {
					result[0] = linInt(uAnodMax, result27[0], uAnodMed,
							result17[0], anodAngle);
					result[1] = linInt(uAnodMax, result27[1], uAnodMed,
							result17[1], anodAngle);
					result[2] = linInt(uAnodMax, result27[2], uAnodMed,
							result17[2], anodAngle);
					result[3] = linInt(uAnodMax, result27[3], uAnodMed,
							result17[3], anodAngle);
					result[4] = linInt(uAnodMax, result27[4], uAnodMed,
							result17[4], anodAngle);
					return result;
				}
			}
		}

		return result;
	}

	// calculeaza m si n din cele 2 puncte p1(x1,y1) si p2(x2,y2)
	// si returneaza y(x).
	// pentru actualul mod nu se poate da eroare
	private static double linInt(double x1, double y1, double x2, double y2,
			double x) {
		double result = -1.0;
		double[] mn = new double[2];
		// insucces
		mn[0] = -1.0;// m
		mn[1] = -1.0;// n
		double num = x1 - x2;
		if (num != 0.0) {
			mn[0] = (y1 - y2) / num;
			mn[1] = (x1 * y2 - y1 * x2) / num;
			result = mn[0] * x + mn[1];
		}
		return result;
	}
}
