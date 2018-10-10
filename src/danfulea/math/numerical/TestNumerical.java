package danfulea.math.numerical;

//import jdf.math.numerical.Function;
//import jdf.math.numerical.SpecFunc;
//import jdf.math.numerical.TwoPointBoundaryValue;

/*
import jdf.math.numerical.Complex; 
import jdf.math.numerical.LinAEq;
import jdf.math.numerical.Interpolator;
import jdf.math.numerical.Integrator;

import jdf.math.numerical.EvalFunc;


import jdf.math.numerical.RandomNrs;
import jdf.math.numerical.Sorting;
import jdf.math.numerical.Graph;
import jdf.math.numerical.RootFind;
import jdf.math.numerical.MinMaxFunc;
import jdf.math.numerical.Eigensystem;
import jdf.math.numerical.FastFourierTransform;
import jdf.math.numerical.FFTAnalysis;
import jdf.math.numerical.Stats;
import jdf.math.numerical.ModelingData;
import jdf.math.numerical.OrdinaryDiffEq;
import jdf.math.numerical.TwoPointBoundaryValue;
import jdf.math.numerical.IntegralEq;
import jdf.math.numerical.PartialDiffEq;
*/
public class TestNumerical implements Function{
	
	public TestNumerical()	
	{
		System.out.println("asinh of 5 = " +SpecFunc.asinh(5));
		System.out.println("asinh of 5 = " +SpecFunc.asinh_ok(5));
//Eq. systms
/*
		int nn=3;int mm=1;
		double[][] acof=new double[nn][nn];
		acof[0][0]=5.0;acof[0][1]=3.0;acof[0][2]=2.0;
		acof[1][0]=2.0;acof[1][1]=-6.0;acof[1][2]=3.0;
		acof[2][0]=3.0;acof[2][1]=2.0;acof[2][2]=-4.0;

		double[][] bcof=new double[nn][mm];
		bcof[0][0]=17.0;bcof[1][0]=-1.0;bcof[2][0]=-5.0;
		//5x + 3y + 2z = 17 //x=1,y=2,z=3
		//2x - 6y + 3z = -1
		//3x + 2y - 4z = -5
//gaussj
		//LinAEq.gaussj(acof,nn,bcof,mm);//OK!!@@@@@@@@@@@1
//gauss
		//double detg = LinAEq.gauss(acof,nn,bcof,mm);//VERY OK!!@@@@@@@@@@@2
		//System.out.println(" detgauss= "+detg);//185
//LU
		int[] idx=new int[nn];
		double[] bcofv=new double[nn];
		bcofv[0]=17.0;bcofv[1]=-1.0;bcofv[2]=-5.0;
//		LinAEq.ludcmp(acof,nn,idx);//OK!!@@@@@@@3a
		//for (int i=0; i<nn; i++)
		//{
		//	System.out.println("i= "+i+" a[i][0]= "+acof[i][0]+"; a[i][1]= "+acof[i][1]+"; a[i][2]= "+acof[i][2]);//OK!
		//}

		//To summarize, this is the preferred way to solve the linear set of equations
		//A � x = b:
		//float **a,*b,d;
		//int n,*indx;
		//...
		//ludcmp(a,n,indx,&d);
		//lubksb(a,n,indx,b);
		//The answer x will be given back in b. Your original matrix A will have
		//been destroyed.
		//If you subsequently want to solve a set of equations with the same A but a
		//different right-hand side b, you repeat only
		//lubksb(a,n,indx,b);
		//not, of course, with the original matrix A, but with a and indx as were already
		//set by ludcmp.

//		LinAEq.lubksb(acof,nn,idx,bcofv);//OK!!@@@@@@@@@@3b
//		System.out.println("x= "+bcofv[0]+"; y= "+bcofv[1]+"; z= "+bcofv[2]);//OK!

		//System.out.println("x= "+bcof[0][0]+"; y= "+bcof[1][0]+"; z= "+bcof[2][0]);//OK!


		//double[][] dim = LinAEq.inverseMatrix(acof,nn,idx);//OK!!@@@@@@@4
		//for (int i=0; i<nn; i++)
		//{
		//	System.out.println("i= "+i+" dim[i][0]= "+dim[i][0]+"; dim[i][1]= "+dim[i][1]+"; dim[i][2]= "+dim[i][2]);//OK!
		//}

		//double det=LinAEq.determinantMatrix(acof,nn,idx);//OK!!@@@@@@@5
		//System.out.println("det= "+det);//=185
//improve LU iterativ odata!!@@@@@3c
		//double[] sol=LinAEq.improveLU(acof,nn,idx,bcofv);
		//System.out.println("improve x= "+sol[0]+"; y= "+sol[1]+"; z= "+sol[2]);//SUPER OK!
//Tridiagonal
		//tridiag:x=1,y=2,z=3,t=4
		//3x+2y+0z+0t=7
		//2x-6y+3z+0t=-1
		//0x+2y-4z+2t=0
		//0x+0y+3z-4t=-7
		int num=4;
		double[] a=new double[num];double[] b=new double[num];double[] c=new double[num];
		double[] r=new double[num];double[] u=new double[num];
		a[0]=0.0;//formaly!!
		a[1]=2.0;a[2]=2.0;a[3]=3.0;
		b[0]=3.0;b[1]=-6.0;b[2]=-4.0;b[3]=-4.0;
		c[0]=2.0;c[1]=3.0;c[2]=2.0;
		c[3]=0.0;//formaly
		r[0]=7.0;r[1]=-1.0;r[2]=0.0;r[3]=-7.0;
		//LinAEq.tridag(a,b,c,r,u,num);//OK!@@@6
		//System.out.println("orig x= "+u[0]+ " y= "+u[1]+ " z= "+u[2]+ " t= "+u[3]);
		//=>x= 2.2162162162162167 y= 2.6756756756756754 z= 3.5405405405405403 t= 4.405405405405405
		//LinAEq.triDiag(a,b,c,r,num);//OK!@@@6'
		//System.out.println(" x= "+r[0]+ " y= "+r[1]+ " z= "+r[2]+ " t= "+r[3]);
		//=>x= 2.2162162162162167 y= 2.6756756756756754 z= 3.5405405405405403 t= 4.405405405405405

//double d1=2.0;double d2=1.0;double d3=-2.0;
///d1=d2=d3=0.0;
//System.out.println(" d1= "+d1+ " d2= "+d2+ " d3= "+d3);
//double ss=15487.0;int i=0;
//for (ss=0.0,i=1;i<=5;i++)
//	ss+=i;
//System.out.println(" @ "+ss);
		double[] x = new double[nn];
		LinAEq.svd(acof, nn, bcofv, x);//OK!@@@7
		System.out.println("x= "+x[0]+"; y= "+x[1]+"; z= "+x[2]);//OK!
*/

//Interpolation
		/*int n=5;//5 points
		double[] xx={1.0,2.0,3.0,4.0,5.0};
		double[] yy={10.0,22.0,35.0,48.0,61.0};
		double xinput=4.5;

		double[] aa=new double[n];
		double[] bb=new double[n];
		double[] cc=new double[n];
		double[] dd=new double[n];

		Interpolator.set_spline(xx,yy,aa,bb,cc,dd,n);
		double rr=Interpolator.spline(xinput,xx,aa,bb,cc,dd,n);
		System.out.println(" spline= "+rr);

		Interpolator.polint(xx, yy, n, xinput);
		System.out.println(" poli: y(x)= "+Interpolator.ypoli+" +/- "+Interpolator.dypoli);
*/
//Neville ~ Lagrange
/*        int n=5;//5 points
		double[] xx={1.0,2.0,3.0,4.0,5.0};
		double[] yy={10.0,22.0,35.0,48.0,61.0};
		double xinput=4.5;

		//Interpolator.polint(xx, yy, n, xinput);
		//System.out.println(" poli: y(x)= "+Interpolator.ypoli+" +/- "+Interpolator.dypoli);
		//Interpolator.ratint(xx, yy, n, xinput);
		//System.out.println(" rat: y(x)= "+Interpolator.yrat+" +/- "+Interpolator.dyrat);

		//double[] yy2=new double[n];
		//Interpolator.spline_default(xx, yy, n, yy2);
		//Interpolator.splint(xx, yy, yy2, n, xinput);
		//System.out.println(" slpine: y(x)= "+Interpolator.yspline);

		//double x=5.5;
		//int j=Interpolator.locate(xx,n,x);
		//System.out.println(" locate= "+j);
		//int jj=Interpolator.hunt(xx,n,x,3);
		//System.out.println(" hunt= "+jj);
		//double[] cf = new double[n];
		//Interpolator.polcoe(xx, yy, n-1, cf);

		//for (int i=0; i<n; i++)
		//	System.out.println("i= "+i+" cf[i]= "+cf[i]);//OK!

		//double[] cf1 = new double[n];
		//Interpolator.polcof(xx, yy, n-1, cf1);

		//for (int i=0; i<n; i++)
		//	System.out.println("i= "+i+"cf[i], cf1[i]= "+cf[i]+" ; "+cf1[i]);//OK!

		double[] x1={10.5,12.5,13.5};//3=m
		double[] x2={1.5,2.5,3.5,4.5};//4=n;| m <-> - n
		double[][] ya={
						{12.5, 13.6, 14.8, 16.9},
						{22.5, 23.6, 24.8, 26.9},
						{32.5, 33.6, 34.8, 36.9},
					  };
		double xin1=11.5;
		double xin2=2.7;
		Interpolator.polin2(x1,x2,ya,3,4,xin1,xin2);
		System.out.println(" poli2: y(x)= "+Interpolator.ypoli2+" +/- "+Interpolator.dypoli2);

		double[][] y2a=new double[3][4];
		Interpolator.splie2(x1,x2,ya,3,4,y2a);
		Interpolator.splin2(x1,x2,ya,y2a,3,4,xin1,xin2);
		System.out.println(" spline2: y(x)= "+Interpolator.yspline2);
*/
//============================================================================
		/*double a=0.0;double b=1.0;
		//double a=1.0;double b=2.0;//=>midinf
		//double a=-0.3300;double b=-0.2300;//=>midexp
		int steps=25;
		Integrator itgr=new Integrator(this);
		double dd=0.0;
		//dd=itgr.trapzd(a,b,steps);
		//System.out.println(" integral trapz= "+dd);
		dd=itgr.qtrap(a,b);
		System.out.println(" integral qtrapz= "+dd);
		dd=itgr.qsimp(a,b);
		System.out.println(" integral qsimp= "+dd);
		dd=itgr.qromb(a,b);
		System.out.println(" integral qromb= "+dd);

		int steps2=10;
		//dd=itgr.midpnt(a,b,steps2);
		//System.out.println(" integral midpnt= "+dd);

		dd=itgr.qromo(a,b,"trapzd");
		System.out.println(" integral qromo trapzd= "+dd);
		dd=itgr.qromo(a,b,"midpnt");
		System.out.println(" integral qromo midpnt= "+dd);
		dd=itgr.qromo(a,b,"midinf");
		System.out.println(" integral qromo midinf= "+dd);
		dd=itgr.qromo(a,b,"midsql");
		System.out.println(" integral qromo midsql= "+dd);
		dd=itgr.qromo(a,b,"midsqu");
		System.out.println(" integral qromo midsqu= "+dd);
		dd=itgr.qromo(a,b,"midexp");
		System.out.println(" integral qromo midexp= "+dd);

		dd=itgr.qgaus(a,b);
		System.out.println(" integral qgauss= "+dd);

		//dd=itgr.gammln(8.98);
		//System.out.println(" gama= "+dd);

		dd=itgr.QD(a,b);
		System.out.println(" DCADRE= "+dd);
*/
//=================eval func
//P(x)=6.0+4x+2x^2+3x^3
 /*   	double[] c={6.0,4.0,2.0,3.0};
    	int nc=3;//degree 3
    	double x=1.0;//eval at 1.0
    	int nd=3;//three derivates!!
    	double[] pd=new double[nd+1];
		//EvalFunc.ddpoly(c,nc,x,pd,nd);
		//for (int i=0;i<=nd;i++)
		//{
		//	System.out.println("i= "+i+" ; value= "+pd[i]);
		//}
		double aa=1.0;
		double bb=-3.0;
		double cc=2.0;
		double[] r = EvalFunc.quadratic(aa,bb,cc);
		System.out.println("x1= "+r[0]+" x2= "+r[1]);//2,1 ok

		aa=1.0;
		bb=-6.0;
		cc=11.0;
		double dd=-6.0;
		r = EvalFunc.cubic(aa,bb,cc,dd);
		System.out.println("x1= "+r[0]+" x2= "+r[1]+" x3= "+r[2]);//1,3,2 OK!!

		EvalFunc ef=new EvalFunc(this);
		//double deri=ef.dfridr(1.0, 0.5);
		//System.out.println("deri = "+deri);//ok
		int ne=30;
		double[] ce=new double[ne];
		double ae=1.0;double be=6.0;
		ef.chebft(ae, be, ce, ne);
		double xe=3.0;

		//double res=ef.chebev(ae, be, ce, ne, xe);
		//System.out.println("cebi eval = "+res);//ok=-5

		double[] cder=new double[ne];
		ef.chder(ae, be, ce, cder, ne);
		res=ef.chebev(ae, be, cder, ne, xe);
		System.out.println("cebi der eval = "+res);//ok=-6!!

		double[] cint=new double[ne];
		ef.chint(ae, be, ce, cint, ne);
		res=ef.chebev(ae, be, cint, ne, xe);
		System.out.println("cebi int eval = "+res);//ok=-6!!
        double[] de = new double[ne];
		ef.chebpc(ce, de, ne);
		ef.pcshft(ae, be, de, ne);
		int ndeg=ne-1;
		EvalFunc.ddpoly(de,ndeg,xe,pd,nd);
		System.out.println(" ; value= "+pd[0]);//ok!!=-5!!*/
//================SPEC Function
		/*double dd=SpecFunc.gammln(15.00);
		System.out.println(" gama= "+dd);

		dd=SpecFunc.factrl(43);
		System.out.println(" fact= "+dd);

		dd=SpecFunc.factln(14);
		System.out.println(" factln= "+dd);

		dd=SpecFunc.bico(5,3);//C5,3=10
		System.out.println(" comb= "+dd);

		dd=SpecFunc.gammp(2,3);
		System.out.println(" p(a,x)= "+dd);
		dd=SpecFunc.gammq(2,3);
		System.out.println(" q(a,x)= "+dd);

		dd=SpecFunc.erff(2.3);
		System.out.println(" erff= "+dd);
		dd=SpecFunc.erffc(2.3);
		System.out.println(" erffc= "+dd);
		dd=SpecFunc.ERF1(2.3);
		System.out.println(" ERF1= "+dd);
		dd=SpecFunc.erfc1(2.3);
		System.out.println(" erfc1= "+dd);
		dd=SpecFunc.erfcc(2.3);
		System.out.println(" erfc= "+dd);

		dd=SpecFunc.expint(2,4.0);
		System.out.println(" en= "+dd);
		dd=SpecFunc.ei(4.0);
		System.out.println(" ei= "+dd);
		dd=SpecFunc.betai(0.2,1.0,0.5);
		System.out.println(" betai= "+dd);*/

		//double ddd=0.0;
		/*ddd=SpecFunc.bessj0(5.5);
		System.out.println(" bessj0= "+ddd);
		ddd=SpecFunc.bessy0(5.5);
		System.out.println(" bessy0= "+ddd);
		ddd=SpecFunc.bessj1(5.5);
		System.out.println(" bessj1= "+ddd);
		ddd=SpecFunc.bessy1(5.5);
		System.out.println(" bessy1= "+ddd);
		ddd=SpecFunc.bessj(5,5.5);
		System.out.println(" bessj= "+ddd);
		ddd=SpecFunc.bessy(5,5.5);
		System.out.println(" bessy= "+ddd);

		ddd=SpecFunc.bessi0(5.5);
		System.out.println(" bessi0= "+ddd);
		ddd=SpecFunc.bessk0(5.5);
		System.out.println(" bessk0= "+ddd);
		ddd=SpecFunc.bessi1(5.5);
		System.out.println(" bessi1= "+ddd);
		ddd=SpecFunc.bessk1(5.5);
		System.out.println(" bessk1= "+ddd);
		ddd=SpecFunc.bessi(5,5.5);
		System.out.println(" bessi= "+ddd);
		ddd=SpecFunc.bessk(5,5.5);
		System.out.println(" bessk= "+ddd);

		SpecFunc.bessjy(5.5,5);//bessjy(double x, double xnu)=>OK!!
		System.out.println(" bessjy=> j= "+SpecFunc.rjb+" y= "+SpecFunc.ryb+" true? "+SpecFunc.failB);
		SpecFunc.bessik(5.5,5);//bessjy(double x, double xnu)=>OK!!
		System.out.println(" bessik=> i= "+SpecFunc.rib+" k= "+SpecFunc.rkb+" true? "+SpecFunc.failB);

		SpecFunc.sphbes(5,5.5);//bessjy(double x, double xnu)=>OK!!
		System.out.println(" sphbes=> j= "+SpecFunc.sjb+" y= "+SpecFunc.syb+" true? "+SpecFunc.failB);
*/
		//ddd=SpecFunc.plgndr(3, 1, 0.5);
		//System.out.println(" legendre= "+ddd);

		/*Complex a=new Complex(1.,2.);
		Complex b=new Complex(3.,4.);
		Complex c=Complex.Cadd(a,b);
		System.out.println(" Re= "+c.r+" im= "+c.i);
		c=Complex.Csub(a,b);
		System.out.println(" Re= "+c.r+" im= "+c.i);
		c=Complex.Cmul(a,b);
		System.out.println(" Re= "+c.r+" im= "+c.i);*/

		/*SpecFunc.frenel(5.5);
		System.out.println(" frenel s= "+SpecFunc.s_fresnel+" c= "+SpecFunc.c_fresnel+" true? "+SpecFunc.failB);

		SpecFunc.cisi(5.5);
		System.out.println(" integral sin= "+SpecFunc.si_integral+" cos= "+SpecFunc.ci_integral+" true? "+SpecFunc.failB);

		ddd=SpecFunc.dawson(5.5);
		System.out.println(" dawson= "+ddd);*/
/*
		ddd=SpecFunc.rf(1.0, 2.0, 3.0);
		System.out.println(" rf= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.rd(1.0, 2.0, 3.0);
		System.out.println(" rd= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.rj(1.0, 2.0, 3.0,4.0);
		System.out.println(" rj= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.rc(1.0, 2.0);
		System.out.println(" rc= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.ellf(0.5, 1.2);
		System.out.println(" ellf= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.elle(0.5, 1.2);
		System.out.println(" elle= "+ddd+" true? "+SpecFunc.failB);
		ddd=SpecFunc.ellpi(0.5, 1.2,1.2);
		System.out.println(" ellpi= "+ddd+" true? "+SpecFunc.failB);
		SpecFunc.sncndn(1.0, 2.0);
		System.out.println(" ellipt jacob sn= "+SpecFunc.sn_jacob+
		" cn= "+SpecFunc.cn_jacob+
		" dn= "+SpecFunc.dn_jacob+" true? "+SpecFunc.failB);*/

		/*RandomNrs.idum_ran0=15;
		for (int i=1; i<=20; i++)
		{
			ddd=RandomNrs.ran0();
			System.out.println(" i= "+i+ " ran0= "+ddd);
		}
		RandomNrs.idum_ran1=15;
		for (int i=1; i<=20; i++)
		{
			ddd=RandomNrs.ran1();
			System.out.println(" i= "+i+ " ran1= "+ddd);
		}
		RandomNrs.idum_ran2=15;
		for (int i=1; i<=20; i++)
		{
			ddd=RandomNrs.ran2();
			System.out.println(" i= "+i+ " ran2= "+ddd);
		}

		RandomNrs.idum_ran3=15;
		for (int i=1; i<=20; i++)
		{
			ddd=RandomNrs.ran3();
			System.out.println(" i= "+i+ " ran3= "+ddd);
		}
		for (int i=1; i<=20; i++)
		{
			ddd=RandomNrs.RANDOMSET();
			System.out.println(" i= "+i+ " egs= "+ddd);
		}*/
		/*RandomNrs.idum_ran1=15;
		for (int i=1; i<=20; i++)
		{
			//ddd=RandomNrs.expdev();
			//System.out.println(" i= "+i+ " expdev= "+ddd);
			//ddd=RandomNrs.gasdev();
			//System.out.println(" i= "+i+ " gasdev= "+ddd);
			//ddd=RandomNrs.gamdev(2);
			//System.out.println(" i= "+i+ " gamdev= "+ddd);
			//ddd=RandomNrs.poidev(20.0);
			//System.out.println(" i= "+i+ " poidev= "+ddd);
			ddd=RandomNrs.bnldev(0.12,40);
			System.out.println(" i= "+i+ " bnldev= "+ddd);

		}*/
		/*double a=1.0;double b=2.0;
		Integrator itgr=new Integrator(this);
		double dd=0.0;
		//dd=itgr.qtrap(a,b);
		//System.out.println(" integral trapz= "+dd);
		dd=itgr.quad3d( 0.0, 1.0);
		System.out.println(" integral 3d= "+dd);*/

		/*int Ndim=1;
		double[] regn=new double[2*Ndim];
		regn[0]=1.0;regn[1]=2.0;

		RandomNrs.func=this;

		RandomNrs.vegas(regn, Ndim, 0, 1000,5, 0);
		RandomNrs.vegas(regn, Ndim, 1, 100000,1, 0);
		System.out.println(" integral vegas= "+RandomNrs.tgral_vegas+
		"; stdev= "+RandomNrs.sd_vegas+";  chi2= "+RandomNrs.chi2a_vegas);//OK!!!!

		//averageeul functiei pe intervalul [a,b]=>ok ~1.5!!
		RandomNrs.miser(regn, Ndim, 50, 0.0);
		System.out.println(" eval miser average= "+RandomNrs.ave_miser+
		"  stdev= "+RandomNrs.var_miser);//OK!!
		*/
//==============SORTING================
        /*int n=7;
		double[] sir={2.0,6.0,1.0,9.0,4.0,5.0,7.0};
		double[] sir2={21.0,16.0,11.0,19.0,41.0,15.0,71.0};
		int[] index=new int[7];
		int[] irk=new int[7];
		double[] hp=new double[7];
		//Sorting.piksrt(n, sir);
		//Sorting.shell(n, sir);
		//Sorting.sort(n, sir);
		//Sorting.sort2(n, sir, sir2);
		//Sorting.hpsort(n, sir);
		//Sorting.indexx(n, sir, index);Sorting.rank(n, index, irk);
		//double sel=Sorting.select(7 ,n, sir);//al doilea e 2!! OK!
		//for (int i=1; i<=n; i++)
		//{
		//	//System.out.println("i= "+i+" sir= "+sir[i-1]+" true? "+Sorting.failB);//OK, 1,2,...9
		//	System.out.println("i= "+i+" sir= "+sir[i-1]+" indx= "+index[i-1]+" arr[indx[j-1]-1]= "+sir[index[i-1]-1]+" irank "+irk[i-1]+" true? "+Sorting.failB);//OK, 1,2,...9
		//}
		//System.out.println(" sel ="+sel);

		//double sel=Sorting.selip(7 ,n, sir);
		//System.out.println(" sel ="+sel);
		Sorting.hpsel(3 ,n, sir, hp);//primele 3 maxime cu primul element= al 3-lea garantat!
		System.out.println(" 1 ="+hp[0]+ " 2= "+hp[1]+" 3= "+hp[2]);*/

		//new Graph(-4.0,4.0,this);
//ROOOOOOOOOOOOOOOOOOOT
	//Bracketing...
		/*RootFind rf=new RootFind(this);
		int itg=rf.zbrac(0.0,1.0);
		//System.out.println("i= "+itg+" x1= "+rf.x1_zbrac+" x2= "+rf.x2_zbrac);
		double deltaAccuracy=1.e-06;

		rf.x2_zbrac=rf.x2_zbrac+1.0;
	//Bracketing...
		if (itg==1)
		{
			ddd=rf.rtbis(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//bisect
			System.out.println("bisect sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.rtflsp(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//false position method
			System.out.println("false pos sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.rtsec(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//secant
			System.out.println("secant sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.zriddr(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//Ridders method
			System.out.println("ridder sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.zbrent(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//BRENT method
			System.out.println("zbrent sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.rtnewt(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//Newton method
			System.out.println("rtnewt sol= "+ddd+" true? "+Sorting.failB);
			ddd=rf.rtsafe(rf.x1_zbrac, rf.x2_zbrac, deltaAccuracy);//Newton - bisectmethod
			System.out.println("rtsafe sol= "+ddd+" true? "+Sorting.failB);

		}

		//a0+a1x+a2x^2;=>4+0*x+(-1)x^2=0

		int poligrd=2;
		Complex[] roots=new Complex[poligrd];
		Complex[] a={new Complex(4.0,0.0),new Complex(0.0,0.0),new Complex(-1.0,0.0)};

		rf.zroots(a, poligrd,  roots, 1);
		System.out.println("Laguerre poli sol1= "+roots[0].r+"; sol2= "+roots[1].r+" true? "+Sorting.failB);

		double[] cof={4.0,0.0,-1.0};
		double[] cof1={-1.0,0.0,4.0};//birge=>a0x^n+....
		double[] rtr=new double[poligrd];double[] rti=new double[poligrd];
		rf.zrhqr(cof,poligrd,rtr,rti);
		System.out.println("Eigenvalue Methods poli sol1= "+rtr[0]+"; sol2= "+rtr[1]+" true? "+Sorting.failB);
		double[] rtrr=rf.birge(cof1);
		System.out.println("Birge poli sol1= "+rtrr[0]+"; sol2= "+rtrr[1]+" true? "+Sorting.failB);
*/
		/*RootFind rf=new RootFind(this);
		double[] x0={0.8,1.2,2.2};//n=3
		double[] x01={0.8,1.2,2.2};//n=3
		double[] x02={0.8,1.2,2.2};//n=3
		rf.mnewt(100, x0, 3, 1.0e-06, 1.0e-06);
		System.out.println("nonlin mnewt x0= "+x0[0]+"; x1= "+x0[1]+"; x2= "+x0[2]);//OK

		rf.newt(x01, 3);
		System.out.println("nonlin newt x0= "+x01[0]+"; x1= "+x01[1]+"; x2= "+x01[2]);//OK

		rf.broydn(x02, 3);
		System.out.println("nonlin broydn x0= "+x02[0]+"; x1= "+x02[1]+"; x2= "+x02[2]);//OK
*/

		//MinMaxFunc mmf=new MinMaxFunc(this);
		/*mmf.mnbrak(-2.0, -1.0);
		System.out.println("mnbrak minim bracket ax= "+
				 MinMaxFunc.ax_mnbrak+" ;bx= "+MinMaxFunc.bx_mnbrak+
		"; cx= "+MinMaxFunc.cx_mnbrak+"; fa= "+MinMaxFunc.fa_mnbrak+
		"; fb= "+MinMaxFunc.fb_mnbrak+"; fc= "+MinMaxFunc.fc_mnbrak);

		ddd=mmf.golden(-2.0, -1.0, 4.0, 1.0e-6);
		System.out.println("golden minim bracket xmin= "+MinMaxFunc.xmin_golden+" fmin= "+ddd);
		ddd=mmf.brent(-2.0, -1.0, 4.0, 1.0e-6);
		System.out.println("brent minim bracket xmin= "+MinMaxFunc.xmin_brent+" fmin= "+ddd+" true? "+MinMaxFunc.failB);
		ddd=mmf.dbrent(-2.0, -1.0, 4.0, 1.0e-6);
		System.out.println("dbrent minim bracket xmin= "+MinMaxFunc.xmin_dbrent+" fmin= "+ddd+" true? "+MinMaxFunc.failB);
		*/
//===========================================================================
//===========================================================================
		/*
		int ndim=2;

		double x01=-2.0;double x11=-3.0;
		double[] x1=new double[ndim];
		x1[0]=x01;x1[1]=x11;
		double yy1=MF(x1);System.out.println(" yy= "+yy1);

		double x02=-1.0;double x12=-2.0;
		double[] x2=new double[ndim];
		x2[0]=x02;x2[1]=x12;
		double yy2=MF(x2);System.out.println(" yy2= "+yy2);

		double x03=4.0;double x13=3.0;
		double[] x3=new double[ndim];
		x3[0]=x03;x3[1]=x13;
		double yy3=MF(x3);System.out.println(" yy3= "+yy3);

		double[][] pmtx=new double[ndim+1][ndim];
		pmtx[0][0]=x1[0];pmtx[0][1]=x1[1];
		pmtx[1][0]=x2[0];pmtx[1][1]=x2[1];
		pmtx[2][0]=x3[0];pmtx[2][1]=x3[1];

		double[] yy=new double[ndim+1];
		yy[0]=yy1;yy[1]=yy2;yy[2]=yy3;

		double ftol=1.0e-06;

		mmf.amoeba(pmtx,yy,ndim,ftol);
		//System.out.println("2dim amoeba xmin= "+" true? "+MinMaxFunc.failB);
		for (int i=1;i<=ndim+1;i++)
		{
			for (int j=1;j<=ndim;j++)
				System.out.println("x min= "+i+" "+j+" = "+pmtx[i-1][j-1]);
			System.out.println("y min= "+i+" = "+yy[i-1]);
		}
		*/
		//int ndim=1;

		/*double x01=-2.0;
		double[] x1=new double[ndim];
		x1[0]=x01;
		double yy1=MF(x1);System.out.println(" yy= "+yy1);

		double x02=-1.0;
		double[] x2=new double[ndim];
		x2[0]=x02;
		double yy2=MF(x2);System.out.println(" yy2= "+yy2);

		double[][] pmtx=new double[ndim+1][ndim];
		pmtx[0][0]=x1[0];
		pmtx[1][0]=x2[0];

		double[] yy=new double[ndim+1];
		yy[0]=yy1;yy[1]=yy2;

		double ftol=1.0e-06;

		mmf.amoeba(pmtx,yy,ndim,ftol);
		//System.out.println("2dim amoeba xmin= "+" true? "+MinMaxFunc.failB);
		for (int i=1;i<=ndim+1;i++)
		{
			for (int j=1;j<=ndim;j++)
				System.out.println("x min= "+i+" "+j+" = "+pmtx[i-1][j-1]);
			System.out.println("y min= "+i+" = "+yy[i-1]);
		}*/
		/*double ftol=1.0e-06;
		double[] xinitialpoint=new double[ndim] ;
		xinitialpoint[0]=-2.0;
		double[][] directions=new double[ndim][ndim];
		directions[0][0]=1.0;*/

/*		mmf.powell(xinitialpoint,directions,ndim,ftol);
		for (int i=1;i<=ndim;i++)
		{
			System.out.println("powel x min= "+i+" = "+xinitialpoint[i-1]);
		}
		System.out.println("powel y min= "+MinMaxFunc.fret);
*/		//===OK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		/*double[] xinitialpoint1=new double[ndim] ;
		xinitialpoint1[0]=-2.0;
		mmf.frprmn(xinitialpoint1,ndim,ftol);
		for (int i=1;i<=ndim;i++)
		{
			System.out.println("frprmn x min= "+i+" = "+xinitialpoint1[i-1]);
		}
		System.out.println("frprmn y min= "+MinMaxFunc.fret_fr);*/

//		double[] direction=new double[ndim];
//		direction[0]=1.0;
//		mmf.dlinmin(xinitialpoint,direction,ndim);
//		for (int i=1;i<=ndim;i++)
//		{
//			System.out.println("dlinmin x min= "+i+" = "+xinitialpoint[i-1]);
//		}
//		System.out.println("dlinmin y min= "+MinMaxFunc.fret_dl);

		/*mmf.dfpmin(xinitialpoint,ndim,1.0e-05);//gtol
		for (int i=1;i<=ndim;i++)
		{
			System.out.println("dfpmin x min= "+i+" = "+xinitialpoint[i-1]);
		}
		System.out.println("dfpmin y min= "+MinMaxFunc.fret_df);*/

		//aij = aji==sy,metric
		/*double[][] aj=new double[2][2];
		aj[0][0]=1.0;aj[0][1]=3.0;
		aj[1][0]=3.0;aj[1][1]=2.0;
		int nj=2;
		double[] dj=new double[2];
		double[][] vj=new double[2][2];
		Eigensystem.jacobi(aj, nj, dj, vj);

		for (int i=1; i<=nj;i++)
		{
			System.out.println("eigvalue jacobi= "+dj[i-1]+" true? "+Eigensystem.failB);
			for (int j=1; j<=nj;j++)
			{
				System.out.println(i+" & "+j+"  ,eigvec jacobi= "+vj[i-1][j-1]);
			}
		}

		Eigensystem.eigsrt(dj, vj, nj);
		for (int i=1; i<=nj;i++)
		{
			System.out.println("eigvalue eigsrt= "+dj[i-1]+" true? "+Eigensystem.failB);
			for (int j=1; j<=nj;j++)
			{
				System.out.println(i+" & "+j+"  ,eigvec eigsrt= "+vj[i-1][j-1]);
			}
		}*/

/*
eigvalue jacobi= -1.5413812651491101 true? false
1 & 1  ,eigvec jacobi= 0.7630199824727256
1 & 2  ,eigvec jacobi= 0.6463748961301958==>OK verif 1.0*0.76-3.0*0.64=-1.54*0.76
eigvalue jacobi= 4.54138126514911 true? false
2 & 1  ,eigvec jacobi= -0.6463748961301958
2 & 2  ,eigvec jacobi= 0.7630199824727256==>OK verif 3.0*0.64+2.0*0.76=4.54*0.76
*/

//=========================================================================
		//double[] dat={12.5,12.8,13.5,14.1,15.6,11.5};
		//FFTAnalysis.avevar(dat,6);
		//System.out.println("ave= "+FFTAnalysis.ave_avevar+" var= "+FFTAnalysis.var_avevar);

		//int n1=6;int n2=6;
		//double[] dat={12.5,12.8,13.5,14.1,15.6,11.5};
		//double[] dat2={18.5,12.8,13.5,14.5,15.6,18.6};
		/*int n1=21;int n2=23;
		double[] dat={24.0,43.0,58.0,71.0,43.0,
					  49.0,61.0,44.0,67.0,49.0,
					  53.0,56.0,59.0,52.0,62.0,
					  54.0,57.0,33.0,46.0,43.0,
					  57.0};
		double[] dat2={42.0,43.0,55.0,26.0,62.0,
						37.0,33.0,41.0,19.0,54.0,
						20.0,85.0,46.0,10.0,17.0,
						60.0,53.0,42.0,37.0,42.0,
						55.0,28.0,48.0};
		*//*double std=0.0;double var1=0.0;double var2=0.0;
		Stats.avevar(dat,n1);
		std=Math.sqrt(Stats.var_avevar);var1=std*std;
		System.out.println("ave= "+Stats.ave_avevar+" var= "+Stats.var_avevar+" stdev= "+std);
		Stats.moment(dat2,n2);
		std=Math.sqrt(Stats.var_moment);var2=std*std;
		System.out.println("mom ave2= "+Stats.ave_moment+" mom var2= "+Stats.var_moment+" stdev= "+std);
		Stats.ttest(dat,n1,dat2,n2);
		System.out.println("eq.var. tcalc = "+Stats.t_ttest+" alpha probab = "+Stats.prob_ttest);
		//+"\n probab one tailed for signif. difference (%) = "+(1.0-Stats.prob_ttest/2.0));

		Stats.tutest(dat,n1,dat2,n2);
		System.out.println("non.eq.var tcalc = "+Stats.t_tutest+" alpha probab = "+Stats.prob_tutest);
		Stats.tptest(dat,dat2,n1);
		System.out.println("pair tcalc = "+Stats.t_tptest+" alpha probab = "+Stats.prob_tptest);
		Stats.ftest(dat,n1,dat2,n2);
		System.out.println("FTEST fcalc = "+Stats.f_ftest+" alpha probab = "+Stats.prob_ftest);

		//if(1.0-Stats.prob_ttest/2.0>0.95)
		//	System.out.println("Significant difference");
		//else
		//	System.out.println("NOT difference");

		System.out.println("================================================================ ");
		System.out.println(" 2tailed probab= "+Stats.prob_2tail_95fordiff_ttest+
							" 1 tailed probab= "+Stats.prob_1tail_95fordiff_ttest);
		System.out.println("	ttest	significance: "+Stats.significance_ttest);

		System.out.println(" 2tailed probab= "+Stats.prob_2tail_95fordiff_tutest+
							" 1 tailed probab= "+Stats.prob_1tail_95fordiff_tutest);
		System.out.println("	tutest	significance: "+Stats.significance_tutest);

		System.out.println(" 2tailed probab= "+Stats.prob_2tail_95fordiff_tptest+
							" 1 tailed probab= "+Stats.prob_1tail_95fordiff_tptest);
		System.out.println("	tptest	significance: "+Stats.significance_tptest);

		System.out.println(" 2tailed probab= "+Stats.prob_2tail_95fordiff_ftest+
							" 1 tailed probab= "+Stats.prob_1tail_95fordiff_ftest);
		System.out.println("	ftest	significance: "+Stats.significance_ftest);

		System.out.println("================================================================ ");

		double df=10;
		double tint=1.812;//1.372;//1.812;//2.634;//2.228;//2.228;//0.700;//3.16927;//1.812;
		System.out.println("================================================================ ");
		double dd=SpecFunc.betai(0.5*df,0.5,df/(df+(tint)*(tint)));
		System.out.println("alpha= "+dd+" twotailed p%= "+(1.0-dd)+" onetailed %= "+(1.0-dd/2.0));//two tailed!!


		double f_ftest=0.0;
		double df1=0.0;double df2=0.0;
		if (var1 > var2)
		{// Make F the ratio of the larger variance to the smaller one.
			f_ftest=var1/var2;
			df1=n1-1;
			df2=n2-1;
		}
		else
		{
			f_ftest=var2/var1;
			df1=n2-1;
			df2=n1-1;
		}
		df1=df2=12;
		f_ftest=2.69;//4.16;//10.97;//5.05;
		dd = 2.0*SpecFunc.betai(0.5*df2,0.5*df1,df2/(df2+df1*(f_ftest)));
		if (dd > 1.0) dd=2.0-dd;
		dd=dd/2;
		System.out.println("alpha= "+dd+"  "+f_ftest+" p% "+(1.0-dd));
		*/
		/*double[] xx={1.0,2.0,3.0,4.0,5.0};
		double[] yy={10.0,20.0,30.0,40.0,50.0};
		int nn=5;
		double[] sig={0.01,0.02,0.015,0.014,0.05};

		double[] xx1={1.01,2.01,3.01,4.01,5.01};
		double[] yy1={10.0,20.0,30.0,40.0,50.0};
		int nn1=5;
		double[] sig1={0.01,0.02,0.015,0.014,0.05};

		ModelingData.fit(xx,yy,nn,sig,1);//1);
		System.out.println(" a= "+ModelingData.a_fit+" b= "+ModelingData.b_fit);
		System.out.println(" sa= "+ModelingData.siga_fit+" sb= "+ModelingData.sigb_fit);
		System.out.println(" chi2= "+ModelingData.chi2_fit+" q goodness of fit= "+ModelingData.q_fit);

		ModelingData.medfit(xx1,yy1,nn1);//not always works!!!
		System.out.println(" amf= "+ModelingData.a_medfit+" bmf= "+ModelingData.b_medfit);
		System.out.println(" abdev= "+ModelingData.abdev_medfit);*/

		/*ModelingData.fitexy(xx,yy,nn,sig,sig);
		System.out.println(" a= "+ModelingData.a_fitexy+" b= "+ModelingData.b_fitexy);
		System.out.println(" sa= "+ModelingData.siga_fitexy+" sb= "+ModelingData.sigb_fitexy);
		System.out.println(" chi2= "+ModelingData.chi2_fitexy+" q goodness of fit= "+ModelingData.q_fitexy);
		*///poli fit
		/*ModelingData.func=this;
		//f(x)=1.0+2.0x+3.0x^2=>
		//f(x)=c-bx;=lny=4.0-2.0*x=>y=54.59exp(-2.0x);
		//f(x)=c-bx^2;=lny=4.0-2.0*x*x=>y=54.59exp(-2.0x^2);
//=2*EXP(-(B1-1)*(B1-1)/(3*3))+3*EXP(-(B1-2)*(B1-2)/(4*4))+4*EXP(-(B1-3)*(B1-3)/(5*5))
		int ndat=6;
		double[] xxx={1.0,2.0,3.0,4.0,5.0,6.0};
		//double[] yyy={6.0,17.0,34.0,57.0,86.0,121.0};//OK
		//double[] yyy={2.0,0.0,-2.0,-4.0,-6.0,-8.0};//OK
		double[] yyy={8.23,8.63,8.10,6.91,5.45,4.02};
		double[] ssig={0.001,0.002,0.0015,0.0014,0.005,0.004};//sigma
		//int mma=3;//1.0,2.0,3.0 =1+ordin poli( 2 )
		int mma=3*3;
		double[] acof={2.0,1.0,3.0,3.0,2.0,4.0,4.0,3.0,5.0};
		//double[] acof=new double[mma];
		//double[] acof={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
		int[] iia={1,1,1,1,1,1,1,1,1};
		double[][] ccovar=new double[mma][mma];
		double[][] alph=new double[mma][mma];*/
		/*double[] acof=new double[mma];//null
		//double[] acof={0.0,2.0,0.0};//null
		int[] iia={1,1,1};//{1,0,1};//all fit=>OK!!
		double[][] ccovar=new double[mma][mma];//covar matrix
		/*ModelingData.lfit(xxx, yyy, ssig, ndat, acof, iia, mma, ccovar);
		for (int i=1; i<=mma; i++)
		{
			System.out.println(" poli i= "+i+" a[i-1]= "+acof[i-1]);
		}//OK!!!
		System.out.println(" chi2=  "+ModelingData.chisq_lfit);
		double[][] uu=new double[ndat][mma];
		double[][] vv=new double[mma][mma];
		double[] ww=new double[mma];
		ModelingData.svdfit(xxx, yyy, ssig, ndat, acof, mma,uu,vv,ww);
		for (int i=1; i<=mma; i++)
		{
			System.out.println(" poli i= "+i+" a[i-1]= "+acof[i-1]);
		}//OK!!!
		System.out.println(" chi2=  "+ModelingData.chisq_svdfit);*/
	//ModelingData.mrqmin(xxx, yyy, ssig, ndat, acof,iia, mma, ccovar, alph,-1.0);
	//ModelingData.mrqmin(xxx, yyy, ssig, ndat, acof,iia, mma, ccovar, alph,ModelingData.alamda_mrqmin);
	//ModelingData.mrqmin(xxx, yyy, ssig, ndat, acof,iia, mma, ccovar, alph,ModelingData.alamda_mrqmin);
	//	for (int i=1; i<=mma; i++)
	//	{
	//		System.out.println(" poli i= "+i+" a[i-1]= "+acof[i-1]);
	//	}//OK!!!
	//	System.out.println(" chi2=  "+ModelingData.chisq_mrqmin+" conv? "+ModelingData.convB);

//y(x; a) is the sum of na/3 Gaussians (15.5.16). The amplitude, center, and width of the
	//Gaussians are stored in consecutive locations of a: a[i] = Bk, a[i+1] = Ek, a[i+2] = Gk,
	//k = 1, ..., na/3. The dimensions of the arrays are a[1..na], dyda[1..na].



		/*double[][] f=new double[NMAX][2];
		double[] t=new double[2*NMAX];
		 f[0][0]=1;f[0][1]=0; f[1][0]=2;f[1][1]=0;
		 f[2][0]=3;f[2][1]=0; f[3][0]=4;f[3][1]=0; f[4][0]=5;f[4][1]=0;
		 f[5][0]=6;f[5][1]=0; f[6][0]=7;f[6][1]=0; f[7][0]=8;f[7][1]=0;
		 int r=3;
		 int n=(int)Math.pow(2,r);
		System.out.println("Transformata Forier :\n");
		 TFR(3,-1,f);
		 for (int i=1;i<=n;i++)
		  {
		//  //printf("F%d = ( %12.6e ) + ( %12.6e )*i\n",i,f[i-1][0],f[i-1][1]);
		  System.out.println("i= "+i+"  [0]= "+f[i-1][0]+"  [1]= "+f[i-1][1]);
 		 }
		 t[0]=1;t[1]=0; t[2]=2;t[3]=0;
		 t[4]=3;t[5]=0; t[6]=4;t[7]=0; t[8]=5;t[9]=0;
		 t[10]=6;t[11]=0; t[12]=7;t[13]=0; t[14]=8;t[15]=0;
		 FastFourierTransform.realft(t,2*NMAX,-1);
		 int k=1;
		 for (int i=1;i<=2*NMAX;i++)
		  {
			 //printf("F%d = ( %12.6e ) + ( %12.6e )*i\n",i,f[i-1][0],f[i-1][1]);
			 if(i % 2 ==0)
			 {
		  		System.out.println("i= "+k+"  Re= "+t[i-2]+"  Im= "+t[i-1]);
				k++;
		  	 }
		  	// System.out.println("i= "+k+"  ReIM= "+t[i-1]);
 		 }*/

/*
 A difference between two means is significant (at the given probability level)
 if the calculated t value is greater than the value given in this table.
 A probability of p = 0.05 (95% probability of making a correct statement)=TWOtailed=
 is usually acceptable for biological work, but p = 0.1 can be used for a "one-tailed" t-test.

 OR:
 Remember that P <0.05 is the arbitrary value that is generally accepted to be significant
(There must be less than a 5% possibility that the difference between means is due to chance.)
*/
//=========================difeq=ordinary
		/*OrdinaryDiffEq.func=this;
		int nvr=2;//1;//2;
		int nsteps=10000;
		double[] initCond_y0={1.0,1.0};//nvar
		double x0_init=0.0;
		double xf_propag=Math.PI/4.0;
		OrdinaryDiffEq.rkdumb_init(nvr,nsteps);
		OrdinaryDiffEq.rkdumb(initCond_y0,nvr,x0_init,xf_propag,nsteps);
		//System.out.println(" sol [0] at x=xinit = "+OrdinaryDiffEq.y_rkdumb[0][0]);
		//System.out.println(" der [1] at x=xinit = "+OrdinaryDiffEq.y_rkdumb[1][0]);
		System.out.println(" sol [0] at x=xfinal = "+OrdinaryDiffEq.y_rkdumb[0][nsteps]);
		System.out.println(" der [1] at x=xfinal = "+OrdinaryDiffEq.y_rkdumb[1][nsteps]);
		//System.out.println(" fail? = "+OrdinaryDiffEq.failB);

		//0.01=first guess of stepzize!
		//kmax_odeint=0.0 defaulrt and dxsav_odeint
		//OrdinaryDiffEq.odeint(initCond_y0, nvr, x0_init, xf_propag, 1.e-05, 0.01, 0.0);
		//System.out.println(" sol [0] at x=xfinal = "+initCond_y0[0]);
		//System.out.println(" der [1] at x=xfinal = "+initCond_y0[1]);

		double[] y_at_xo={1.0};//y=f(x0)
		double[] yp_at_xo={0.0};//yp=f'(x0)
		double[] yot=new double[nvr];
		double htot=xf_propag;
		//OrdinaryDiffEq.rk4(y_at_xo,yp_at_xo,nvr, x0_init, htot, yot);
		//System.out.println(" sol rk4 [0] at x=xfinal = "+yot[0]);
		//OrdinaryDiffEq.mmid(y_at_xo,yp_at_xo,nvr, x0_init, htot, nsteps, yot);
		//System.out.println(" sol mmid [0] at x=xfinal = "+yot[0]);
		//double[] yscal={100.0};//y=f(x0)
		//OrdinaryDiffEq.bsstep(y_at_xo,yp_at_xo,nvr,x0_init,xf_propag,1.e-06,yscal);
		//System.out.println(" sol bsstep [0] at x=xfinal = "+y_at_xo[0]+"  "+OrdinaryDiffEq.xx_bsstep+" fail? "+OrdinaryDiffEq.failB);

//y''=2x
		double x00=0.0;
		double xfin=Math.PI/4.0;
		double[] yy={1.0,0.0};//f and f' at xo; f=x^3/3+Const
		double[] ddy={0.0};//=funx(x,y)=2x
		double htott=xfin-x00;
		double[] yott=new double[2];
//		OrdinaryDiffEq.stoerm(yy,ddy,2,x00,htott,10000,yott);
//		System.out.println(" sol stoerm [0] at xfin = "+yott[0]);
//		System.out.println(" derivata 1 at stoerm [1] at xfin = "+yott[1]);
		//double[] y3={1.0,1.0,0.0};
		double[] y3={1.0,1.0};double[] y33={1.0,1.0};
		double[] dy3={1.0,-1.0};//new double[3];//{1.0,1.0,0.0};
		double[] dy33={1.0,-1.0};
		//dy3[0]=-.013*y3[0] - 1000.0*y3[0]*y3[2];
		//dy3[1]=- 2500.0*y3[1]*y3[2];
		//dy3[2]=-.013*y3[0] - 1000.0*y3[0]*y3[2]- 2500.0*y3[1]*y3[2];
		double[] ysc3={1.0,1.0};//{1.0,1.0,1.0};
		double[] ysc33={1.0,1.0};
		double x3_0=0.0;
		double htry3=Math.PI/4.0;//0.20;//xfin;//50.0;
		//OrdinaryDiffEq.stiff(y3, dy3, 3, x3_0, htry3, 1.0e-04, ysc3,0);
		OrdinaryDiffEq.stiff(y3, dy3, 2, x3_0, htry3, 1.0e-04, ysc3,1);
		System.out.println(" sol stiff [0] at xfin = "+y3[0]);
		System.out.println(" sol stiff [1] at xfin = "+y3[1]+"  xfin= "+OrdinaryDiffEq.x_stiff+"  pi/4 "+xfin);
//DO NOT REACH TO SOLUTION BREAK =>MUST BE CALLED MULTIPLE TIME!!

		OrdinaryDiffEq.stifbs(y33, dy33, 2, x3_0, htry3, 1.0e-04,ysc33, 1);
		System.out.println(" sol stifbs [0] at xfin = "+y33[0]+"  xfin= "+OrdinaryDiffEq.xx_stifbs+"  pi/4 "+xfin);
		System.out.println(" deriv stifbs [1] at xfin = "+y33[1]);*/
		//OK WORKS!!!!

//=========================================================

		TwoPointBoundaryValue.func=this;
		int NE =3;
		int M =41;
		int NB =1;
		int NSI =NE;
		int NYJ =NE;
		int NYK =M;
		int NCI =NE;
		int NCJ =(NE-NB+1);
		int NCK =(M+1);
		int NSJ =(2*NE+1);

		mpt=M;
		x=new double[M+1];
//Sample program using solvde. Computes eigenvalues of spheroidal harmonics Smn(x; c) for
//m ? 0 and n ? m. In the program, m is mm, c2 is c2, and ? of equation (17.4.20) is anorm.
	int i,itmax,k;
	int[] indexv=new int[NE+1];
	double conv,deriv,fac1,fac2,q1,slowc;
	double[] scalv=new double[NE+1];
	//float **y,**s,***c;
	double[][] y=new double[NYJ][NYK];//matrix(1,NYJ,1,NYK);
	double[][]s=new double[NSI][NSJ];//matrix(1,NSI,1,NSJ);
	double[][][]c=new double[NCI][NCJ][NCK];//f3tensor(1,NCI,1,NCJ,1,NCK);
	itmax=100;
	conv=5.0e-6;
	slowc=1.0;
	h=1.0/(M-1);
	//printf("\nenter m n\n");
	//scanf("%d %d",&mm,&n);
	mm=2;
	n=2;
	c2=4.0;
	/*
m n 	c2 Lambdaexact Lambdasfroid
2 2 	0.1 6.01427 6.01427
		1.0 6.14095 6.14095
		4.0 6.54250 6.54253
2 5 	1.0 30.4361 30.4372
		16.0 36.9963 37.0135
4 11	-1.0 131.560 131.554
	*/
	if ((n+mm & 1)!=0)//if (n+mm & 1)
	{// No interchanges necessary.
		indexv[0]=1;
		indexv[1]=2;
		indexv[2]=3;
	}
	else
	{// Interchange y1 and y2.
		indexv[0]=2;
		indexv[1]=1;
		indexv[2]=3;
	}
	anorm=1.0; //Compute ?.
	if (mm!=0)//if (mm)
	{
		q1=n;
		for (i=1;i<=mm;i++) anorm = -0.5*anorm*(n+i)*(q1--/i);
	}
	for (k=1;k<=(M-1);k++)
	{// Initial guess.
		x[k-1]=(k-1)*h;//x[k]=(k-1)*h;
		fac1=1.0-x[k-1]*x[k-1];//fac1=1.0-x[k]*x[k];
		fac2=Math.exp((-mm/2.0)*Math.log(fac1));
		//y[1][k]=plgndr(n,mm,x[k])*fac2; Pmn from �6.8.
		y[0][k-1]=SpecFunc.plgndr(n,mm,x[k-1])*fac2;
		//deriv = -((n-mm+1)*plgndr(n+1,mm,x[k])- Derivative of Pmn from a recurrencerelation.
		//(n+1)*x[k]*plgndr(n,mm,x[k]))/fac1;
		deriv = -((n-mm+1)*SpecFunc.plgndr(n+1,mm,x[k-1])-
				(n+1)*x[k]*SpecFunc.plgndr(n,mm,x[k-1]))/fac1;
		//y[2][k]=mm*x[k]*y[1][k]/fac1+deriv*fac2;
		y[1][k-1]=mm*x[k-1]*y[0][k-1]/fac1+deriv*fac2;
		//y[3][k]=n*(n+1)-mm*(mm+1);
		y[2][k-1]=n*(n+1)-mm*(mm+1);
	}
	x[M-1]=1.0;//x[M]=1.0; Initial guess at x = 1 done separately.
	y[0][M-1]=anorm;//y[1][M]=anorm;
	y[2][M-1]=n*(n+1)-mm*(mm+1);//y[3][M]=n*(n+1)-mm*(mm+1);
	//y[2][M]=(y[3][M]-c2)*y[1][M]/(2.0*(mm+1.0));
	y[1][M-1]=(y[2][M-1]-c2)*y[0][M-1]/(2.0*(mm+1.0));
	scalv[0]=Math.abs(anorm);//scalv[1]=fabs(anorm);
	//scalv[2]=(y[2][M] > scalv[1] ? y[2][M] : scalv[1]);
	scalv[1]=(y[1][M-1] > scalv[0] ? y[1][M-1] : scalv[0]);
	//scalv[3]=(y[3][M] > 1.0 ? y[3][M] : 1.0);
	scalv[2]=(y[2][M-1] > 1.0 ? y[2][M-1] : 1.0);
	//for (;;)
	//{

		//printf("\nEnter c**2 or 999 to end.\n");
		//scanf("%f",&c2);
		//if (c2 == 999)
		//{
			//free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
			//free_matrix(s,1,NSI,1,NSJ);
			//free_matrix(y,1,NYJ,1,NYK);
			//return;// 0;
		//	break;
		//}
		TwoPointBoundaryValue.solvde(itmax,conv,slowc,scalv,indexv,NE,NB,M,y,c,s);
		System.out.println("========results RELAXATION=================");
		System.out.println(" m= "+mm+"; n= "+n+"; c**= "+c2+"; lambda= "+(y[2][1]+mm*(mm+1)));
	//	printf("\n %s %2d %s %2d %s %7.3f %s %10.6f\n",
	//	"m =",mm," n =",n," c**2 =",c2,
	//	" lamda =",y[3][1]+mm*(mm+1));
	//} //Return for another value of c2.
	//}
	/////////////WORKS++++++++++++++++++++++++++++
	shootfB=false;
	//int check;//,i;
	//double q1;//,*v;
	double[] v=new double[N2];//(1,N2);//N2=1
	dx=1.0e-4; //Avoid evaluating derivatives exactly at x = -1.
	nvar=3;// Number of equations.
	TwoPointBoundaryValue.nvar=nvar;
	//for (;;)
	//{
	//printf("input m,n,c-squared\n");
	//if (scanf("%d %d %f",&m,&n,&c2) == EOF) break;
	m=mm;
	nn=n;
	c22=c2;
	//if (nn < m || m < 0) continue;
	gmma=1.0; //Compute ? of equation (17.4.20).
	q1=nn;
	for (i=1;i<=m;i++) gmma *= -0.5*(nn+i)*(q1--/i);
	//v[1]=n*(n+1)-m*(m+1)+c2/2.0; Initial guess for eigenvalue.
	v[0]=nn*(nn+1)-m*(m+1)+c22/2.0;
	x1 = -1.0+dx; //Set range of integration.
	x2=0.0;
	TwoPointBoundaryValue.x1=x1;
	TwoPointBoundaryValue.x2=x2;
	TwoPointBoundaryValue.newt(v,N2,"shoot");//,&check,TwoPointBoundaryValue.shoot);// Find v that zeros function f in score.
System.out.println("===========Results Shoot===========");
	if (TwoPointBoundaryValue.check_ln!=0)
	{
		//printf("shoot failed; bad initial guess\n");
		System.out.println("shoot failed; bad initial guess\n");
	}
	else
	{
	//printf("\tmu(m,n)\n");
	//printf("%12.6f\n",v[1]);
		System.out.println("tmu(m,n)= m: "+m+" n:= "+nn+" c2= "+c22+" lambda= "+(v[0]+m*(m+1)));
	}
	//WORKS============================!!!!!!!!!!!!!!!!!!!!!!!!!!
	shootfB=true;
System.out.println("===========Resluts SHOOTF===========");
	//int check,i;
	//float q1,*v1,*v2,*v;
	v=new double[NTOT];//vector(1,NTOT);
	v1=new double[NTOT];//[1..N22]//N22=1
	v2=new double[NTOT];//[1..N11]//N11=2

	nvarr=NTOT; //Number of equations.
	TwoPointBoundaryValue.nvar_shootf=nvarr;
	mmm=mm;//=m
	nnn=n;
	c22=c2;


	nn2=N22;TwoPointBoundaryValue.nn2_shootf=nn2;
	dx=DXX; //Avoid evaluating derivatives exactly at x =
	//�1. for (;;) {
	//printf("input m,n,c-squared\n");
	//if (scanf("%d %d %f",&m,&n,&c2) == EOF) break;
	//if (n < m || m < 0) continue;
	gmma=1.0; //Compute ? of equation (17.4.20).
	q1=nnn;
	for (i=1;i<=mmm;i++) gmma *= -0.5*(nnn+i)*(q1--/i);
	v1[0]=nnn*(nnn+1)-mmm*(mmm+1)+c22/2.0; //Initial guess for eigenvalue and function value.
	v2[1]=v1[0];
	v2[0]=gmma*(1.0-(v2[1]-c22)*dx/(2*(mmm+1)));
	//construct V array!!!!!!!!!!
	for(i=0;i<NTOT;i++)	v[i]=v1[i];
	int jjj=0;
	for(i=N22;i<=NTOT;i++)
	{
		v[i-1]=v2[jjj];
		jjj++;
	}
	//end construct V array!!!!!!!!!!
	x11 = -1.0+dx; //Set range of integration.
	x22=1.0-dx;
	xff=0.0; //Fitting point.
	TwoPointBoundaryValue.x1_shootf=x11;
	TwoPointBoundaryValue.x2_shootf=x22;
	TwoPointBoundaryValue.xf_shootf=xff;

	TwoPointBoundaryValue.newt(v,NTOT,"shootf");//,&check,shootf); //Find v that zeros function f in score.
	if (TwoPointBoundaryValue.check_ln!=0)
	{
	//printf("shootf failed; bad initial guess\n");
	System.out.println("shootf failed; bad initial guess\n");
	}
	else
	{
	//printf("\tmu(m,n)\n");
	//printf("%12.6f\n",v[1]);
	System.out.println("shootf: tmu(m,n) m: "+mmm+", nnn: "+nnn+", c2: "+c22+" lambda= "
	+(v[0]+mmm*(mmm+1)));
	}
	//}
	//free_vector(v,1,NTOT);
	//return 0;
	//}
	//}
	//free_vector(v,1,N2);
//========================================================



		//yot=derivF(2.0,y_at_xo);//works
		//System.out.println(" test = "+yot[0]);
//============TEST=====================================
		//test(y_at_xo,yot);
		//System.out.println(" test = "+yot[0]);
		//============OK================================================================
	}
	////////////
	int mm=0;
	int n=0;
	int mpt=0;//M;
	double h=0;//,
	double c2=0.0;
	double anorm=0.0;
	double[] x;//[M+1];
	boolean shootfB=false;

	int N2 =1;
	int m,nn; //Communicates with load, score, and derivs.
	double c22,dx,gmma;
	int nvar;// Communicates with shoot.
	double x1,x2;

	int N11 =2;
	int N22 =1;
	int NTOT =(N11+N22);
	double DXX =1.0e-4;
	int mmm,nnn; //Communicates with load1, load2, score,
	//and derivs.
	//float c2,dx,gmma;
	int nn2,nvarr;// Communicates with shootf.
	double x11,x22,xff;
	double[] v1;
	double[] v2;
//===========================
	public void test(double[]yy ,double[] yout)
	{
		//yout=dervF(2.0);
		yout=dervF(2.0,yout);
		yout[0]=yout[0]+89.89;
		System.out.println(" inside = "+yout[0]);
	}
	//public double[] dervF(double x)
	public double[] dervF(double x, double[] f)
    {
		//double[] f=new double[1];
//??NOTE: this new double destroy the object=>on return may be compromised in the anove situation!!
		f[0]=2.*x;
		return f;
	}
//===========================end TEST=================================================
    public void load(double x1, double[] v, double[] y)
    {
		double y1 = ((nn-m & 1)!=0 ? -gmma : gmma);
		y[2]=v[0];//y[3]=v[1];
		y[1] = -(y[2]-c22)*y1/(2*(m+1));//y[2] = -(y[3]-c2)*y1/(2*(m+1));
		y[0]=y1+y[1]*dx;//y[1]=y1+y[2]*dx;
	}

    public void load1(double x1, double[] v, double[] y)
    {
		double y1 = ((nnn-mmm & 1)!=0 ? -gmma : gmma);
		y[2]=v[0];//y[2]=v1[0];
		y[1] = -(y[2]-c22)*y1/(2*(mmm+1));
		y[0]=y1+y[1]*dx;
	}

    public void load2(double x2, double[] v, int nn2, double[] y)
    {
		//y[2]=v2[1];
		//y[0]=v2[0];
		y[2]=v[N22+1];//y[3]=v2[2];
		y[0]=v[N22];//y[1]=v2[1];//y[2]=(y[3]-c2)*y[1]/(2*(m+1));
		//prima valoare a lui v2 este v[N22]!!!!
		y[1]=(y[2]-c22)*y[0]/(2*(mmm+1));
	}

    //public double[] score(double x2, double[] y)
    public void score(double x2, double[] y, double[] f)
    {
		//double[] res=new double[10];//nvar=> to be implemented
		//f[1]=(nn-m & 1 ? y[1] : y[2]);
		if(!shootfB)
			f[0]=((nn-m & 1)!=0 ? y[0] : y[1]);
		else
			for (int i=1;i<=3;i++) f[i-1]=y[i-1];
		//return res;
	}

	public void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
		int indexv[], int ne, double[][] s, double[][] y)
	{
		//Returns matrix s for solvde.
//==================================
		double temp,temp1,temp2;
		if (k == k1)
		{// Boundary condition at first point.
			if ((n+mm & 1)!=0)
			{
				//s[3][3+indexv[1]]=1.0; Equation (17.4.32).
				s[2][2+indexv[0]]=1.0;
				s[2][2+indexv[1]]=0.0;//s[3][3+indexv[2]]=0.0;
				s[2][2+indexv[2]]=0.0;//s[3][3+indexv[3]]=0.0;
				s[2][jsf-1]=y[0][0];//s[3][jsf]=y[1][1]; Equation (17.4.31).
			}
			else
			{
				//s[3][3+indexv[1]]=0.0; Equation (17.4.32).
				s[2][2+indexv[0]]=0.0;
				s[2][2+indexv[1]]=1.0;//s[3][3+indexv[2]]=1.0;
				s[2][2+indexv[2]]=0.0;//s[3][3+indexv[3]]=0.0;
				s[2][jsf-1]=y[1][0];//s[3][jsf]=y[2][1]; Equation (17.4.31).
			}
		}
		else if (k > k2)
		{// Boundary conditions at last point.
			//s[1][3+indexv[1]] = -(y[3][mpt]-c2)/(2.0*(mm+1.0)); (17.4.35).
			s[0][2+indexv[0]] = -(y[2][mpt-1]-c2)/(2.0*(mm+1.0));
			s[0][2+indexv[1]]=1.0;//s[1][3+indexv[2]]=1.0;
			s[0][2+indexv[2]] = -y[0][mpt-1]/(2.0*(mm+1.0));//s[1][3+indexv[3]] = -y[1][mpt]/(2.0*(mm+1.0));
			s[0][jsf-1]=y[1][mpt-1]-(y[2][mpt-1]-c2)*y[0][mpt-1]/(2.0*(mm+1.0));//s[1][jsf]=y[2][mpt]-(y[3][mpt]-c2)*y[1][mpt]/(2.0*(mm+1.0)); (17.4.33).
			s[1][2+indexv[0]]=1.0;//s[2][3+indexv[1]]=1.0; Equation (17.4.36).
			s[1][2+indexv[1]]=0.0;//s[2][3+indexv[2]]=0.0;
			s[1][2+indexv[2]]=0.0;//s[2][3+indexv[3]]=0.0;
			s[1][jsf-1]=y[0][mpt-1]-anorm;//s[2][jsf]=y[1][mpt]-anorm; Equation (17.4.34).
		}
		else
		{// Interior point.
			s[0][indexv[0]-1] = -1.0;//s[1][indexv[1]] = -1.0; Equation (17.4.28).
			s[0][indexv[1]-1] = -0.5*h;//s[1][indexv[2]] = -0.5*h;
			s[0][indexv[2]-1]=0.0;//s[1][indexv[3]]=0.0;
			s[0][2+indexv[0]]=1.0;//s[1][3+indexv[1]]=1.0;
			s[0][2+indexv[1]] = -0.5*h;//s[1][3+indexv[2]] = -0.5*h;
			s[0][2+indexv[2]]=0.0;//s[1][3+indexv[3]]=0.0;
			temp1=x[k-1]+x[k-2];//temp1=x[k]+x[k-1];
			temp=h/(1.0-temp1*temp1*0.25);//temp=h/(1.0-temp1*temp1*0.25);
			temp2=0.5*(y[2][k-1]+y[2][k-2])-c2*0.25*temp1*temp1;//temp2=0.5*(y[3][k]+y[3][k-1])-c2*0.25*temp1*temp1;
			s[1][indexv[0]-1]=temp*temp2*0.5;//s[2][indexv[1]]=temp*temp2*0.5; Equation (17.4.29).
			s[1][indexv[1]-1] = -1.0-0.5*temp*(mm+1.0)*temp1;//s[2][indexv[2]] = -1.0-0.5*temp*(mm+1.0)*temp1;
			s[1][indexv[2]-1]=0.25*temp*(y[0][k-1]+y[0][k-2]);//s[2][indexv[3]]=0.25*temp*(y[1][k]+y[1][k-1]);
			s[1][2+indexv[0]]=s[1][indexv[0]-1];//s[2][3+indexv[1]]=s[2][indexv[1]];
			s[1][2+indexv[1]]=2.0+s[1][indexv[1]-1];//s[2][3+indexv[2]]=2.0+s[2][indexv[2]];
			s[1][2+indexv[2]]=s[1][indexv[2]-1];//s[2][3+indexv[3]]=s[2][indexv[3]];
			s[2][indexv[0]-1]=0.0;//s[3][indexv[1]]=0.0; Equation (17.4.30).
			s[2][indexv[1]-1]=0.0;//s[3][indexv[2]]=0.0;
			s[2][indexv[2]-1] = -1.0;//s[3][indexv[3]] = -1.0;
			s[2][2+indexv[0]]=0.0;//s[3][3+indexv[1]]=0.0;
			s[2][2+indexv[1]]=0.0;//s[3][3+indexv[2]]=0.0;
			s[2][2+indexv[2]]=1.0;//s[3][3+indexv[3]]=1.0;
			//s[1][jsf]=y[1][k]-y[1][k-1]-0.5*h*(y[2][k]+y[2][k-1]); (17.4.23).
			s[0][jsf-1]=y[0][k-1]-y[0][k-2]-0.5*h*(y[1][k-1]+y[1][k-2]);
			//s[2][jsf]=y[2][k]-y[2][k-1]-temp*((x[k]+x[k-1]) (17.4.24).
			s[1][jsf-1]=y[1][k-1]-y[1][k-2]-temp*((x[k-1]+x[k-2])
			//*0.5*(mm+1.0)*(y[2][k]+y[2][k-1])-temp2
			*0.5*(mm+1.0)*(y[1][k-1]+y[1][k-2])-temp2
			//*0.5*(y[1][k]+y[1][k-1]));
			*0.5*(y[0][k-1]+y[0][k-2]));
			//s[3][jsf]=y[3][k]-y[3][k-1]; Equation (17.4.27).
			s[2][jsf-1]=y[2][k-1]-y[2][k-2];
		}
//==================================
	}

	public double[] derivF(double x, double[] y)
    {
		double[] f=new double[y.length];

		//a single differential equation: y'=xy;
		//f[0]=x*y[0];

		//a single differential equation: y'=2*x;
		//f[0]=2.*x;

		//a single differential equation: y''+y=0;
		//it is a eq. system: y1'=y2; y2'=-y1;
		//f[0]=y[1];
		//f[1]=-y[0];
/*
		dydx[1] = -0.013*y[1]-1000.0*y[1]*y[3];
		dydx[2] = -2500.0*y[2]*y[3];
		dydx[3] = -0.013*y[1]-1000.0*y[1]*y[3]-2500.0*y[2]*y[3];
*/
		//f[0] = -0.013*y[0]-1000.0*y[0]*y[2];
		//f[1] = -2500.0*y[1]*y[2];
		//f[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];

		f[0]=y[1];//dydx[1]=y[2];
		//dydx[2]=(2.0*x*(m+1.0)*y[2]-(y[3]-c2*x*x)*y[1])/(1.0-x*x);
		f[1]=(2.0*x*(m+1.0)*y[1]-(y[2]-c22*x*x)*y[0])/(1.0-x*x);
		//dydx[3]=0.0;
		f[2]=0.0;

		return f;
	}

//y=f(x)
	public double F(double x)
	{
		//double result = 2*x*x+5*x+9-Math.sin(x);
		//double result = 4.*x*x;
		//double result = 4.-x*x;
		double result = -4.+x*x;
		//double result = 1./1.+x*x;
		//double result = 1./Math.sqrt(1.-x*x);
		//double result = Math.exp(-1.0*x);
		//double result = Math.sin(x);
		//double result = 3.*Math.sqrt(x)-x*x+4.*x-6.;
		//double result = 2.*(Math.sqrt(100.-x*x)+x*x/36.-9.);
		//double result = Math.PI/4.*(Math.log(x)*Math.log(x)+2.*Math.log(x)+1.0);
		return result;
	}
//y=f(x) and yy=f'(x)
	public double[] FD(double x)
	{
		double[] res=new double[2];

		//res[0]=4.-x*x;
		//res[1]=-2.0*x;
		res[0]=-4.+x*x;
		res[1]=2.0*x;

		return res;
	}
//y=f(x1,x2,...)
	public double MF(double[] x)
	{
		double res=0.0;

		//res=-8.0+x[0]*x[0]+x[1]*x[1];
		res=-4.0+x[0]*x[0];

		return res;
	}
//y'1=f'(x1,x2,...)/dx1,...the vector gradient df[1..n] evaluated at the input point x
	public double[] DMF(double[] x)
	{
		double[] res=new double[1];//0.0;

		//res=-8.0+x[0]*x[0]+x[1]*x[1];
		res[0]=2.0*x[0];

		return res;
	}

//3d function for integration:
	public double F3D(double x, double y, double z)
	{
		double res=0.0;
		res=x*y*z;
		return res;
	}

	public double yy1(double x)
	{
		return 0.0;
	}

	public double yy2(double x)
	{
		return 2.0;
	}

	public double z1(double x,double y)
	{
		return 0.0;
	}

	public double z2(double x,double y)
	{
		return 3.0;
	}
//==end 3d function integration
	//The user supplies a routine funcs(x,afunc,ma) that
    //returns the ma basis functions evaluated at x = x in the array afunc[1..ma].
	public double[] aF(double x, int ma)
	{
		double[] res=new double[ma];
		//y(x) = a1 + a2x + a3x^2 + � � � + aMx^M-1;
		//p[1]=1.0;
		//for (j=2;j<=np;j++) p[j]=p[j-1]*x;=>p
//POLINOME!!
/*
		for (int i=1;i<=ma;i++)
			res[i-1]=Math.pow(x,i-1);
*/
    //Fitting routine for an expansion with nl Legendre polynomials pl,
    //evaluated using the recurrence
//LEGENDRE
/*		int j=0;
		double twox=0.0;double f2=0.0;double f1=0.0;double d=0.0;
		res[0]=1.0;
		res[1]=x;
		if (ma > 2)
		{
			twox=2.0*x;
			f2=x;
			d=1.0;
			for (j=3;j<=ma;j++)
			{
				f1=d++;
				f2 += twox;
				res[j-1]=(f2*res[j-2]-f1*res[j-3])/d;
			}
		}
*/	//========legendre
//GAUSS,EXPONENTIAL y=a*exp(-bx)=>lny=loga - bx=c-bx;
		//res[0]=1.0;
		//res[1]=-x;
		//GAussian:
		//y(x)=a*exp(-bx^2); general form. miu=0
		//lny=lna-bx*x;
		res[0]=1.0;
		res[1]=-x*x;
		return res;
	}

	//public void fgauss(double x, double[] a, double y, double[] dyda, int na)
	public double fdf(double x, double[] a, double[] dyda, int na)
	//y(x; a) is the sum of na/3 Gaussians (15.5.16). The amplitude, center, and width of the
	//Gaussians are stored in consecutive locations of a: a[i] = Bk, a[i+1] = Ek, a[i+2] = Gk,
	//k = 1, ..., na/3. The dimensions of the arrays are a[1..na], dyda[1..na].
	{
		int i=0;
		double fac=0.0;double ex=0.0;double arg=0.0;
		double y=0.0;
		for (i=1;i<=na-1;i+=3)
		{
			arg=(x-a[i])/a[i+1];//arg=(x-a[i+1])/a[i+2];
			ex=Math.exp(-arg*arg);
			fac=a[i-1]*ex*2.0*arg;//fac=a[i]*ex*2.0*arg;
			y += a[i-1]*ex;//a[i]*ex;
			dyda[i-1]=ex;//dyda[i]=ex;
			dyda[i]=fac/a[i+1];//dyda[i+1]=fac/a[i+2];
			dyda[i+1]=fac*arg/a[i+1];//dyda[i+2]=fac*arg/a[i+2];
		}
		return y;
	}
//eq system nonlinear:
/*
ex: ax1+bx2-c=0 and
	dx1+ex2-f=0//=>functions
	//5x + 3y + 2z = 17 //x=1,y=2,z=3
	//2x - 6y + 3z = -1
	//3x + 2y - 4z = -5
*/
	public double[] vecfunc(int n, double[] x)
	{
		double[] res =new double[n];

		res[0]=5.0*x[0]+3.0*x[1]+2.0*x[2]-17.0;
		res[1]=2.0*x[0]-6.0*x[1]+3.0*x[2]+1.0;
		res[2]=3.0*x[0]+2.0*x[1]-4.0*x[2]+5.0;
		return res;
	}

	public void printSequence(String s)
	{
		System.out.println(s);
	}

	public double g(double t)
	{
		return 12.0;
	}

	public double ak(double t, double s)
	{
		return 21.0;
	}

	public double g(int k, double t)
	{
		return 12.0;
	}

	public double ak(int k, int l, double t, double s)
	{
		return 21.0;
	}

	public void kermom(double[] w, double y,int m)
	{

	/*
		double d,df,clog,x2,x3,x4,y2;
		if (y >= x) {
		d=y-x;
		df=2.0*sqrt(d)*d;
		w[1]=df/3.0;
		w[2]=df*(x/3.0+d/5.0);
		w[3]=df*((x/3.0 + 0.4*d)*x + d*d/7.0);
		w[4]=df*(((x/3.0 + 0.6*d)*x + 3.0*d*d/7.0)*x+d*d*d/9.0);
		} else {
		x3=(x2=x*x)*x;
		x4=x2*x2;
		y2=y*y;
		d=x-y;
		w[1]=d*((clog=log(d))-1.0);
		w[2] = -0.25*(3.0*x+y-2.0*clog*(x+y))*d;
		w[3]=(-11.0*x3+y*(6.0*x2+y*(3.0*x+2.0*y))
		+6.0*clog*(x3-y*y2))/18.0;
		w[4]=(-25.0*x4+y*(12.0*x3+y*(6.0*x2+y*
		(4.0*x+3.0*y)))+12.0*clog*(x4-(y2*y2)))/48.0;

	*/
	}

	public static void main (String[] args)
	{
		new TestNumerical();
	}
//===============================================================================
	public int NMAX =8;
	public void TFR(int m9,int q,double[][] v)//double v[NMAX][2])
	{
	// *************************************************************
	//  ===   FUNCTIA PENTRU TRANSFORMAREA FOURIER RAPIDA   ====/
	//   CALCULEAZA TRANSFORMAREA FOURIER DUPA ALGORITMUL COOLEY-
	//   TUKEY PENTRU  UN SINGUR SET DE DATE COMPLEXE, [V[][0],V[][1]]
	//       - TRANSFORMAREA DIRECTA SE OBTINE PENTRU  Q=-1.
	//       - TRANSFORMAREA INVERSA SE OBTINE PENTRU  Q=1.
	//       - VECTORUL V[] CONTINE PARTEA REALA.
	//       - VECTORUL W[] CONTINE PARTEA IMAGINARA.
	//       - m9 ESTE PUTEREA LUI 2 CARE DA LUNGIMEA LUI V[]
	// *************************************************************
	 	int i,j,k,l,j0,j1,m,i1=0;
	 	int j2,m1,m2,n;
	 	double u,c,c1,c2;
	 	double m0, s,v9,r1,r2,ee;
	 	boolean l1B=false;

	 	int[] e=new int[12];//int e[12];
	 	double PI=4*Math.atan(1.0);
	  	n=(int)Math.pow(2.0,m9);
	  	u=2*PI*q/n;
	  	for(i=1;i<=m9;i++)
	    {
	  		ee=Math.pow(2,i);//Math.pow(2,m9+i);
	  		e[i-1]=(int)Math.ceil(ee);//16,32,64
	  		//System.out.println(" "+e[i-1]);
	    }
	    for (i=1;i<=m9;i++)
	    {
	    	m0=Math.pow(2,i-1);
	   		m1=n/(int)m0;
	   		m2=m1/2;
	   		l=0;
	   		for(j=1;j<=m0;j++)
	     	{
	  			v9=u*l;
	  			s=Math.sin(v9);
	  			c=Math.cos(v9);
	  			j0=m1*(j-1);
	    		for(k=1;k<=m2;k++)
	  			{
	  			    j1=j0+k;
	  			    j2=j1+m2;
	  			    c1=v[j2-1][0]*c-v[j2-1][1]*s;
	  			    c2=v[j2-1][0]*s+v[j2-1][1]*c;
	  			    v[j2-1][0]=v[j1-1][0]-c1;
	  			    v[j2-1][1]=v[j1-1][1]-c2;
	  			    v[j1-1][0]=v[j1-1][0]+c1;
	  			    v[j1-1][1]=v[j1-1][1]+c2;
	  			}
	  		    for(m=2;m<=m9;m++)
	     		{
	     		  i1=m;
	       		  if (l<e[m-1]) break;// goto l3;
	       		  l=l-e[m-1];
	      		}
	  //l3:       l=l+e[i1-1];
	  			l=l+e[i1-1];
	     	}
	 	 }
	     l=0;//System.out.println(" trece");
	     for(i=1;i<=n;i++)
	     {
			 l1B=false;
	   		if (l>i) l1B=true;//goto l1;
	   		if(!l1B)
	   		{
	       		r1=v[i-1][0];
	       		r2=v[i-1][1];
	       		v[i-1][0]=v[l][0];
	       		v[i-1][1]=v[l][1];
	       		v[l][0]=r1;
	       		v[l][1]=r2;
			}
	//l1:    for(m=1;m<=m9;m++)
			for(m=1;m<=m9;m++)
	       {
	       		//i=m;
	       		if (l<e[m-1])  break;//goto l2;
	       		l=l-e[m-1];
	       }
	//l2:       l=l+e[i1-1];
		   l=l+e[i1-1];//System.out.println(" i "+i);
	    }//System.out.println(" trece");
	}
}

