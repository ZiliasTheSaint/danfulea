package danfulea.math.numerical;

//import jdf.math.numerical.Function;
//import jdf.math.numerical.Integrator;

public class TestIntegral implements Function{
	public double limInf=0.;
	public double limSup=1.;
	//------------------------
	public double x0=0.;
	public double xf=0.;
	public int npoints=10;
	public double step=0.5;
	public double[] y0;
	public int order=1;
/*
 * import Jad.mathemathics.Function;
import Jad.mathemathics.DifferentialFunction;
import Jad.mathemathics.AdvancedIntegralSolver;
import Jad.mathemathics.AdvancedDifferentialEquation;
import Jad.mathemathics.AdvancedEquations;

public class TestIntegral implements Function,DifferentialFunction
{
	

	
	
    //diferential conditions:
    //x=initial x0 at which we have y=y0 conditions (y0,y0',...,y0'(order-1))
    //h=step is distance from x0 to evaluate function and its derivates
    //n=order is the differential equation order
	private void setDifferentialConditions(double x,double h, double[] y, int n)
	{
		x0=x;
		step=h;
		y0=new double[y.length];
		for (int i=0; i<y.length; i++)
		{
			y0[i]=y[i];
		}
		order=n;
	}
    //diferential conditions:
    //x=initial x0 at which we have y=y0 conditions (y0,y0',...,y0'(order-1))
    //xx=xf is distance from x0 to evaluate function and its derivates
    //nstep=npoints is number of the substeps in which the interval [x0,xf] is divided
    //n=order is the differential equation order
    //@For each steps the RK4 is computed and new condition are updated!!
	private void setDifferentialConditionsForInterval(double x,double xx,int nstep, double[] y, int n)
	{
		x0=x;
		xf=xx;
		npoints=nstep;
		y0=new double[y.length];
		for (int i=0; i<y.length; i++)
		{
			y0[i]=y[i];
		}
		order=n;
	}

    //solve the differential equations
	private void evaluateDifferentialEquation()
	{
		AdvancedDifferentialEquation adf=new AdvancedDifferentialEquation(this);
		//adf.RK4(x0,step,y0,order);
		adf.allRK4(x0,xf,npoints,y0,order);//more accurate!!
	}

    //setup the interval in which the solution is supposed to be found!
	private void setLimits(double inf, double sup)
	{
		limInf=inf;
		limSup=sup;
	}
    //solve the equations
	private void evaluateEquation()
	{
		AdvancedEquations ae=new AdvancedEquations(this);
		double x0=ae.Secant(limInf,limSup);
		ae.Iter(x0);

		double[] coef={-1.0,0.0,4.0};//-x*x +0*x+4=0
		double[] sol=AdvancedEquations.birge(coef);

		for(int i=0;i<sol.length;i++)
		{
			System.out.println("birge sol= "+sol[i]);
		}
	}

//====================================================================================
	public double[] func(double x, double[] y)
    {
		double[] f=new double[y.length];
		//a single differential equation: y'=xy;
		//f[0]=x*y[0];

		//a single differential equation: y'=2*x;
		//f[0]=2.*x;

		//a single differential equation: y''+y=0;
		//it is a eq. system: y1'=y2; y2'=-y1;
		f[0]=y[1];
		f[1]=-y[0];

		return f;
	}

	public double F(double x)
	{
		//double result = 2*x*x+5*x+9-Math.sin(x);
		//double result = 4.*x*x;
		double result = 4.-x*x;
		//double result = 1./1.+x*x;
		//double result = 1./Math.sqrt(1.-x*x);
		//double result = Math.exp(-1.0*x);
		//double result = Math.sin(x);
		//double result = 3.*Math.sqrt(x)-x*x+4.*x-6.;
		//double result = 2.*(Math.sqrt(100.-x*x)+x*x/36.-9.);
		//double result = Math.PI/4.*(Math.log(x)*Math.log(x)+2.*Math.log(x)+1.0);
		return result;
	}

	public void printSequence(String s)
	{
		System.out.println(s);
	}
//====================================================================================
///    public static void main(String[] args)/
	//{
	//	new TestIntegral();
	//}
//
}
 */
	public TestIntegral()
	{
//Integrale=================
		//setIntegrationLimits(0.,1.);
		//setIntegrationLimits(0.,3.);//=>36, ok!!
		//setIntegrationLimits(0.,1./2.);//=>0.458, ok!!
		//setIntegrationLimits(0.,1./2.);//=>0.542 vs real arctg=0.464!!
		//setIntegrationLimits(0.,0.99);//=>1.429 vs real PI/2!!
		//setIntegrationLimits(0.,1.e+2);//=>1.0 ok!!
		//setIntegrationLimits(0.,Math.PI);//=>2.0 ok!!
		//setIntegrationLimits(1.,4.);//=>5.0 ok!!
		//setIntegrationLimits(0.,6.);//=>8.350 vs 8.360! ok!!
		setIntegrationLimits(1.,10.);//=>48.710 vs 48.73! ok!!
	    evaluateIntegral();
//======================================================================================
//DiffEq=================

	 //@@@@@@@@@@@@@@@Note if order=1 <=> init.length must be equal with order, e.g. 1

	    //double[] init={1.0};//at x0 we have y0=y(x0)=1.0
	    //setDifferentialConditions(0.,4,init,1);
	    //setDifferentialConditionsForInterval(0.,4.,10, init, 1);

	   // double[] init={1.0,1.0};//at x0 we have y0=y(x0)=1.0; y0'=1.0
	    //setDifferentialConditions(0.,4,init,2);
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@increse npoints until we cannot gain more precision!!
	  //  setDifferentialConditionsForInterval(0.,Math.PI/4.,100000, init, 2);

	  //  evaluateDifferentialEquation();
	//System.out.println("e1*1/2= "+Math.sqrt(Math.E));
	//System.out.println("e2*2/2= "+Math.E*Math.E);
	//System.out.println("e4*4/2= "+Math.pow(Math.E,4*4./2.));
	//System.out.println("cospi/4+sinpi/4 "+Math.cos(Math.PI/4.)+Math.sin(Math.PI/4.));
//======================================================================================
//EQ evaluation==============================
	//	setLimits(0.,5.);//=>48.710 vs 48.73! ok!!
	//    evaluateEquation();
//======================================================================================

	}
	
	
	 public static void main(String[] args)
		{
		new TestIntegral();
		}
	
	
	//integration limits!
		private void setIntegrationLimits(double inf, double sup)
		{
			limInf=inf;
			limSup=sup;
		}
	    //solve the integral
		private void evaluateIntegral()
		{
			Integrator ais=new Integrator(this);
			double intg = ais.QD(limInf,limSup);
			System.out.println("DCADRE= "+intg);
			ais.qTrapez(limInf,limSup);
			ais.qSimpson(limInf,limSup);
			ais.qGauss5(limInf,limSup);
		}
	
	
	
	
	public void printSequence(String s) {// for printing
		System.out.println(s);
	}

	// -------------------------------------------------------
	
	public double F(double x) {// y=f(x)
		double result = 4.+2*x+x*x;

		return result;
		// return 0.0;
	}

	// --------------------------------------------------------
	public double[] FD(double x) {// 0=->y=f(x) and 1-> dy=f'(x)
		return null;
	}

	public double MF(double[] x) {// y=f(x1,x2,...)
		return 0.0;
	}

	public double[] DMF(double[] x) {// y'1=df(x1,x2,...)/dx1;...the vector
										// gradient
		// df[1..n] evaluated at the input point x
		return null;
	}

	// ==============3D func====================
	public double F3D(double x, double y, double z) {
		return 0.0;
	}

	public double yy1(double x) {
		return 0.0;
	}

	public double yy2(double x) {
		return 0.0;
	}

	public double z1(double x, double y) {
		return 0.0;
	}

	public double z2(double x, double y) {
		return 0.0;
	}

	// end 3d====================================
	// non linear equation systems: root finding
	public double[] vecfunc(int n, double[] x) {
		return null;
	}

	// ==========================================
	public double[] aF(double x, int ma) {// poli fit
		double[] y = new double[ma];

		for (int i = 0; i < ma; i++) {
			y[i] = Math.pow(x, i);
		}
		// y[0]=1.0;
		// y[1]=x;
		// y[2]=x*x;
		// y[3]=x*x*x;
		return y;
	}

	// fdf used in mrqmin-mrqcof routine!!!!!!!!!!!!!!!!!!
	public double fdf(double x, double[] a, double[] dyda, int na) {
		// gauss fit etc. nonliniar
		// double y = 0.0;

		// double expt = Math.exp(-(x - a[1]) * (x - a[1]) / (2.0 * a[2] *
		// a[2]));
		// y = a[0] * expt;
		// dyda[0] = expt;
		// dyda[1] = y * ((x - a[1]) / (a[2] * a[2]));
		// dyda[2] = y * ((x - a[1]) * (x - a[1]) / (a[2] * a[2] * a[2]));
		return 0.0;// y;
	}

	// ===========================================================
	public double[] derivF(double x, double[] y) {
		return null;
	}

	// ============2point
	public void load(double x1, double[] v, double[] y) {
		return;
	}

	public void load1(double x1, double[] v, double[] y) {
		return;
	}

	public void load2(double x2, double[] v, int nn2, double[] y) {
		return;
	}

	public void score(double x2, double[] y, double[] f) {
		return;
	}

	public void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
			int indexv[], int ne, double[][] s, double[][] y) {
		return;
	}

	// ================================================
	public double g(double t) {// g(t)=FREDHOLM
		return 0.0;
	}

	public double ak(double t, double s) {// KERNEL
		return 0.0;
	}

	public double g(int k, double t) {// voltera
		return 0.0;
	}

	public double ak(int k, int l, double t, double s) {// voltera
		return 0.0;
	}

	public void kermom(double[] w, double y, int m) {
		return;
	}

	// =========END IMPLEMENTATION OF FUNCTION INTERFACE===================
}
