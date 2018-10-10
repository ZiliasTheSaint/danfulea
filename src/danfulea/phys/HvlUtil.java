package danfulea.phys;

import danfulea.math.numerical.Interpolator;

/**
 * This class computes HVL from experimental data. Redundant methods since all you need 
 * is really a single function mmAl = f(Exposure) and compute mmAl for HVL (Exposure0/2). 
 * This function can be either a polynomial or better a spline function. Therefore, 
 * this class is deprecated. 
 * @author Dan Fulea, 2006
 *
 */
@Deprecated
public class HvlUtil {
	private static boolean zero = false;
	private static final double eps = 1e-10;

	public static boolean validate()
	{
		return zero;
	}
	
	//calculeaza valoarea unei functii pe baza coeficientilor polinomiali
	//e vorba de hvl function-->y=a0+...+anx^n-->y=a0/2(hvl1) ||y=a0/4(hvl2)
	//1 pentru hvl1 si orice altceva pentru hvl2!!!
	//zero==valoarea maxima a functiei tabelate--la zero mmAl!!
	////valoarea de zero poate fi diferita de coef[0]!!!->statistica!!
	public static double hvlFunc(double x,double[] coef,int hvl,double zero)
	{
	    double c = 0.0;

		for(int i =0; i<coef.length; i++)
		    c=c+coef[i]*Math.pow(x,i);
		if (hvl==1)
		    return c-zero/2;
		else
		    return c-zero/4;
	}

	//calculeaza derivata in punctul x!!a functiei data de coeficientii polinomiali
	//f(x)=a0+a1x+....+anx^n<->coef =a0,...,an
	public static double dFunc(double x,double[] coef)
	{
		double c = 0.0;
		for(int i=1; i<coef.length; i++)
		   c=c+i*coef[i]*Math.pow(x,i-1);
		return c;
	}

	//double[][] coef->pentru legatura cu sysEq!!!!a0,a1,....,an
	//converteste in sir matricea data de indexul n!!->pentru sisteme simple n=0!!!
	public static double[] convertMatrix(double[][] tl, int n)
	{
		double[] d = new double[tl.length];
		for(int i=0; i<tl.length; i++)
		   d[i]=tl[i][n];
		return d;
	}

//rezolvarea ecutatiilor prin metoda secantei
//a,b=capetele intervalului functiei,restul is parametrii functiei de evaluat(hvl)
	public static double hvlSecant(double a, double b,double[] coef,int hvl,double zeroVal)
	{
        zero=true;
        double x=a;
        double dx=1;
        double fx =1;
        //---particularizare
        double fa=hvlFunc(a,coef,hvl,zeroVal);//functia de prelucrat
        //------------------
        if (Math.abs(fa)<=eps)
           return x;
        x=b;
        double fb=hvlFunc(b,coef,hvl,zeroVal);
        if (Math.abs(fb)<=eps)
           return x;

		if (fa*fb<0)
        {
		   while(zero)
		   {

              x=(a*fb-b*fa)/(fb-fa);
	          fx=hvlFunc(x,coef,hvl,zeroVal);
	          if (fa*fx>0)
	          {
	            a=x;
	            fa=fx;
		      }
	          else
	          {
	            b=x;
	            fb=fx;
		      }
	          dx=b-a;
	          if (x!=0)
                dx=dx/x;

              if(Math.abs(dx)<=eps || Math.abs(fx)<=eps)
                 break;
		   }
	    }
		else
	    {
	 	  zero=false;
	    }

		if(zero)
		{
		  return x;
	    }
	    else
	      return 0.0;

	}
//pentru evaluarea maximului pe baza derivatei->folosit la MAB la mKcl vs effic
    public static double maxSecant(double a, double b,double[] coef)
	{
        zero=true;
        double x=a;
        double dx=1;
        double fx =1;
        //---particularizare
        double fa=dFunc(a,coef);//functia de prelucrat
        //------------------
        if (Math.abs(fa)<=eps)
           return x;
        x=b;
        double fb=dFunc(b,coef);
        if (Math.abs(fb)<=eps)
           return x;

		if (fa*fb<0)
        {
		   while(zero)
		   {

              x=(a*fb-b*fa)/(fb-fa);
	          fx=dFunc(x,coef);
	          if (fa*fx>0)
	          {
	            a=x;
	            fa=fx;
		      }
	          else
	          {
	            b=x;
	            fb=fx;
		      }
	          dx=b-a;
	          if (x!=0)
                dx=dx/x;

              if(Math.abs(dx)<=eps || Math.abs(fx)<=eps)
                 break;
		   }
	    }
		else
	    {
	 	  zero=false;
	    }

		if(zero)
		{
		  return x;
	    }
	    else
	      return 0.0;

	}


//rezolvarea ecutatiilor prin metoda secantei
//x=solutia grosiera=>se incearca o solutie rafinata
	public static double hvlIter(double x, double[] coef,int hvl,double zeroVal)
	{
         int imax = 100;
         zero =true;
         int i=0;
         double test,dx = 1.0;

         while(zero)
         {
			 i++;
			 dx=hvlFunc(x,coef,hvl,zeroVal);//functia de prelucrat
			 test=1/dx;//de control al neconvergentei----+/-infinit
			 if (test==0.0)
		     {
				 //cand sirul e neconvergent
				 zero=false;
				 break;
			 }

			 x=x-dx;
			 if(x!=0)
			   dx=dx/x;

			 if(Math.abs(dx)<=eps || i>=imax)
                 break;
		 }

		 if (i>imax)
		    zero=false;

		  if (zero)
		     return x;
		  else
		     return 0.0;
	}

    public static double maxIter(double x, double[] coef)
	{
         int imax = 100;
         zero =true;
         int i=0;
         double test,dx = 1.0;

         while(zero)
         {
			 i++;
			 dx=dFunc(x,coef);//functia de prelucrat
			 test=1/dx;//de control al neconvergentei----+/-infinit
			 if (test==0.0)
		     {
				 //cand sirul e neconvergent
				 zero=false;
				 break;
			 }

			 x=x-dx;
			 if(x!=0)
			   dx=dx/x;

			 if(Math.abs(dx)<=eps || i>=imax)
                 break;
		 }

		 if (i>imax)
		    zero=false;

		  if (zero)
		     return x;
		  else
		     return 0.0;
	}
    
    public static double[] invert(double[] coef)
	{
        double[] c = new double[coef.length];
        for(int i=0; i<coef.length; i++)
        {
			c[i]=coef[coef.length-1-i];
		}
		return c;
	}

	public static double[] updateToHvlBirge(double[] coef, int hvl, double zeroVal)
	{
		double[] c = new double[coef.length];
		if (hvl==1)
		  c[0]=coef[0]-zeroVal/2;
		else
		  c[0]=coef[0]-zeroVal/4;
		for (int i=1; i<coef.length; i++)
		{
		     c[i]=coef[i];
		}
		return c;
	}

    public static double[] updateToMaxBirge(double[] coef)
	{
		//derivata pentru acesti coeficienti
		double[] c = new double[coef.length-1];

		for (int i=1; i<coef.length; i++)
		{
		     c[i-1]=coef[i]*i;
		}
		return c;
	}

	public static double revertCoordonatesForHvl(double[] x,double[] y,int n,int hvl,double zeroVal,double[][] a)
	{
         Interpolator.polynomial(y,x,n,a);
         double[] d = convertMatrix(a,0);

         double result = 0.0;
         for (int i=0; i<d.length; i++)
            if (hvl==1)
               result = result + d[i]*Math.pow(zeroVal/2,i);
            else
               result = result + d[i]*Math.pow(zeroVal/4,i);
         return result;
	}
	//=============SOLUTIONS================
	
}
