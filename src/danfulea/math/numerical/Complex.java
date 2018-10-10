package danfulea.math.numerical;

/**
 * Complex number class. 
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea,  09 OCT. 2006.
 */
public class Complex {
	public double r;//real
	public double i;//imaginar

	/**
	 * Constructor
	 * @param re real part
	 * @param im imaginary part
	 */
	public Complex(double re, double im)
	{
		this.r=re;
		this.i=im;
	}

	/**
	 * Adds two complex numbers
	 * @param a first number
	 * @param b second number
	 * @return the result
	 */
	public static Complex Cadd(Complex a, Complex b)
	{
		double cr=a.r+b.r;
		double ci=a.i+b.i;
		Complex c=new Complex(cr,ci);
		return c;
	}
	
	/**
	 * Subtracts two complex numbers
	 * @param a first number
	 * @param b second number
	 * @return the result
	 */
	public static Complex Csub(Complex a, Complex b)
	{
		double cr=a.r-b.r;
		double ci=a.i-b.i;
		Complex c=new Complex(cr,ci);
		return c;
	}
	
	/**
	 * Multiplies two complex numbers
	 * @param a first number
	 * @param b second number
	 * @return the result
	 */
	public static Complex Cmul(Complex a, Complex b)
	{
		double cr=a.r*b.r-a.i*b.i;
		double ci=a.i*b.r+a.r*b.i;
		Complex c=new Complex(cr,ci);
		return c;
	}

	/**
	 * Complex conjugate of a given complex number.
	 * @param z the number
	 * @return the result
	 */
	public static Complex Conjg(Complex z)
	{
		double cr=z.r;
		double ci = -z.i;
		Complex c=new Complex(cr,ci);
		return c;
	}

	/**
	 * Division of 2 complex numbers
	 * @param a first number
	 * @param b second number
	 * @return the result
	 */
	public static Complex Cdiv(Complex a, Complex b)
	{
		double r=0.0;double den=0.0;double cr=0.0;double ci=0.0;
		if (Math.abs(b.r) >= Math.abs(b.i))
		{
			r=b.i/b.r;
			den=b.r+r*b.i;
			cr=(a.r+r*a.i)/den;
			ci=(a.i-r*a.r)/den;
		}
		else
		{
			r=b.r/b.i;
			den=b.i+r*b.r;
			cr=(a.r*r+a.i)/den;
			ci=(a.i*r-a.r)/den;
		}
		Complex c=new Complex(cr,ci);
		return c;
	}

	/**
	 * The absolute value, the module, of a complex number.
	 * @param z the number
	 * @return the result
	 */
	public static double Cabs(Complex z)
	{
		double x=0.0;double y=0.0;double ans=0.0;double temp=0.0;
		x=Math.abs(z.r);
		y=Math.abs(z.i);
		if (x == 0.0)
			ans=y;
		else if (y == 0.0)
			ans=x;
		else if (x > y)
		{
			temp=y/x;
			ans=x*Math.sqrt(1.0+temp*temp);
		}
		else
		{
			temp=x/y;
			ans=y*Math.sqrt(1.0+temp*temp);
		}
		return ans;
	}

	/**
	 * Formal complex square root
	 * @param z the number
	 * @return the result
	 */
	public static Complex Csqrt(Complex z)
	{
		double x=0.0;double y=0.0;double w=0.0;double r=0.0;
		double cr=0.0;double ci=0.0;

		if ((z.r == 0.0) && (z.i == 0.0))
		{
			cr=0.0;
			ci=0.0;
			Complex c=new Complex(cr,ci);
			return c;
		}
		else
		{
			x=Math.abs(z.r);
			y=Math.abs(z.i);
			if (x >= y)
			{
				r=y/x;
				w=Math.sqrt(x)*Math.sqrt(0.5*(1.0+Math.sqrt(1.0+r*r)));
			}
			else
			{
				r=x/y;
				w=Math.sqrt(y)*Math.sqrt(0.5*(r+Math.sqrt(1.0+r*r)));
			}
			if (z.r >= 0.0)
			{
				cr=w;
				ci=z.i/(2.0*w);
			}
			else
			{
				ci=(z.i >= 0) ? w : -w;
				cr=z.i/(2.0*ci);
			}
			Complex c=new Complex(cr,ci);
			return c;
		}
	}

	/**
	 * Multiplication by a real number
	 * @param x the real number
	 * @param a the complex number
	 * @return the result
	 */
	public static Complex RCmul(double x, Complex a)
	{
		double cr=x*a.r;
		double ci=x*a.i;
		Complex c=new Complex(cr,ci);
		return c;
	}
}
