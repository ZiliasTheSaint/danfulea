package danfulea.phys.egs;

/**
 * The interface for implementing the sensitive routines for Monte Carlo simulation of radiation transport based on 
 * EGS4/EGSnrc simulation toolkit.
 * @author Dan Fulea, 07 NOV. 2005
 *
 */
public interface EgsQuestion {
	void printSequence(String s);
    void AUSGAB(int IARG);
    void HOWNEAR();
    void HOWFAR();
}
