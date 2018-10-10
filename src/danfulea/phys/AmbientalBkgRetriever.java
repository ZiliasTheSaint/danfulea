package danfulea.phys;

/**
 * A simple interface for ambiental BKG fetching data!
 * 
 * @author Dan Fulea, 15 APR. 2011
 * 
 */
public interface AmbientalBkgRetriever {
	
	double getAmbientalBKGPulsesAtChannel(double channel);
	double[] getAmbientalNetAreaAndUnc(double startChannel, double endChannel);
}
