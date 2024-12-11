package mcmc;

import mutations.DCJ;
/**
 * This class is a container. It is used for returning with a DCJ operation and its proposal probability.
 * @author miklosi
 *
 */
public class DCJMetropolisContainer {
	
	public double proposalProbability;
	public DCJ mutation;
	
	public String print(){
		return "probability: "+proposalProbability+"\n"+(mutation != null ? mutation.print() : "no DCJ mutation in this container");
	}

}
