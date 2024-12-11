package mcmc;

import io.ReadGenome;

import java.awt.print.Printable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import mutations.DCJ;

import genome.Genome;
import phylogeny.Path;
import utils.Utils;

/**
 * This class is for testing our idea on two genomes. It can be also used if we would like to analyze only two genomes
 * @author miklosi
 *
 */
public class PathMCMC {

	double invBoltzmannTemperature;
	double invDcjTemperature;
	Path path;
	double time = 0.01;

	public PathMCMC(Path path, double invBoltzmannTemperature, double invDcjTemperature){
		this.path = path;
		this.invBoltzmannTemperature = invBoltzmannTemperature;
		this.invDcjTemperature = invDcjTemperature;
	}

	void mcmcstep(){
		if(Utils.generator.nextDouble() < 0.5){
			int length = path.length();
			double oldloglikelihood = path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature);
			/*
		int pos1;
		if(path.length() > 1){
			pos1 = Utils.generator.nextInt(length-1);
			//path.swap(pos1);
		}

		//here we re-sample a sub-path
		int pos2;
		pos1 = 1; pos2 = 0;
		while(pos1 > pos2){
			pos1 = (length == 0 ? 0 : Utils.generator.nextInt(length));
			pos2 = (length == 0 ? 0 : Utils.generator.nextInt(length));
		}
		pos2++;
			 */
			int[] interval = Utils.proposeInterval(length+1);
			int pos1 = interval[0];
			int pos2 = interval[1];
			PathMetropolisContainer pmc = path.proposeNewSubpath(pos1, pos2);
			//System.out.println("Starting with genome\n"+pmc.path.startGenome.print()+"\nthe newly proposed path: \n"+pmc.path.print()+"\n----------------");
		//	System.out.println("The old path was:\n"+path.print());
			ArrayList<DCJ> oldDcjs = path.dcjs;
			path.dcjs = new ArrayList<DCJ>();
			for(int i = 0; i < pos1; i++){
				path.dcjs.add(oldDcjs.get(i));
			}
			for(int i = 0; i < pmc.path.dcjs.size(); i++){
				path.dcjs.add(pmc.path.dcjs.get(i));
			}
			for(int i = pos2; i < oldDcjs.size(); i++){
				path.dcjs.add(oldDcjs.get(i));
			}
		//	System.out.println(path.startGenome.print());
			//System.out.println("Resampling path between "+pos1+" and "+pos2+"\n");
			path.check();
			double intervalProposalRatio = Utils.proposalProbability(path.length(), new int[] {pos1,pos1+pmc.path.length()})/Utils.proposalProbability(oldDcjs.size(), interval);
			double newloglikelihood = path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature);
			if(Utils.generator.nextDouble() < intervalProposalRatio*Math.exp(newloglikelihood-oldloglikelihood)*pmc.proposalProbability/pmc.path.backproposalProbability()){
				//do nothing
				System.out.println(" accepting a new path"/*+path.print()*/);
			}
			else{
				System.out.println(" rejecting a new path"/*+path.print()*/);
				path.dcjs = oldDcjs;
				
			}
		}else{
			
			double diff = (Utils.generator.nextDouble() - 0.5) * Path.timespan;
			if(time+diff > 0){
				double oldLoglikelihood = loglikelihood();
				time += diff;
				double newLoglikelihood = loglikelihood();
				if(Utils.generator.nextDouble() < Math.exp(newLoglikelihood - oldLoglikelihood)){
					//accept the change
					System.out.println(" Accepting the time change. Time is now "+time);
				}else{
					time -= diff;
					System.out.println(" Rejecting the time change");
				}
				
			}else{
				System.out.println(" rejection of time change due to proposing negative value");
			}
		}
	}

	File doMCMC(int burnin, int samples, int period){
		File output = new File("output.txt");//at the moment, we discard this file...
		for(int i = 0; i < burnin; i++){
			mcmcstep();
			//System.out.println("loglikelihood: "+path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature));
		}
		for(int i = 0; i < samples; i++){
			for(int j = 0; j < period; j++){
				mcmcstep();
			}
			System.out.println("loglikelihood: "+path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature)+
					"\ttime: "+time+
					"\taverage chromosome number: "+((double)path.chromosomeSum()/(path.length()+1))+"\tcircular: "+path.circularChromosomeSum()+"\tpathlength: "+path.length());
			System.out.println(path.print());
		}
		return output;
	}

	double loglikelihood(){
		return path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature);
	}

	public static void main(String[] args) throws IOException{
		Genome g2 = new Genome(new int[][] {{0,1},{2,3}, {4,5}, {6,7}, {8,9}, {10,11}, {12,13}, {14, 15}, {16,0}});
		Genome g1 = new Genome(new int[][] {{0,1}, {2,9}, {10,7}, {8,3}, {4,14}, {13,5}, {6,11}, {12,15}, {16,0}});
		//	Genome g1 = new Genome(new int[][] {{0,11}, {12,8}, {7,1}, {2,13}, {14,10}, {9,3}, {4,15}, {16,17}, {18,5}, {6,0}});
	//	Genome g1 = ReadGenome.read("../data/halved_ancestor.genome");
	//	System.out.println(g1.print());
		//Genome g2 = new Genome(new int[][] {{0,1}, {2,3}, {4,5}, {6,7}, {8,9}, {10,11}, {12,13}, {14,15}, {16,17}, {18,0}});
	//	Genome g2 = ReadGenome.read("../data/Skulyveri.genome");
		Path path = new Path(g1,g2);

		PathMCMC mcmc = new PathMCMC(path, 1.0, 0.00);
		mcmc.doMCMC(0, 1000, 1000);

	}

}
