package mcmc;

import genome.Genome;

import io.ReadDuplicatedGenome;
import io.ReadGenome;

import java.io.File;
import java.io.IOException;

import mutations.DCJ;

import phylogeny.Path;
import utils.Utils;

public class DCJCoupledPathMCMC {
	PathMCMC[] chain;
	
	public DCJCoupledPathMCMC(Path path, int numberOfChains, double scale){
		chain = new PathMCMC[numberOfChains];
		for(int i = 0; i < chain.length; i++){
			chain[i] = new PathMCMC(new Path(path), 1.0, (numberOfChains-1.0-i)*scale);
		}
		
	}
	
	public DCJCoupledPathMCMC(Genome g1, Genome g2, int numberOfChains, double scale){
		chain = new PathMCMC[numberOfChains];
		DCJ.chooseSorting = 1.0;
		DCJ.chooseNeutral = 0.0;
		DCJ.chooseBad = 0.0;
	//	DCJ.chooseSorting = 0.8;
	//	DCJ.chooseNeutral = 0.15;
	//	DCJ.chooseBad = 0.05;
		for(int i = 0; i < chain.length; i++){
			chain[i] = new PathMCMC(new Path(g1,g2), 1.0, (numberOfChains-1.0-i)*(numberOfChains-1.0-i)*scale);
			//System.out.println("chain "+i+" constructed");
		}
		DCJ.chooseSorting = 0.995;
		DCJ.chooseNeutral = 0.004;
		DCJ.chooseBad = 0.001;
	}
	
	void mcmcstep(){
		for(int i = 0; i < chain.length; i++){
			System.out.print("#"+i);
			chain[i].mcmcstep();
		}
		if(chain.length == 1){
			return;
		}
		int i = Utils.generator.nextInt(chain.length-1);
		double oldloglikelihood = chain[i].loglikelihood() + chain[i+1].loglikelihood();
		Path temp = chain[i].path;
		chain[i].path = chain[i+1].path;
		chain[i+1].path = temp;
		if(Utils.generator.nextDouble() < Math.exp(chain[i].loglikelihood() + chain[i+1].loglikelihood() - oldloglikelihood)){
			//accept, do nothing
			System.out.println("#accepting swap between chains "+i+" and "+(i+1));
		}
		else{
			temp = chain[i].path;
			chain[i].path = chain[i+1].path;
			chain[i+1].path = temp;
			System.out.println("#rejecting swap between chains "+i+" and "+(i+1));
	}
	}
	
	File doMCMC(int burnin, int samples, int period){
		for(int i = 0; i < burnin; i++){
			mcmcstep();
			//System.out.println("loglikelihood: "+path.logCountLikelihood(time, invBoltzmannTemperature, invDcjTemperature));
		}
		for(int i = 0; i < samples; i++){
			for(int j = 0; j < period; j++){
				mcmcstep();
			}
			for(int k = 0; k < chain.length; k++){
				System.out.println("loglikelihood of chain "+k+": "+chain[k].loglikelihood()+
						"\ttime: "+chain[k].time+
						"\taverage chromosome number: "+((double)chain[k].path.chromosomeSum()/(chain[k].path.length()+1))+
						"\tcircular: "+chain[k].path.circularChromosomeSum()+"\tpathlength: "+chain[k].path.length());
				System.out.println("\nPath "+k+"\n"+chain[k].path.print());
			}
		}
		
		return null;
	}

	/**
	 * only for testing reasons
	 * @param args unused
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException{
	//	Genome g1 = new Genome(new int[][] {{0,11}, {12,8}, {7,1}, {2,13}, {14,10}, {9,3}, {4,15}, {16,17}, {18,5}, {6,0}});
	//	Genome g2 = new Genome(new int[][] {{0,1}, {2,3}, {4,5}, {6,7}, {8,9}, {10,11}, {12,13}, {14,15}, {16,17}, {18,0}});
		//Genome g2 = new Genome(new int[][] {{0,1}, {2,3}, {4,5}, {6,7}, {8,9}, {10,11}, {12,13}, {14,15}, {16,17}, {18,0}});

		//Genome g1 = ReadGenome.read("data/euarchontoglires.genome");
		//Genome g2 = ReadGenome.read("data/Homo_sapiens.genome");
		Genome g1 = null;
		Genome g2 = null;
		
		if(args[0].equals("s")){
			//System.out.println("simple");
			g1 = ReadGenome.read(args[1]);
			g2 = ReadGenome.read(args[2]);
		}else if(args[0].equals("d")){
			//System.out.println("double");
			g1 = ReadDuplicatedGenome.read(args[1]);
			g2 = ReadDuplicatedGenome.read(args[2]);
		}else{
			System.out.println(args[0]);
			throw new Error("usage:\njava -Xmx512M -server -cp.:.. mcmc/DCJCoupledPathMCMC [s|d] file1 file2 burnin_steps samples repeats_between_two_samples number_of_parallel_chains s");
		}
		
		//Path path = new Path(g1,g2);

		//DCJCoupledPathMCMC mcmc = new DCJCoupledPathMCMC(g1,g2,10,0.1);
		//mcmc.doMCMC(0, 10, 10);
	
		DCJCoupledPathMCMC mcmc = new DCJCoupledPathMCMC(g1,g2,Integer.parseInt(args[6]),Double.parseDouble(args[7]));
		//System.out.println("mcmc constructed...");
		mcmc.doMCMC(Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]));
	}
	
}
