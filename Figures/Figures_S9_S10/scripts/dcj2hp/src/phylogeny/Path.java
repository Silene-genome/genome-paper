package phylogeny;

import java.util.ArrayList;

import utils.Utils;

import mcmc.DCJMetropolisContainer;
import mcmc.PathMetropolisContainer;
import mutations.DCJ;

import genome.Genome;

public class Path {
	
	public ArrayList<DCJ> dcjs;
	public Genome startGenome;
	Genome endGenome;
	public static final double continueProbability = 0.05;
	public static final double timespan = 0.002;

	public Path(Genome g1, Genome g2){
		dcjs = new ArrayList<DCJ>();
		startGenome = new Genome(g1);
		endGenome = new Genome(g2);
		Genome current = new Genome(g1);
		while(!current.equals(g2)){
			DCJ dcj = DCJ.randomMutation(current,g2).mutation;
		//	current = current.mutate(dcj);
			current.mutate_this(dcj);
			dcjs.add(dcj);
		}
	}
	
	/**
	 * trivial constructor for packing genomes outside
	 */
	public Path(){
		dcjs = new ArrayList<DCJ>();
	}

	/**
	 * Copy constructor
	 * @param p this path is copied
	 */
	public Path(Path p){
		startGenome = new Genome(p.startGenome);
		endGenome = new Genome(p.endGenome);
		dcjs = new ArrayList<DCJ>();
		for(int i = 0; i < p.dcjs.size(); i++){
			dcjs.add(new DCJ(p.dcjs.get(i)));
		}
	}
	
	public String print(){
		String s = "";
		for(DCJ dcj : dcjs){
			s += dcj.print()+"\n";
		}
		return s;
	}
	
	/**
	 * This function swaps two dcjs
	 */
	public void swap(int i){
		if(dcjs.size() < 2) return;
		
		int[] m1 = dcjs.get(i).mutationInfo;
		int[] m2 = dcjs.get(i+1).mutationInfo;
		//System.out.println(dcjs.get(i).print());
		//System.out.println(dcjs.get(i+1).print());
		int count = 0;
		for(int x : m1){
			if(x == 0){
				count++;
			}
		}
		for(int x : m2){
			if(x == 0){
				count++;
			}
		}
		if(count > 1){
			return;
		}
		
		ArrayList<int[]> same = new ArrayList<int[]>();
		for(int k = 0; k < m1.length; k++){
			for(int j = 0; j < m2.length; j++){
				if(m1[k] != 0 && m1[k] == m2[j]){
					same.add(new int[]{k,j});
				}
			}
		}
		switch(same.size()){
			case 0:{
				int temp;
				for(int l = 0; l < m1.length; l++){
					temp = m1[l];
					m1[l] = m2[l];
					m2[l] = temp;
				}	
				break;
			}
			case 1:break;
			case 2:{
				if((m1[0] == 0 && m1[1] == 0)||
				   (m1[2] == 0 && m1[3] == 0)||
				   (m2[0] == 0 && m2[1] == 0)||
				   (m2[2] == 0 && m2[3] == 0)){
					break;
				}
				int temp;
				int[] s1 = same.get(0);
				int[] s2 = same.get(1);
				//System.out.println("s1[0]: "+s1[0]+" s1[1]: "+s1[1]);
				//System.out.println("s2[0]: "+s2[0]+" s2[1]: "+s2[1]);
				if(s1[0] == 0){
					temp = m1[0];
					m1[0] = m1[3];
					m1[3] = temp;
					temp = m1[1];
					m1[1] = m1[2];
					m1[2] = temp;
					temp = s1[1];
					s1[1] = s2[1];
					s2[1] = temp;
				}
				if(s1[1] > 1){
					temp = m2[0];
					m2[0] = m2[3];
					m2[3] = temp;
					temp = m2[1];
					m2[1] = m2[2];
					m2[2] = temp;
					s1[1] = 3 - s1[1];
					s2[1] = 3 - s2[1];
				}
				if(s1[1] != 0){
					temp = m2[0];
					m2[0] = m2[1];
					m2[1] = temp;
					temp = m2[2];
					m2[2] = m2[3];
					m2[3] = temp;
				}
				int a = m1[0], b = m1[1], c = m1[2], d = m1[3], e = m2[2], f = m2[3];
				if(Utils.generator.nextInt(2) == 0){
					//System.out.println("case 1");
					m1[0] = b; m1[1] = a; m1[2] = e; m1[3] = f;
					m2[0] = a; m2[1] = f; m2[2] = c; m2[3] = d;
				} else{
					//System.out.println("case 2");
					m1[0] = c; m1[1] = d; m1[2] = e; m1[3] = f;
					m2[0] = c; m2[1] = e; m2[2] = a; m2[3] = b;
				}
				break;
			}
			case 4:	break;
			default: break;
		}
//		System.out.println("After change:");
//		System.out.println(dcjs.get(i).print());
	//	System.out.println(dcjs.get(i+1).print());
	//	System.out.println("-------");
	
	}
	
	public PathMetropolisContainer proposeNewSubpath(int startstep, int endstep){
		PathMetropolisContainer pmc = new PathMetropolisContainer();
		pmc.path = new Path();
		Genome current = new Genome(startGenome);
		for(int i = 0; i < startstep; i++){
			current.mutate_this(dcjs.get(i));
		}
		pmc.path.startGenome = new Genome(current);
		for(int i = startstep; i < endstep; i++){
			current.mutate_this(dcjs.get(i));
		}
		pmc.path.endGenome = new Genome(current);
		/////////////
		current = new Genome(pmc.path.startGenome);
		//System.out.println("starting proposing a new path from genome:\n"+current.print());
		pmc.proposalProbability = 1.0;
		while(!current.equals(pmc.path.endGenome) || Utils.generator.nextDouble() < continueProbability){
			DCJMetropolisContainer dcjmc = DCJ.randomMutation(current,pmc.path.endGenome);
		//	current = current.mutate(dcj);
			//System.out.println(dcjmc.proposalProbability);
			pmc.proposalProbability *= dcjmc.proposalProbability;
			if(current.equals(pmc.path.endGenome)){
				pmc.proposalProbability *= continueProbability;
			}
			current.mutate_this(dcjmc.mutation);
			pmc.path.dcjs.add(dcjmc.mutation);
			//System.out.println("current genome: "+current.print());
		}
		pmc.proposalProbability *= 1.0 - continueProbability;
		////////////
		//System.out.println("The target genome: "+pmc.path.endGenome.print());
		return pmc;
	}
	
	/**
	 * This function calculates the proposal probability of this path
	 * @return
	 */
	public double backproposalProbability(){
		double pr = 1.0;
		Genome current = new Genome(startGenome);
		for(int i = 0; i < dcjs.size(); i++){
			DCJ dcj = dcjs.get(i);
			pr *= DCJ.backProposalProbability(current, endGenome, dcj);
			if(current.equals(endGenome)){
				pr *= continueProbability;
			}
			current.mutate_this(dcj);
		}
		pr *= (1.0 - continueProbability);
		
		return pr;
	}
	
	/**
	 * This function calculates the sum of the chromosomes along the path
	 * @return sum of the chromosomes
	 */
	public int chromosomeSum(){
		int c = 0;
		Genome current = new Genome(startGenome);
		int[] str = current.chromosomeStructure();
		c += str[0] + str[1];
		for(int i = 0; i < dcjs.size(); i++){
			current.mutate_this(dcjs.get(i));
			str = current.chromosomeStructure();
			c += str[0] + str[1];
		}
		return c;
	}
	
	/**
	 * This function calculates the variance of chromosomes along the path
	 * @return
	 */
	public double chromosomeVariance(){
		double v = (double)chromosomeSum()/(double)(length()+1);
		int c = 0;
		Genome current = new Genome(startGenome);
		int[] str = current.chromosomeStructure();
		c += (str[0] + str[1])*(str[0] + str[1]);
		for(int i = 0; i < dcjs.size(); i++){
			current.mutate_this(dcjs.get(i));
			str = current.chromosomeStructure();
			c += (str[0] + str[1])*(str[0] + str[1]);
		}
		return (double)c/(double)(length()+1) - v*v;
	}
	
	/**
	 * Tells the length of the Path, namely, the number of mutations applied on it.
	 * @return
	 */
	public int length(){
		return dcjs.size();
	}
	
	/**
	 * 
	 */
	public int numberOfExtremities(){
		int e = 0;
		for(int i = 0; i < startGenome.adjacencies.length; i++){
			if(startGenome.adjacencies[i][0] != 0 ){
				e++;
			}
			if(startGenome.adjacencies[i][0] != 0 ){
				e++;
			}
		}
		return e;
	}
	
	/**
	 * This function calculates the likelihood in a simple model where the probability depends on the
	 * number of steps and the sum of the chromosomes
	 * @param time evolutionary time
	 * @param invBoltzmannTemperature inverse (!!!) of the boltzmann temperature
	 * @param invDcjTemperature inverse (!!!) of the DCJ temperature setting it to 0 makes the likelihood being 
	 * independent from the type of DCJ (Hannenhalli-Pevzner type mutation or not)
	 * @return the logarithm of the likelihood in the simple count model
	 */
	public double logCountLikelihood(double time, double invBoltzmannTemperature, double invDcjTemperature){
		//int chr = chromosomeSum();
		//double var = chromosomeVariance();
		double[] exitRate = new double[length()];
		for(int i = 0; i < exitRate.length; i++){
			exitRate[i] = time;
		}
		double[] totalExitRates = new double[length()+1];
		Genome current = new Genome(startGenome);
		totalExitRates[0] = current.numberofHPs()*time;
		for(int i = 0; i < dcjs.size(); i++){
			DCJ dcj = dcjs.get(i);
			current.mutate_this(dcj);
			totalExitRates[i+1] = current.numberofHPs()*time;
		}
		double loglikelihood;// = Utils.logTrajectoryLikelihood(time,exitRate, totalExitRates);
		
		double var = circularChromosomeSum();
		
		int length = length();
		int e = numberOfExtremities();
		e = e/2;
		double x = -time*(e+1.0)*e/2.0 + length * Math.log(time);
		for(int i = 2; i <= length; i++){
			x -= Math.log((double)i);
		}
		//System.out.println("Comparing likelihoods: "+x+"\t"+loglikelihood);
		loglikelihood = x;
		loglikelihood *= invBoltzmannTemperature;
		//x -= (double)chr/(length+1) * invDcjTemperature;
		loglikelihood -= var * invDcjTemperature;
		return loglikelihood - time;
	}
	
	public double circularChromosomeSum() {
		int c = 0;
		Genome current = new Genome(startGenome);
		int[] str = current.chromosomeStructure();
		c += str[0];
		for(int i = 0; i < dcjs.size(); i++){
			current.mutate_this(dcjs.get(i));
			str = current.chromosomeStructure();
			c += str[0];
		}
		return c;
	}

	public void check(){
		Genome current = new Genome(startGenome);
		for(int i = 0; i < dcjs.size(); i++){
			current.mutate_this(dcjs.get(i));
		}
	}
	
	/**
	 * for testing purposes only
	 * @param args
	 */
	public static void main(String[] args){
		Genome g1 = new Genome(new int[][] {{4,0},{1,0}, {2,3}, {5,6}, {7,8}, {0,9}, {10,11}, {0, 12}, {13,14}, {15,16}, {17, 18}});
		Genome g2 = new Genome(new int[][] {{1,2}, {3,4}, {0,5}, {6,7}, {0,8}, {9,10}, {11,0}, {12,13}, {14,0}, {15,18}, {16,17}});
		Path path = new Path(g1,g2);
		System.out.println(path.print());
		path.swap(3);
		System.out.println(path.print());
		for(int i = 0; i < 1000; i++){
			PathMetropolisContainer pmc = path.proposeNewSubpath(0, path.dcjs.size());
			System.out.println("Proposal: "+pmc.proposalProbability);
			System.out.println("BackProposal: "+pmc.path.backproposalProbability());
			System.out.println("Sum of chromosomes: "+pmc.path.chromosomeSum());
			System.out.println("Loglikelihood: "+pmc.path.logCountLikelihood(1.0, 1.0, 0.0));
		}
		
	}
}
