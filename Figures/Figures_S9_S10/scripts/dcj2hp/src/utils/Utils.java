package utils;

import java.util.Random;

public class Utils {

	public static Random generator = new Random(841);
	//public static Random generator = new Random(1012);
	
	public static double[] lengthProbabilities;
	public static double[] cumulativeLP;
	
	static{
		lengthProbabilities = new double [10];
		cumulativeLP = new double[11];
		setLengthProbabilities();
	}
	
	static void setLengthProbabilities(){
		lengthProbabilities[0] = 0.1;
		lengthProbabilities[1] = 0;
		cumulativeLP[0] = 0.1;
		cumulativeLP[1] = 0.1;
		for(int i = 2; i < lengthProbabilities.length; i++){
			lengthProbabilities[i] = 1.0/(double)i;
			cumulativeLP[i] = cumulativeLP[i-1]+lengthProbabilities[i-1];
		}
		cumulativeLP[cumulativeLP.length-1] = cumulativeLP[cumulativeLP.length-2] + lengthProbabilities[cumulativeLP.length-2];
	}
	
	public static int[] proposeInterval(int totalLength){
		while(totalLength >= lengthProbabilities.length){
			lengthProbabilities = new double[lengthProbabilities.length*2];
			cumulativeLP = new double[lengthProbabilities.length+1];
			setLengthProbabilities();
			//System.out.println("Now totellength: "+totalLength+" lengthProbabilities length: "+lengthProbabilities.length);
		}
		double x = generator.nextDouble()*cumulativeLP[totalLength];
		int a = 0;
		int b = totalLength+1;
		while(b-a != 1){
			if(cumulativeLP[(a+b)/2] < x){
				a = (a+b)/2;
			}
			else{
				b = (a+b)/2;
			}
			//System.out.println("cumulative["+a+"]: "+cumulativeLP[a]+", cumulative["+b+"]: "+cumulativeLP[b]+" x: "+x);
		}
		int start = totalLength - a == 0 ? 0 : generator.nextInt(totalLength-a);
		return new int[] {start, start+a};
	}
	
	public static double proposalProbability(int totalLength, int[] interval){
		while(totalLength >= lengthProbabilities.length){
			lengthProbabilities = new double[lengthProbabilities.length*2];
			cumulativeLP = new double[lengthProbabilities.length+1];
			setLengthProbabilities();
		}
		int span = interval[1]-interval[0];
		return lengthProbabilities[span]/((totalLength - span == 0 ? 1 : totalLength - span)*cumulativeLP[totalLength]);
	}
	
	public static double logTrajectoryLikelihood(double timeunit, double[] exitRate, double[]totalExitRates){
		double value = 0.0;
		for (int i = 0; i < exitRate.length; i++) {
			value += Math.log(exitRate[i]);
		}

		//sorting of total exit rates
		boolean sorted = false;
		while (!sorted) {
			sorted = true;
			for (int j = 1; j < totalExitRates.length; j++) {
				if(Math.abs(totalExitRates[j] - totalExitRates[j-1]) < 0.1){
					totalExitRates[j] += timeunit/totalExitRates.length;
					sorted = false;
				}
				if (totalExitRates[j - 1] > totalExitRates[j]) {
					double temp = totalExitRates[j - 1];
					totalExitRates[j - 1] = totalExitRates[j];
					totalExitRates[j] = temp;
					
					sorted = false;
				}
			}	
		}
		/*
		 System.out.print("\nTotal exit rates after sorting:\t");
		 for(int i = 0; i < totalExitRates.length; i++){
		 System.out.print(totalExitRates[i]+", ");
		 }
		 System.out.println();
		 */
		 double[][] table = new double[totalExitRates.length][totalExitRates.length];
		 for(int i = 0; i <totalExitRates.length; i++){
			 table[i][i] = -totalExitRates[i];
		 }
		 for(int k =1; k < totalExitRates.length; k++){
			 //check for numerical errors. If exist, take exponential interpolation
			 for(int i = 0; i < totalExitRates.length - k; i++){
				 if(table[i][i+k-1] < table[i+1][i+k]){
					 int j = 1;
					 while(i+k+j < totalExitRates.length && table[i][i+k-1] < table[i+j+1][i+k+j]){
						 j++;
					 }
					 if(i+j+k == totalExitRates.length){
					//	 System.out.println("correcting at the end, ("+i+","+(i+k-1)+")");
					//	 System.out.println("table["+i+"]["+(i+k-1)+"]: "+table[i][i+k-1]);
						 for(j = 0; i+j+k < totalExitRates.length; j++){
							 table[i+j+1][i+k+j] = table[i][i+k-1] - (totalExitRates[i+k+j]+totalExitRates[i+j+1] - 
									 								totalExitRates[i]-totalExitRates[i+k-1])/2;
						//	 System.out.println("table["+(i+j+1)+"]["+(i+k+j)+"]: "+table[i+j+1][i+k+j]);
							}
					 }else{
						// System.out.println("correcting in the middle, ("+i+","+(i+k-1)+")");
						 //System.out.println("table["+i+"]["+(i+k-1)+"]: "+table[i][i+k-1]);
						 for(int l = 0; l < j; l++){
							 table[i+l+1][i+k+l] = (totalExitRates[i+l+1]+totalExitRates[i+k+l] - totalExitRates[i] - totalExitRates[i+k-1])*
							 (table[i][i+k-1] - table[i+j+1][i+k+j])/(totalExitRates[i]+totalExitRates[i+k-1]-totalExitRates[i+j+1]-totalExitRates[i+k+j])+
							 table[i][i+k-1];
							// System.out.println("table["+(i+l+1)+"]["+(i+k+l)+"]: "+table[i+l+1][i+k+l]);
						 }
					 }
				 }
			 }
			for(int i = 0; i < totalExitRates.length - k; i++){
				if(table[i][i+k-1] < table[i+1][i+k]){
					System.out.println("problem: ("+i+","+(i+k)+") "+table[i][i+k-1]+", "+table[i+1][i+k]);
					table[i][i+k] = Double.NEGATIVE_INFINITY;
				}else{
					table[i][i+k] = logsubstract(table[i][i+k-1], table[i+1][i+k]);
					table[i][i+k] -= Math.log(totalExitRates[i+k]-totalExitRates[i]);
				}
			}
		 }
		
		return value + table[0][totalExitRates.length-1];
	}
	
	
	private static double logadd(double a, double b){
		if(a < b){
			double temp = a;
			a = b;
			b = temp;
		}
		return a + Math.log1p(Math.exp(b-a));
	}
	
	private static double logsubstract(double a, double b){
		if(a < b){
			return b + Math.log(Math.exp(a-b) - 1.0);
		}
		return a + Math.log1p(-Math.exp(b-a));
	}
	
	public static void main(String[] args){
		//testing length proposal
		System.out.println("\n"+logTrajectoryLikelihood(0.01, new double[] {1}, new double[] {12.23551980215192, 13.272695964382525, 13.82024074038244, 14.346726101920817, 14.369374638074657, 14.38463949169004, 14.40463949169004, 14.425698906151574, 14.546790539305402, 14.86268175622843, 14.925859999613035, 15.020627364689943, 15.11539472976685}));
		System.out.println("\n"+logTrajectoryLikelihood(0.01, new double[] {1}, new double[] {0,
				0.1,
				0.2894,
				0.3894,
				0.5785,
				0.7773,
				0.9611,
				1.0629,
				1.1641,
				1.2668,
				1.4605,
				1.5668,
				1.7398,
				1.8668,
				2.0471}));
		/*
		int totalLength = 10;
		int[] interval;
		do{
			interval = proposeInterval(totalLength);
			System.out.println("Propsed interval: ["+interval[0]+","+interval[1]+"]");
			System.out.println("Its probability: "+proposalProbability(totalLength, interval));
		}while(interval[1] != 10);
		*/
	}
}
