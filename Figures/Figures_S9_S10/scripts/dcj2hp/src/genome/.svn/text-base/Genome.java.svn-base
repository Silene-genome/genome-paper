package genome;

import java.util.ArrayList;

import representation.*;
import utils.Utils;

import mcmc.DCJMetropolisContainer;
import mutations.DCJ;

public class Genome {
	
	public int[][] adjacencies;
	private boolean[] visited;
	/**
	 * maps from the gene codes to the adjacency indexes;
	 */
	private int[] map;
	
	/**
	 * constructor for deep copy
	 */
	public Genome(Genome g){
		adjacencies = new int[g.adjacencies.length][2];
		for(int i = 0; i < adjacencies.length; i++){
			adjacencies[i][0] = g.adjacencies[i][0];
			adjacencies[i][1] = g.adjacencies[i][1];
		}
		visited = new boolean[adjacencies.length];
		map = new int[g.map.length];
		for(int i = 0; i < map.length; i++){
			map[i] = g.map[i];
		}
		
	}
	
	public Genome(int[][] adjacencies){
		this.adjacencies = adjacencies;
		visited = new boolean[adjacencies.length];
		mapIt();
	}
	
	private void mapIt(){
		int geneends = 0;
		for(int i = 0; i < adjacencies.length; i++){
			if(adjacencies[i][0] != 0){
				geneends++;
			}
			if(adjacencies[i][1] != 0){
				geneends++;
			}
		}
	//	System.out.println("geneends: "+geneends);
		map = new int[geneends+1];
		for(int i = 0; i < adjacencies.length; i++){
			if(adjacencies[i][0] != 0){
				map[adjacencies[i][0]] = i;
			}
			if(adjacencies[i][1] != 0){
				map[adjacencies[i][1]] = i;
			}
		}
	}
	
	public boolean equals(Genome g){
		if(g.adjacencies.length != adjacencies.length){
			return false;
		}
		boolean answer = true;
		for(int i = 0; answer && i < adjacencies.length; i++){
			boolean found = false;
			for(int j = 0; j < g.adjacencies.length && !found; j++){
				if((adjacencies[i][0] == g.adjacencies[j][0] && adjacencies[i][1] == g.adjacencies[j][1]) ||
				   (adjacencies[i][1] == g.adjacencies[j][0] && adjacencies[i][0] == g.adjacencies[j][1])){
					found = true;
				}
			}
			answer = found;
		}
		return answer;
	}
		
	/**
	 * This function generates a new genome which is obtained by applying DCJ m to this genome.
	 * @param m DCJ that is pplied to the genome
	 * @return the generated genome
	 */
	public Genome mutate(DCJ m){
		Genome g = new Genome(this);
		
		if(m.mutationInfo[2] == 0 && m.mutationInfo[3] == 0){
			// it is a split or linearization
			int x = -1;
			for(int i = 0; x == -1 && i < g.adjacencies.length; i++){
				if((adjacencies[i][0] == m.mutationInfo[0] && adjacencies[i][1] == m.mutationInfo[1]) ||
				   (adjacencies[i][1] == m.mutationInfo[0] && adjacencies[i][0] == m.mutationInfo[1])){
					x = i;
				}
			}
			if(x == -1){
				throw new Error("Cannot apply this mutation:"+m.print()+" to this genome. x: "+x);
			}
			int[][] a = new int[g.adjacencies.length+1][2];
			for(int k = 0; k < x; k++){
				a[k][0] = g.adjacencies[k][0];
				a[k][1] = g.adjacencies[k][1];
			}
			for(int k = x+1; k <g.adjacencies.length; k++){
				a[k-1][0] = g.adjacencies[k][0];
				a[k-1][1] = g.adjacencies[k][1];
			}
			a[a.length-2][0] = m.mutationInfo[0];
			a[a.length-2][1] = 0;
			a[a.length-1][0] = m.mutationInfo[1];
			a[a.length-1][1] = 0;
			g.adjacencies = a;
			g.visited = new boolean[g.adjacencies.length];
			g.mapIt();
			return g;
		}
				
		int x = -1;
		int y = -1;
		for(int i = 0; (x == -1 || y == -1) && i < g.adjacencies.length; i++){
			if((adjacencies[i][0] == m.mutationInfo[0] && adjacencies[i][1] == m.mutationInfo[1]) ||
			   (adjacencies[i][1] == m.mutationInfo[0] && adjacencies[i][0] == m.mutationInfo[1])){
				x = i;
			}
			else if((adjacencies[i][0] == m.mutationInfo[2] && adjacencies[i][1] == m.mutationInfo[3]) ||
				    (adjacencies[i][1] == m.mutationInfo[2] && adjacencies[i][0] == m.mutationInfo[3])){
				y = i;
			}
		}
		if(x == -1 || y == -1){
			throw new Error("Cannot apply this mutation:"+m.print()+" to this genome. x: "+x+", y: "+y);
		}
		int i = (adjacencies[x][0] == m.mutationInfo[0] ? 1 : 0);
		int j = (adjacencies[y][0] == m.mutationInfo[3] ? 1 : 0);
		g.map[adjacencies[x][i]] = y;
		g.map[adjacencies[y][j]] = x;
		int temp = adjacencies[x][i];
		g.adjacencies[x][i] = g.adjacencies[y][j];
		g.adjacencies[y][j] = temp;
		if((g.adjacencies[x][0] == 0 && g.adjacencies[x][1] == 0) ||
		    g.adjacencies[y][0] == 0 && g.adjacencies[y][1] == 0){
			int z = (g.adjacencies[x][0] == 0 && g.adjacencies[x][1] == 0 ? x : y);
			int[][] a = new int[g.adjacencies.length-1][2];
			for(int k = 0; k < z; k++){
				a[k][0] = g.adjacencies[k][0];
				a[k][1] = g.adjacencies[k][1];
			}
			for(int k = z+1; k <g.adjacencies.length; k++){
				a[k-1][0] = g.adjacencies[k][0];
				a[k-1][1] = g.adjacencies[k][1];
			}
			g.adjacencies = a;
			g.visited = new boolean[g.adjacencies.length];
			g.mapIt();
		}
		return g;
	}

	/**
	 * this function modifies the genome itself applying the DCJ operation m on it.
	 * @param m DCJ that is applied to this genome
	 */
	public void mutate_this(DCJ m){
		
		if(m.mutationInfo[2] == 0 && m.mutationInfo[3] == 0){
			// it is a split or linearization
			int x = -1;
			for(int i = 0; x == -1 && i < adjacencies.length; i++){
				if((adjacencies[i][0] == m.mutationInfo[0] && adjacencies[i][1] == m.mutationInfo[1]) ||
				   (adjacencies[i][1] == m.mutationInfo[0] && adjacencies[i][0] == m.mutationInfo[1])){
					x = i;
				}
			}
			if(x == -1){
				throw new Error("Cannot apply this mutation:"+m.print()+" to this genome. x: "+x);
			}
			int[][] a = new int[adjacencies.length+1][2];
			for(int k = 0; k < x; k++){
				a[k][0] = adjacencies[k][0];
				a[k][1] = adjacencies[k][1];
			}
			for(int k = x+1; k < adjacencies.length; k++){
				a[k-1][0] = adjacencies[k][0];
				a[k-1][1] = adjacencies[k][1];
			}
			a[a.length-2][0] = m.mutationInfo[0];
			a[a.length-2][1] = 0;
			a[a.length-1][0] = m.mutationInfo[1];
			a[a.length-1][1] = 0;
			adjacencies = a;
			if(adjacencies.length == visited.length){
				for(int i = 0; i < visited.length; i++){
					visited[i] = false;
				}
			}
			else{
				visited = new boolean[adjacencies.length];
			}
			mapIt();
			return;
		}
				
		int x = -1;
		int y = -1;
		for(int i = 0; (x == -1 || y == -1) && i < adjacencies.length; i++){
			if((adjacencies[i][0] == m.mutationInfo[0] && adjacencies[i][1] == m.mutationInfo[1]) ||
			   (adjacencies[i][1] == m.mutationInfo[0] && adjacencies[i][0] == m.mutationInfo[1])){
				x = i;
			}
			else if((adjacencies[i][0] == m.mutationInfo[2] && adjacencies[i][1] == m.mutationInfo[3]) ||
				    (adjacencies[i][1] == m.mutationInfo[2] && adjacencies[i][0] == m.mutationInfo[3])){
				y = i;
			}
		}
		if(x == -1 || y == -1){
			throw new Error("Cannot apply this mutation:"+m.print()+" to this genome: "+this.print()+" x: "+x+", y: "+y);
		}
		int i = (adjacencies[x][0] == m.mutationInfo[0] ? 1 : 0);
		int j = (adjacencies[y][0] == m.mutationInfo[3] ? 1 : 0);
		map[adjacencies[x][i]] = y;
		map[adjacencies[y][j]] = x;
		int temp = adjacencies[x][i];
		adjacencies[x][i] = adjacencies[y][j];
		adjacencies[y][j] = temp;
		if((adjacencies[x][0] == 0 && adjacencies[x][1] == 0) ||
		    adjacencies[y][0] == 0 && adjacencies[y][1] == 0){
			int z = (adjacencies[x][0] == 0 && adjacencies[x][1] == 0 ? x : y);
			int[][] a = new int[adjacencies.length-1][2];
			for(int k = 0; k < z; k++){
				a[k][0] = adjacencies[k][0];
				a[k][1] = adjacencies[k][1];
			}
			for(int k = z+1; k <adjacencies.length; k++){
				a[k-1][0] = adjacencies[k][0];
				a[k-1][1] = adjacencies[k][1];
			}
			adjacencies = a;
			if(adjacencies.length == visited.length){
				for(int k = 0; k < visited.length; k++){
					visited[k] = false;
				}
			}
			else{
				visited = new boolean[adjacencies.length];
			}
			mapIt();
		}
		return;
	}

	
	/**
	 * 
	 * @return the first value is the number of circular chromosomes, the second value is the number of
	 * linear chromosomes
	 */
	public int[] chromosomeStructure(){
		//System.out.println(print());
		int[] numbers = new int[2];
		//for(int i = 0; i < map.length; i++){
		//	System.out.print(map[i]+" ");
		//}
		//System.out.println();
		int current;
		for(int i = 0; i < visited.length; i++){
			visited[i] = false;
		}
		//first search for linear chromosomes
		for(int i = 0; i < adjacencies.length; i++){
			if(!visited[i] && (adjacencies[i][0] == 0 || adjacencies[i][1] == 0)){
				current = 0;
				int j = i;
				while(current != -1){
					visited[j] = true;
					current = (current == adjacencies[j][0] ? adjacencies[j][1] : adjacencies[j][0]);
					current += (current % 2 == 0 ? -1 : 1);
					if(current != -1){
						j = map[current];
					}
				}
				numbers[1]++;
			}
		}
		//then the number of circular chromosomes
		for(int i = 0; i < adjacencies.length; i++){
			if(!visited[i]){
				current = adjacencies[i][0];
				int j = i;
				while(!visited[j]){
					visited[j] = true;
					//System.out.println("Visiting ("+adjacencies[j][0]+","+adjacencies[j][1]+")");
					current = (current == adjacencies[j][0] ? adjacencies[j][1] : adjacencies[j][0]);
					current += (current % 2 == 0 ? -1 : 1);
					j = map[current];
				}
				//System.out.println("end of chromosome, now at ("+adjacencies[j][0]+","+adjacencies[j][1]+")");
				numbers[0]++;
			}
		}
		return numbers;
	}
	
	public ArrayList<Chromosome> chromosomes(){
		ArrayList<Chromosome> list = new ArrayList<Chromosome>();
		int current;
		for(int i = 0; i < visited.length; i++){
			visited[i] = false;
		}
		//first search for linear chromosomes
		for(int i = 0; i < adjacencies.length; i++){
			if(!visited[i] && (adjacencies[i][0] == 0 || adjacencies[i][1] == 0)){
				current = 0;
				int j = i;
				ArrayList<int[]> a = new ArrayList<int[]>();
				while(current != -1){
					a.add(new int[] {adjacencies[j][0],adjacencies[j][1]});
					visited[j] = true;
					current = (current == adjacencies[j][0] ? adjacencies[j][1] : adjacencies[j][0]);
					current += (current % 2 == 0 ? -1 : 1);
					if(current != -1){
						j = map[current];
					}
				}
				list.add(new Chromosome(a));
			}
		}
		//then the number of circular chromosomes
		for(int i = 0; i < adjacencies.length; i++){
			if(!visited[i]){
				current = adjacencies[i][0];
				int j = i;
				ArrayList<int[]> a = new ArrayList<int[]>();
				while(!visited[j]){
					visited[j] = true;
					a.add(new int[] {adjacencies[j][0],adjacencies[j][1]});
					//System.out.println("Visiting ("+adjacencies[j][0]+","+adjacencies[j][1]+")");
					current = (current == adjacencies[j][0] ? adjacencies[j][1] : adjacencies[j][0]);
					current += (current % 2 == 0 ? -1 : 1);
					j = map[current];
				}
				//System.out.println("end of chromosome, now at ("+adjacencies[j][0]+","+adjacencies[j][1]+")");
				list.add(new Chromosome(a));
			}
		}
		return list;
	}
	
	
	public String print(){
		String s = "Adjacencies: \n";
		for(int i = 0; i < adjacencies.length; i++){
			s += "("+adjacencies[i][0]+","+adjacencies[i][1]+"), ";
		}
		s+="\n";
		return s;
	}
	
	public String printChromosomes(){
		String s = "";
		ArrayList<Chromosome> list = chromosomes();
		for(Chromosome ch : list){
			s += ch.print()+"\n";
		}
		
		return s;
	}
	
	public int dcjDistance(Genome g){
		int d = 0;
		for(int i = 0; i < adjacencies.length; i++){
			if(adjacencies[i][0] != 0){d++;}
			if(adjacencies[i][1] != 0){d++;}
		}
		AdjacencyGraph ag = new AdjacencyGraph(this, g);
		Decomposition deco = ag.getDecomposition();
		return d/2 - deco.circles.size() - deco.oddpaths.size()/2;
	}
	
	/**
	 * if returns with 0, it acts on a circular chromosome
	 * if returns with 1, it acts on a linear chromosome
	 * @param dcj
	 * @return
	 */
	public int typeOfChromosomeOnWhichActs(int extremity){
		int start = map[extremity];
		int actual = extremity + (extremity %2 == 0 ? -1 : +1);
		int actualposition = map[actual];
		while(actualposition != start){
			actual = adjacencies[actualposition][0] == actual ? adjacencies[actualposition][1] : adjacencies[actualposition][0];
			if(actual == 0){
				return 1;
			}
			actual += actual % 2 == 0 ? -1 : +1;
			actualposition = map[actual];
		}
		return 0;
	}
	
	public boolean actsOnOneChromosome(DCJ dcj){
		int start = map[dcj.mutationInfo[0]];
		int end = map[dcj.mutationInfo[2]];
		if(end == -1){
			end = map[dcj.mutationInfo[3]];
		}
		int actual = dcj.mutationInfo[0] + (dcj.mutationInfo[0] %2 == 0 ? -1 : +1);
		if(actual != -1){
			int actualposition = map[actual];
			while(actualposition != start){
				actual = adjacencies[actualposition][0] == actual ? adjacencies[actualposition][1] : adjacencies[actualposition][0];
				if(actual == 0){
					break;
				}
				if(actualposition == end){
					return true;
				}
				actual += actual % 2 == 0 ? -1 : +1;
				actualposition = map[actual];
			}
		}
		actual = dcj.mutationInfo[1] + (dcj.mutationInfo[1] %2 == 0 ? -1 : +1);
		if(actual != -1){
			int actualposition = map[actual];
			while(actualposition != start){
				actual = adjacencies[actualposition][0] == actual ? adjacencies[actualposition][1] : adjacencies[actualposition][0];
				if(actual == 0){
					break;
				}
				if(actualposition == end){
					return true;
				}
				actual += actual % 2 == 0 ? -1 : +1;
				actualposition = map[actual];
			}
		}
		return false;
	}
	
	public int numberofHPs(){
		ArrayList<Chromosome> chromosomes = chromosomes();
		int count = 0;
		int cum = 0;
		int linearchr = 0;
		for(Chromosome c : chromosomes){
		//	System.out.println("This is a " + (c.isLinear() ? "linear" : "circular")+" genome with "+c.adjacencies.length+" adjacencies");
			count += (c.isLinear() ? (c.adjacencies.length*(c.adjacencies.length-1))/2 - 1 + c.adjacencies.length-2 :(c.adjacencies.length*(c.adjacencies.length-1))/2);
			if(c.isLinear()){
				count += 2 * cum * c.adjacencies.length - 4 * linearchr; 
				cum += c.adjacencies.length;
				linearchr++;
			}
		}
		return count;
	}
	
	/**
	 * For testing purposes
	 * @param args
	 */
	public static void main(String[] args){
		//Genome g1 = new Genome(new int[][] {{4,0},{1,0}, {2,3}, {5,6}, {7,8}, {0,9}, {10,11}, {0, 12}, {13,14}, {15,16}, {17, 18}});
		Genome g2 = new Genome(new int[][] {{1,2}, {3,4}, {0,5}, {6,7}, {0,8}, {9,10}, {11,0}, {12,13}, {14,0}, {15,18}, {16,17}});
		Genome g1 = new Genome(new int[][] {{2,3}, {4,5}, {6,0}, {0,7}, {8,9}, {10,11}, {12,13}, {14,0}, {0,15}, {16,17}, {18,19}, {20,1}});
		//Genome g2 = new Genome(new int[][] {{1,4}, {3,6}, {5,8}, {7,10}, {9,12}, {11,14}, {13,16}, {15,18}, {17,20}, {19,2}});
		System.out.println(g1.actsOnOneChromosome(new DCJ(new int[] {4,5,12,13})));
		System.out.println(g2.numberofHPs());
		/*
		for(int i = 0; i < 1000; i++){
			DCJMetropolisContainer mc = DCJ.randomMutation(g1,g2);
			System.out.println(mc.print());
			g1 = g1.mutate(mc.mutation);
		}
		*/
		
		/*
		Genome g2 = new Genome(new int[][] {{1,2}, {0,3}, {4,5}, {6,7}, {7,8}, {8,9}});
		System.out.println(g1.equals(g2));
		int[] structure = g1.chromosomeStructure();
		System.out.println(structure[0]+" "+structure[1]);
		Genome g3 = g1.mutate(new DCJ(new int[]{0,1,0,10}));
		System.out.println(g3.printChromosomes());
		g3 = g3.mutate(new DCJ(new int[] {0,17,0,18}));
		System.out.println(g3.printChromosomes());
		g3 = g3.mutate(new DCJ(new int[] {12,13,2,3}));
		System.out.println(g3.printChromosomes());
		g3 = g3.mutate(new DCJ(new int[] {12, 2, 0, 0}));
		System.out.println(g3.print());
		structure = g3.chromosomeStructure();
		System.out.println(structure[0]+" "+structure[1]);
		System.out.println(g3.printChromosomes());
		*/
	}

}
