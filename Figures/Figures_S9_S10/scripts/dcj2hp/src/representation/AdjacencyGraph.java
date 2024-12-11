package representation;

import java.util.ArrayList;

import genome.Genome;
import mutations.DCJ;

public class AdjacencyGraph {
	
	
	public AdjacencyGraphNode[] set1;
	public AdjacencyGraphNode[] set2;
	public int[] location1;
	public int[] location2;
	Decomposition myDecomposition = null;
	
	public AdjacencyGraph(Genome g1, Genome g2){
		//copying the genomes into the nodes
		set1 = new AdjacencyGraphNode[g1.adjacencies.length];
		for(int i = 0; i < set1.length; i++){
			set1[i] = new AdjacencyGraphNode();
			set1[i].labels[0] = g1.adjacencies[i][0];
			set1[i].labels[1] = g1.adjacencies[i][1];
		}
		set2 = new AdjacencyGraphNode[g2.adjacencies.length];
		for(int i = 0; i < set2.length; i++){
			set2[i] = new AdjacencyGraphNode();
			set2[i].labels[0] = g2.adjacencies[i][0];
			set2[i].labels[1] = g2.adjacencies[i][1];
		}	
		//finding the largest endpoint
		int max = 0;
		for(int i = 0; i < g1.adjacencies.length; i++){
			max = Math.max(max, g1.adjacencies[i][0]);
			max = Math.max(max, g1.adjacencies[i][1]);			
		}
		max++;
		//setting up the mapping for fast constructions
		location1 = new int[max];
		location2 = new int[max];
		for(int i = 0; i < g1.adjacencies.length; i++){
			location1[g1.adjacencies[i][0]] = i;
			location1[g1.adjacencies[i][1]] = i;
		}
		for(int i = 0; i < g2.adjacencies.length; i++){
			location2[g2.adjacencies[i][0]] = i;
			location2[g2.adjacencies[i][1]] = i;
		}
		location1[0] = -1;
		location2[0] = -1;
	//and now construction of the graph
		for(int i = 0; i < set1.length; i++){
			set1[i].neighbours[0] = location2[set1[i].labels[0]];
			set1[i].neighbours[1] = location2[set1[i].labels[1]];
		}
		for(int i = 0; i < set2.length; i++){
			set2[i].neighbours[0] = location1[set2[i].labels[0]];
			set2[i].neighbours[1] = location1[set2[i].labels[1]];
		}
		
	}
	
	public Decomposition getDecomposition(){
		if(myDecomposition != null){
			return myDecomposition;
		}
		else{
			myDecomposition = new Decomposition(this);
			return myDecomposition;
		}
	}
	
	public String print(){
		String s = "";
		for(int i = 0; i < set1.length; i++){
			int n1 = set1[i].neighbours[0];
			int n2 = set1[i].neighbours[1];
			s += "("+set1[i].labels[0]+","+set1[i].labels[1]+") connected to ";
			if(n1 != -1){
				s += "("+set2[n1].labels[0]+","+set2[n1].labels[1]+")";
			}
			if(n1 != -1 && n2 != -1){
				s += " and ";
			}
			if(n2 != -1){
				s += "("+set2[n2].labels[0]+","+set2[n2].labels[1]+")";
			}
			s += "\n";
		}
		s += "-------------\n";
		for(int i = 0; i < set2.length; i++){
			int n1 = set2[i].neighbours[0];
			int n2 = set2[i].neighbours[1];
			s += "("+set2[i].labels[0]+","+set2[i].labels[1]+") connected to ";
			if(n1 != -1){
				s += "("+set1[n1].labels[0]+","+set1[n1].labels[1]+")";
			}
			if(n1 != -1 && n2 != -1){
				s += " and ";
			}
			if(n2 != -1){
				s += "("+set1[n2].labels[0]+","+set1[n2].labels[1]+")";
			}
			s += "\n";
		}
		
		return s;
	}

	public static void main(String[] args){
		Genome g1 = new Genome(new int[][] {{0,1}, {12,13}, {14,15}, {0, 17}, {0, 18}, {16,11}, {2,3}, {4,5}, {6,7}, {8,9}, {10,0}});
		g1 = g1.mutate(new DCJ(new int[] {12,13,6,7}));
		AdjacencyGraph ag = new AdjacencyGraph(
				g1,
				new Genome(new int[][] {{0,1}, {12,13}, {14,15}, {0, 17}, {0, 18}, {16,11}, {2,3}, {4,5}, {6,7}, {8,9}, {10,0}}));
		System.out.println(ag.print());
	}
	
}
