package io;

import genome.Genome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReadDuplicatedGenome {

	private static int transform(String s){
		int x0 = Integer.valueOf(s.substring(0, s.length()-1));
		int par = s.charAt(s.length()-1) == 'a' ? 1 :0;
		x0 *= 2;
		x0 -= x0 < 0 ? -par : par;
		//System.out.println("transformed value: "+x0);
		return x0;
	}
	
	public static Genome read(String filename) throws IOException{
		BufferedReader bf = new BufferedReader(new FileReader(filename));
		String s = "";
		ArrayList<int[]> adjacencylist = new ArrayList<int[]>();
		while((s = bf.readLine()) != null){
		//	System.out.println("s now: "+s);
			String[] ss = s.split(" ");
		//	System.out.println("Size of ss: "+ss.length);
		//	System.out.println("adjacencylist: "+adjacencylist);
		//	System.out.println("its size: "+adjacencylist.size());
		//	System.out.println(ss[0]);
			int x0 = transform(ss[0]);
			adjacencylist.add(new int[]{0, x0 > 0 ? 2*x0 - 1 : -2*x0});
			for(int i = 1; i < ss.length; i++){
				x0 = transform(ss[i-1]);
				int x1 = transform(ss[i]);
				adjacencylist.add(new int[] 
				       {x0 > 0 ? 2*x0 : -2*x0-1,
						x1 > 0 ? 2*x1 - 1 : -2*x1});
			}
			x0 = transform(ss[ss.length-1]);
			adjacencylist.add(new int[] 
			       {x0 > 0 ? 2*x0 : -2*x0-1,
					0});
		}
		int[][] adjacencies = new int[adjacencylist.size()][2];
		for(int i = 0; i < adjacencies.length; i++){
			int[] x = adjacencylist.get(i);
			adjacencies[i][0] = x[0];
			adjacencies[i][1] = x[1];
		}
		//System.out.println("adjacencies size: "+adjacencies.length);
		//for(int[] x : adjacencies){
			//System.out.println("("+x[0]+","+x[1]+")");
	//	}
		return new Genome(adjacencies);
		//	return null;
	}

	public static void main(String[] args) throws IOException{
		Genome g = read("data/duplicated_ancestor.genome");
	}
	
}
