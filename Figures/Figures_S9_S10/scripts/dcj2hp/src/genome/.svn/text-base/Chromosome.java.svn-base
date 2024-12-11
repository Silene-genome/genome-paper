package genome;

import java.util.ArrayList;

public class Chromosome {
	
	int[][] adjacencies;
	//boolean isCircular;
	
	public Chromosome(int[][] a){
		adjacencies = a;
	}
	
	public String print(){
		String s = "";
		for(int i = 0; i < adjacencies.length; i++){
			s += adjacencies[i][0] != 0 ? adjacencies[i][0] : "";
			s += (adjacencies[i][0] != 0 && adjacencies[i][1] != 0 ? " - " : "");
			s += adjacencies[i][1] != 0 ? adjacencies[i][1]+" " : " ";
		}
		return (adjacencies[0][0] == 0 ? "|" : "(")+s.substring(0, s.length()-1)+(adjacencies[0][0] == 0 ? "|" : ")");
	}
	
	public Chromosome(ArrayList<int[]> list){
		adjacencies = new int[list.size()][2];
		int[] a;
		////////
		a = list.get(0);
		adjacencies[0][0] = Math.min(a[0], a[1]);
		adjacencies[0][1] = Math.max(a[0], a[1]);
		if(a[0] != 0 && a[1] != 0 && list.size() > 1){
			int[] b = list.get(1);
			if(Math.abs(adjacencies[0][0] - b[0]) == 1 || Math.abs(adjacencies[0][0] - b[1]) == 1){
				int temp = adjacencies[0][0];
				adjacencies[0][0] = adjacencies[0][1];
				adjacencies[0][1] = temp;
			}
		//	else{
		//		System.out.println("Igz jo, "+a[0]+" "+a[1]+" "+b[0]+" "+b[1]);
		//	}
		}
		////////
		for(int i = 1; i < adjacencies.length; i++){
			a = list.get(i);
			adjacencies[i][0] = (Math.abs(adjacencies[i-1][1] - a[0]) == 1 ? a[0] : a[1]);
			adjacencies[i][1] = (adjacencies[i][0] == a[0] ? a[1] : a[0]);
		}
		
	}
	
	public boolean isLinear(){
		return adjacencies[0][0] == 0;
	}

}
