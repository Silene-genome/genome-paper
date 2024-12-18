package utils;

import java.io.*;

public class SwapAnalysis {
	
	public static void main(String[] args) throws IOException{
		int[] permutation = new int[10];
		//for(int i = 0; i < permutation.length; i++){
		//	permutation[i] = i;
		//}
		permutation[permutation.length-1] = 1;
		BufferedReader br = new BufferedReader(new FileReader("../data1/mouse_10_0.01_swap.txt"));
		String s;
		while((s = br.readLine()) != null){
			int pos = Integer.parseInt(s.split(" ")[4]);
			if(pos >= 0){
				if(s.contains("accept")){
				
					int temp = permutation[pos];
					permutation[pos] = permutation[pos+1];
					permutation[pos+1] = temp;
					permutation[permutation.length-1] = 1;
				}
			
				int sum = 0;
				for(int i = 0; i < permutation.length; i++){
					//System.out.print(permutation[i]+" ");
					sum += permutation[i];
				}
				System.out.println(sum);
			}
			/*
			if(permutation[0] == permutation.length-1){
				for(int i = 0; i < permutation.length; i++){
					System.out.print(permutation[i]+" ");
				}
				System.out.println();
				for(int i = 0; i < permutation.length; i++){
					permutation[i] = i;
				}
			}
			*/
		}
	}

}
