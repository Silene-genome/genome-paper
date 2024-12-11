package utils;

import java.io.*;

import mutations.DCJ;

import genome.Genome;
import io.ReadDuplicatedGenome;
import io.ReadGenome;

public class MCMCAnalysis {
	
	public static void main(String[] args) throws IOException{
		//Genome g1 = ReadGenome.read("../data/euarchontoglires.genome");
		Genome g1 = ReadDuplicatedGenome.read("../data1/duplicated_ancestor.genome");

		BufferedReader br = new BufferedReader(new FileReader("../data1/duplicated_samplesAugust2.mcmc"));
		String s;
		while((s = br.readLine()) != null){
			String[] ss = s.split("\t");
			//if(ss[2].equals("variance: 0.0")){
			if(ss[3].equals("circular: 0.0")){
				
				//System.out.println("found!");
				br.readLine();
				br.readLine();
				Genome current = new Genome(g1);
				while((s = br.readLine()).contains("(")){
					DCJ dcj = new DCJ(s);
					System.out.println(dcj.typeOn(current));
					current.mutate_this(dcj);
				}
				System.out.println("-");
			}
			else{
				//System.out.println("not found, "+ss[2]);
				br.readLine();
				br.readLine();
				while((s = br.readLine()).contains("("));
			}
		}
	}

}
