package mutations;

import representation.AdjacencyGraph;
import representation.Circle;
import representation.Decomposition;
import representation.MPath;
import representation.OddPath;
import representation.WPath;
import utils.Utils;
import genome.Genome;
import mcmc.DCJMetropolisContainer;

public class DCJ {
	
	public int[] mutationInfo;

	public DCJ(int[] m){
		mutationInfo = m;
	}

	/**
	 * copy constructor
	 * @param d this DCJ is copied
	 */
	public DCJ(DCJ d){
		mutationInfo = new int[4];
		for(int i = 0; i < mutationInfo.length; i++){
			mutationInfo[i] = d.mutationInfo[i];
		}
	}
	
	public DCJ(String descriptor){
		String[] ss = descriptor.split("[(,)|]");
		mutationInfo = new int[4];
		//for(int i = 0; i < ss.length; i++){
		//	System.out.println(ss[i]);
		//}
		for(int i = 0; i < mutationInfo.length; i++){
			mutationInfo[i] = Integer.valueOf(ss[i+1]);
		}
	}
	
	public String print() {
		return "("+mutationInfo[0]+","+mutationInfo[1]+"|"+mutationInfo[2]+","+mutationInfo[3]+")";
	}
	
	public static double chooseSorting = 0.95;
	public static double chooseNeutral = 0.04;
	public static double chooseBad = 0.01;
	
	public static double backProposalProbability(Genome startGenome, Genome towardsGenome, DCJ dcj){
		double pr = 1.0;
		AdjacencyGraph ag = new AdjacencyGraph(startGenome, towardsGenome);
		Decomposition deco = ag.getDecomposition();
		int sorting = deco.c2+ // cycle increasing in cycles
						deco.w2+ // cycle increasing in w-paths
						deco.m2+ // cycle increasing in m-paths
						deco.odd2+ // cycle increasing in odd paths
						deco.m*(2*deco.w+1); // cycle increasing by combining an m and w path or splitting an m path
		int neutral = deco.c2+ //acting on a cycle without splitting
		  			  deco.w2-deco.wpaths.size()+ //acting on a W path without splitting or circularization
		  			  deco.m2+ //acting on an M path without splitting
		  			  deco.odd2+ //acting on an odd path without splitting
		  			  2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1))+ //acting on or joining two W paths
		  			  2*deco.cumm+ //acting on two M paths
		  			  deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2+ // acting on two odd paths
		  			  2*(deco.m+deco.w)*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size()+ // acting on or merging a W and an odd path or acting on an M and an odd path
		  			  deco.w-2*deco.wpaths.size()+ // splitting a W path
		  			  deco.odd-deco.oddpaths.size(); // splitting an odd path
		int bad = 2*deco.cumc+ // merging two cycles
				  deco.cumodd+ // combining two odd paths into an M and a W path
				  2*(deco.m+deco.w+deco.odd)*deco.c+ //merging a cycle and an odd, W or M path
				  deco.c; // opening a cycle
		int oldDistance = startGenome.dcjDistance(towardsGenome);
		Genome newGenome = startGenome.mutate(dcj);
		int newDistance = newGenome.dcjDistance(towardsGenome);
		if(oldDistance - newDistance == 1){
			pr = (chooseSorting + (neutral + bad == 0 ? 1.0 - chooseSorting : 0.0))/sorting; 
			//System.out.println("sorting: "+sorting+"\tneutral: "+neutral+"\tbad: "+bad+"\tSorting!");
			//System.out.println(pr);
		}
		else if(oldDistance - newDistance == 0){
			pr = ((sorting == 0 ? chooseSorting : 0.0)+ chooseNeutral + (bad == 0? chooseBad : 0.0))/neutral;
			//System.out.println("sorting: "+sorting+"\tneutral: "+neutral+"\tbad: "+bad+"\tNeutral!");
			//System.out.println(pr);
		}
		else if(oldDistance - newDistance == -1){
			pr = ((sorting == 0 ? chooseSorting : 0.0)+ (neutral == 0? chooseNeutral : (sorting == 0 ? -chooseSorting : 0.0)) + chooseBad)/bad;		
			//System.out.println("sorting: "+sorting+"\tneutral: "+neutral+"\tbad: "+bad+"\tBad!");
			//System.out.println(pr);
		}
		else{
			throw new Error("DCJ distance changed by more than 1, or less than -1!!! \n"+dcj.print()+"\n start"+startGenome.print()+"\ntowards"+towardsGenome.print());
		}
		return pr;
	}
	
	public static DCJMetropolisContainer randomMutation(Genome startGenome, Genome towardsGenome){
		DCJMetropolisContainer mc = new DCJMetropolisContainer();
		AdjacencyGraph ag = new AdjacencyGraph(startGenome, towardsGenome);
		Decomposition deco = ag.getDecomposition();
	//	System.out.println("Decomposition: \n"+deco.print());
		mc.proposalProbability = 1.0;
		double x = Utils.generator.nextDouble();
		int sorting = deco.c2+ // cycle increasing in cycles
					  deco.w2+ // cycle increasing in w-paths
					  deco.m2+ // cycle increasing in m-paths
					  deco.odd2+ // cycle increasing in odd paths
					  deco.m*(2*deco.w+1); // cycle increasing by combining an m and w path or splitting an m path
		int neutral = deco.c2+ //acting on a cycle without splitting
					  deco.w2-deco.wpaths.size()+ //acting on a W path without splitting or circularization
					  deco.m2+ //acting on an M path without splitting
					  deco.odd2+ //acting on an odd path without splitting
					  2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1))+ //acting on or joining two W paths
					  2*deco.cumm+ //acting on two M paths
					  deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2+ // acting on two odd paths
					  2*(deco.m+deco.w)*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size()+ // acting on or merging a W and an odd path or acting on an M and an odd path
					  deco.w-2*deco.wpaths.size()+ // splitting a W path
					  deco.odd-deco.oddpaths.size(); // splitting an odd path
		int bad = 2*deco.cumc+ // merging two cycles
				  deco.cumodd+ // combining two odd paths into an M and a W path
				  2*(deco.m+deco.w+deco.odd)*deco.c+ //merging a cycle and an odd, W or M path
				  deco.c; // opening a cycle
	//	System.out.println("sorting: "+sorting+"\tneutral: "+neutral+"\tbad: "+bad);
		if((sorting > 0 && x < chooseSorting) || (neutral + bad == 0)){
			//Sorting mutations
			//System.out.println("\t Sorting! c2: "+deco.c2);
			mc.proposalProbability *= (chooseSorting + (neutral + bad == 0 ? (1.0-chooseSorting) : 0.0))/sorting;
			int t = Utils.generator.nextInt(sorting);
			int count = 0;
		//	System.out.println("\t Sorting! c2: "+deco.c2+" t: "+t);
			//trying to find a cycle increasing in cycles
			if(t < deco.c2){
				int length;
				Circle c;
				for(int i = 0; i < deco.circles.size(); i++){
					c = deco.circles.get(i);
					length = c.labels.length/2;
					count += length*(length-1)/2;
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{c.labels[2*a+1],c.labels[2*a],c.labels[2*b],c.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.c2;
			// or cycle increasing in m-paths
			if(t < count + deco.m2){
				int length;
				MPath m;
				for(int i = 0; i < deco.mpaths.size(); i++){
					m = deco.mpaths.get(i);
					length = m.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{m.labels[2*a+1],m.labels[2*a],m.labels[2*b],m.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.m2;
			// or cycle increasing in w-paths
			if(t < count + deco.w2){
				int length;
				WPath w;
				for(int i = 0; i < deco.wpaths.size(); i++){
					w = deco.wpaths.get(i);
					length = w.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{w.labels[2*a+1],w.labels[2*a],w.labels[2*b],w.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.w2;
			// or cycle increasing in odd paths
			if(t < count + deco.odd2){
				int length;
				OddPath odd;
				for(int i = 0; i < deco.oddpaths.size(); i++){
					odd = deco.oddpaths.get(i);
					length = odd.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{odd.labels[2*a+1],odd.labels[2*a],odd.labels[2*b],odd.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.odd2;
			// or combining m and w paths to two odd paths
			if(t < count + 2*deco.m*deco.w){
				int length;
				WPath w;
				int t1 = Utils.generator.nextInt(deco.w);
				int count1 = 0;
				for(int i = 0; i < deco.wpaths.size(); i++){
					w = deco.wpaths.get(i);
					length = w.labels.length/2;
					count1 += length;
					if(count1 > t1){
						int a = Utils.generator.nextInt(length);
						MPath m;
						t1 = Utils.generator.nextInt(deco.m);
						count1 = 0;
						for(int j = 0; j < deco.mpaths.size(); j++){
							m = deco.mpaths.get(j);
							length = m.labels.length/2;
							count1 += length;
							if(count1 > t1){
								int b = Utils.generator.nextInt(length);
								int parite = Utils.generator.nextInt(2);
								mc.mutation = new DCJ(new int[]{m.labels[2*b+parite],m.labels[2*b+1-parite],w.labels[2*a],w.labels[2*a+1]});
								return mc;
							}
						}
						
					}
				}
			}
			// or splitting an m path into two odd paths
			count += 2*deco.m*deco.w;
			if(t < count + deco.m){
				MPath m;
				int t1 = Utils.generator.nextInt(deco.m);
				int count1 = 0;
				int length;
				for(int j = 0; j < deco.mpaths.size(); j++){
					m = deco.mpaths.get(j);
					length = m.labels.length/2;
					count1 += length;
					if(count1 > t1){
						int b = Utils.generator.nextInt(length);
						mc.mutation = new DCJ(new int[]{m.labels[2*b],m.labels[2*b+1],0,0});
						return mc;
					}
				}
			}
			//System.out.println(deco.c2+" "+ // cycle increasing in cycles
			//		  deco.w2+" "+ // cycle increasing in w-paths
			//		  deco.m2+" "+ // cycle increasing in m-paths
			//		  deco.odd2+" "+ // cycle increasing in odd paths
			//		  deco.m*(2*deco.w+1)+" count: "+count+" t:"+t);
				
			
			throw new Error("Couldn't generate a new sorting DCJ mutation!!!");
			
		}
		else if((neutral > 0 && x < chooseSorting + chooseNeutral) || bad == 0){
			//Neutral mutations
			//System.out.println("\t Neutral!");
			mc.proposalProbability *= ((sorting == 0 ? chooseSorting : 0.0) + chooseNeutral + (bad == 0 ? chooseBad : 0.0))/neutral;
			int t = Utils.generator.nextInt(neutral);
			int count = 0;
		//	deco.c2+ //acting on a cycle without splitting
			if(t < deco.c2){
				//System.out.println("deco.c2");
				int length;
				Circle c;
				for(int i = 0; i < deco.circles.size(); i++){
					c = deco.circles.get(i);
					length = c.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{c.labels[2*a],c.labels[2*a+1],c.labels[2*b],c.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.c2;
		//	  deco.w2-deco.wpaths.size()+ //acting on a W path without splitting or circularization
			if(t < count + deco.w2 - deco.wpaths.size()){
				//System.out.println("deco.w2 - deco.wpaths.size()");
				int length;
				WPath w;
				for(int i = 0; i < deco.wpaths.size(); i++){
					w = deco.wpaths.get(i);
					length = w.labels.length/2;
					count += length*(length-1)/2-1;
					if(count > t){
						int a = 0,b = length - 1;
						while(a == 0 && b == length - 1){
							a = Utils.generator.nextInt(length);
							b = Utils.generator.nextInt(length-1);
							if(b >= a){
								b++;
							}else{
								int temp = a;
								a = b;
								b = temp;
							}
						}
						mc.mutation = new DCJ(new int[]{w.labels[2*a],w.labels[2*a+1],w.labels[2*b],w.labels[2*b+1]});
						return mc;
					}
				}
			}
			count+= deco.w2 - deco.wpaths.size();
		//	  deco.m2+ //acting on an M path without splitting
			if(t < count + deco.m2){
				//System.out.println("deco.m2");
				int length;
				MPath m;
				for(int i = 0; i < deco.mpaths.size(); i++){
					m = deco.mpaths.get(i);
					length = m.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{m.labels[2*a],m.labels[2*a+1],m.labels[2*b],m.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.m2;
		//	  deco.odd2+ //acting on an odd path without splitting
			if(t < count + deco.odd2){
				//System.out.println("deco.odd2");
				int length;
				OddPath odd;
				for(int i = 0; i < deco.oddpaths.size(); i++){
					odd = deco.oddpaths.get(i);
					length = odd.labels.length/2;
					count += length*(length-1);
					if(count > t){
						int a = Utils.generator.nextInt(length);
						int b = Utils.generator.nextInt(length-1);
						if(b >= a){
							b++;
						}else{
							int temp = a;
							a = b;
							b = temp;
						}
						mc.mutation = new DCJ(new int[]{odd.labels[2*a],odd.labels[2*a+1],odd.labels[2*b],odd.labels[2*b+1]});
						return mc;
					}
				}
			}
			count += deco.odd2;
		//	  2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1))+ //acting on or joining two W paths
			if(t < count + 2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1))){
				//System.out.println("2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1)");
				int i = 0, j = 0;
				int a = 0, b = 0;
				int count1;
				int length;
				boolean accept = false;
				while(!accept){
					i = j;
					while(i == j){
						count1 = Utils.generator.nextInt(deco.w);
						length = 0;
						for(int k = 0; k < deco.wpaths.size(); k++){
							length += deco.wpaths.get(k).labels.length/2;
							if(length > count1){
								i = k;
								break;
							}
						}
						count1 = Utils.generator.nextInt(deco.w);
						length = 0;
						for(int k = 0; k < deco.wpaths.size(); k++){
							length += deco.wpaths.get(k).labels.length/2;
							if(length > count1){
								j = k;
								break;
							}
						}
					}
					a = Utils.generator.nextInt(deco.wpaths.get(i).labels.length/2);
					b = Utils.generator.nextInt(deco.wpaths.get(j).labels.length/2);
					if((a == 0 || a == deco.wpaths.get(i).labels.length/2-1) && (b == 0 || b == deco.wpaths.get(j).labels.length/2-1)){
						accept = Utils.generator.nextDouble() < 0.5;
					} else{
						accept = true;
					}
				}
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.wpaths.get(i).labels[2*a+parite], deco.wpaths.get(i).labels[2*a+1-parite],
									deco.wpaths.get(j).labels[2*b], deco.wpaths.get(j).labels[2*b+1]};
				if((m[0] == 0 && m[3] == 0) ||(m[1] == 0 && m[2] == 0)){
					int temp = m[0];
					m[0] = m[1];
					m[1] = temp;
				}
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1));
		//	  2*deco.cumm+ //acting on two M paths
			if(t < count + 2*deco.cumm){
				//System.out.println("2*deco.cumm");
				int i = 0, j = 0;
				int a = 0, b = 0;
				int count1;
				int length;
				
				i = j;
				while(i == j){
					count1 = Utils.generator.nextInt(deco.m);
					length = 0;
					for(int k = 0; k < deco.mpaths.size(); k++){
						length += deco.mpaths.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}
					}
					count1 = Utils.generator.nextInt(deco.m);
					length = 0;
					for(int k = 0; k < deco.mpaths.size(); k++){
						length += deco.mpaths.get(k).labels.length/2;
						if(length > count1){
							j = k;
							break;
						}
					}
				}
				a = Utils.generator.nextInt(deco.mpaths.get(i).labels.length/2);
				b = Utils.generator.nextInt(deco.mpaths.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.mpaths.get(i).labels[2*a+parite], deco.mpaths.get(i).labels[2*a+1-parite],
									deco.mpaths.get(j).labels[2*b], deco.mpaths.get(j).labels[2*b+1]};
				mc.mutation = new DCJ(m);
				return mc;
				
			}
			count += 2*deco.cumm;
		//	  deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2+ // acting on two odd paths
			if(t < count + deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2){
				//System.out.println("deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2");
				int i = 0, j = 0;
				int a = 0, b = 0;
				int count1;
				int length;
				while(a == 0 && b == 0){
					i = j;
					while(i == j){
						count1 = Utils.generator.nextInt(deco.odd);
						length = 0;
						for(int k = 0; k < deco.oddpaths.size(); k++){
							length += deco.oddpaths.get(k).labels.length/2;
							if(length > count1){
								i = k;
								break;
							}
						}
						count1 = Utils.generator.nextInt(deco.odd);
						length = 0;
						for(int k = 0; k < deco.oddpaths.size(); k++){
							length += deco.oddpaths.get(k).labels.length/2;
							if(length > count1){
								j = k;
								break;
							}
						}
					}
					a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
					b = Utils.generator.nextInt(deco.oddpaths.get(j).labels.length/2);
					
					
				}
				int[] m = new int[]{deco.oddpaths.get(i).labels[2*a], deco.oddpaths.get(i).labels[2*a+1],
									deco.oddpaths.get(j).labels[2*b+1], deco.oddpaths.get(j).labels[2*b]};
				mc.mutation = new DCJ(m);
				return mc;
				
			}
			count += deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2;
		//	  2*(deco.m+deco.w)*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size()+ // acting on or merging a W and an odd path or acting on an M and an odd path
			if(t < count + 2*deco.m*deco.odd){
				//System.out.println("2*deco.m*deco.odd");
				int i = 0;
				int j = 0;
				int a = 0, b = 0;
				int count1 = Utils.generator.nextInt(deco.odd);
				int length = 0;
				for(int k = 0; k < deco.oddpaths.size(); k++){
					length += deco.oddpaths.get(k).labels.length/2;
					if(length > count1){
						i = k;
						break;
					}
				}
				count1 = Utils.generator.nextInt(deco.m);
				length = 0;
				for(int k = 0; k < deco.mpaths.size(); k++){
					length += deco.mpaths.get(k).labels.length/2;
					if(length > count1){
						j = k;
						break;
					}
				}
				a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
				b = Utils.generator.nextInt(deco.mpaths.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.oddpaths.get(i).labels[2*a], deco.oddpaths.get(i).labels[2*a+1],
						deco.mpaths.get(j).labels[2*b+parite], deco.mpaths.get(j).labels[2*b+1-parite]};
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2*deco.m*deco.odd;
			if(t < count + 2*deco.w*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size()){
				//System.out.println("2*deco.w*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size()");
				int i = 0;
				int j = 0;
				int a = 0, b = 0;
				boolean accept = false;
				while(!accept){
					int count1 = Utils.generator.nextInt(deco.odd);
					int length = 0;
					for(int k = 0; k < deco.oddpaths.size(); k++){
						length += deco.oddpaths.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}
					}
					count1 = Utils.generator.nextInt(deco.w);
					length = 0;
					for(int k = 0; k < deco.wpaths.size(); k++){
						length += deco.wpaths.get(k).labels.length/2;
						if(length > count1){
							j = k;
							break;
						}
					}
					a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
					b = Utils.generator.nextInt(deco.wpaths.get(j).labels.length/2);
					if(a == 0 && (b == 0 || b == deco.wpaths.get(j).labels.length/2-1)){
						accept = Utils.generator.nextDouble() < 0.5;
					}else{
						accept = true;
					}
				}
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.oddpaths.get(i).labels[2*a+parite], deco.oddpaths.get(i).labels[2*a+1-parite],
									deco.wpaths.get(j).labels[2*b], deco.wpaths.get(j).labels[2*b+1]};
				if((m[0] == 0 && m[3] == 0) ||(m[1] == 0 && m[2] == 0)){
					int temp = m[0];
					m[0] = m[1];
					m[1] = temp;
				}
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2*deco.w*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size();
		//	  deco.w-2*deco.wpaths.size()+ // splitting a W path
			if(t < count + deco.w-2*deco.wpaths.size()){
				//System.out.println("deco.w-2*deco.wpaths.size()");
				int i = 0;
				int a = 0;
				boolean accept = false;
				while(!accept){
					int count1 = Utils.generator.nextInt(deco.w);
					int length = 0;
					for(int k = 0; k < deco.wpaths.size(); k++){
						length += deco.wpaths.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}						
					}
					a = Utils.generator.nextInt(deco.wpaths.get(i).labels.length/2);
					accept = (a !=0 && a != deco.wpaths.get(i).labels.length/2-1);
				}
				mc.mutation = new DCJ(new int[] {deco.wpaths.get(i).labels[2*a],deco.wpaths.get(i).labels[2*a+1],0,0});
				return mc;
			}
			count += deco.w-2*deco.wpaths.size();
		//	  deco.odd-deco.oddpaths.size(); // splitting an odd path
			if(t < count + deco.odd-deco.oddpaths.size()){
				//System.out.println("deco.odd-deco.oddpaths.size()");
				int i = 0;
				int a = 0;
				while(a == 0){
					int count1 = Utils.generator.nextInt(deco.odd);
					int length = 0;
					for(int k = 0; k < deco.oddpaths.size(); k++){
						length += deco.oddpaths.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}						
					}
					a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
				}
				mc.mutation = new DCJ(new int[] {deco.oddpaths.get(i).labels[2*a],deco.oddpaths.get(i).labels[2*a+1],0,0});
				return mc;
			}
			//System.out.println(deco.c2+" "+ //acting on a cycle without splitting
			//		  (deco.w2-deco.wpaths.size())+" "+ //acting on a W path without splitting or circularization
			//		  deco.m2+" "+ //acting on an M path without splitting
			//		  deco.odd2+" "+ //acting on an odd path without splitting
			//		  (2*deco.cumw-2*(deco.wpaths.size() * (deco.wpaths.size()-1)))+" "+ //acting on or joining two W paths
			//		  (2*deco.cumm)+" "+ //acting on two M paths
			//		  (deco.cumodd-(deco.oddpaths.size()*(deco.oddpaths.size()-1))/2)+" "+ // acting on two odd paths
			//		  (2*(deco.m+deco.w)*deco.odd-2*deco.wpaths.size()*deco.oddpaths.size())+" "+ // acting on or merging a W and an odd path or acting on an M and an odd path
			//		  (deco.w-2*deco.wpaths.size())+" "+ // splitting a W path
			//		  (deco.odd-deco.oddpaths.size())+" count: "+count+" t:"+t); // splitting an odd path);
			throw new Error("Cannot generate a neutral dcj!");
		} 
		else{
			//Bad mutations
			//System.out.println("\t Bad!");
			mc.proposalProbability *= ((sorting == 0 ? chooseSorting : 0.0) + (neutral == 0 ? chooseNeutral : (sorting == 0 ? -chooseSorting : 0.0)) + chooseBad)/bad;
			int t = Utils.generator.nextInt(bad);
			int count = 0;
			//2*deco.cumc+ // merging two cycles
			if(t < 2*deco.cumc){
				int i = 0, j = 0;
				while(i == j){
					int count1 = Utils.generator.nextInt(deco.c);
					int length = 0;
					for(int k = 0; k < deco.circles.size(); k++){
						length += deco.circles.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}
					}
					count1 = Utils.generator.nextInt(deco.c);
					length = 0;
					for(int k = 0; k < deco.circles.size(); k++){
						length += deco.circles.get(k).labels.length/2;
						if(length > count1){
							j = k;
							break;
						}
					}
				}
				int a = Utils.generator.nextInt(deco.circles.get(i).labels.length/2);
				int b = Utils.generator.nextInt(deco.circles.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.circles.get(i).labels[2*a], deco.circles.get(i).labels[2*a+1],
						deco.circles.get(j).labels[2*b+parite], deco.circles.get(j).labels[2*b+1-parite]};
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2*deco.cumc;
			 // deco.cumodd+ // combining two odd paths into an M and a W path
			if(t < count + deco.cumodd){
				int i = 0, j = 0;
				while(i == j){
					int count1 = Utils.generator.nextInt(deco.odd);
					int length = 0;
					for(int k = 0; k < deco.oddpaths.size(); k++){
						length += deco.oddpaths.get(k).labels.length/2;
						if(length > count1){
							i = k;
							break;
						}
					}
					count1 = Utils.generator.nextInt(deco.odd);
					length = 0;
					for(int k = 0; k < deco.oddpaths.size(); k++){
						length += deco.oddpaths.get(k).labels.length/2;
						if(length > count1){
							j = k;
							break;
						}
					}
				}
				int a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
				int b = Utils.generator.nextInt(deco.oddpaths.get(j).labels.length/2);
				int[] m = new int[]{deco.oddpaths.get(i).labels[2*a], deco.oddpaths.get(i).labels[2*a+1],
						deco.oddpaths.get(j).labels[2*b], deco.oddpaths.get(j).labels[2*b+1]};
				mc.mutation = new DCJ(m);
				return mc;
				
			}
			count += deco.cumodd;
			 // 2*(deco.m+deco.w+deco.odd)*deco.c+ //merging a cycle and an odd, W or M path
			if(t < count + 2* deco.m*deco.c){
				int i = 0;
				int j = 0;
				int a = 0, b = 0;
				int count1 = Utils.generator.nextInt(deco.m);
				int length = 0;
				for(int k = 0; k < deco.mpaths.size(); k++){
					length += deco.mpaths.get(k).labels.length/2;
					if(length > count1){
						i = k;
						break;
					}
				}
				count1 = Utils.generator.nextInt(deco.c);
				length = 0;
				for(int k = 0; k < deco.circles.size(); k++){
					length += deco.circles.get(k).labels.length/2;
					if(length > count1){
						j = k;
						break;
					}
				}
				a = Utils.generator.nextInt(deco.mpaths.get(i).labels.length/2);
				b = Utils.generator.nextInt(deco.circles.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.mpaths.get(i).labels[2*a], deco.mpaths.get(i).labels[2*a+1],
						deco.circles.get(j).labels[2*b+parite], deco.circles.get(j).labels[2*b+1-parite]};
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2* deco.m*deco.c;
			if(t < count + 2* deco.w*deco.c){
				int i = 0;
				int j = 0;
				int a = 0, b = 0;
				int count1 = Utils.generator.nextInt(deco.w);
				int length = 0;
				for(int k = 0; k < deco.wpaths.size(); k++){
					length += deco.wpaths.get(k).labels.length/2;
					if(length > count1){
						i = k;
						break;
					}
				}
				count1 = Utils.generator.nextInt(deco.c);
				length = 0;
				for(int k = 0; k < deco.circles.size(); k++){
					length += deco.circles.get(k).labels.length/2;
					if(length > count1){
						j = k;
						break;
					}
				}
				a = Utils.generator.nextInt(deco.wpaths.get(i).labels.length/2);
				b = Utils.generator.nextInt(deco.circles.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.wpaths.get(i).labels[2*a], deco.wpaths.get(i).labels[2*a+1],
						deco.circles.get(j).labels[2*b+parite], deco.circles.get(j).labels[2*b+1-parite]};
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2* deco.w*deco.c;
			if(t < count + 2* deco.odd*deco.c){
				int i = 0;
				int j = 0;
				int a = 0, b = 0;
				int count1 = Utils.generator.nextInt(deco.odd);
				int length = 0;
				for(int k = 0; k < deco.oddpaths.size(); k++){
					length += deco.oddpaths.get(k).labels.length/2;
					if(length > count1){
						i = k;
						break;
					}
				}
				count1 = Utils.generator.nextInt(deco.c);
				length = 0;
				for(int k = 0; k < deco.circles.size(); k++){
					length += deco.circles.get(k).labels.length/2;
					if(length > count1){
						j = k;
						break;
					}
				}
				a = Utils.generator.nextInt(deco.oddpaths.get(i).labels.length/2);
				b = Utils.generator.nextInt(deco.circles.get(j).labels.length/2);
				int parite = Utils.generator.nextInt(2);
				int[] m = new int[]{deco.oddpaths.get(i).labels[2*a], deco.oddpaths.get(i).labels[2*a+1],
						deco.circles.get(j).labels[2*b+parite], deco.circles.get(j).labels[2*b+1-parite]};
				mc.mutation = new DCJ(m);
				return mc;
			}
			count += 2* deco.odd*deco.c;
			 //deco.c; // opening a cycle
			if(t < count + deco.c){
				int i = 0;
				int count1 = Utils.generator.nextInt(deco.c);
				int length = 0;
				for(int k = 0; k < deco.circles.size(); k++){
					length += deco.circles.get(k).labels.length/2;
					if(length > count1){
						i = k;
						break;
					}
				}
				int a = Utils.generator.nextInt(deco.circles.get(i).labels.length/2);
				mc.mutation = new DCJ(new int[] {deco.circles.get(i).labels[2*a],deco.circles.get(i).labels[2*a+1],0,0});
				return mc;
			}
		}
		System.out.println("Total: "+(sorting+neutral+bad)+", of which sorting: "+sorting+", neutral: "+neutral+", bad: "+bad);
		System.out.println("c: "+deco.c+" c2: "+deco.c2+" cumc: "+deco.cumc);
		System.out.println("m: "+deco.m+" m2: "+deco.m2+" cumm: "+deco.cumm);
		System.out.println("w: "+deco.w+" w2: "+deco.w2+" cumw: "+deco.cumw);
		System.out.println("odd: "+deco.odd+" odd2: "+deco.odd2+" cumodd: "+deco.cumodd);
		throw new Error("Cannot generate a mutation!\nsorting: "+sorting+"\tneutral: "+neutral+"\tbad: "+bad+" x: "+x);
		//return mc;

	}
	
	public String typeOn(Genome g){
		Genome g1 = g.mutate(this);
		int[] c = g.chromosomeStructure();
		int[] c1 = g1.chromosomeStructure();
		if(c[0] == c1[0] && c[1] < c1[1]){
			if(g.typeOfChromosomeOnWhichActs(mutationInfo[0]) == 0){
				return "circular chromosome fission";
			}
			else{
				return "fission of circular chromosome from a linear";
			}
		}
		if(c[0] == c1[0] && c[1] > c1[1]){
			if(g.typeOfChromosomeOnWhichActs(mutationInfo[0]) == 0 &&
		       g.typeOfChromosomeOnWhichActs(mutationInfo[3]) == 0){
				return "circular chromosome fusion";
			}else{
				return "fusion of a circular chromosome into a linear";
			}
		}
		if(c[0] < c1[0] && c[1] == c1[1]){
			return "linear chromosome break";
		}
		if(c[0] > c1[0] && c[1] == c1[1]){
			return "fusion of two linear chromosome";
		}
		if(c[0] < c1[0] && c[1] > c1[1]){
			return "circular chromosome opemning";
		}
		if(c[0] > c1[0] && c[1] < c1[1]){
			return "circularization of a linear chromosome";
		}
		if(c[0] == c1[0] && c[1] == c1[1]){
			if(g.actsOnOneChromosome(this)){
				return "reversal";
			}else{
				if(mutationInfo[0] == 0 || mutationInfo[1] == 0 || mutationInfo [2] == 0 || mutationInfo[3] == 0){
					return "translocation";
				}else{
					return "reciprocal translocation";
				}
			}
		}

		return null;
	}
	
	/**
	 * For testing reasons only
	 * @param args
	 */
	public static void main(String[] args){
		DCJ dcj = new DCJ("(54,651|439,436)");
		System.out.println(dcj.print());
		
	}
}
