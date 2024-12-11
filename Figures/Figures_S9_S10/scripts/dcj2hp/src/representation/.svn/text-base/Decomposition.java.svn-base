package representation;

import java.util.ArrayList;

public class Decomposition {
	
	public ArrayList<MPath> mpaths = new ArrayList<MPath>();
	public ArrayList<WPath> wpaths = new ArrayList<WPath>();
	public ArrayList<OddPath> oddpaths = new ArrayList<OddPath>();
	public ArrayList<Circle> circles = new ArrayList<Circle>();
	public int c2 = 0; // total possible sorting dcjs on cycles
	public int w2 = 0; // total possible sorting dcjs on w paths
	public int m2 = 0; // total possible sorting dcjs on m paths
	public int w = 0; // total length of w paths
	public int m = 0; // total length of m paths;
	public int c = 0; // total length of cycles;
	public int odd2 = 0; // total possible sorting dcjs on odd paths
	public int odd = 0; // total length of odd paths
	public int cumc = 0;
	public int cumm = 0;
	public int cumw = 0;
	public int cumodd = 0;
	
	public Decomposition(AdjacencyGraph ag){
		mpaths = new ArrayList<MPath>();
		wpaths = new ArrayList<WPath>();
		oddpaths = new ArrayList<OddPath>();
		circles = new ArrayList<Circle>();
		AdjacencyGraphNode[] set1 = ag.set1;
		AdjacencyGraphNode[] set2 = ag.set2;
		boolean[] used1 = new boolean[set1.length];
		for(int i = 0; i < used1.length; i++){used1[i]= false;}
		boolean[] used2 = new boolean[set2.length];
		for(int i = 0; i < used2.length; i++){used2[i]= false;}
		//searching for odd and W paths
		for(int i = 0; i<used1.length; i++){
			if(!used1[i] && (set1[i].labels[0] == 0 || set1[i].labels[1] == 0)){
				int j = i;
				boolean isW = true;
				int count = 2;
				int currentLabel = set1[i].labels[0] + set1[i].labels[1]; //ugly but probably a bit faster
				while(currentLabel != 0){
					used1[j] = true;
					int k = ag.location2[currentLabel];
					used2[k] = true;
					currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
					count++;
					if(currentLabel != 0){
						j = ag.location1[currentLabel];
						used1[j] = true;
						currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
						count++;
					}
					else{
						isW = false;
					}
				}
				j = i;
				if(isW){
					WPath wpath = new WPath(count);
					cumw += w * count/2;
					w += count/2;
					w2 += (count*(count-2))/8;
					
					count = 1;
					currentLabel = set1[i].labels[0] + set1[i].labels[1];
					wpath.labels[0] = 0;
					wpath.labels[1] = currentLabel; 
					while(currentLabel != 0){
						int k = ag.location2[currentLabel];
						currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
						count++;
						wpath.labels[count] = currentLabel;
						if(currentLabel != 0){
							j = ag.location1[currentLabel];
							currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
							count++;
							wpath.labels[count] = currentLabel;
						}
					}
					wpaths.add(wpath);
				}
				else{
					OddPath oddpath = new OddPath(count-1);
					cumodd += odd * (count-1)/2;
					odd += (count-1)/2;
					odd2 += (count-1)*(count-3)/8;
					count = 1;
					currentLabel = set1[i].labels[0] + set1[i].labels[1];
					oddpath.labels[0] = 0;
					oddpath.labels[1] = currentLabel; 
					while(currentLabel != 0){
						int k = ag.location2[currentLabel];
						currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
						count++;
						if(currentLabel != 0){
							oddpath.labels[count] = currentLabel;
							j = ag.location1[currentLabel];
							currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
							count++;
							oddpath.labels[count] = currentLabel;
						}
					}
					oddpaths.add(oddpath);
				}
			}
		}
		//and now M paths
		for(int i = 0; i<used2.length; i++){
			//System.out.println("i: "+i);
			if(!used2[i] && (set2[i].labels[0] == 0 || set2[i].labels[1] == 0)){
				int k = i;
				int count = 1;
				int currentLabel = set2[i].labels[0] + set2[i].labels[1]; //ugly but probably a bit faster
				while(currentLabel != 0){
					used2[k] = true;
					int j = ag.location1[currentLabel];
					used1[j] = true;
					currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
					count++;
					k = ag.location2[currentLabel];
					currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
					count++;
					used2[k] = true;
				}
				MPath mpath = new MPath(count-1);
				cumm += m * (count - 1)/2;
				m += (count-1)/2;
				m2 += ((count-1)*(count-3))/8;
				k = i; count = 0;
				currentLabel = set2[i].labels[0] + set2[i].labels[1];
				while(currentLabel != 0){
					mpath.labels[count] = currentLabel;
					int j = ag.location1[currentLabel];
					currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
					count++;
					mpath.labels[count] = currentLabel;
					k = ag.location2[currentLabel];
					currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
					count++;
					//mpath.labels[count] = currentLabel;
				}
				mpaths.add(mpath);
			}
		}
		//and eventually the cycles
		for(int i = 0; i<used1.length; i++){
			if(!used1[i]){
				int j = i;
				int count = 1;
				int currentLabel = set1[i].labels[0]; //ugly but probably a bit faster
				do{
					used1[j] = true;
					int k = ag.location2[currentLabel];
					used2[k] = true;
					currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
					count++;
					j = ag.location1[currentLabel];
					currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
					count++;
					
				}while(j!=i);
				Circle circle = new Circle(count-1);
				cumc += c * (count-1)/2;
				c += (count-1)/2;
				c2 += ((count-1)*(count-3))/8;
				currentLabel = set1[i].labels[0];
				count = 0;
				do{
					int k = ag.location2[currentLabel];
					currentLabel = set2[k].labels[0] == currentLabel ? set2[k].labels[1] : set2[k].labels[0];
					circle.labels[count] = currentLabel;
					count++;
					j = ag.location1[currentLabel];
					currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
					circle.labels[count] = currentLabel;
					count++;
					//circle.labels[count] = currentLabel;
				}while(j!=i);
			//	currentLabel = set1[j].labels[0] == currentLabel ? set1[j].labels[1] : set1[j].labels[0];
			//	circle.labels[count] = currentLabel;
			
				circles.add(circle);
			}
		}
		///////////////////////////
		/*
		System.out.println("Odd paths:\n^^^^^^^^^^\n");
		for(OddPath oddpath : oddpaths){
			for(int i = 0; i < oddpath.labels.length; i++){
				System.out.print(oddpath.labels[i]+" ");
			}
			System.out.println();
		}
		System.out.println("W paths:\n^^^^^^^^^^\n");
		for(WPath wpath : wpaths){
			for(int i = 0; i < wpath.labels.length; i++){
				System.out.print(wpath.labels[i]+" ");
			}
			System.out.println();
		}
		System.out.println("M paths:\n^^^^^^^^^^\n");
		for(MPath mpath : mpaths){
			for(int i = 0; i < mpath.labels.length; i++){
				System.out.print(mpath.labels[i]+" ");
			}
			System.out.println();
		}
		System.out.println("Circles:\n^^^^^^^^^^\n");
		for(Circle circle : circles){
			for(int i = 0; i < circle.labels.length; i++){
				System.out.print(circle.labels[i]+" ");
			}
			System.out.println();
		}

*/
	}
	
	public String print(){
		String s = "";
		s += "M paths:\n";
		for(MPath m : mpaths){
			s += m.print()+"\n";
		}
		s += "W paths:\n";
		for(WPath w : wpaths){
			s += w.print()+"\n";
		}
		s += "Odd paths:\n";
		for(OddPath o : oddpaths){
			s += o.print()+"\n";
		}
		s += "Circles:\n";
		for(Circle c : circles){
			s += c.print()+"\n";
		}
		return s;
	}

}
