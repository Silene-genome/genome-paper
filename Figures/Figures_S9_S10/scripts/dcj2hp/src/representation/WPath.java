package representation;

public class WPath {
public int[] labels;
	
	public WPath(int n){
		labels = new int[n];
	}

	public String print(){
		String s = "";
		for(int i = 0; i < labels.length; i++){
			s += labels[i] + " ";
		}
		
		return s;
	}
	
}