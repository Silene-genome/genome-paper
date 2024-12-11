package representation;

public class MPath {
public int[] labels;
	
	public MPath(int n){
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
