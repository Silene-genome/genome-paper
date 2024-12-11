package representation;

public class AdjacencyGraphNode{
	public int[] labels;
	int[] neighbours;
	
	public AdjacencyGraphNode() {
		neighbours = new int[2];
		labels = new int[2];
	}
}