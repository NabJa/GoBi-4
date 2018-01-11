package src.genomicUtils;

public class Tuple<X, Y> {
	public final X first;
	public final Y second;

	
	public Tuple(X first, Y second) {
		this.first = first;
		this.second = second;
	}
	
//	public void setTuple(X first, Y second) {
//		this.first = first;
//		this.second = second;
//	}
	
	public X getFirst() {
		return first;
	}

	public Y getSecond() {
		return second;
	}
	
	@Override
	public String toString() {
		return "" + first + " " + second;
	}
}