package src.genomicUtils;

public class Triplet<X, Y, Z> {

	private final X first;
	private final Y second;
	private final Z third;

	public Triplet(X first, Y second, Z third) {
		this.first = first;
		this.second = second;
		this.third = third;
	}

	public X getFirst() {
		return first;
	}

	public Y getSecond() {
		return second;
	}

	public Z getThird() {
		return third;
	}

	@Override
	public String toString() {
		return "" + first + " " + second + " " + third;
	}

}