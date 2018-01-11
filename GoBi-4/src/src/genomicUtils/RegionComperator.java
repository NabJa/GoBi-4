package src.genomicUtils;

public class RegionComperator implements Comperator<Region> {

	@Override
	public int compare(Region r1, Region r2) {
		if(r1.getX1() == r2.getX1() && r1.getX2() == r2.getX2()) {
			return 0;
		} else {
			return 1;			
		}
	}
}
