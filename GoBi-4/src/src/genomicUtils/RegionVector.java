package src.genomicUtils;

import java.util.ArrayList;
import java.util.TreeSet;

public class RegionVector {

	public String id; // transcript ID
	public int x1 = Integer.MAX_VALUE; // transcript start
	public int x2 = 0; // transcript end
	public int strandness = 0;

	public ArrayList<Region> regions = new ArrayList<Region>();
	public TreeSet<Region> uniqueRegions = new TreeSet<Region>();

	public RegionVector() {
	}

	public RegionVector(String id) {
		this.id = id;
	}

	public RegionVector(String id, int x1, int x2) {
		this.id = id;
		this.x1 = x1;
		this.x2 = x2;
	}

	public String getID() {
		return id;
	}

	public int getX1() {
		return x1;
	}

	public int getX2() {
		return x2;
	}

	public int getLength() {
		return x2 - x1;
	}

	public int getRegionLength() {
		int len = 0;
		for (Region r : regions) {
			len += r.getLength();
		}
		return len;
	}

	/**
	 * adds the Region in a sorted way. Smallest start first.
	 * 
	 * @param region
	 */
	public void addRegion(Region region) {
		if (!regions.contains(region)) {
			if (regions.size() == 0) {
				regions.add(region);
			} else {
				for (int i = 0; i < regions.size(); i++) {
					if (regions.get(i).getX1() < region.getX1())
						continue;
					regions.add(i, region);
					return;
				}
				regions.add(region);
			}

			if (region.getX1() < this.x1) {
				this.x1 = region.getX1();
			}
			if (region.getX2() > this.x2) {
				this.x2 = region.getX2();
			}
		}
	}

	public void addNewRegions(RegionVector rv) {
		for (Region r : rv.regions) {
			if (!regions.contains(r)) {
				regions.add(r);
			}
		}
	}

	public void addUniqueRegion(Region region) {

		this.uniqueRegions.add(region);

		if (region.getX1() < this.x1) {
			this.x1 = region.getX1();
		}
		if (region.getX2() > this.x2) {
			this.x2 = region.getX2();
		}
	}

	public int[] getTransLoc() {
		int[] region = { x1, x2 };
		return region;
	}

	public boolean contains(RegionVector rv) {

		int containedCount = 0;

		for (int i = 0; i < rv.regions.size(); i++) {
			for (int j = 0; j < regions.size(); j++) {
				if (rv.regions.get(i).getX1() >= regions.get(j).getX1()
						&& rv.regions.get(i).getX2() <= regions.get(j).getX2()) {
					containedCount++;
				}
			}
		}

		if (containedCount == rv.regions.size()) {
			return true;
		} else {
			return false;
		}
	}
	
	public boolean overlapsAll(RegionVector rv) {

		int overlapCount = 0;

		for (int i = 0; i < rv.regions.size(); i++) {
			for (int j = 0; j < regions.size(); j++) {
				if(rv.regions.get(i).overlaps(regions.get(j))) {
					overlapCount++;
					break;
				}
			}
		}

		if (overlapCount == rv.regions.size()) {
			return true;
		} else {
			return false;
		}
	}

	public boolean containsRegion(Region r) {
		for (Region thisR : this.regions) {
			if (thisR.contains(r)) {
				return true;
			}
		}
		return false;
	}

	public boolean intersects(RegionVector rv) {

		int containedCount = 0;

		for (int i = 0; i < rv.regions.size(); i++) {
			for (int j = 0; j < regions.size(); j++) {
				if (rv.regions.get(i).getX1() >= regions.get(j).getX1()
						&& rv.regions.get(i).getX2() <= regions.get(j).getX2()) {
					containedCount++;
				}
			}
		}

		if (containedCount == rv.regions.size()) {
			return true;
		} else {
			return false;
		}
	}

	public boolean isEqual(RegionVector pair) {
		if (pair.regions.size() == regions.size()) {
			for (int i = 0; i < pair.regions.size(); i++) {
				if (pair.regions.get(i).compareTo(regions.get(i)) != 0) {
					return false;
				}
			}
		} else {
			for (int i = 0; i < regions.size() - 1; i++) {
				if (regions.get(i).getX2() + 1 == regions.get(i + 1).getX1()) {
					regions.get(i).setRegion(regions.get(i).getX1(), regions.get(i + 1).getX2());
					regions.remove(i + 1);
					if (this.isEqual(pair)) {
						return true;
					}
				}
			}
			for (int i = 0; i < pair.regions.size() - 1; i++) {
				if (pair.regions.get(i).getX2() + 1 == pair.regions.get(i + 1).getX1()) {
					pair.regions.get(i).setRegion(pair.regions.get(i).getX1(), pair.regions.get(i + 1).getX2());
					pair.regions.remove(i + 1);
					if (this.isEqual(pair)) {
						return true;
					}
				}
			}
			return false;
		}
		return true;
	}

	public RegionVector glue() {
		if (regions.size() < 2) {
			return this;
		}
		for (int i = 1; i < regions.size(); i++) {
			if (regions.get(i).getX1() - 1 == regions.get(i - 1).getX2()) {
				Region r = new Region(regions.get(i - 1).getX1(), regions.get(i).getX2());
				regions.remove(i);
				regions.remove(i - 1);
				this.addRegion(r);
			}
		}
		return this;
	}

	public RegionVector cutRV(int a, int b) {
		RegionVector cutRV = new RegionVector();

		for (Region r : regions) {

			if (r.getX1() < a && r.getX2() > b) {
				Region nr = new Region(a, b);
				cutRV.addRegion(nr);

			} else if (r.getX1() < a) {
				Region nr = new Region(a, r.getX2());
				if (nr.getX1() <= nr.getX2()) {
					cutRV.addRegion(nr);
				}

			} else if (r.getX2() > b) {
				Region nr = new Region(r.getX1(), b);
				if (nr.getX1() <= nr.getX2()) {
					cutRV.addRegion(nr);
				}
				break;
			} else {
				cutRV.addRegion(r);
			}
		}
		return cutRV;
	}

	public RegionVector mergeRV(RegionVector rv) {

		for (Region r : regions) {
			rv.addRegion(r);
		}

		RegionVector mergedRV = new RegionVector();

		if (rv.regions.size() == 0)
			return mergedRV;

		Region pre = rv.regions.get(0);

		for (int i = 0; i < rv.regions.size(); i++) {
			Region curr = rv.regions.get(i);
			if (curr.getX1() > pre.getX2()) {
				mergedRV.addRegion(pre);
				pre = curr;
			} else {
				Region merged = new Region(pre.getX1(), Math.max(pre.getX2(), curr.getX2()));
				pre = merged;
			}
		}
		mergedRV.addRegion(pre);

		return mergedRV;
	}

	/**
	 * Returns an ArrayList<Region> that are the inverse of given RegionVector.
	 * Example: Input{1,2; 3,4; 5,6} -> Output{2,3; 4,5}
	 * 
	 * @return
	 */
	public ArrayList<Region> arrayInverse() {

		ArrayList<Region> introns = new ArrayList<Region>();

		for (int i = 0; i < regions.size() - 1; i++) {
			Region intron = new Region(regions.get(i).getX2(), regions.get(i + 1).getX1());
			introns.add(intron);
		}

		return introns;
	}

	/**
	 * Returns an RegionVector that are the inverse of given RegionVector. Has to be
	 * save in new Variable!!! Example: Input{1,2; 3,4; 5,6} -> Output{2,3; 4,5}
	 * 
	 * @return
	 */
	public RegionVector inverse() {

		RegionVector introns = new RegionVector();

		for (int i = 0; i < regions.size() - 1; i++) {
			Region intron = new Region(regions.get(i).getX2(), regions.get(i + 1).getX1());
			if (intron.getX1() < intron.getX2()) {
				introns.addRegion(intron);
			}
		}

		introns.id = id;

		return introns;
	}

	public TreeSet<Region> inverseToSet() {

		TreeSet<Region> introns = new TreeSet<Region>();

		for (int i = 0; i < regions.size() - 1; i++) {
			Region intron = new Region(regions.get(i).getX2() + 1, regions.get(i + 1).getX1());
			introns.add(intron);
		}
		return introns;

	}

	public void printRegions() {
		System.out.println("ID: " + id);
		for (Region r : regions) {
			System.out.println(r.getX1() + " " + r.getX2() + " ");
		}
		System.out.println();
	}

	public int getSize() {
		int i = 0;
		for (Region r : regions) {
			i++;
		}
		return i;
	}

	@Override
	public String toString() {
		String result = id;
		for (int i = 0; i < regions.size(); i++) {
			result += " " + regions.get(i);
		}
		return result;
	}

	@Override
	public int hashCode() {
		int result = regions.hashCode();
		result = 101 * result * (strandness + 7);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof RegionVector))
			return false;
		if (obj == this)
			return true;

		RegionVector rv = (RegionVector) obj;

		if (!rv.regions.containsAll(this.regions)) {
			return false;
		}

		if (this.x1 == rv.getX1() && this.x2 == rv.getX2() && this.strandness == rv.strandness) {
			return true;
		} else {
			return false;
		}
	}
}