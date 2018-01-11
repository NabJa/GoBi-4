package src.genomicUtils;

import augmentedTree.IntervalTree;

public class testMain {

	public static void main(String[] args) {

//		testIntervallTree();
//		testRVCut();
//		testMergeRV();
//		testContains();
//		testEquals();
		testGlue();
	}

	public static void testIntervallTree() {

		RegionVector testVec = new RegionVector();

		Region r21 = new Region(110, 130);
		Region r31 = new Region(111, 112);
		testVec.addRegion(r21);
		testVec.addRegion(r31);

		Region r1 = new Region(110, 130);
		Region r2 = new Region(110, 140);
		Region r3 = new Region(111, 112);
		RegionVector list = new RegionVector();
		list.addRegion(r1);
		list.addRegion(r2);
		list.addRegion(r3);

		IntervalTree<Region> testTree = new IntervalTree<Region>();
		testTree.addAll(testVec.regions);

		boolean t = testTree.containsAll(list.regions);

		System.out.println(t);

	}

	public static void testMergeRV() {

		RegionVector testVec = new RegionVector();
		RegionVector testVec1 = new RegionVector();

		Region r1 = new Region(5, 10);
		Region r2 = new Region(10, 15);
		Region r3 = new Region(20, 21);
		Region r21 = new Region(10, 17);
		Region r31 = new Region(20, 21);

		testVec.addRegion(r1);
		testVec.addRegion(r2);
		testVec.addRegion(r3);

		testVec1.addRegion(r21);
		testVec1.addRegion(r31);

		testVec = testVec.mergeRV(testVec1);

		System.out.println(testVec);
	}

	public static void testRVCut() {

		RegionVector testVec = new RegionVector();
		RegionVector testVec1 = new RegionVector();

		Region r1 = new Region(5, 10);
		Region r2 = new Region(12, 15);
		Region r3 = new Region(20, 21);

		Region r21 = new Region(13, 15);
		Region r31 = new Region(20,21);
		testVec1.addRegion(r21);
		testVec1.addRegion(r31);
		
		testVec.addRegion(r1);
		testVec.addRegion(r2);
		testVec.addRegion(r3);

		testVec = testVec.cutRV(1, 21);

		System.out.println(testVec);
		System.out.println(testVec.equals(testVec1));
	}

	public static void testContains() {
		
		RegionVector testVec = new RegionVector();
		RegionVector testVec1 = new RegionVector();

		Region r1 = new Region(5, 10);
		Region r2 = new Region(10, 15);
		Region r3 = new Region(20, 21);
		
		Region r21 = new Region(9, 15);
		Region r31 = new Region(20,21);
		
		testVec.addRegion(r1);
		testVec.addRegion(r2);
		testVec.addRegion(r3);

		testVec1.addRegion(r21);
		testVec1.addRegion(r31);

		System.out.println(testVec.contains(testVec1));
		
		
	}
	
	public static void testEquals() {
		
		RegionVector testVec = new RegionVector();
		RegionVector testVec1 = new RegionVector();

		Region r1 = new Region(9, 15);
//		Region r2 = new Region(10, 15);
		Region r3 = new Region(20, 21);
		
		Region r21 = new Region(9, 15);
		Region r31 = new Region(20,21);
		
		testVec.addRegion(r1);
//		testVec.addRegion(r2);
		testVec.addRegion(r3);

		testVec1.addRegion(r21);
		testVec1.addRegion(r31);

		System.out.println(testVec.isEqual(testVec1));
		
		
	}

public static void testGlue() {
		
		RegionVector testVec = new RegionVector();

		Region r1 = new Region(9, 15);
		Region r2 = new Region(16, 32);
		
		testVec.addRegion(r1);
		testVec.addRegion(r2);


		testVec.glue();
		
		System.out.println(testVec);
		
		
	}
	
}
