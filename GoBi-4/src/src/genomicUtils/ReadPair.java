package src.genomicUtils;

import java.util.ArrayList;

import augmentedTree.IntervalTree;
import src.feature_extraction.GeneAnnotation;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

public class ReadPair {

	public SAMRecord first;
	public SAMRecord second;
	public String readName;
	public boolean firstIsNeg;
	public int start;
	public int end;
	public int gcount = 0;
	public boolean intronic = false;
	public boolean mergedTrans = false;
	public boolean trans = false;
	public boolean antisense = false;
	public RegionVector genomicRegions = new RegionVector();
	public RegionVector firstRegions = new RegionVector();
	public RegionVector secondRegions = new RegionVector();
	public IntervalTree<Region> regionTree = new IntervalTree<Region>();
	public GeneAnnotation geneAnno = new GeneAnnotation();
	public int gdist = Integer.MAX_VALUE;

	public ReadPair(SAMRecord sr1, SAMRecord sr2) {
		if (sr1.getFirstOfPairFlag()) {
			this.first = sr1;
			this.second = sr2;
			this.firstIsNeg = sr1.getReadNegativeStrandFlag();
		} else {
			this.first = sr2;
			this.second = sr1;
			this.firstIsNeg = sr2.getReadNegativeStrandFlag();
		}
		this.readName = sr1.getReadName();
		this.start = getMin(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());
		this.end = getMax(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());
	}

	private int getMin(int v, int x, int y, int z) {
		int minVX = Math.min(v, x);
		int minYZ = Math.min(y, z);
		int gloabalMin = Math.min(minVX, minYZ);
		return gloabalMin;
	}

	private int getMax(int v, int x, int y, int z) {
		int maxVX = Math.max(v, x);
		int maxYZ = Math.max(y, z);
		int globalMax = Math.max(maxVX, maxYZ);
		return globalMax;
	}

	public int getFirstStart() {
		return firstRegions.regions.get(0).getX1();
	}

	public int getFirstEnd() {
		return firstRegions.regions.get(firstRegions.regions.size() - 1).getX2();
	}

	public int getSecondStart() {
		return secondRegions.regions.get(0).getX1();
	}

	public int getSecondEnd() {
		return secondRegions.regions.get(secondRegions.regions.size() - 1).getX2();
	}

	public void calcGenomicRegions(Integer frstrand) {
		RegionVector genRegions = new RegionVector();
		// RegionVector frRegions = new RegionVector();
		// RegionVector scRegions = new RegionVector();

		for (AlignmentBlock ab : first.getAlignmentBlocks()) {
			Region abRegion = new Region(ab.getReferenceStart(), (ab.getReferenceStart() + ab.getLength() - 1));
			genRegions.addRegion(abRegion);
			this.firstRegions.addRegion(abRegion);
		}
		for (AlignmentBlock ab : second.getAlignmentBlocks()) {
			Region abRegion = new Region(ab.getReferenceStart(), (ab.getReferenceStart() + ab.getLength() - 1));
			genRegions.addRegion(abRegion);
			this.secondRegions.addRegion(abRegion);
		}

		genRegions = genRegions.mergeRV(genomicRegions);
		// firstRegions = firstRegions.mergeRV(frRegions);
		// secondRegions = secondRegions.mergeRV(scRegions);

		genRegions.id = first.getReadName();
		if (firstIsNeg && frstrand != 0) {
			genRegions.strandness = -1;
		} else {
			genRegions.strandness = 1;
		}

		this.genomicRegions = genRegions;

		firstRegions.glue();
		secondRegions.glue();
		genomicRegions.glue();
		genomicRegions.glue();

	}

	public RegionVector getMergedRegions(SAMRecord sr1, SAMRecord sr2) {
		RegionVector genRegions = new RegionVector();
		RegionVector genRegions1 = new RegionVector();

		for (AlignmentBlock ab : sr1.getAlignmentBlocks()) {
			Region abRegion = new Region(ab.getReadStart(), (ab.getReadStart() + ab.getLength()));
			genRegions.addRegion(abRegion);
		}
		for (AlignmentBlock ab : sr2.getAlignmentBlocks()) {
			Region abRegion = new Region(ab.getReadStart(), (ab.getReadStart() + ab.getLength()));
			genRegions1.addRegion(abRegion);
		}
		genRegions = genRegions.mergeRV(genRegions1);
		return genRegions;
	}

	public void getGeneDistance(ArrayList<Gene> neighbour) {
		int dist;
		for (int i = 0; i < neighbour.size(); i++) {
			if (firstRegions.regions.get(0).getX2() < neighbour.get(i).start
					&& secondRegions.regions.get(secondRegions.regions.size() - 1).getX1() > neighbour.get(i).end) {
				int firstDist = neighbour.get(i).start - firstRegions.regions.get(i).getX2();
				int secDist = secondRegions.regions.get(secondRegions.regions.size()).getX1() - neighbour.get(i).end;
				dist = Math.min(firstDist, secDist);
				if (dist < this.gdist) {
					this.gdist = dist;
				}
			} else if (neighbour.get(i).end < start) {
				dist = start - neighbour.get(i).end - 1;
				if (dist < this.gdist) {
					this.gdist = dist;
				}
			} else if (neighbour.get(i).start > end) {
				dist = neighbour.get(i).start - end - 1;
				if (dist < this.gdist) {
					this.gdist = dist;
				}
			} else {
				this.gdist = 0;
			}
		}
	}

	@Override
	public String toString() {
		return "" + genomicRegions.toString();

	}

}
