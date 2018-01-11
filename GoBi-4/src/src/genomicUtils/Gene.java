package src.genomicUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import augmentedTree.*;

public class Gene implements Interval {

	public String geneID;
	public int start = Integer.MAX_VALUE;
	public int end = 0;
	public String strand;
	public String source;
	public String type;
	public String geneChr;
	public String geneName;

	public HashMap<String, RegionVector> transcripts = new HashMap<String, RegionVector>();
	
	public Collection<String> wtStarts = new ArrayList<String>();
	public Set<Integer> wtEnds = new HashSet<Integer>();
	public Set<Integer> wts = new HashSet<Integer>();

	public Gene() {
	}

	public Gene(int start, int end) {
		this.start = start;
		this.end = end;
	}

	public void updatePos(int start, int end) {
		if (start < this.start) {
			this.start = start;
		}
		if (end > this.end) {
			this.end = end;
		}
	}

	public void setGene(String id, String geneChr, String strand, int start, int end) {
		this.geneID = id;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.geneChr = geneChr;
		this.transcripts = new HashMap<String, RegionVector>();
	}

	public void setGene(String id, String geneChr, String strand, int start, int end, String biotype) {
		this.geneID = id;

		if (start < this.start) {
			this.start = start;
		}
		if (end > this.end) {
			this.end = end;
		}
		this.strand = strand;
		this.geneChr = geneChr;
		this.type = biotype;
		this.transcripts = new HashMap<String, RegionVector>();
	}

	public void setGene(String id, String geneChr, String geneName, String strand) {
		this.geneID = id;
		this.geneChr = geneChr;
		this.geneName = geneName;
		this.strand = strand;
		this.transcripts = new HashMap<String, RegionVector>();
	}

	public void setGene(String id, int start, int end, String strand, String geneChr, String source, String type,
			String geneName) {
		this.geneID = id;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.source = source;
		this.type = type;
		this.geneChr = geneChr;
		this.geneName = geneName;
		this.transcripts = new HashMap<String, RegionVector>();
	}

	/**
	 * inserts RegionVector in HashMap transcripts
	 * 
	 * @param rv
	 */
	public void insertRV(RegionVector rv) {
		transcripts.put(rv.getID(), rv);
	}

	public String getID() {
		return geneID;
	}

	public int nTrans() {
		int i = 0;
		for (RegionVector r : transcripts.values()) {
			i += r.getSize();
		}
		return i;
	}

	public void getStarts2() { // Collection<Integer>

		transcripts.forEach((k, v) -> {
			System.out.println(v.getX1() + " " + v.getX2() + "\t" + "\t" + k);
		});
	}

	public void printWtStarts() {
		for (String e : wtStarts) {
			System.out.println(e);
		}
	}

	public void printTranscriptsInverse() {
		transcripts.forEach((k, v) -> {
			v.inverse().printRegions();
		});

	}

	public String getChr() {
		return geneChr;
	}

	public boolean getGeneNegativeStrandFlag() {
		if (this.strand.equals("-")) {
			return true;
		} else {
			return false;
		}

	}

	@Override
	public int getStart() {
		return this.start;
	}

	@Override
	public int getStop() {
		return this.end;
	}

	@Override
	public String toString() {
		return "" + geneID + " " + start + " " + end;
	}

}