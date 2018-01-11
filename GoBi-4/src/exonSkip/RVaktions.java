package exonSkip;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import src.genomicUtils.*;

public class RVaktions {

	public void printRV(RegionVector rv) {
		for (Region r : rv.regions) {
			System.out.println(r.getX1() + " " + r.getX2() + " ");
		}
	}

	/**
	 * Takes 2 RegionVectors and returns a set with all overlapping starts
	 * 
	 * @param o1
	 * @param o2
	 * @return
	 */
	public Set<Integer> compareStarts(RegionVector o1, RegionVector o2) {

		Set<Integer> sameStarts = new TreeSet<Integer>();
		for (int i = 0; i < o1.regions.size(); i++) {
			if (o1.regions.get(i).getX1() == o2.regions.get(i).getX1()) { // o1.inverse().regions.get(i).getX1() ==
																			// o2.inverse().regions.get(i).getX1()
				sameStarts.add(i);
			}
		}
		return sameStarts;
	}

	/**
	 * Takes 2 RegionVectors and returns a set with all overlapping ends
	 * 
	 * @param o1
	 * @param o2
	 * @return
	 */
	public Set<Integer> compareEnds(RegionVector o1, RegionVector o2) {
		Set<Integer> sameEnds = new TreeSet<Integer>();
		for (int i = 0; i < o1.regions.size(); i++) {
			if (o1.regions.get(i).getX2() == o2.regions.get(i).getX2()) { // o1.inverse().regions.get(i).getX2() ==
																			// o2.inverse().regions.get(i).getX2())
				sameEnds.add(i);
			}
		}
		return sameEnds;
	}

	/**
	 * If position == 0 it returns true if RV contains same start as input intron.
	 * If position != 0 it returns true if RV contains same end as input intron.
	 * 
	 * @param intron
	 * @param rv
	 * @param position
	 * @return
	 */
	public boolean containsSamePosition(Region intron, RegionVector rv, int position) {
		if (position == 0) {
			int intronStart = intron.getX1();
			for (int i = 0; i < rv.regions.size(); i++) {
				Region RVintrons = rv.regions.get(i);
				int RVstarts = RVintrons.getX1();
				if (intronStart == RVstarts) {
					return true;
				}
			}
			return false;
		} else {
			int intronEnd = intron.getX2();
			for (int i = 0; i < rv.regions.size(); i++) {
				Region RVintrons = rv.regions.get(i);
				int RVends = RVintrons.getX2();
				if (intronEnd == RVends) {
					return true;
				}
			}
			return false;
		}
	}

	/**
	 * It will return all IDs of transcripts with same start as input start
	 * 
	 * @param start
	 * @param transcripts
	 * @return
	 */
	public HashSet<String> getStarts(int start, HashMap<String, RegionVector> transcripts) {

		HashSet<String> sameStarts = new HashSet<String>();

		transcripts.keySet().forEach((id) -> {
			RegionVector rv = transcripts.get(id);

			for (Region intron : rv.inverse().regions) {
				if (intron.getX1() == start) {
					sameStarts.add(rv.id);
				}
			}
		});

		return sameStarts;
	}

	/**
	 * Returns all IDs of transcripts with same ends as input start
	 * 
	 * @param start
	 * @param transcripts
	 * @return
	 */
	public HashSet<String> getEnds(int end, HashMap<String, RegionVector> transcripts) {

		HashSet<String> sameEnds = new HashSet<String>();

		transcripts.keySet().forEach((id) -> {
			RegionVector rv = transcripts.get(id);
			for (Region intron : rv.inverse().regions) {
				if (intron.getX2() == end) {
					sameEnds.add(rv.id);
				}
			}
		});
		return sameEnds;
	}

	/**
	 * Retains 2 sets. You get all elements that are in BOTH sets.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public HashSet<String> getSameIntron(HashSet<String> setA, HashSet<String> setB) {
		setA.retainAll(setB);
		return setA;
	}

	public RegionVector subtract(Region intron, RegionVector rv, OutputMap outMap, Output output, RegionVector rv1,
			HashSet<String> sameIntrons) {

		RegionVector introns = new RegionVector();

		int i = 0;

		if (rv.regions.size() > 0) {
			while (rv.regions.get(i).getX1() != intron.getX1() && i < rv.regions.size() - 1) {
				i++;
			}

			if (regionInRV(intron, rv) == false) {

				while (rv.regions.get(i).getX2() < intron.getX2()) {
					int start = rv.regions.get(i).getX1();
					int end = rv.regions.get(i).getX2();
					Region skippedIntron = new Region(start, end);
					introns.addRegion(skippedIntron);
					i++;
				}
				int start = rv.regions.get(i).getX1();
				int end = rv.regions.get(i).getX2();
				Region skippedIntron = new Region(start, end);
				introns.addRegion(skippedIntron);
			}
		}
		return introns;
	}

	public boolean regionInRV(Region r, RegionVector rv) {
		boolean bol = false;
		for (Region a : rv.regions) {
			if (a.getX1() == r.getX1() && a.getX2() == r.getX2()) {
				bol = true;
			}
		}
		return bol;
	}

	/**
	 * Returns a HashSet with strings of the ids of RVs with same start AND same end
	 * in any intron combination
	 * 
	 * @param start
	 * @param end
	 * @param transcripts
	 * @return
	 */
	public HashSet<String> getOverlappingIntrons(int start, int end, HashMap<String, RegionVector> transcripts) {

		HashSet<String> sameStarts = new HashSet<String>();
		HashSet<String> sameEnds = new HashSet<String>();
		HashSet<String> sameIntrons = new HashSet<String>();

		sameStarts = getStarts(start, transcripts);
		sameEnds = getEnds(end, transcripts);

		sameIntrons = getSameIntron(sameStarts, sameEnds);

		return sameIntrons;
	}

	public HashSet<String> getSVs(Region intron, HashMap<String, RegionVector> transcripts, Output out,
			RegionVector splivedV) {

		HashSet<String> svs = new HashSet<String>();
		HashSet<String> svProts = new HashSet<String>();

		for (String id : transcripts.keySet()) { // For every transcript ID

			RegionVector rv = transcripts.get(id); // Get transcript values

			RegionVector inverseRV = rv.inverse();
			
			for (Region r : inverseRV.regions) // For every Intron in transcript:
			{
				if (r.getX1() == intron.getX1() && r.getX2() == intron.getX2()) // If intron == intron of any other
																				// transcript
				{
					svs.add(rv.id);
					svProts.add(rv.regions.get(0).regionID);
					out.exonToSV.put(r, rv);
				}
			}
		}
		out.sv_prots = svProts;
		return svs;
	}

	public HashSet<String> getWTs(Region intron, String id, HashMap<String, RegionVector> transcripts, Output out) {

		HashSet<String> wts = new HashSet<String>();
		Set<String> skippedExons = new HashSet<String>();
		RegionVector rv = transcripts.get(id); // Get transcript values

		for (Region r : rv.regions) // For every Region of transcript
		{
			if (intron.getX1() < r.getX1() && intron.getX2() > r.getX2()) // if region is in intron:
			{
				wts.add(r.regionID);
				skippedExons.add(r.regionID);
				out.wt_prots.addAll(skippedExons);
				out.exonToWT.put(r, rv);
			}
		}
		
		return wts;
	}

	public HashSet<String> removeSV(Region intron, HashSet<String> sameIntrons, Gene gene) {
		HashSet<String> wildtypes = new HashSet<String>();
		for (String id : sameIntrons) {
			RegionVector rv = gene.transcripts.get(id).inverse();
			for (Region r : rv.regions) {
				if (r.getX1() != intron.getX1() || r.getX2() != intron.getX2()) {
					wildtypes.add(id);
				}
			}
		}
		return wildtypes;
	}

	/**
	 * 
	 * @param intron
	 * @param gene
	 * @param transID
	 * @return
	 */
	public void getSkippedExon(Region intron, Gene gene, OutputMap outMap, RegionVector rv) {

		intron.geneID = gene.geneID;
		
		Output output = new Output();
		HashSet<String> sameIntrons = new HashSet<String>();

		sameIntrons = getOverlappingIntrons(intron.getX1(), intron.getX2(), gene.transcripts);
		// sameIntrons.remove(rv.id);

		HashSet<String> wildtypes = removeSV(intron, sameIntrons, gene);

		RegionVector skippedIntrons = new RegionVector();
		RegionVector skipped = new RegionVector();

		if (!outMap.isInResultMap(intron)) {
			for (String id : wildtypes) {

				HashSet<String> sv_prots = getSVs(intron, gene.transcripts, output, rv);
				HashSet<String> wt_prots = getWTs(intron, id, gene.transcripts, output);
				
//				if(wt_prots.size() > 0 && sv_prots.size() > 0 && gene.transcripts.get(sv_prots.toArray()[0]).getSize() > 0 && gene.transcripts.get(wt_prots.toArray()[0]).getSize() > 0) {
//					System.out.println("sv" + " " + gene.transcripts.get(sv_prots.toArray()[0]));
//					System.out.println("wt" + " " + gene.transcripts.get(wt_prots.toArray()[0]));
//					
//					output.svTrans.addNewRegions(gene.transcripts.get(sv_prots.toArray()[0]));
//					output.wtTrans.addNewRegions(gene.transcripts.get(wt_prots.toArray()[0]));
//				}

				RegionVector queryRV1 = gene.transcripts.get(id);

				RegionVector queryRV = queryRV1.inverse();

				skippedIntrons = subtract(intron, queryRV, outMap, output, queryRV1, wildtypes);
				skipped.addNewRegions(skippedIntrons);
				
				
				if (skipped.regions.size() > 1) {
					output.setOutput(gene, intron, skipped);
					outMap.resultMap.put(intron, output);
				}
			}
		}

	}

	public void getSkippedExonFromGen(Gene gene, OutputMap outMap) {
		for (RegionVector rv : gene.transcripts.values()) {
			RegionVector inverseRV = rv.inverse();
			for (Region r : inverseRV.regions) {
				getSkippedExon(r, gene, outMap, rv);
			}
		}
	}

}
