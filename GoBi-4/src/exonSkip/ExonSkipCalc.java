package exonSkip;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import src.genomicUtils.Gene;
import src.genomicUtils.Region;
import src.genomicUtils.RegionVector;

public class ExonSkipCalc {

	public void getExonFromGen(Gene gene, OutputMap outMap) {

		for (RegionVector rv : gene.transcripts.values()) {
			for (Region intron : rv.inverse().regions) {
				if (!outMap.isInResultMap(intron)) {
					getResults(intron, gene, outMap, rv);
				}
			}
		}
	}

	public void getResults(Region intron, Gene gene, OutputMap outMap, RegionVector transcript) {
		
		Output output = new Output();
		HashSet<String> wtsv = getOverlappingIntrons(intron, gene.transcripts);
		HashSet<String> wt = removeSV(intron, wtsv, gene);
		RegionVector finalskipped = new RegionVector();

		for (String id : wt) {
			RegionVector rv = gene.transcripts.get(id);
			RegionVector skipped = subtract(intron, rv);
			finalskipped.addNewRegions(skipped);
			HashSet<String> sv_prots = getSVs(intron, gene.transcripts, output, transcript);
			HashSet<String> wt_prots = getWTs(intron, id, gene.transcripts, output);
		}
				
		if (finalskipped.regions.size() > 1) {
			output.setOutput(gene, intron, finalskipped);
			outMap.resultMap.put(intron, output);
		}
		
	}

	public RegionVector subtract(Region intron, RegionVector rv) {

		RegionVector introns = new RegionVector();

		int i = 0;

		if (rv.regions.size() > 0) {
			
			while (rv.regions.get(i).getX1() != intron.getX1() && i < rv.regions.size() - 1) {
				i++;
			}

			// if (regionInRV(intron, rv) == false) {}
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
		return introns;
	}

	public HashSet<String> getSVs(Region intron, HashMap<String, RegionVector> transcripts, Output out,
			RegionVector splivedV) {

		HashSet<String> svs = new HashSet<String>();
		HashSet<String> svProts = new HashSet<String>();

		for (String id : transcripts.keySet()) { // For every transcript ID

			RegionVector rv = transcripts.get(id); // Get transcript values

			for (Region r : rv.inverse().regions) // For every Intron in transcript:
			{
				if (r.getX1() == intron.getX1() && r.getX2() == intron.getX2()) // If intron == intron of any other
																				// transcript
				{
					svs.add(rv.id);
					svProts.add(rv.regions.get(0).regionID);
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

		int skippedExonsInTranscript = 1;
		for (Region r : rv.regions) // For every Region of transcript
		{
			if (intron.getX1() < r.getX1() && intron.getX2() > r.getX2()) // if region is in intron:
			{
				wts.add(r.regionID);

				out.insertMaxSkippedExons(skippedExonsInTranscript);
				out.insertMinSkippedExons(skippedExonsInTranscript);
				skippedExonsInTranscript++;

				skippedExons.add(r.regionID);
				
				out.insertMinSkippedBases(r.getLength());
				out.insertMaxSkippedBases(r.getLength());
				out.wt_prots.addAll(skippedExons);
			}
		}

		return wts;
	}

	public HashSet<String> removeSV(Region intron, HashSet<String> sameIntrons, Gene gene) {
		HashSet<String> wildtypes = new HashSet<String>();

		System.out.println("SVs: " + sameIntrons.size());

		for (String id : sameIntrons) {
			RegionVector rv = gene.transcripts.get(id).inverse();
			for (Region r : rv.regions) {
				if (r.getX1() != intron.getX1() && r.getX2() != intron.getX2()) {
					wildtypes.add(id);
				}
			}
		}

		System.out.println("WTs: " + wildtypes.size());
		System.out.println();

		return wildtypes;
	}

	public HashSet<String> getOverlappingIntrons(Region intron, HashMap<String, RegionVector> transcripts) {

		int start = intron.getX1();
		int end = intron.getX2();

		HashSet<String> sameStarts = new HashSet<String>();
		HashSet<String> sameEnds = new HashSet<String>();
		HashSet<String> sameIntrons = new HashSet<String>();

		sameStarts = getStarts(start, transcripts);
		sameEnds = getEnds(end, transcripts);

		sameIntrons = getSameIntron(sameStarts, sameEnds);

		return sameIntrons;
	}

	public HashSet<String> getSameIntron(HashSet<String> setA, HashSet<String> setB) {
		setA.retainAll(setB);
		return setA;
	}

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

	public RegionVector HashSetToRegionVector(HashSet<Region> set) {
		RegionVector sv = new RegionVector();
		for(Region r : set) {
			sv.addRegion(r);
		}
		return sv;
	}
	
}
