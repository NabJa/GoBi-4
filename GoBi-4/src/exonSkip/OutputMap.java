package exonSkip;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import src.genomicUtils.Region;

public class OutputMap {

	FileWriter file;
	BufferedWriter writer;
	String outputDestination;
	public HashMap<Region, Output> resultMap = new HashMap<Region, Output>();
	public HashMap<Region, Output> resultMapExons = new HashMap<Region, Output>();
	public ArrayList<Region> skippedExons = new ArrayList<Region>();
	public HashMap<String, ArrayList<Region>> geneToExon = new HashMap<String, ArrayList<Region>>();

	ArrayList<Output> resultList = new ArrayList<Output>();

	public OutputMap() {

	} 

	public OutputMap(String outputDestination) {
		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination);
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while writing output from OutputMap.", e);
		}
	}

	public boolean isInResultMap(Region r) {
		if (resultMap.containsKey(r)) {
			return true;
		} else {
			return false;
		}
	}

	public boolean isNotInResults(Region r) {
		if (resultMap.containsKey(r)) {
			return true;
		} else {
			return false;
		}

	}

	public void addToResultIfNew(Region region, Output output) {
		if (!resultMap.containsKey(region)) {
			addToResultMap(region, output);
		}
	}

	public void addToResultMap(Region r, Output out) {
		resultMap.put(r, out);
	}

	public void printOutput() {
		try {
			writer.write("id" + "\t" + "symbol" + "\t" + "chr" + "\t" + "strand" + "\t" + "nprots" + "\t" + "ntrans"
					+ "\t" + "SV" + "\t" + "WT" + "\t" + "WT_prots" + "\t" + "SV_prots" + "\t" + "min_skipped_exon"
					+ "\t" + "max_skipped_exon" + "\t" + "min_skipped_bases" + "\t" + "max_skipped_bases");
			writer.newLine();

		} catch (Exception e) {
			throw new RuntimeException("got error while printing Headline!", e);
		}

		for (Region r : resultMap.keySet()) {
			resultMap.get(r).printOutTxt(file, writer);
		}

		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("got error while closing writer!", e);
		}
	}

	public void printOutput1() {
		try {
			writer.write("id" + "\t" + "symbol" + "\t" + "chr" + "\t" + "strand" + "\t" + "nprots" + "\t" + "ntrans"
					+ "\t" + "SV" + "\t" + "WT" + "\t" + "WT_prots" + "\t" + "SV_prots" + "\t" + "min_skipped_exon"
					+ "\t" + "max_skipped_exon" + "\t" + "min_skipped_bases" + "\t" + "max_skipped_bases");
			writer.newLine();

		} catch (Exception e) {
			throw new RuntimeException("got error while printing Headline!", e);
		}

		for (Output out : resultMap.values()) {
			out.printOutTxt(file, writer);
		}
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("got error while printing Output!", e);
		}
	}

	public void printExons() {
		try {
			writer.write("id" + "\t" + "SV" + "\t" + "WT");
			writer.newLine();

		} catch (Exception e) {
			throw new RuntimeException("got error while printing Headline!", e);
		}

		for (Output out : resultMap.values()) {
			for (Region wt : out.wt.inverse().regions) {
				out.printSkippedExons(file, writer, wt);
				skippedExons.add(wt);
				wt.geneID = out.geneID;
				resultMapExons.put(wt, out);

				ArrayList<Region> exonOfGene = geneToExon.get(out.geneID);
				if (exonOfGene == null) {
					ArrayList<Region> newExonOfGene = new ArrayList<Region>();
					newExonOfGene.add(wt);
					geneToExon.put(out.geneID, newExonOfGene);

				} else {
					exonOfGene.add(wt);
					geneToExon.put(out.geneID, exonOfGene);
				}

			}
		}
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("got error while printing Output!", e);
		}
	}

}
