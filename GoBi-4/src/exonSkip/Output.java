package exonSkip;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import src.genomicUtils.*;

public class Output {

	String geneID; //
	String geneName; //
	String chr; // chromosome
	String strand; //
	int nprots; // number of annotated CDS in the gene
	int ntrans; // number of annotated transcripts in the gene
	Region sv; // the SV intron as start:end
	RegionVector wt; // the WT introns within the SV intron seperated by | as start:end
	
	HashSet<String> sv_prots = new HashSet<String>(); // IDs of the SV CDSs, seperated by |
	HashSet<String> wt_prots = new HashSet<String>(); // IDs of the WT CDSs, seperated by |
	
	public HashMap<Region, RegionVector> exonToWT = new HashMap<Region, RegionVector>();
	public HashMap<Region, RegionVector> exonToSV = new HashMap<Region, RegionVector>();
	
	public RegionVector svTrans = new RegionVector();
	public RegionVector wtTrans = new RegionVector();
	
	int minSkippedExons = Integer.MAX_VALUE -1;
	int maxSkippedExons = 0; 
	int minSkippedBases = Integer.MAX_VALUE - 1;
	int maxSkippedBases = 0;
    
	public void addSV_prots(String id) {
		sv_prots.add(id);
	}

	public void getAllSVProtIDs(RegionVector rv) {
		for(Region r : rv.regions) {
			sv_prots.add(r.getID());
		}
	}

	public void addWT_prots(String id) {
		wt_prots.add(id);
	}

	public void getAllWTProtIDs(RegionVector rv) {
		for (Region r : rv.regions) {
			if (!wt_prots.contains(r.regionID)) {
				wt_prots.add(r.regionID);
			}
		}
	}

	public void addWT(Region r) {
		wt.addRegion(r);
	}

	public void insertMinSkippedExons(int mse) {
		if (mse < minSkippedExons) {
			this.minSkippedExons = mse;
			}
	}

	public void insertMaxSkippedExons(int mxse) {
		if (mxse > maxSkippedExons) {
			this.maxSkippedExons = mxse;
		}
	}

	public void insertMinSkippedBases(int msb) {
		if (msb < minSkippedBases) {
			this.minSkippedBases = msb;
		}
	}

	public void insertMaxSkippedBases(int mxsb) {
		if (mxsb > maxSkippedBases) {
			this.maxSkippedBases = mxsb;
		}
	}
	
	public void printSkippedExons(FileWriter file, BufferedWriter writer, Region r) {
		try {
			
			writer.write(geneID + "\t" + sv.getX1() + ":" + sv.getX2() + "\t" + r.getX1() + ":" + r.getX2());
			writer.newLine();
			writer.write("WT: " + exonToWT.get(r).toString());
			writer.newLine();
			writer.write("SV: " + exonToSV.values().toString());
			writer.newLine();
			writer.newLine();
			
		} catch (IOException e) {

			e.printStackTrace();
		}
	}
	
	public void printExon(FileWriter file, BufferedWriter writer){
		
		try {
			
			writer.write(geneID + "\t" + sv.getX1() + ":" + sv.getX2() + "\t");

			String str = "";

			for (Region r : wt.regions) {
				str += ("|" + r.getX1()) + ":" + r.getX2();
			}
			if (str.length() > 1) { // maybe could be removed. only important if RegionVector is empty. Shouldnt be
									// the case
				str = str.substring(1, str.length());
				writer.write(str + "\t");
			}
			writer.newLine();
			
		} catch (Exception e) {
			throw new RuntimeException("got error while printing exon output.", e);
		}
		
	}

	/**
	 * prints a tsv-file of everything in object Output
	 * 
	 * @throws Exception
	 */
	public void printOutTxt(FileWriter file, BufferedWriter writer) {

		try {

			writer.write(geneID + "\t" + geneName + "\t" + chr + "\t" + strand + "\t" + nprots + "\t" + ntrans + "\t");

			writer.write(sv.getX1() + ":" + sv.getX2() + "\t");

			String str = "";
			for (Region r : wt.regions) {
				str += ("|" + Integer.toString(r.getX1()) + ":" + Integer.toString(r.getX2()));
			}
			if (str.length() > 1) { // maybe could be removed. only important if RegionVector is empty. Shouldnt be
									// the case
				str = str.substring(1, str.length());
				writer.write(str + "\t");
			}

			String str1 = "";
			for (String id : wt_prots) {
				str1 += ("|" + id);
			}
			str1 = str1.substring(1, str1.length());
			writer.write(str1 + "\t");

			String str2 = "";
			for (String id : sv_prots) {
				str2 += ("|" + id);
			}
			str2 = str2.substring(1, str2.length());
			writer.write(str2 + "\t");
			
			writer.write(minSkippedExons + "\t" + maxSkippedExons + "\t" + minSkippedBases + "\t" + maxSkippedBases);
			writer.newLine();

		} catch (Exception e) {
			throw new RuntimeException("got error while printing output in txt.", e);
		}
	}

	public void setOutput(Gene gene, Region intron, RegionVector skippedExons) {
		this.geneID = gene.geneID;
		this.geneName = gene.geneName;
		this.chr = gene.geneChr;
		this.strand = gene.strand;
		this.nprots = gene.transcripts.size();
		this.ntrans = gene.nTrans();
		this.sv = intron;
		this.wt = skippedExons;
	}
	
}
