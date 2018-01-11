package src.feature_extraction;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import src.genomicUtils.Gene;
import src.genomicUtils.ReadPair;
import src.genomicUtils.Region;
import src.genomicUtils.Tuple;

public class Feature_writer {

	String outputDestination;
	FileWriter file;
	BufferedWriter writer;

	public Feature_writer(String outputDestination) {
		this.outputDestination = outputDestination;

		try {
			this.file = new FileWriter(outputDestination);
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating BufferedWriter for Feature_writer.", e);
		}
	}

	public double calcPSI(double incl, double excl) {
		return incl / (incl + excl);
	}

	public void writeAltSplicing(HashMap<Region, Tuple<Integer, Integer>> countMap) {

		try {
			writer.write("gene" + "\t" + "exon" + "\t" + "num_incl_reads" + "\t" + "num_excl_reads" + "\t"
					+ "num_total_reads" + "\t" + "psi");
			writer.newLine();
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		for (Region exon : countMap.keySet()) {
			Tuple<Integer, Integer> count = countMap.get(exon);
			try {
				writer.write(String.format("%s\t%s-%s\t%s\t%s\t%s\t%.3f", exon.geneID, exon.getX1(), exon.getX2()+1,
						count.getFirst(), count.getSecond(), (count.getFirst() + count.getSecond()),
						calcPSI(count.getFirst(), count.getSecond())));
				writer.newLine();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	public void writeBAMFeatures(ReadPair pair, int mm, int clipping, int nsplit, int pcridx) {

		String readid = pair.readName;
		String gene = writeGene(pair);

		try {

			if (pair.trans) {
				writer.write(readid + "\t" + "mm:" + mm + "\t" + "clipping:" + clipping + "\t" + "nsplit:" + nsplit
						+ "\t" + "gcount:" + pair.geneAnno.annotGenes.size() + "\t" + gene + "\t" + "pcrindex: "
						+ pcridx);

			} else if (pair.mergedTrans) {
				writer.write(readid + "\t" + "mm:" + mm + "\t" + "clipping:" + clipping + "\t" + "nsplit:" + nsplit
						+ "\t" + "gcount:" + pair.geneAnno.genes.size() + "\t" + gene + "\t" + "pcrindex: " + pcridx);

			} else if (pair.intronic) {
				writer.write(readid + "\t" + "mm:" + mm + "\t" + "clipping:" + clipping + "\t" + "nsplit:" + nsplit
						+ "\t" + "gcount:" + pair.geneAnno.genes.size() + "\t" + gene + "\t" + "pcrindex: " + pcridx);

			} else {
				writer.write(readid + "\t" + "mm:" + mm + "\t" + "clipping:" + clipping + "\t" + "nsplit:" + nsplit
						+ "\t" + "gcount:" + 0 + "\t" + "gdist:" + pair.gdist + "\t" + "antisense:" + pair.antisense
						+ "\t" + "pcrindex: " + pcridx);
			}

			writer.newLine();
		} catch (Exception e) {
			throw new RuntimeException("Got error while writing Feature writer." + e);
		}
	}

	public void writeBAMFeatureSplitIncon(String readid) {
		try {

			writer.write(readid + "\t" + "split-inconsistent:true");
			writer.newLine();

		} catch (Exception e) {
			throw new RuntimeException("Got error while writing Feature writer." + e);
		}
	}

	public void writeRPKMout(ReadPair pair, int pcridx) {

		String gene = writeGeneOnly(pair);

		try {
			writer.write(pair.readName + "\t" + gene + "\t" + pcridx);
			writer.newLine();
		} catch (Exception e) {
			// TODO: handle exception
		}

	}

	public void closeBAMFeatures() {
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("Got error while closing Feature writer", e);
		}
	}

	public String writeGeneOnly(ReadPair pair) {

		String gene = "";
		if (pair.intronic) {
			String preGene = "";
			for (Gene g : pair.geneAnno.genes) {
				preGene += "|" + g.geneID;
			}
			gene = preGene.substring(1, preGene.length());

		}
		if (pair.trans) {
			String preGene = "";
			for (String g : pair.geneAnno.annotGenes.keySet()) {
				preGene += "|" + g;
			}
			gene = preGene.substring(1, preGene.length());
			return gene;
		} else if (pair.mergedTrans) {
			String preGene = "";
			for (Gene g : pair.geneAnno.genes) {
				preGene += "|" + g.geneID;
			}
			gene = preGene.substring(1, preGene.length());
		} else {
			for (String g : pair.geneAnno.annotGenes.keySet()) {
				gene += "|" + g;
			}
			if (gene.length() > 0) {
				gene.substring(1, gene.length());

			}
		}
		return gene;
	}

	public String writeGene(ReadPair pair) {

		String gene = "";
		if (pair.intronic) {
			String preGene = "";
			for (Gene g : pair.geneAnno.genes) {
				preGene += "|" + g.geneID + "," + g.type + ":INTRON";
			}
			gene = preGene.substring(1, preGene.length());

		}
		if (pair.trans) {
			String preGene = "";
			for (String g : pair.geneAnno.annotGenes.keySet()) {
				String preTrans = "";
				preGene += "|" + g + "," + pair.geneAnno.transGenes.get(g).type + ":";
				for (String t : pair.geneAnno.annotGenes.get(g)) {
					preTrans += "," + t;
				}
				preGene += preTrans.substring(1, preTrans.length());
			}
			gene = preGene.substring(1, preGene.length());
			return gene;
		} else if (pair.mergedTrans) {
			String preGene = "";
			for (Gene g : pair.geneAnno.genes) {
				preGene += "|" + g.geneID + "," + g.type + ":MERGED";
			}
			gene = preGene.substring(1, preGene.length());
		} else {
			for (String g : pair.geneAnno.annotGenes.keySet()) {
				gene += g + "," + pair.geneAnno.transGenes.get(g).type + ":INTRON" + "|";
			}
		}
		return gene;
	}
}
