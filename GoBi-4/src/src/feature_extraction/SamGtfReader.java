package src.feature_extraction;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import src.genomicUtils.*;

public class SamGtfReader {

	// HashMap<String, Gene> genesPos = new HashMap<String, Gene>();
	// HashMap<String, Gene> genesNeg = new HashMap<String, Gene>();

	HashMap<String, Chromosom> allChromosomes = new HashMap<String, Chromosom>();
	public HashMap<String, Gene> genes = new HashMap<String, Gene>();

	public void readExon(String path) {

		BufferedReader reader = null;

		try {
			File file = new File(path);
			reader = new BufferedReader(new FileReader(file));

			String rline = "";

			int count = 0;

			while ((rline = reader.readLine()) != null) {

				if (rline.indexOf('#') != 0) {

					int firstTab = rline.indexOf('\t');
					int secondTab = rline.indexOf('\t', firstTab + 1);

					if (rline.substring(secondTab + 1, secondTab + 4).toLowerCase().equals("cds") || rline.substring(secondTab + 1, secondTab + 4).toLowerCase().equals("exon")) {

						String[] line = rline.split("(\t)|(;)");

						String chr = line[0];
						int start = Integer.parseInt(line[3]);
						int end = Integer.parseInt(line[4]);
						String strand = line[6];

						String bioType = biotypeSearch(line);
						String transID = transcriptIDSearch(line);
						String geneID = geneIDSearch(line);
						// String geneName = geneNameSearch(line);
						String proteinID = proteinIDSearch(line);

						Gene gene = new Gene();

						gene.setGene(geneID, chr, strand, start, end, bioType);

						fillGenes(gene, geneID, proteinID, transID, start, end);

						Chromosom chromosom = allChromosomes.get(chr);

						HashMap<String, Gene> genesPos;
						HashMap<String, Gene> genesNeg;

						if (chromosom != null) {
							genesPos = chromosom.genesPos();
							genesNeg = chromosom.genesNeg();
						} else {
							genesPos = new HashMap<String, Gene>();
							genesNeg = new HashMap<String, Gene>();
							Chromosom newChrom = new Chromosom(genesPos, genesNeg);
							allChromosomes.put(chr, newChrom);
							// System.out.println(count + " " + chr);
							count++;
						}

						Gene newGene;
						if (gene.strand.equals("-")) {
							newGene = genesNeg.putIfAbsent(geneID, gene);
						} else {
							newGene = genesPos.putIfAbsent(geneID, gene);
						}

						Region cds = new Region(start, end, bioType);
						RegionVector transcript = new RegionVector(transID);

						if (newGene == null) // means this gene is new
						{
							transcript.addRegion(cds);
							gene.transcripts.put(transID, transcript);

						} else // gene already exists
						{
							Gene correspondingGene;
							if (gene.strand.equals("-")) {
								allChromosomes.get(gene.getChr()).genesNeg().get(geneID).updatePos(start, end);
								correspondingGene = allChromosomes.get(gene.getChr()).genesNeg().get(geneID);
							} else {
								allChromosomes.get(gene.getChr()).genesPos().get(geneID).updatePos(start, end);
								correspondingGene = allChromosomes.get(gene.getChr()).genesPos().get(geneID);
							}

							RegionVector newTrans = correspondingGene.transcripts.putIfAbsent(transID, transcript);

							if (newTrans == null) // means this Trans is new
							{
								transcript.addRegion(cds);
								correspondingGene.transcripts.put(transID, transcript);
							} else // Trans already exists
							{
								if (gene.strand.equals("-")) {
									allChromosomes.get(gene.getChr()).genesNeg().get(geneID).transcripts.get(transID)
											.addRegion(cds);
								} else {
									allChromosomes.get(gene.getChr()).genesPos().get(geneID).transcripts.get(transID)
											.addRegion(cds);
								}
							}
						}
					}
				}
			}

		} catch (Exception e) {
			throw new RuntimeException("got error while reading gtf.", e);
		} finally {

			try {
				reader.close();
			} catch (Exception e) {
				throw new RuntimeException("got error while closing gtf.", e);
			}

		}
	}

	public void fillGenes(Gene gene, String geneID, String proteinID, String transID, int start, int end) {
		Gene newGene = genes.putIfAbsent(geneID, gene);

		Region cds = new Region(start, end, proteinID);
		RegionVector transcript = new RegionVector(transID);

		if (newGene == null) // means this gene is new
		{
			transcript.addRegion(cds);
			gene.transcripts.put(transID, transcript);

		} else // gene already exists
		{
			Gene correspondingGene = genes.get(geneID);
			RegionVector newTrans = correspondingGene.transcripts.putIfAbsent(transID, transcript);

			if (newTrans == null) // means this Trans is new
			{
				transcript.addRegion(cds);
				correspondingGene.transcripts.put(transID, transcript);
			} else // Trans already exists
			{
				genes.get(geneID).transcripts.get(transID).addRegion(cds);
			}
		}
	}

	public String transcriptIDSearch(String[] line) {
		String targetID = "";
		for (int i = 8; i < line.length; i++) {
			if (line[i].indexOf('t') == 0 && line[i].indexOf('d') == 12) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String proteinIDSearch(String[] line) {
		String targetID = "";
		for (int i = 8; i < line.length; i++) {
			if (line[i].indexOf('p') == 1 && line[i].indexOf('d') == 10) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String geneNameSearch(String[] line) {
		String targetID = "";
		for (int i = 8; i < line.length; i++) {
			if (line[i].indexOf('g') == 1 && line[i].indexOf('m') == 8) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String geneIDSearch(String[] line) {
		String targetID = "";
		for (int i = 8; i < line.length; i++) {
			line[i] = line[i].trim();
			if (line[i].indexOf('g') == 0 && line[i].indexOf('d') == 6) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String biotypeSearch(String[] line) {
		String targetTy = "";
		for (int i = 8; i < line.length; i++) {
			line[i] = line[i].trim();
			if (line[i].indexOf('g') == 0 && line[i].indexOf('y') == 9) {
				targetTy = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetTy;
	}

}
