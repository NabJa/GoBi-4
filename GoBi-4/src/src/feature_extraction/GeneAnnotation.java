package src.feature_extraction;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import src.genomicUtils.Gene;

public class GeneAnnotation {

	// HashMap<Tuple<String, String>, ArrayList<String>> annotGenes = new
	// HashMap<Tuple<String, String>, ArrayList<String>>();
	public HashMap<String, ArrayList<String>> annotGenes = new HashMap<String, ArrayList<String>>();
	public HashMap<String, Gene> transGenes = new HashMap<String, Gene>();
	public ArrayList<Gene> genes = new ArrayList<Gene>();
	public Gene mergGene;

	public void annotateGene(Gene gene, String trans) {
		if (annotGenes.containsKey(gene.geneID)) {
			ArrayList<String> transcripts = annotGenes.get(gene.geneID);
			transcripts.add(trans);
		} else {
			ArrayList<String> transcripts = new ArrayList<String>();
			transcripts.add(trans);
			annotGenes.put(gene.geneID, transcripts);
		}
		transGenes.put(gene.geneID, gene);
	}

	public void annotateMerged(Gene gene) {
		this.genes.add(gene);
	}

	public void annotateIntronic(ArrayList<Gene> genes) {
		this.genes = genes;
	}

}
