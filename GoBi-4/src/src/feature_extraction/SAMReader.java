package src.feature_extraction;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import augmentedTree.IntervalTree;
import exonSkip.Output;
import exonSkip.OutputMap;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import src.genomicUtils.Chromosom;
import src.genomicUtils.Gene;
import src.genomicUtils.ReadPair;
import src.genomicUtils.Region;
import src.genomicUtils.RegionVector;
import src.genomicUtils.Tuple;

public class SAMReader {

	public static HashMap<String, SAMRecord> lookup = null;

	HashMap<String, Chromosom> allChromosomes;

	IntervalTree<Gene> interTreeNeg = new IntervalTree<Gene>();
	IntervalTree<Gene> interTreePos = new IntervalTree<Gene>();
	HashMap<RegionVector, Integer> pcrIndex = new HashMap<RegionVector, Integer>();
	HashMap<Region, Tuple<Integer, Integer>> countMap = new HashMap<Region, Tuple<Integer, Integer>>(); // Exon -> Incl
																										// Excl
	OutputMap outMap = new OutputMap();
	public Feature_writer writer;

	Integer frstrand;

	int flagCount = 0;
	int firstRegionCheck = 0;
	int secondRegionCheck = 0;
	int readPairs = 0;

	boolean transcriptomic = false;
	boolean merged = false;
	boolean intronic = false;

	public SAMReader(SamGtfReader samGTF, int frstrand, String outputDestination, OutputMap outMap) {
		this.writer = new Feature_writer(outputDestination);
		this.allChromosomes = samGTF.allChromosomes;
		this.frstrand = frstrand;
		this.outMap = outMap;

		for (Region r : outMap.skippedExons) {
			countMap.put(r, new Tuple<Integer, Integer>(0, 0));
		}
	}

	public SAMReader(SamGtfReader samGTF, int frstrand, String outputDestination) {
		this.writer = new Feature_writer(outputDestination);
		// this.interTreeNeg.addAll(samGTF.genesNeg.values());
		// this.interTreePos.addAll(samGTF.genesPos.values());
		this.allChromosomes = samGTF.allChromosomes;
		this.frstrand = frstrand;
	}

	public void readSAM(String bamPath) {
		System.out.println("All chroms: " + allChromosomes.size());

		File bamF = new File(bamPath);
		SAMFileReader sam_reader = new SAMFileReader(bamF, false);
		sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		Iterator<SAMRecord> it = sam_reader.iterator();

		String chromosom = "";

		while (it.hasNext()) {

			SAMRecord sr = it.next();

			if (!sr.getReferenceName().equals(chromosom)) // On new chromosome
			{
				lookup = new HashMap<String, SAMRecord>();

				IntervalTree<Gene> intrTreePos = new IntervalTree<Gene>();
				IntervalTree<Gene> intrTreeNeg = new IntervalTree<Gene>();

				intrTreePos.addAll(allChromosomes.get(sr.getReferenceName()).genesPos().values());
				intrTreeNeg.addAll(allChromosomes.get(sr.getReferenceName()).genesNeg().values());

				this.interTreePos = intrTreePos;
				this.interTreeNeg = intrTreeNeg;

				pcrIndex = new HashMap<RegionVector, Integer>();
			}

			chromosom = sr.getReferenceName();
			SAMRecord other_seen = lookup.get(sr.getReadName());

			
			if (other_seen != null) // Found Read Pair
			{
				if(other_seen.getReadName().equals("75928") || sr.getReadName().equals("75928")) {
					System.out.println();
				}
				readPairs++;

				if (!isSuperIntergenic(other_seen, sr)) // Pair shouldnt be filtered
				{

					ReadPair pair = new ReadPair(sr, other_seen);

					int nsplit = getNSplit(sr, other_seen);

					if (nsplit != -1) {

						pair.calcGenomicRegions(frstrand);
						pair = searchGenes(pair);

						if (pair.trans == true) {

							for (String geneID : pair.geneAnno.annotGenes.keySet()) {
								ArrayList<Region> exonsOfGene = outMap.geneToExon.get(geneID);

								if (exonsOfGene != null) {

									for (Region exon : exonsOfGene) {
										Output out;

										if (outMap.resultMapExons.containsKey(exon)) {
											out = outMap.resultMapExons.get(exon);
											RegionVector sv = new RegionVector();
											RegionVector wt;
											
											if (out.exonToSV.values().size() > 0 && out.exonToWT.containsKey(exon)) {
												for (Region intron : out.exonToSV.keySet()) {
													sv = out.exonToSV.get(intron);
												}
												wt = out.exonToWT.get(exon);

												Region read1 = new Region(pair.firstRegions.x1, pair.firstRegions.x2);
												Region read2 = new Region(pair.secondRegions.x1, pair.secondRegions.x2);

												Region readPair = new Region(pair.genomicRegions.x1,
														pair.genomicRegions.x2);
												
												Region exclusivExon = new Region(exon.x1, exon.x2 + 1);
																								
												if (readPair.overlaps(exclusivExon)) {
													updateCountMap(exon, exclusivExon, read1, read2);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			// Test if read-pair can be ignored:
			if (!can_ignore(sr)) {
				lookup.put(sr.getReadName(), sr);
			}
		}

		System.out.println("Filtered by Flag: " + flagCount);
		System.out.println("Filtered first Region check: " + firstRegionCheck);
		System.out.println("Filtered second Region check: " + secondRegionCheck);
		System.out.println("Annotated read Pairs: " + (readPairs - secondRegionCheck));

		try {
			sam_reader.close();
		} catch (Exception e) {
			throw new RuntimeException("Error while closing BAM reader!", e);
		}
	}

	public void updateCountMap(Region exon, Region exclusivExon, Region read1, Region read2) {
		if (read1.overlaps(exclusivExon) || read2.overlaps(exclusivExon)) {
			Tuple<Integer, Integer> count = countMap.get(exon);
			int incl = count.getFirst() + 1;
			int excl = count.getSecond();
			Tuple<Integer, Integer> newCounts = new Tuple<Integer, Integer>(incl, excl);
			countMap.put(exon, newCounts);
		} else {
			Tuple<Integer, Integer> count = countMap.get(exon);
			int incl = count.getFirst();
			int excl = count.getSecond() + 1;
			Tuple<Integer, Integer> newCounts = new Tuple<Integer, Integer>(incl, excl);
			countMap.put(exon, newCounts);
		}
	}

	/**
	 * Test if there is at least one gene between read pair and no gene spanning
	 * them.
	 * 
	 * @param sr1
	 * @param sr2
	 * @return
	 */
	public boolean isSuperIntergenic(SAMRecord sr1, SAMRecord sr2) {

		boolean sr1_IsFirst = sr1.getFirstOfPairFlag();
		boolean sr1Neg = sr1.getReadNegativeStrandFlag();
		boolean sr2Neg = sr2.getReadNegativeStrandFlag();

		switch (frstrand) {
		case -1:
			if (sr1_IsFirst) {
				if (sr1Neg) {
					if (searchPosTree(sr1, sr2))
						return true;
				} else {
					if (searchNegTree(sr1, sr2))
						return true;
				}
			} else {
				if (sr2Neg) {
					if (searchPosTree(sr1, sr2))
						return true;
				} else {
					if (searchNegTree(sr1, sr2))
						return true;
				}
			}
			break;

		case 1:
			if (sr1_IsFirst) {
				if (sr1Neg) {
					if (searchNegTree(sr1, sr2))
						return true;
				} else {
					if (searchPosTree(sr1, sr2))
						return true;
				}
			} else {
				if (sr2Neg) {
					if (searchNegTree(sr1, sr2))
						return true;
				} else {
					if (searchPosTree(sr1, sr2))
						return true;
				}
			}
			break;

		case 0:
			if (searchNegAndPosTree(sr1, sr2))
				return true;
			break;

		default:
			System.out.println("No frstrand given!");
			if (searchNegAndPosTree(sr1, sr2))
				return true;
			break;
		}
		return false;
	}

	public void writeAltSplicing() {
		writer.writeAltSplicing(countMap);
	}

	public boolean searchNegAndPosTree(SAMRecord sr1, SAMRecord sr2) {

		int start = getMin(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());
		int end = getMax(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());

		ArrayList<Gene> betweenReadsNeg = interTreeNeg.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());
		ArrayList<Gene> betweenReadsPos = interTreePos.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());

		int genesBetweenReadsNeg = 0;
		int genesBetweenReadsPos = 0;

		for (Gene g : betweenReadsNeg) {
			if (g.getChr().equals(sr1.getReferenceName())) {
				genesBetweenReadsNeg++;
			}
		}

		for (Gene g : betweenReadsPos) {
			if (g.getChr().equals(sr1.getReferenceName())) {
				genesBetweenReadsPos++;
			}
		}

		// Is there a gene between the reads?
		if (genesBetweenReadsNeg > 0 || genesBetweenReadsPos > 0) {

			ArrayList<Gene> genesSpanningNeg = interTreeNeg.getIntervalsSpanning(start, end, new ArrayList<Gene>());
			ArrayList<Gene> genesSpanningPos = interTreePos.getIntervalsSpanning(start, end, new ArrayList<Gene>());

			int readsSpanningNeg = 0;
			int readsSpanningPos = 0;

			for (Gene g : genesSpanningNeg) {
				if (g.getChr().equals(sr1.getReferenceName())) {
					readsSpanningNeg++;
				}
			}

			for (Gene g : genesSpanningPos) {
				if (g.getChr().equals(sr1.getReferenceName())) {
					readsSpanningPos++;
				}
			}

			// Are the reads inside of a gene?
			if (readsSpanningNeg == 0 && readsSpanningPos == 0) {
				secondRegionCheck++;
				return true;
			}
		}
		return false;

	}

	public boolean searchNegTree(SAMRecord sr1, SAMRecord sr2) {

		int start = getMin(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());
		int end = getMax(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());

		ArrayList<Gene> betweenReadsNeg = interTreeNeg.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());

		int genesBetweenReadsNeg = 0;

		for (Gene g : betweenReadsNeg) {
			if (g.getChr().equals(sr1.getReferenceName())) {
				genesBetweenReadsNeg++;
			}
		}

		// Is there a gene between the reads?
		if (genesBetweenReadsNeg > 0) {

			ArrayList<Gene> genesSpanningNeg = interTreeNeg.getIntervalsSpanning(start, end, new ArrayList<Gene>());
			int readsSpanningNeg = 0;

			for (Gene g : genesSpanningNeg) {
				if (g.getChr().equals(sr1.getReferenceName())) {
					readsSpanningNeg++;
				}
			}

			// Are the reads inside of a gene?
			if (readsSpanningNeg == 0) {
				secondRegionCheck++;
				return true;
			}
		}
		return false;
	}

	public boolean searchPosTree(SAMRecord sr1, SAMRecord sr2) {

		int start = getMin(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());
		int end = getMax(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(),
				sr2.getAlignmentEnd());

		ArrayList<Gene> betweenReadsPos = interTreePos.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());

		int genesBetweenReadsPos = 0;

		for (Gene g : betweenReadsPos) {
			if (g.getChr().equals(sr1.getReferenceName())) {
				genesBetweenReadsPos++;
			}
		}

		// Is there a gene between the reads?
		if (genesBetweenReadsPos > 0) {

			ArrayList<Gene> genesSpanningPos = interTreePos.getIntervalsSpanning(start, end, new ArrayList<Gene>());
			int readsSpanningPos = 0;

			for (Gene g : genesSpanningPos) {
				if (g.getChr().equals(sr1.getReferenceName())) {
					readsSpanningPos++;
				}
			}

			// Are the reads inside of a gene?
			if (readsSpanningPos == 0) {
				secondRegionCheck++;
				return true;
			}
		}
		return false;
	}

	public boolean isIntergenic(SAMRecord sr) {
		int start1 = Math.min(sr.getAlignmentStart(), sr.getMateAlignmentStart());
		int end1 = Math.max(sr.getAlignmentStart(), sr.getMateAlignmentStart());
		int start = Math.min(start1, sr.getAlignmentEnd());
		int end = Math.max(end1, sr.getAlignmentEnd());

		ArrayList<Gene> betweenReadsNeg = interTreeNeg.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());
		ArrayList<Gene> betweenReadsPos = interTreePos.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());

		if (betweenReadsNeg.size() > 0 || betweenReadsPos.size() > 0) {
			ArrayList<Gene> genesSpanningNeg = interTreeNeg.getIntervalsSpanning(start, end, new ArrayList<Gene>());
			ArrayList<Gene> genesSpanningPos = interTreePos.getIntervalsSpanning(start, end, new ArrayList<Gene>());
			if (genesSpanningNeg.size() == 0 && genesSpanningPos.size() == 0) {
				firstRegionCheck++;
				return true;
			}
		}
		return false;
	}

	public boolean can_ignore(SAMRecord sr) {

		if (sr.getMateUnmappedFlag() || sr.getReadUnmappedFlag() || // is one of the reads unmapped?
				!sr.getReferenceName().equals(sr.getMateReferenceName()) || // are reads on same Chromosome?
				(sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag()) || // are reads on same Strand?
				sr.getMateAlignmentStart() == 0 || sr.getAlignmentStart() == 0 || sr.getNotPrimaryAlignmentFlag()) {
			flagCount++;
			return true;
		}

		// if (isIntergenic(sr))
		// return true;

		return false;
	}

	public int getMissmatches(SAMRecord sr) {
		Integer nm = (Integer) sr.getAttribute("NM");
		nm = (nm != null) ? nm : (Integer) sr.getAttribute("nM");
		nm = (nm != null) ? nm : (Integer) sr.getAttribute("XM");
		return nm;
	}

	public int getClipping(SAMRecord sr) {
		int startClips = sr.getAlignmentStart() - sr.getUnclippedStart();
		int endClips = sr.getUnclippedEnd() - sr.getAlignmentEnd();
		return startClips + endClips;
	}

	public int getNSplit(SAMRecord sr, SAMRecord other_seen) {

		RegionVector readRegions = new RegionVector();
		RegionVector readRegionsRw = new RegionVector();

		for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
			int ref_s = ab.getReferenceStart();
			int ref_end = ref_s + ab.getLength();

			Region readRegion = new Region(ref_s, ref_end);
			readRegions.addRegion(readRegion);
		}
		for (AlignmentBlock ab : other_seen.getAlignmentBlocks()) {

			int ref_s = ab.getReferenceStart();
			int ref_end = ref_s + ab.getLength();

			Region readRegion = new Region(ref_s, ref_end);
			readRegionsRw.addRegion(readRegion);
		}

		RegionVector readIntrons = readRegions.inverse();
		RegionVector readRwIntrons = readRegionsRw.inverse();

		IntervalTree<Region> intronTree = new IntervalTree<Region>();
		IntervalTree<Region> intronRwTree = new IntervalTree<Region>();

		intronTree.addAll(readRegions.regions);
		intronRwTree.addAll(readRegionsRw.regions);

		for (Region intron : readIntrons.regions) {
			if (intronRwTree.getIntervalsIntersecting(intron.getX1() + 1, intron.getX2() - 1, new ArrayList<Region>())
					.size() > 0) {
				return -1;
			}
		}

		for (Region intron : readRwIntrons.regions) {
			if (intronTree.getIntervalsIntersecting(intron.getX1() + 1, intron.getX2() - 1, new ArrayList<Region>())
					.size() > 0) {
				return -1;
			}
		}

		Set<Region> unique = new HashSet<Region>();

		unique.addAll(readIntrons.regions);
		unique.addAll(readRwIntrons.regions);

		return unique.size();
	}

	
	public ReadPair searchGenes(ReadPair pair) {
		switch (frstrand) {
		case -1:
			if (pair.firstIsNeg) // First of pair is negativ and negativ experiment: Search positiv Tree:
			{
				ArrayList<Gene> spanningGenes = interTreePos.getIntervalsSpanning(pair.start, pair.end,
						new ArrayList<Gene>());
				pair = spanningTrans(spanningGenes, pair);

				if (pair.trans) {
					return pair;
				} else {
					pair = searchMergedTrans(spanningGenes, pair);
				}
				if (pair.mergedTrans) {
					return pair;
				} else {
					if (spanningGenes.size() > 0) {
						pair.intronic = true;
						pair.geneAnno.genes.addAll(spanningGenes);
						return pair;
					} else {
						ArrayList<Gene> spanningGenesAntisense = interTreeNeg.getIntervalsSpanning(pair.start, pair.end,
								new ArrayList<Gene>());
						ArrayList<Gene> neighbour = interTreePos.getIntervalsNeighbor(pair.start, pair.end,
								new ArrayList<Gene>());
						pair.getGeneDistance(neighbour);
						// pair = spanningTrans(spanningGenesAntisense, pair);

						if (spanningGenesAntisense.size() > 0) {
							pair.antisense = true;
							pair.trans = false;
						} else {
							pair.antisense = false;
							pair.trans = false;
						}
					}
				}
			} else {
				ArrayList<Gene> spanningGenes = interTreeNeg.getIntervalsSpanning(pair.start, pair.end,
						new ArrayList<Gene>());
				pair = spanningTrans(spanningGenes, pair);
				if (pair.trans) {
					return pair;
				} else {
					pair = searchMergedTrans(spanningGenes, pair);
				}
				if (pair.mergedTrans) {
					return pair;
				} else {
					if (spanningGenes.size() > 0) {
						pair.intronic = true;
						pair.geneAnno.genes.addAll(spanningGenes);
						return pair;
					} else {
						ArrayList<Gene> spanningGenesAntisense = interTreePos.getIntervalsSpanning(pair.start, pair.end,
								new ArrayList<Gene>());
						ArrayList<Gene> neighbour = interTreeNeg.getIntervalsNeighbor(pair.start, pair.end,
								new ArrayList<Gene>());
						pair.getGeneDistance(neighbour);
						// pair = spanningTrans(spanningGenesAntisense, pair);

						if (spanningGenesAntisense.size() > 0) {
							pair.antisense = true;
							pair.trans = false;
						} else {
							pair.antisense = false;
							pair.trans = false;
						}
					}
				}
			}
			return pair;

		case 1:
			if (!pair.firstIsNeg) // First of pair is negativ and negativ experiment: Search positiv Tree:
			{
				ArrayList<Gene> spanningGenes = interTreePos.getIntervalsSpanning(pair.start, pair.end,
						new ArrayList<Gene>());
				pair = spanningTrans(spanningGenes, pair);

				if (pair.trans) {
					return pair;
				} else {
					pair = searchMergedTrans(spanningGenes, pair);
				}
				if (pair.mergedTrans) {
					return pair;
				} else {
					if (spanningGenes.size() > 0) {
						pair.intronic = true;
						pair.geneAnno.genes.addAll(spanningGenes);
						return pair;
					} else {
						ArrayList<Gene> spanningGenesAntisense = interTreeNeg.getIntervalsSpanning(pair.start, pair.end,
								new ArrayList<Gene>());
						ArrayList<Gene> neighbour = interTreePos.getIntervalsNeighbor(pair.start, pair.end,
								new ArrayList<Gene>());
						pair.getGeneDistance(neighbour);
						// pair = spanningTrans(spanningGenesAntisense, pair);

						if (spanningGenesAntisense.size() > 0) {
							pair.antisense = true;
							pair.trans = false;
						} else {
							pair.antisense = false;
							pair.trans = false;
						}
					}
				}
			} else {
				ArrayList<Gene> spanningGenes = interTreeNeg.getIntervalsSpanning(pair.start, pair.end,
						new ArrayList<Gene>());
				pair = spanningTrans(spanningGenes, pair);
				if (pair.trans) {
					return pair;
				} else {
					pair = searchMergedTrans(spanningGenes, pair);
				}
				if (pair.mergedTrans) {
					return pair;
				} else {
					if (spanningGenes.size() > 0) {
						pair.intronic = true;
						pair.geneAnno.genes.addAll(spanningGenes);
						return pair;
					} else {
						ArrayList<Gene> spanningGenesAntisense = interTreePos.getIntervalsSpanning(pair.start, pair.end,
								new ArrayList<Gene>());
						ArrayList<Gene> neighbour = interTreeNeg.getIntervalsNeighbor(pair.start, pair.end,
								new ArrayList<Gene>());
						pair.getGeneDistance(neighbour);
						// pair = spanningTrans(spanningGenesAntisense, pair);

						if (spanningGenesAntisense.size() > 0) {
							pair.antisense = true;
							pair.trans = false;
						} else {
							pair.antisense = false;
							pair.trans = false;
						}
					}
				}
			}
			return pair;

		case 0:
			ArrayList<Gene> spanningGenesN = interTreeNeg.getIntervalsSpanning(pair.start, pair.end,
					new ArrayList<Gene>());
			ArrayList<Gene> spanningGenesP = interTreePos.getIntervalsSpanning(pair.start, pair.end,
					new ArrayList<Gene>());
			spanningGenesN.addAll(spanningGenesP);
			pair = spanningTrans(spanningGenesN, pair);

			if (pair.trans) {
				return pair;
			} else {
				pair = searchMergedTrans(spanningGenesN, pair);
			}
			if (pair.mergedTrans) {
				return pair;
			} else {
				if (spanningGenesN.size() > 0) {
					pair.intronic = true;
					pair.geneAnno.genes.addAll(spanningGenesN);
					return pair;
				} else {
					ArrayList<Gene> neighbourN = interTreeNeg.getIntervalsNeighbor(pair.start, pair.end,
							new ArrayList<Gene>());
					ArrayList<Gene> neighbourP = interTreePos.getIntervalsNeighbor(pair.start, pair.end,
							new ArrayList<Gene>());
					pair.getGeneDistance(neighbourN);
					pair.getGeneDistance(neighbourP);

					// pair = spanningTrans(spanningGenesAntisense, pair);
					pair.trans = false;
				}
			}
			return pair;

		default:
			System.out.println("No frstrand given!");
			return pair;
		}
	}

	public ReadPair spanningTrans(ArrayList<Gene> spanningGenes, ReadPair pair) {

		for (Gene g : spanningGenes) {
			for (RegionVector rv : g.transcripts.values()) {

				RegionVector cutRV1 = rv.cutRV(pair.getFirstStart(), pair.getFirstEnd());
				RegionVector cutRV2 = rv.cutRV(pair.getSecondStart(), pair.getSecondEnd());

				cutRV1.id = rv.id;
				cutRV2.id = rv.id;

				if (cutRV1.isEqual(pair.firstRegions) && cutRV2.isEqual(pair.secondRegions)) {
					pair.geneAnno.annotateGene(g, rv.id);
					pair.trans = true;
				}
			}
		}
		return pair;

	}

	public ReadPair searchMergedTrans(ArrayList<Gene> spanningGenes, ReadPair pair) {
		for (Gene g : spanningGenes) {
			RegionVector mergedTrans = new RegionVector();
			for (RegionVector rv : g.transcripts.values()) {
				mergedTrans = rv.mergeRV(mergedTrans);
			}
			if (mergedTrans.contains(pair.genomicRegions)) {
				pair.geneAnno.annotateMerged(g);
				pair.mergedTrans = true;
			}
		}
		return pair;
	}

	public int getMin(int v, int x, int y, int z) {
		int minVX = Math.min(v, x);
		int minYZ = Math.min(y, z);
		int gloabalMin = Math.min(minVX, minYZ);
		return gloabalMin;
	}

	public int getMax(int v, int x, int y, int z) {
		int maxVX = Math.max(v, x);
		int maxYZ = Math.max(y, z);
		int globalMax = Math.max(maxVX, maxYZ);
		return globalMax;
	}

	public void updatePCRidx(RegionVector rv) {

		Integer queryIDX = pcrIndex.get(rv);

		if (!pcrIndex.containsKey(rv)) {
			pcrIndex.put(rv, 0);
		} else {
			queryIDX++;
			pcrIndex.put(rv, queryIDX);
		}
	}
}
