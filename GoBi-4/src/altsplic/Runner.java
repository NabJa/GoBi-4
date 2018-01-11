package altsplic;

import exonSkip.*;
import src.feature_extraction.*;
import src.genomicUtils.*;

public class Runner {

	public static String gtf;
	public static String bam;
	public static String out;

	public static void main(String[] args) {

		readArgs(args);

		// CDSReader cdsReader = new CDSReader();
		// cdsReader.readCDS(gtf);

		SamGtfReader gtfFile = new SamGtfReader();
		gtfFile.readExon(gtf);

		OutputMap outMap = new OutputMap(out);
		RVaktions compr = new RVaktions();

		for (Gene gene : gtfFile.genes.values()) {
			compr.getSkippedExonFromGen(gene, outMap);
		}
		
		outMap.printExons();  //prints Gene \t Exon with corresponding WT and SV

		SAMReader samReader = new SAMReader(gtfFile, 0, out, outMap);
		System.out.println("Start reading SAM");
		samReader.readSAM(bam);
		samReader.writeAltSplicing();
		samReader.writer.closeBAMFeatures();

		
//		Region r1 = new Region(10, 15);
//		Region sr = new Region(170,180);
//		Region r2 = new Region(20, 25);
//		
//		Region r3 = new Region(8, 13);
//		Region r4 = new Region(12, 15);
//		Region r5 = new Region(21, 25);
//		
//		RegionVector rv1 = new RegionVector();
//		rv1.addRegion(r1);
//		rv1.addRegion(r2);
//		rv1.addRegion(sr);
//		
//		RegionVector rv2 = new RegionVector();
//		rv2.addRegion(r3);
//		rv2.addRegion(r4);
//		rv2.addRegion(r5);
//		
//		System.out.println(rv1.overlapsAll(rv2));
//		System.out.println(rv2.overlapsAll(rv1));

	}

	public static void readArgs(String[] args) {
		String noInp = "Pls enter -o fllowed by output path and -gtf followed by input path !!";

		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-gtf":
				gtf = args[i + 1];
				i++;
				break;
			case "-bam":
				bam = args[i + 1];
				i++;
				break;
			case "-o":
				out = args[i + 1];
				i++;
				break;
			default:
				System.out.println(noInp);
			}
		}

	}

}
