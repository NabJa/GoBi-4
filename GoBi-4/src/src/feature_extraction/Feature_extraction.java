package src.feature_extraction;

import java.util.ArrayList;

import src.genomicUtils.*;

public class Feature_extraction {

	ArrayList<Triplet<String, String, String>> bams = new ArrayList<Triplet<String, String, String>>();

	static String gtf;
	static String bam;
	static String outDest;
	static int frstrand = 0;
	
	public static void main(String args[]) {

		readArguments(args);
		
		SamGtfReader samGTF = new SamGtfReader();
		System.out.println("Start reading GTF");
		samGTF.readExon(gtf);
		
		SAMReader samReader = new SAMReader(samGTF, frstrand, outDest);
		
		System.out.println("Start reading SAM");
		samReader.readSAM(bam);
		samReader.writer.closeBAMFeatures();

	}

	public static void readArguments(String args[]) {

		String noInp = "Error while reading input instructions!!";

		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-gtf":
				gtf = args[i + 1];
				i++;
				System.out.println("GTF: " + gtf);
				break;

			case "-bam":
				bam = args[i + 1];
				i++;
				System.out.println("BAM: " + bam);
				break;
				
			case "-o":
				outDest = args[i + 1];
				i++;
				System.out.println("Output Destination: " + outDest);
				break;
				
			case "-frstrand":
				Boolean frstrandBol = Boolean.parseBoolean(args[i + 1]);
				if(frstrandBol == true) {
					frstrand = 1;
				} else if(frstrandBol == false) {
					frstrand = -1;
				}
				i++;
				System.out.println("Fragment strand: " + frstrand);
				break;
			
			default:
				System.out.println(noInp);
			}
		}	
		System.out.println();
	}
}
