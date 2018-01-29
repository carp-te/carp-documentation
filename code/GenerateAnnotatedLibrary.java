import java.io.*;
import java.util.*;

/*************************
GenerateAnnotatedLibrary
	Generates an annotated library of the non-protein coding, non SSR repeating patterns in Pogona (PV)
	Author: Joy Raison
	Data: October 16, 2014
	Updated: October 17, 2014 to remove the insertion of the satellite
	Updated: March 3, 2016 to strip member size and length information from consensus sequence names
	   when reading their fasta file. 
	Updated: July 12, 2016 to allow differential coverage limits for SINE elements. 
	Inputs: ConsensusSequences.fa (The consensus sequences (fasta format))
		ConsensusSequences.fa.map (The map file from the censor run)
		known.txt (list of censor IR "identified" sequences with name of the library sequence they match)
		notKnown.fa.tewb.gff (gff file of the GB_TE matched sequences)
		GB_TE.01092014.fa (The GB_TE library)
		notKnown.fa.ervwb.gff (gff file of the all_retrovirus matched sequences)
		all_retrovirus.fasta (The all_retrovirus library)
		SSR.txt (list of sequences identified as SSRs, and the SSR they matched)
		protein.txt (list of the sequences identified as proteins and the protein they matched)
		//LA4v2-satellite.fa (a satellite sequence for which no consensus sequence was found)
		/home/a1635743/RepBase20.04.fasta/*rep.ref (RepBase libraries to base 
			classification on)
	Outputs: wantedCSHeaders.txt (for checking individual sequence headers)
		R4_Library.fasta (The annotated library)
***************************/


public class GenerateAnnotatedLibrary {
	private static class Hit {
		private String target;
		private int start;
		private int end;
		
		public Hit (String target, int start, int end) {
			this.target = target;
			this.start = start;
			this.end = end;
			}

		public int getEnd () {return end;}

		public int getLength() {return end - start + 1;}
		
		public int getStart() {return start;}
		
		public String getTarget() {return target;}
		}

		
	private static class HitComparator implements Comparator {
		public int compare (Object o1, Object o2) throws ClassCastException {
			try {
				Hit h1 = (Hit) o1;
				Hit h2 = (Hit) o2;
				int result = h1.getStart() - h2.getStart();
				if (result == 0) result = h2.getEnd() - h1.getEnd();
				return result;
				}
			catch (ClassCastException ce) {
				throw new ClassCastException("Only Hits can be compared in HitComparator");
				}
			}
				
		public boolean equals (Object o) {return o instanceof HitComparator;}
		}
			

	private static class CS {
		private static HitComparator hitComparator = new HitComparator();
		private static double maxOverlap = 0.8;
		private static double minCoverageBP = 0.05;
		private String name;
		private String classification;
		private String annotation;
		private int length;
		private List<Hit> hits = new ArrayList<Hit>();
		
		public CS (String name, Integer len) {
			this.name = name;
			length = len.intValue();
			}

		public void addHit (Hit hit) {hits.add(hit);}

		public void addHits (List<Hit> hits) {this.hits.addAll(hits);}
		
		public void annotateAndClassify(Hashtable<String, String> retroAnnots, Set<String> irs, 
			RBClassifier classifier, double sineMinCoverage, double restMinCoverage) {
			if (hits.size() == 0) classifyAndAnnotateNone();
			else {
				removeSubHits();
				Hashtable<String, Integer> typeCoverage = new Hashtable<String, Integer>();
				for (Iterator<Hit> iter=hits.iterator(); iter.hasNext();) {
					Hit hit = iter.next();
					String target = hit.getTarget();
					int length = hit.getLength();
					if (typeCoverage.containsKey(target)) length += typeCoverage.get(target).intValue();
					typeCoverage.put(target, new Integer(length));
					}
				double minCoverage = minCoverageBP * length;
				Set<String> tooLittleCoverage = new HashSet<String>();
				for (Enumeration<String> e=typeCoverage.keys(); e.hasMoreElements();) {
					String target = e.nextElement();
					if (typeCoverage.get(target).intValue() < minCoverage) tooLittleCoverage.add(target);
					}
				for (Iterator<String> iter=tooLittleCoverage.iterator(); iter.hasNext();) 
					typeCoverage.remove(iter.next());
				int nTargets = typeCoverage.size();
				if (nTargets == 0) classifyAndAnnotateNone();
				else {
					annotation = getTargetAnnotation(typeCoverage);
					int coverage = getCoverage();
					if (((double) coverage)/length < (allSINE(classifier)?sineMinCoverage:restMinCoverage)) 
						classification = "#PartialAnnotation";
					else if (hits.size() > 1 && !irs.contains(name) && !isSatellite()) 
						classification = "#Chimeric";
					else {
						String target = annotation.substring(0, annotation.indexOf(" "));
						if (target.indexOf("|") > 0) {
							classification = "#Retrovirus_like";
							annotation += " " + retroAnnots.get(target);
							}
						else classification = ":" + target;
						}						
					}
				}
			}

		public String getFastaHeader() {return ">" + name + classification + " " + annotation + "\n";}

		private boolean allSINE(RBClassifier classifier) {// check if all targets are SINEs
			boolean allSINE = true;
			for (Iterator<Hit> iter=hits.iterator(); allSINE && iter.hasNext();)
				allSINE = classifier.classify(iter.next().getTarget()).startsWith("SINE");
			return allSINE;
			}
						
		private void classifyAndAnnotateNone() {
			classification = "#Unclassified";
			annotation = "Matches no similar sequence";
			}

		private int getCoverage() {
			int coverage = 0;
			Hit lastHit = null;
			for (Iterator<Hit> iter=hits.iterator(); iter.hasNext();) {
				Hit hit = iter.next();
				if (lastHit == null) coverage = hit.getLength();
				else {
					if (hit.getStart() > lastHit.getEnd()) coverage += hit.getLength();
					else coverage += hit.getEnd() - lastHit.getEnd();
					}
				lastHit = hit;
				}
			return coverage;
			}
			
		private String getTargetAnnotation (Hashtable<String, Integer> typeCoverage) {
			int n = typeCoverage.size();
			String[] targets = new String[n];
			double[] coveragePC = new double[n];
			double pcDivisor = length / 100.;
			for (Enumeration<String> e=typeCoverage.keys(); e.hasMoreElements();) {
				String target = e.nextElement();
				double pcCoverage = typeCoverage.get(target).intValue() / pcDivisor;
				boolean inserted = false;
				for (int i=0; !inserted && i<n; i++)	
					if (pcCoverage > coveragePC[i]) {
						for (int j=n-2; j>=i; j--) {
							coveragePC[j+1] =coveragePC[j];
							targets[j+1] = targets[j];
							}
						coveragePC[i] = pcCoverage;
						targets[i] = target;
						inserted = true;
						}	
				}
			String anno = "";
			for (int i=0; i<n; i++) anno += (i>0?"; ":"") + targetAnnotation (targets[i], coveragePC[i]);
			return anno;
			}

		private boolean isSatellite () {
			String[] fields = annotation.split(" ");
			return fields.length == 2 && fields[0].indexOf("SAT") >= 0;
			}
																
		private void removeSubHits() {
			Collections.sort(hits, hitComparator);
			List<Hit> subHits = new ArrayList<Hit>();
			Hit lastHit = null;
			for (Iterator<Hit> iter=hits.iterator(); iter.hasNext();) {
				Hit hit = iter.next();
				if (lastHit == null) lastHit = hit;
				else 
					if (hit.getEnd() < lastHit.getEnd()) subHits.add(hit);
					else if (hit.getStart() > lastHit.getEnd()) lastHit = hit;
					else {
						double minLengthLimit = (lastHit.getEnd() - hit.getStart() + 1)/maxOverlap;
						int lastLength = lastHit.getLength();
						int hitLength = hit.getLength();
						if (lastLength < hitLength) {
							if (lastLength < minLengthLimit) subHits.add(lastHit);
							lastHit = hit;
							}
						else if (hitLength < minLengthLimit) subHits.add(hit);	
						else lastHit = hit;						
						}					
				}
			hits.removeAll(subHits);
			}
			
		private String targetAnnotation (String target, double pc) {
			return target + " (" + ((int) (pc + .5)) + ")"; 
			}
		}

	private static class RBClassifier {
  	private String libDir = getLibraryDirectory();
		private Hashtable<String, String> classification = getClassifications(libDir);

  
		public String classify (String name) {
			if (name.startsWith("AFROSINE")) return "SINE2/AFROSINE";
			if (name.startsWith("BTLTR1")) return "LTR/BTLTR1";
			if (name.startsWith("Bov-tA")) return "SINE2/BOV.tA";
			if (name.startsWith("BovB")) return "LINE/RTE_BovB";
			if (name.startsWith("ERE"))
				switch (name.charAt(3)) { 
					case '1': return "SINE2/ERE1";
					case '2': return "SINE2/ERE2";
					case '3': return "SINE2/ERE3";
					case '4': return "SINE2/ERE4";
					}
			if (name.startsWith("MIR") && !name.startsWith("MIRAGE")) return "SINE2/MIR";
			if (name.startsWith("SINEC")) return "SINE2/CanSINE";
			if (name.startsWith("THER")) return "SINE2/MIR";
			if (name.startsWith("ERV-1_PM")) return "ERV/ERV1";
			if (name.startsWith("ERV-2_PM")) return "ERV";
			if (name.startsWith("AluY")) return "SINE1/7SL";
			if (name.startsWith("Dada")) return "DNA/Dada";
			if (name.startsWith("DNA-TTAA0-")) return "DNA";		
			if (name.startsWith("Kolobok")) return "DNA/Kolobok";
			if (name.startsWith("Sola")) return "DNA/Sola";
			String classif = classification.get(name);
			if (classif == null) {if (name.startsWith("ERV")) return "ERV";}
			else {
				if (classif.equals("SINE")) return classif + "/Unclassified";
				if (classif.indexOf("CR1") >= 0) {
		  		if (name.indexOf("LINE") >= 0) classif = (name.indexOf("LINE2") >= 0)?classif + "_L2":"LINE";
		  		else if (name.startsWith("L2") || name.startsWith("CR1-L2")) {
		    		int nameLength = name.indexOf("L2") + 2;
			  		if (name.length()==nameLength) classif += "_L2";
			  		else switch (name.charAt(nameLength)) {
							case ('0'):
							case ('1'):
							case ('2'):
							case ('3'):
							case ('4'):
							case ('5'):
							case ('6'):
							case ('7'):
							case ('8'):
							case ('9'): break;
							default: classif += "_L2";
							}
						}
			  	}
				}
			return classif==null?"Unknown":classif;
			}

  	private void reClassify (Hashtable<String, String> classifications) {
    	String[][] seqClass = {
      	{"ERV/ERV3_MaLR", "MLT1J-int", "MLT1M", "MTE-int", "ORR1B1-int", "ORR1D-int", "ORR1E", 
        	"MLT1G1", "MLT1H", "MLT1H1", "MLT1H2", "MLT1I", "MLT1J", "MLT1J1", "MLT1J2", "MLT1K", 
        	"MTE", "ORR1A0", "ORR1B1", "ORR1B2", "ORR1C1", "ORR1C2", "MLT-int", "MLT1A0", "MLT1A1", 
        	"MLT1C", "MLT1C1", "MLT1D", "MLT1E", "MLT1E1", "MLT1E1A", "MLT1E2", "MLT1F", "MLT1F_I", 
        	"MLT1G", "MLT1G2", "MLT1G3", "MLT1L", "MLT1O", "MLT1_I", "MSTA", "MSTA1", "MSTA2", "MSTB1", 
        	"MSTD", "MST_I", "MTAI", "ORR1AI", "ORR1BI", "THE1A", "THE1B", "THE1C", "THE1D", "MTA", 
        	"MTB", "MTC", "MTD", "ORR1A", "ORR1B", "ORR1C", "ORR1D",  "MLT1B", "MLT1F1", "MLT1F2", 
        	"MSTB", "MSTC"},
      	{"DNA", "UCON71_Crp"},
				{"LTR/ALTR2", "ALTR2"},
				{"SINE/BOVA2", "BOVA2"},       
      	{"SINE2/AFRO_LA", "AFRO_LA"}, 
				{"SINE2/BOV.tA", "BOVTA", "Bovc-tA2"},
				{"SINE2/CanSINE", "CAN"},
				{"SINE2/MIR", "MAR1", "MAR1_MD", "MAR1c_Mdo", "MON1", "SINE-1_MD", "SINE-2_MD", "WSINE1"}
				};
    	for (int c=0; c<seqClass.length; c++) 
      	for (int s=1; s<seqClass[c].length; s++) classifications.put(seqClass[c][s], seqClass[c][0]);
    	}
    
  	private void addUnknownClassifications (Hashtable<String, String> classifications) {
    	String[][] seqClass = {
     		{"DNA", "DNA-1-1_DR", "DNA-1-2_DR", "DNA-1-4_DR", "DNA-1-9_DR", "DNA-2-3_DR", "DNA-2-9_DR",
        	"DNA-2-12_DR", "DNA-2-31_DR", "DNA-2-33_DR", "DNA-4-2_DR", "DNA-8-2_DR", "DNA-8-5_DR",
        	"DNA-8-6_DR", "DNA-8-15_DR", "DNA-8-20_DR", "DNA-8-25_DR", "DNA-8-32_DR", "DNA-8-33_DR",
        	"DNA-8-34_DR", "DNA-8-36_DR", "DNA-TTAA-2_DR", "DNA2-5_DR", "DNA25TWA1_DR", "DNA5-10_CGi",
        	"DNA8-1_DR", "DNA8-5_DR", "DNA8-41_AP", "DNA9-7_STu"},
      	{"DNA/hAT", "hAT-N71_DR", "HATN9_DR"},
      	{"DNA/Mariner/Tc1", "PARISa_DPo", "Sagan-1_PMa"},
      	{"ERV", "LTR6_Ami"},
     		{"ERV/ERV2", "IAPLTR1_Mm_LTR"},
				{"LINE/CR1", "CR1-7_DR", "CR1-15_DR", "CR1-19_DR", "CR1-20_DR", "CR1-21_DR", "CR1-23_DR",
			  	"CR1-26_DR", "CR1-28_DR", "CR1-27_DR", "CR1-29_DR", "CR1-30_DR", "CR1-31_DR", "CR1-38_DR",
			  	"CR1-40_DR", "CR1-43_DR"},
				{"LINE/I", "I-2_DR"},
				{"LINE/Jockey", "TART"},
				{"LTR", "LTR-1_Crp", "LTR-11_DR", "LTR-1302_Crp", "LTR-1B_Crp", "LTR-3_Crp", 
			  	"LTR-775_Gav_odd", "LTR_PP", "LTR1_CR", "LTR3_CR", "PtPiedmont_I", "PtPiedmont_LTR"},
				{"NonLTR", "SSSINESAT"},
				{"SINE", "S1_BN", "SINE_LC", "SINE_SO"},
				{"Transposable Element", "TE-447_AMi"},
				{"Unknown", "DRP_EG", "NTLTR1", "SSNHEI", "Tc1N1_DR"},
				{"Unknown/tandem repeat", "Mf3_MF"} 
				};
			for (int c=0; c<seqClass.length; c++) 
		  	for (int s=1; s<seqClass[c].length; s++) 
		    	if (!classifications.containsKey(seqClass[c][s]))
		      	classifications.put(seqClass[c][s], seqClass[c][0]);
    	}
    
		private Hashtable<String, String> getClassifications (String libDir) {
			Hashtable<String, String> classifications = new Hashtable<String, String>();
			Hashtable<String, String> families = getFamilyClassifications();
			String line = null;
			String[] fields = null;
			String id = null;
			BufferedReader in = null;
			try {
		  	File dir = new File(libDir);
		  	if (dir.exists()) {
		    	File[] libs = dir.listFiles();
			  	for (int i=0; i<libs.length; i++) {
			    	if (libs[i].getName().endsWith("rep.ref")) {
				    	in = new BufferedReader(new FileReader(libs[i]));
				   		while ((line = in.readLine()) != null) 
					    	if (line.charAt(0) == '>') {
						    	fields = line.split("\t");
						    	id = fields[0].substring(1);
						    	if (!classifications.containsKey(id)) 
							    	if (fields.length < 2) 
								    	classifications.put(id, fields[0].endsWith("(n)")?"Simple Repeat/" + 
								      	id:"Unknown");
							    	else if (families.containsKey(fields[1])) 
							      	classifications.put(id, families.get(fields[1]));
							    	else classifications.put(id, "Unknown/" + fields[1]);
						    	}
				    	in.close();
				    	}
				  	}
					}
				else System.out.println("The directory " + libDir + " does not exist.");
				}
			catch (IOException ie) {ie.printStackTrace();}
			reClassify(classifications);
			addUnknownClassifications(classifications);
			//System.out.println("There are " + classifications.size() + " classifications");
			return classifications;
			}

		private Hashtable<String, String> getFamilyClassifications() {
			Hashtable<String, String> families = new Hashtable<String, String>();
			String[][] classes = {
		  	{"DNA", "Mariner/Tc1", "hAT", "Repetitive element", "Repeat DNA", "DNA transposon", "AMTAM2", 
		    	"APO1_AP", "APO2_AP", "ARS_TA", "BHIKHARI_I", "BMRP1", "Ginger2/TDD", "Academ", "Zisupton", 
		    	"BREP1", "BS1", "BstUI repeat", "CAM2_GG", "CEREP3", "CERP2", "CERP3", "CERP4", "CEU86951",
		    	"Charlie-Galluhop", "CRTOC1", "CryptonS", "C_OC", "CHD", "Chapaev", "CSP2034", "DMRP1",
		    	"MuDR", "EnSpm", "Eutr1", "Eutr10", "Eutr11", "Eutr12", "Eutr13", "Eutr14", "Eutr15",
					"Eutr16", "Eutr17", "Eutr18", "Eutr2", "Eutr3", "Eutr4", "Eutr5", "Eutr6", "Eutr9", "DMRPR",
					"EUTREP11", "EUTREP12", "EUTREP14", "EUTREP15", "EUTREP16", "EUTREP2", "EUTREP4", "DRB_RN",
					"EUTREP5", "EUTREP6", "EUTREP7", "EUTREP8", "MARE10", "MARE11", "MARE4", "MARE8", "LVU1",
					"MARE9", "BDDF1", "P", "Merlin", "Harbinger", "Transib", "Novosib", "Helitron", "Polinton", 
					"Kolobok", "ISL2EU", "Crypton", "MER35", "OOREP1", "piggyBac", "Zator", "Ginger1", 
					"IS3EU", "2109A repetitive sequence", "ALBAMH1", "ARS406", "AVIXHoI", "CERP1", "CMREP", 
					"AY1 repetitive sequence", "CSP2090", "CSP2111", "CSP2112", "DDTDD", "DMFUSHI", "DMHMR2",
					"DQ524338", "EcoR1 repeat region", "EHINV1", "EHINV2", "EnSpm/CACTA", "ERACSI34_EA", 
					"ERASCI228", "FPREP1", "FR1", "FUGUREP4B", "GISH1_AC", "GPRP1", "GQRP1", "HHAI", "KER",
					"HIND3_MS", "HTE1", "IKIRARA1", "Interspersed repeat", "Inverted repeat", "JH12_XL",
					"Interspersed repetitive sequence", "KRISPIE", "LARP1", "LARP2", "LARRP1", "LDRP2", "LR9A",
					"LGRP1", "LIRP1", "LMRP1", "MCMREP", "MER122", "MER22", "MICROPON-LIKE-1", "MITE", "OARP1",
					"MICROPON-LIKE-2", "miniature inverted repeat", "Minicircle", "MINIME_DN", "MRE1_OL",
					"MSRBMI", "DMHMR1", "NTS_DM", "OFU85403", "Origin of replication", "P-element", "PAT", 
					"PEN1", "PEN2", "PEN4", "PFRP1", "PFRP5", "pSOS family", "R1A_SS", "R1B_DS", "R1B_SS",
					"RC14", "RCS5", "REP-1_Croc", "REP-540_Croc", "Repeat region", "Repetitive element Vi",
					"Repetitive sequence", "RMER1", "RMER1A", "RMER1B", "RP1_GL", "RP2_HV", "RP5S", "SCAI_EH",
					"RS3 repeat region", "SIRE", "STREPB_FA", "STREPE_PF", "SZ23_TC", "Tcn760", 
					"transposon", "TREP_CE", "VADER", "VEGE_DW", "XbaI", "Y\' element", "5S_DM", "ALAD", 
					"AFRP1", "MU4_ZM"},
				{"ERV", "ERV1", "ERV2", "ERV3", "Lentivirus", "ERV4"},
				{"LINE", "R4", "L1", "RTE", "I", "Jockey", "CR1", "RTEX", "L2", "Vingi", "CRE", "NeSL", "R2", 
			  	"Rex1", "RandI", "Tx1", "Crack", "Nimb", "Proto1", "Proto2", "RTETP", "Hero", "Tad1", 
			  	"Loa", "Ingi", "Outcast", "R1", "Daphne", "Ambal", "Kiri", "L2A", "L2B"},
				{"LTR", "Gypsy", "Copia", "BEL", "DIRS"},
				{"Satellite", "SAT", "MSAT"},
				{"Pseudogene", "rRNA", "tRNA", "snRNA"},
				{"Integrated Virus", "DNA Virus", "Caulimoviridae"}};
			String clas = null;
			for (int c=0; c<classes.length; c++) {
				clas = classes[c][0] + "/";
				for (int f=1; f<classes[c].length; f++) families.put(classes[c][f], clas + classes[c][f]);
				}
			String[][] fams = {
		  	{"DNA/Male specific", "Male-specific DNA", "Y chromosome", },
		  	{"ERV", "Endogenous Retrovirus"},
		  	{"LTR", "LTR Retrotransposon", "Long terminal repeat", "LTR-775_Gav_odd"},
		  	{"NonLTR", "Non-LTR Retrotransposon"},
		  	{"Simple Repeat", "AT-rich DNA repeat", "(CA)n related"},
		  	{"SINE", "short interspersed element"},
		  	{"Unknown", "conserved", "Internal sequence of mouse VL30 retro-element", "Nonautonomous",
		    	"TREP16", "TREP17"},
		  	{"Unknown/centromere repeat", "CENSTRIG", "RCH2", "CEN1_SP"},
		  	{"Unknown/centromere-associated repeat", "ATREPTSEQ"},
		  	{"Unknown/direct repeat", "direct repeat", "ISFUN1"},
		  	{"Unknown/dispersed repeat", "ATREP19", "Dispersed repeat", "CEREP4", "CEREP5", "SCAR_MA",
		    	"SCAR_MI"},
		  	{"Unknown/midrepetitive element", "SNAPBACK_TC", "TELREP_AG"},
		  	{"Unknown/Pericentromeric repeat", "IID2-12_AT"},
		  	{"Unknown/tandem repeat", "AlKe1_AL", "ATREP18", "BamHI repetitive sequence", "CARO_CA", 
		    	"D1100 family", "DEC1_DS", "DMHETRP", "CPTAN", "EcoRI family", "ECORI_Hm", "HHA1_BT", 
		    	"OVRP1", "PTR_XL", "SAL_CL", "SAU3A_TR", "STTREP_Mp", "tandem repeat", "TANDREP_TG"}
		  	};
			for (int c=0; c<fams.length; c++)
		  	for (int f=1; f<fams[c].length; f++) families.put(fams[c][f], fams[c][0]);
			//Classes that are subclasses
			String[] fam = {"SINE", "Simple Repeat", "Satellite", "Transposable Element", 
		  	"Integrated Virus", "Pseudogene", "SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE4", "Penelope"};
			for (int f=0; f<fam.length; f++) families.put(fam[f], fam[f]);
			//System.out.println("There are " + families.size() + " family classifications");
			return families;
			}

		private String getLibraryDirectory () {
			return "/home/a1635743/RepBase20.04.fasta";
			}
		}

	private static String iDir = "./";
	private static String oDir = "library/";
	private static String library = oDir + "Denovo_TE_Library.fasta";
	private static String headers = oDir + "wantedCSHeaders.txt";
	//private static String satFile = sDir + "LA4v2-satellite.fa";//File moved to all-repeats
	private static String CSFile = iDir + "ConsensusSequences.fa";
	private static String TEgff = iDir + "notKnown.fa.tewb.gff";
	private static String GBTE = iDir + "GB_TE.21032016.fa";
	private static String ERVgff = iDir + "notKnown.fa.ervwb.gff";
	private static String ALLR = iDir + "all_retrovirus.fasta";
	private static String SSR = iDir + "SSR.txt";
	private static String Proteins = iDir + "protein.txt";
	private static String IRS = iDir + "ConsensusSequences.fa.map"; 
	private static String IRM = iDir + "known.txt";
	private static double restMinCoverage = .9;
	private static double sineMinCoverage = .9;
	private static boolean debug = false;

	public static void main (String[] args) {
		setSineMinCoverage (args);
		//debug = true;
		boolean headersOnly = false;
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(headersOnly?headers:library));
			//This was not needed as it was from all-repeats
			//writeSatellite(out, satFile, "family011387#Satellite", headersOnly);
			writeConsensusSequences(out, headersOnly);
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}


	private static void addIRHits (Hashtable<String, CS> wantedCS) {
		String line = null;
		try {
			StringTokenizer st = null;
			BufferedReader in = new BufferedReader(new FileReader(IRS));
			while ((line = in.readLine()) != null) {
				st = new StringTokenizer(line);
				String seq = st.nextToken();
				if (wantedCS.containsKey(seq)) {
					int start = Integer.parseInt(st.nextToken());
					int end = Integer.parseInt(st.nextToken());
					wantedCS.get(seq).addHit(new Hit(st.nextToken(), start, end));
					}
				}
				in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);}
		}
		
	private static void addRetroHits(Hashtable<String, List<Hit>> hits, String gffFile) {
		String line = null;
		try {
			BufferedReader in = new BufferedReader(new FileReader(gffFile));
			String[] fields = null;
			while ((line = in.readLine()) != null) {
				fields = line.split("\t");
				if (!hits.containsKey(fields[0])) {
					hits.put(fields[0], new ArrayList<Hit>());
					//trace("Adding Sequence: " + fields[0]);
					}
				int i1 = fields[8].indexOf(" ") + 1;
				int i2 = fields[8].indexOf(" ", i1);
				String target = fields[8].substring(i1, i2);
				hits.get(fields[0]).add(new Hit(target, Integer.parseInt(fields[3]), 
					Integer.parseInt(fields[4])));
				}
			in.close();
			trace("There are " + hits.size() + " RetroHits after adding " + gffFile);
			}
		catch (IOException ie) {ie.printStackTrace();}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);}
		}

	private static Hashtable<String, String> getAllRetroAnnotations() {
		String[] files = {GBTE, ALLR};
		Hashtable<String, String> annos = new Hashtable<String, String>();
		String line = null;
		BufferedReader in = null;
		for (int i=0; i<files.length; i++) 
			try {
				in = new BufferedReader(new FileReader(files[i]));
				while ((line = in.readLine()) != null) 
					if (line.length() > 0 && line.charAt(0) == '>') {
						int i1 = line.indexOf("|", 4) + 1;
						int i2 = line.indexOf("|", line.indexOf("|", i1)+1);
						annos.put(line.substring(i1, i2), line.substring(i2+2));
						}
				in.close();
				trace("There are " + annos.size() + " uniquely identified sequences after " + files[i]); 
				}
			catch (IOException ie) {ie.printStackTrace();}
		return annos;
		}

	private static Set<String> getFamilies(String inFile) {
		Set<String> families = new HashSet<String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			in.readLine();
			String line = null;
			while ((line = in.readLine()) != null) families.add(line.substring(0, line.indexOf(" ")));
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		trace("There are " + families.size() + " families");
		return families;
		}
		
	private static Hashtable<String, Integer> getLengths (String inFile) {
		Hashtable<String, Integer> lengths = new Hashtable<String, Integer>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			int length = 0;
			String id = null;
			String line = null;
			while ((line = in.readLine()) != null) 
				if (line.charAt(0) == '>') {
					if (id != null) lengths.put(id, new Integer(length));
					int index = line.indexOf(" ");
					id = index<0?line.substring(1):line.substring(1, index);
					length = 0;
					}
				else length += line.length();
			lengths.put(id, new Integer(length));
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		trace("There are " + lengths.size() + " consensus sequence lengths");
		return lengths;
		}

	private static Hashtable<String, List<Hit>> getRetroHits() {
		Hashtable<String, List<Hit>> retroHits = new Hashtable<String, List<Hit>>();
		String[] retroHitFiles = {TEgff, ERVgff};
		for (int i=0; i<retroHitFiles.length; i++) addRetroHits(retroHits, retroHitFiles[i]);
		trace("There are " + retroHits.size() + " sequences with hits after TEs and all Retrovirus");
		return retroHits;
		}
		
	private static Hashtable<String, CS> getWantedCSs (Hashtable<String, List<Hit>> retroHits) {
		Hashtable<String, Integer> lengths = getLengths(CSFile);
		trace("There are " + lengths.size() + " sequence lengths from " + CSFile);
		removeFamilies(lengths, getFamilies(SSR));
		trace("There are " + lengths.size() + " sequence lengths after removing " + SSR + " sequences");
		Set<String> proteins = getFamilies(Proteins);
		proteins.removeAll(retroHits.keySet());
		removeFamilies(lengths, proteins);
		Hashtable<String, CS> wantedCS = new Hashtable<String, CS>();
		for (Enumeration<String> e=lengths.keys(); e.hasMoreElements();) {
			String key = e.nextElement();
			CS cs = new CS(key, lengths.get(key));
			if (retroHits.containsKey(key)) cs.addHits(retroHits.get(key));
			wantedCS.put(key, cs);
			}
		trace("There are " + wantedCS.size() + " wanted consensus sequences"); 
		return wantedCS;
		}


	private static void removeFamilies (Hashtable<String, Integer> lengths, Set<String> set) {
		for (Iterator iter=set.iterator(); iter.hasNext();) lengths.remove(iter.next());
		}

	private static void setSineMinCoverage (String[] args) {
		if (args != null && args.length > 0)
			try {sineMinCoverage = Double.parseDouble(args[0]);}
			catch (NumberFormatException ne) {System.out.println("Could not parse: " + args[0]);} 
		}
		
	private static void trace (String text) {if (debug) System.out.println(text);}
								
	private static void writeConsensusSequences (BufferedWriter out, boolean headersOnly) {
		//get wanted
		Hashtable<String, CS> wantedCS = getWantedCSs(getRetroHits());
		//get IR annots for wanted
		addIRHits(wantedCS);
		//process wanted and output them to library
		Hashtable<String, String> retroAnno = getAllRetroAnnotations();
		Set<String> irs = getFamilies(IRM);
		RBClassifier classifier = new RBClassifier();
		for (Enumeration<CS> e=wantedCS.elements(); e.hasMoreElements();) 
			e.nextElement().annotateAndClassify(retroAnno, irs, classifier, sineMinCoverage, 
			restMinCoverage);
		if (headersOnly) writeHeaders(out, wantedCS);
		else writeWantedSequences(out, CSFile, wantedCS);
		}

	private static void writeHeaders(BufferedWriter out, Hashtable<String, CS> wantedCS) {
		try {
			for (Enumeration<CS> e=wantedCS.elements(); e.hasMoreElements();) 
				out.write(e.nextElement().getFastaHeader());			
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
				
	private static void writeSatellite(BufferedWriter out, String inFile, String id, boolean headersOnly) 
		{
		try {
			out.write(">" + id + "\n");
			if (!headersOnly) {
				BufferedReader in = new BufferedReader(new FileReader(inFile));
				in.readLine();
				String line = null;
				while ((line = in.readLine()) != null) out.write(line + "\n");
				in.close();
				}
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
		
	private static void writeWantedSequences (BufferedWriter out, String seqFile, 
						  Hashtable<String, CS> wantedCS) {
	    try {
		BufferedReader in = new BufferedReader(new FileReader(seqFile));
			String line = null;
			boolean wanted = false;
			while ((line = in.readLine()) != null) 
			    if (line.charAt(0) == '>') {
				int index = line.indexOf(" ");
				String name = index<0?line.substring(1):line.substring(1,index);
				wanted = wantedCS.containsKey(name);
				if (wanted) out.write(wantedCS.get(name).getFastaHeader());
			    }
			    else if (wanted) out.write(line + "\n");
		}
		catch (IOException ie) {ie.printStackTrace();}
	}		
}

