import java.io.*;
import java.util.*;

/*************************
ClassifyConsensusSequences 
	Classifies sequences based on their censor output using the mam and our_known_reps_20130520.fasta libraries
	Author: Joy Raison
	Data: September 4, 2014
	Updated: September 5, 2014 To include our_known_reps_20130520.fasta library
	Updated:September 9, 2014 To print the notKnown.fa.gff file
	Inputs: The consensus sequences (fasta format)
		The repbase mam library (fasta format)
		The our_known_reps library (fasta format)
		The map file from the censor run
	Outputs: check.txt (for checking individual sequence mappings)
		known.txt (list of "identified" sequences with name of the library sequence they match)
		partial.txt (list of partially matched sequences with a name list of the library sequences they match)
		notKnown.fa (fasta file of the not known sequences, including seqs with hits to CENSOR or not) 
		notKnown.fa.gff (only include family names have hit(s) with CENSOR) 		
***************************/


public class ClassifyConsensusSequences {
	private static String lFile = "./Vertebrate_use.fa";
	private static String bDir = "./";
	private static String dDir = bDir + "results_classify/";
	private static String cFile = bDir + "ConsensusSequences.fa";
	private static String map = bDir + "ConsensusSequences.fa.map";
	private static String olFile = "./our_known_reps_20130520.fasta";
	private static String known = dDir + "known.txt";
	private static String partial = dDir + "partial.txt";
	private static String check = dDir + "check.txt";
	private static String notKnown = dDir + "notKnown.fa";
	private static String gff = dDir + "notKnown.fa.gff";

	public static void main (String[] args) {
		Hashtable<String, Integer> cLengths = getLengths(cFile);
		String[] libs = {lFile, olFile};
		Hashtable<String, Integer> libLengths = getLibraryLengths(libs);
		Set<String> knowns = classify (map, cLengths, libLengths, known, partial, check, gff);
		printNotKnownSequences (cFile, knowns, notKnown);
		}


	private static Set<String> classify (String mapFile, Hashtable<String, Integer> btL, 
		Hashtable<String, Integer> libL, String kOutFile, String pOutFile, String cOutFile, String gffFile) {
		Hashtable<String, String> kMaps = new Hashtable<String, String>();
		Hashtable<String, String> pMaps = new Hashtable<String, String>();
		String line = null;
		try {
			BufferedWriter cOut = new BufferedWriter(new FileWriter(cOutFile));
			BufferedWriter gOut = new BufferedWriter(new FileWriter(gffFile));
			BufferedReader in = new BufferedReader(new FileReader(map));
			while ((line = in.readLine()) != null) processMapping(btL, libL, kMaps, pMaps, line, cOut, gOut);
			in.close();
			cOut.close();
			gOut.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		System.out.println("" + kMaps.size() + " consensus families have known mapping");
		printMappings(kMaps, kOutFile);
		System.out.println("" + pMaps.size() + " consensus families have partial mappings");
		printMappings(pMaps, pOutFile);
		return kMaps.keySet();
		}

	private static Hashtable<String, Integer> getLengths (String file) {
		Hashtable<String, Integer> lengths = new Hashtable<String, Integer>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			int count = 0;
			String line = null;
			String id = null;
			StringTokenizer st = null;
			while ((line = in.readLine()) != null) 
				if (line.length() > 0)
					if (line.charAt(0) == '>') {
						if (count > 0) {
							lengths.put(id, Integer.valueOf(count));
							count = 0;
							}
						st = new StringTokenizer(line);
						id = st.nextToken().substring(1);
						}
					else count += line.length();
			lengths.put(id, Integer.valueOf(count));
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		System.out.println("There are " + lengths.size() + " sequence lengths");
		return lengths;
		}

	private static Hashtable<String, Integer> getLibraryLengths (String[] libs) {
		Hashtable<String, Integer> lengths = new Hashtable<String, Integer>();
		for (int i=0; i<libs.length; i++) {
			Hashtable<String, Integer> fileLengths = getLengths(libs[i]);
			for (Enumeration<String> keys=fileLengths.keys(); keys.hasMoreElements();) {
				String key = keys.nextElement();
				lengths.put(key, fileLengths.get(key));
				}
			}
		System.out.println("There are " + lengths.size() + " library sequence lengths");
		return lengths;
		}
		
	private static void printMappings (Hashtable<String, String> map, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			out.write("Sequence MappedTo\n");
			String key = null;
			for (Enumeration<String> e= map.keys(); e.hasMoreElements();) {
				key = e.nextElement();
				out.write(key + " " + map.get(key) + "\n");
				}
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}


	private static void printNotKnownSequences(String seqFile, Set<String> notWanted, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			BufferedReader in = new BufferedReader(new FileReader(seqFile));
			boolean wanted = false;
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.charAt(0) == '>') wanted = !notWanted.contains(line.substring(1).split(" ")[0]);
				if (wanted) out.write(line + "\n");
				}
			in.close();
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
						
	private static void processMapping (Hashtable<String, Integer> btL, Hashtable<String, Integer> libL,
		Hashtable<String, String> kMaps, Hashtable<String, String> pMaps, String line, BufferedWriter cOut,
		BufferedWriter gOut) {
		StringTokenizer st = new StringTokenizer(line);
		String bt = st.nextToken();
		try {
			int bStart = Integer.parseInt(st.nextToken());
			int bEnd = Integer.parseInt(st.nextToken());
			String lib = st.nextToken();
			int lStart = Integer.parseInt(st.nextToken());
			int lEnd = Integer.parseInt(st.nextToken());
			int bl = 0;
			if (btL.containsKey(bt)) bl = btL.get(bt).intValue();
			else System.out.println("Length of " + bt + " could not be found");
			int ll = 0;
			if (libL.containsKey(lib)) ll = libL.get(lib).intValue();
			else System.out.println("Length of " + lib + " could not be found");
			if (bl > 0 && ll > 0) {
				int bpc = Math.round (100f * (bEnd - bStart + 1) / bl);
				int lspc = Math.round (100f * lStart / ll);
				int lepc = Math.round (100f * lEnd / ll);
				cOut.write(bt + " " + bStart + " " + bEnd + " (" + bpc + ") " + lib + " "  + lStart + " " + lEnd +
					" (" + lspc + "," + lepc + ")\n");
				if (bpc >= 85 && lspc <= 5 && lepc >= 95) kMaps.put(bt, lib); 
				else {
					if (!pMaps.containsKey(bt)) pMaps.put(bt, lib);
					else {
						//System.out.println("bt: " + bt);
						String name = pMaps.get(bt);
						name += ":" + lib;
						pMaps.put(bt, name);		
						}	
					String strand = st.nextToken().charAt(0) == 'c'?"-":"+";
					for (int i=1; i<3; i++) st.nextToken();
					gOut.write(bt + "\tcensor\trepeat\t" + bStart + "\t" + bEnd + "\t" + st.nextToken() + "\t" + 
						strand + "\t.\tRepeat " + lib + " . " + lStart + " " + lEnd + " " + (ll - lEnd) + "\n");
					}		 
				} 
			}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);} 
		catch (IOException ie) {ie.printStackTrace();}
		}

	}
