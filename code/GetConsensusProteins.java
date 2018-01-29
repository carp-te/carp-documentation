import java.io.*;
import java.util.*;

public class GetConsensusProteins {
	private static String dir = "./";
	private static String protein = "protein.txt";
	private static String consensus = "ConsensusSequences.fa";
	private static String outFile = "Proteins.fa";
	
	public static void main (String[] args) {writeConsensus(getProteins());}
		
		private static Hashtable<String, String> getProteins () {
			Hashtable<String, String> proteins = new Hashtable<String, String>();
			try {
				BufferedReader in = new BufferedReader(new FileReader(dir + protein));
				in.readLine();
				String line = null;
				String[] fields = null;
				while ((line = in.readLine()) != null) {
					fields = line.split(" ");
					proteins.put(fields[0], fields[1]);
					}
				in.close();
				}
			catch (IOException ie) {ie.printStackTrace();}
			System.out.println("There are " + proteins.size() + 
				" consensus sequence families identified as proteins");
			return proteins;
			}
			
		private static void writeConsensus (Hashtable<String, String> proteins) {
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(dir + outFile));
				BufferedReader in = new BufferedReader(new FileReader(dir + consensus));
				boolean wanted = false;
				String line = null;
				String id = null;
				while ((line = in.readLine()) != null) {
					if (line.charAt(0) == '>') {
						id = line.substring(1);
						wanted = proteins.containsKey(id);
						if (wanted) line += " " + proteins.get(id);
						}
					if (wanted) out.write(line + "\n");
					}
				in.close();
				out.close();
				}
			catch (IOException ie) {ie.printStackTrace();}
			}
			
		}
