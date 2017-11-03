import java.io.*;
import java.util.*;


public class IdentifySSRs {
	private static String dDir = "/data/rc003/lu/echidna/results_classify/";
	private static String inPrefix = dDir + "notKnownNotProtein";
	private static String pFile = inPrefix + ".phobos";
	
	
	public static void main (String[] args) {
		Hashtable<String, String> ssrs = getSSRs(pFile);
		if (ssrs.size() > 0) {
			outputSSRs (ssrs, dDir + "SSR.txt");
			outputUnknownFasta(ssrs, inPrefix + ".fa", inPrefix + "NotSSR.fa");
			}
		}
		

	private static double addFamilySSR (String name, BufferedReader in, Hashtable<String, String> ssrs) {
		String line = null;
		double maxCovered = 0.0;
		try {
			line = in.readLine();
			int length = Integer.parseInt(line.substring(line.indexOf(":")+2));
			double lowerLimit = length * .1;
			double upperLimit = length * .9;
			in.readLine();
			StringTokenizer st = null;
			while ((line = in.readLine()) != null && line.charAt(0) != '#') {
				st = new StringTokenizer(line);
				st.nextToken();
				int start = Integer.parseInt(st.nextToken());
				if (start <= lowerLimit) { 
					st.nextToken();
					int end = Integer.parseInt(st.nextToken());
					double covered = (end - start + 1) * 100./length;
					if (covered > maxCovered) maxCovered = covered;
					if (end >= upperLimit) {
						for (int i=1; i<15; i++) st.nextToken();
						ssrs.put(name, st.nextToken());
						}
					}
				}
			}
		catch (IOException ie) {ie.printStackTrace();}
		catch (NumberFormatException ne) {System.out.println("Could not parse: " + line);}
		return maxCovered;
		}
		 		
	private static Hashtable<String, String> getSSRs (String inFile) {
		Hashtable<String, String> ssrs = new Hashtable<String, String> ();
		double maxCovered = 0.0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			String line = null;
			double covered = 0.0;
			while ((line = in.readLine()) != null) 
				if (line.charAt(0) == '>') {
					covered = addFamilySSR(line.substring(1), in, ssrs); 
					if (covered > maxCovered) maxCovered = covered;
					}
			in.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		System.out.println("There are " + ssrs.size() + " families that are SSRs, maximum coverage was " + 
			maxCovered + " percent");
		return ssrs;
		}

	private static void outputSSRs (Hashtable<String, String> ssrs, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			out.write("Sequence MappedTo\n");
			for (Enumeration<String> e=ssrs.keys(); e.hasMoreElements();) {
				String key = e.nextElement();
				out.write(key + " " + ssrs.get(key) + "\n");
				}
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		}
		
	private static void outputUnknownFasta (Hashtable<String, String> ssrs, String inFile, String outFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			boolean wanted = false;
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.charAt(0) == '>') wanted = !ssrs.containsKey(line.substring(1));
				if (wanted) out.write(line + "\n");
				}
			in.close();
			out.close();
			}
		catch (IOException ie) {ie.printStackTrace();}
		} 
}
