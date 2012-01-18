package org.gersteinlab.age;

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.io.*;

public class AGEAlign {
	public final static short DEFAULT_MATCH      = 1;
	public final static short DEFAULT_MISMATCH   = -2;
	public final static short DEFAULT_GAP_OPEN   = -2;
	public final static short DEFAULT_GAP_EXTEND = -1;
	
	public final static String AGE_VERSION = "AGE Java 0.0.1";
	
	public final static int MAX_ARGS = 1024;
	
	public static int flag;
	
	public static int start1;
	public static int start2;
	public static int end1;
	public static int end2;
	
	public static boolean revcom1;
	public static boolean revcom2;
	public static boolean both;
	public static boolean error;
	public static boolean version;
	
	public static short match;
	public static short mismatch;
	public static short gapOpen;
	public static short gapExtend;
	
	public static String file1 = null;
	public static String file2 = null;
	
	public static void printUsage()
	{
		System.out.println("Usage:");
		System.out.println("\tagealign");
		System.out.println("\t\t[-version]");
		System.out.println("\t\t[-indel|-tdup|-inv|-invl|-invr]");
		System.out.println("\t\t[-match=value:" + DEFAULT_MATCH + "]");
		System.out.println("\t\t[-mismatch=value:" + DEFAULT_MISMATCH + "]");
		System.out.println("\t\t[-mismatch=value:" + DEFAULT_MISMATCH + "]");
		System.out.println("\t\t[-go=value:" + DEFAULT_GAP_OPEN + "]");
		System.out.println("\t\t[-ge=value:" + DEFAULT_GAP_EXTEND + "]");
		System.out.println("\t\t[-both] [-revcom1] [-revcom2]");
		System.out.println("\t\t[-coor1=start-end] [-coor2=start-end]");
		System.out.println("\t\t file1.fa file2.fa");
	}
	
	public static void processArgs(ArrayList<String> ageArgs)
	{
		int numArgs = ageArgs.size();
		
		for (int i = 0; i < numArgs; i++) {
			String arg = ageArgs.get(i);
			if (arg.charAt(0) == '-') {
				int kvIndex;
				if ((kvIndex = arg.indexOf('=')) > -1) {
					if (arg.length() <= kvIndex + 1)
						continue;
					
					String key = arg.substring(1, kvIndex);
					String value = arg.substring(kvIndex + 1);
					
					if (key.equals("match")) {
						try {
							match = Short.parseShort(value);
						} catch (NumberFormatException nfe) {
							System.err.println("Error in match specification " + value + ".");
							System.err.println("Using default value " + match + ".");
						}
					} else if (key.equals("mismatch")) {
						try {
							mismatch = Short.parseShort(value);
						} catch (NumberFormatException nfe) {
							System.err.println("Error in mismatch specification " + value + ".");
							System.err.println("Using default value " + mismatch + ".");
						}
					} else if (key.equals("go")) {
						try {
							gapOpen = Short.parseShort(value);
						} catch (NumberFormatException nfe) {
							System.err.println("Error in gap open specification "  + value + ".");
							System.err.println("Using default value " + gapOpen + ".");
						}
					} else if (key.equals("ge")) {
						try {
							gapExtend = Short.parseShort(value);
						} catch (NumberFormatException nfe) {
							System.err.println("Error in gap extend specification " + value + ".");
							System.err.println("Using default value " + gapExtend + ".");
						}
					} else if (key.equals("coor1")) {
						int rIndex;
						if ((rIndex = value.indexOf('-')) == -1) {
							System.err.println("Invalid format for coor1: " + value + ".");
							printUsage();
							error = true;
							break;
						}
						start1 = Integer.parseInt(value.substring(0, rIndex - 1));
						end1   = Integer.parseInt(value.substring(rIndex + 1));
						if (start1 == -1 || end1 == -1) {
							System.err.println("Invalid format for coor1: " + value + ".");
							printUsage();
							error = true;
							break;
						}
					} else if (key.equals("coor2")) {
						int rIndex;
						if ((rIndex = value.indexOf('-')) == -1) {
							System.err.println("Invalid format for coor2: " + value + ".");
							printUsage();
							error = true;
							break;
						}
						start2 = Integer.parseInt(value.substring(0, rIndex - 1));
						end2   = Integer.parseInt(value.substring(rIndex + 1));
						if (start1 == -1 || end1 == -1) {
							System.err.println("Invalid format for coor2 " + value + ".");
							printUsage();
							error = true;
							break;
						}
					} else {
						System.err.println("Unknown option " + arg + ":" + key + ".");
						printUsage();
						error = true;
						break;
					}
				} else {
					String key = arg.substring(1);
					
					if (key.equals("version")) {
						version = true;
					} else if (key.equals("indel")) {
						flag |= AGEAligner.INDEL_FLAG;
					} else if (key.equals("inv")) {
						flag |= AGEAligner.INVERSION_FLAG;
					} else if (key.equals("invl")) {
						flag |= AGEAligner.INVL_FLAG;
					} else if (key.equals("invr")) {
						flag |= AGEAligner.INVR_FLAG;
					} else if (key.equals("tdup")) {
						flag |= AGEAligner.TDUPLICATION_FLAG;
					} else if (key.equals("revcom1")) {
						revcom1 = true;
					} else if (key.equals("revcom2")) {
						revcom2 = true;
					} else if (key.equals("both")) {
						both = true;
					} else {
						System.err.println("Unkonwn option " + arg + ".");
						printUsage();
						error = true;
						break;
					}
				}
				continue;
			}

			if (file1 == null)
				file1 = ageArgs.get(i);
			else if (file2 == null)
				file2 = ageArgs.get(i);
			else {
				System.err.println("Too many input files");
				printUsage();
				error = true;
				break;
			}
		}
	}
	
	public static void main(String[] args) throws IOException
	{
		flag = 0;
		
		if (args.length == 1) {
			if (args[0].equals("-version")) {
				System.out.println(AGE_VERSION);
				return;
			}
		}
		
		ArrayList<String> ageArgs = new ArrayList<String>();
		BufferedReader in = null;
		
		if (args.length > 1) {
			for (int i = 1; i < args.length; i++) {
				ageArgs.add(args[i]);
			}
		} else {
			// Reading arguments from input
			try {
				in = new BufferedReader(new InputStreamReader(System.in));
				String line = in.readLine();
				if (line == null) {
					printUsage();
					return;
				}
				do {
					if (line.length() > 0) {
						StringTokenizer st = new StringTokenizer(line, " ");
						while (st.hasMoreTokens()) {
							ageArgs.add(st.nextToken());
						}
						break;
					}
				} while ((line = in.readLine()) != null);
			} catch (IOException e) {
				System.err.println("Cannot open System.in.");
				return;
			}
		}

		String oldFile1 = null;
		String oldFile2 = null;
		Sequence oldSeq1 = null;
		Sequence oldSeq2 = null;
		
		int numArgs = ageArgs.size();
		while (numArgs > 0) {
			start1 = -1;
			start2 = -1;
			end1   = -1;
			end2   = -1;
			
			revcom1 = false;
			revcom2 = false;
			both    = false;
			error   = false;
			version = false;
			
			match     = DEFAULT_MATCH;
			mismatch  = DEFAULT_MISMATCH;
			gapOpen   = DEFAULT_GAP_OPEN;
			gapExtend = DEFAULT_GAP_EXTEND;

			file1 = null;
			file2 = null;
			
			processArgs(ageArgs);
				
			if (version) {
				System.out.println(AGE_VERSION);
				return;
			}

			if (flag == 0)
				flag = AGEAligner.INDEL_FLAG;

			int numModes = 0;
			if ((flag & AGEAligner.INDEL_FLAG) == 1)
				numModes++;
			if ((flag & AGEAligner.TDUPLICATION_FLAG) == 1)
				numModes++;
			if ((flag & AGEAligner.INVR_FLAG) == 1)
				numModes++;
			if ((flag & AGEAligner.INVL_FLAG) == 1)
				numModes++;
			if ((flag & AGEAligner.INVERSION_FLAG) == 1)
				numModes++;

			if (numModes != 1) {
				System.err.print("Error in mode specification. ");
				if (numModes == 0)
					System.err.println("No mode is specified.");
				else if (numModes > 1)
					System.err.println("More than one mode is specified");
				error = true;
			}

			if (file1 != null && file2 == null)
				System.err.println("No second file is given.");
				
			System.out.println("Configuration:");
			System.out.println("\tS1: [" + start1 + ", " + end1 + ", " + revcom1 + "]");
			System.out.println("\tS2: [" + start2 + ", " + end2 + ", " + revcom2 + "]");
			System.out.println("\tFile1: " + file1);
			System.out.println("\tFile2: " + file2);
			System.out.println("\tMatch: " + match + "\n\tMismatch: " + mismatch);
			System.out.println("\tGap open: " + gapOpen + "\n\tGap extend: " + gapExtend);
			System.out.println();

			if (!error && file1 != null && file2 != null) {
				if (file2 != oldFile2) {
					oldSeq2 = Sequence.parseSequences(file2);
					oldFile2 = file2;
				}
				if (file1 != oldFile1) {
					oldSeq1 = Sequence.parseSequences(file1);
					oldFile1 = file1;
				}

				for (Sequence s1 = oldSeq1; s1 != null; s1 = s1.next()) {
					for (Sequence s2 = oldSeq2; s2 != null; s2 = s2.next()) {
						Sequence seq1 = s1.substr(start1, end1);
						Sequence seq2 = s2.substr(start2, end2);

						if (revcom1)
							seq1.revcom();
						if (revcom2)
							seq2.revcom();

						seq1 = new Sequence("ACTGGTGTCAACTG");
						seq2 = new Sequence("ACTGACTG");

						System.out.println("seq1.revcom, seq.revcom completed if necessary");

						Scorer scorer = new Scorer(match, mismatch, gapOpen, gapExtend);
						if (both) {
							Sequence seq3 = seq2.clone();
							seq3.revcom();

							System.out.println("seq3.revcom completed");
							
							AGEAligner aligner1 = new AGEAligner(seq1, seq2);
							AGEAligner aligner2 = new AGEAligner(seq1, seq3);

							System.out.println("Aligners initialized");
							// f00

							boolean res1 = aligner1.align(scorer, flag);
							boolean res2 = aligner2.align(scorer, flag);
							
							if (res1 == false && res2 == false)
								System.err.println("No alignment made.");
							else if (aligner1.score() >= aligner2.score())
								aligner1.printAlignment();
							else
								aligner2.printAlignment();
							seq3 = null;
						} else {
							AGEAligner aligner = new AGEAligner(seq1, seq2);
							if (aligner.align(scorer, flag) == true)
								aligner.printAlignment();
							else
								System.err.println("No alignment made");
						}
					}
				}
			}
			
			numArgs = 0;
			ageArgs = new ArrayList<String>();
			if (in != null) {
				String line;
				try {
					while ((line = in.readLine()) != null) {
						if (line.length() > 0) {
							StringTokenizer st = new StringTokenizer(line, " ");
							while (st.hasMoreTokens()) {
								ageArgs.add(st.nextToken());
							}
							break;
						}
					}
				} catch (IOException e) {
					System.err.println("Cannot read from System.in.");
					return;
				}
				numArgs = ageArgs.size();
			}
		}
	}
}
