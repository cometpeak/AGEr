package org.gersteinlab.age;

import java.util.ArrayList;

public class AGEAlign {
	public final static short DEFAULT_MATCH      = 1;
	public final static short DEFAULT_MISMATCH   = -2;
	public final static short DEFAULT_GAP_OPEN   = -2;
	public final static short DEFAULT_GAP_EXTEND = -1;
	
	public final static String AGE_VERSION = "AGE Java 0.0.1";
	
	public final static int MAX_ARGS = 1024;
	
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
	
	public static void main(String[] args)
	{
		int flag = 0;
		
		if (args.length == 1) {
			if (args[0] == "-version") {
				System.out.println(AGE_VERSION);
				return;
			}
		}
		
		if (args.length < 2) {
			printUsage();
			return;
		}
		
		ArrayList<String> ageArgs = new ArrayList<String>();
		BufferedReader in = null;
		
		// XXX
		if (args.length > 1) {
			for (int i = 1; i < args.length; i++) {
				ageArgs.add(args[i]);
			}
		} else {
			// Reading arguments from input
			try {
				in = new BufferedReader(new InputStreamReader(System.in));
				String s;
				while ((s = in.readLine()) != null) {
					
				}

		}
		
		int numArgs = ageArgs.size();
		while (numArgs > 0) {
			int start1 = -1;
			int start2 = -1;
			int end1   = -1;
			int end2   = -2;
			
			boolean revcom1 = false;
			boolean revcom2 = false;
			boolean both    = false;
			boolean error   = false;
			boolean version = false;
			
			short match     = DEFAULT_MATCH;
			short mismatch  = DEFAULT_MISMATCH;
			short gapOpen   = DEFAULT_GAP_OPEN;
			short gapExtend = DEFAULT_GAP_EXTEND;
			
			for (int i = 0; i < numArgs; i++) {
				String arg = ageArgs.get(i);
				if (arg.charAt(0) == '-') {
					int kvIndex;
					if ((kvIndex = arg.indexOf('=')) > -1) {
						if (arg.length() <= kvIndex + 1)
							continue;
						
						String key = arg.substring(1, kvIndex - 1);
						String value = arg.substring(kvIndex + 1);
						
						if (key.equals("match")) {
							match = Short.parseShort(value);
						} else if (key.equals("mismatch")) {
							mismatch = Short.parseShort(value);
						} else if (key.equals("go")) {
							gapOpen = Short.parseShort(value);
						} else if (key.equals("ge")) {
							gapExtend = Short.parseShort(value);
						} else if (key.equals("coor1")) {
							int rIndex;
							if ((rIndex = value.indexOf('-')) == -1) {
								System.err.println("Invalid format for coor1: " + value + ".");
								printUsage();
								return;
							}
							start1 = Integer.parseInt(value.substring(0, rIndex - 1));
							end1   = Integer.parseInt(value.subString(rIndex + 1));
							if (start1 == -1 || end1 == -1) {
								System.err.println("Invalid format for coor1: " + value + ".");
								printUsage();
								return;
							}
						} else if (key.equals("coor2")) {
							int rIndex;
							if (value.indexOf('-') == -1) {
								System.err.println("Invalid format for coor2: " + value + ".");
								printUsage();
								return;
							}
							start2 = Integer.parseInt(value.substring(0, rIndex - 1));
							end2   = Integer.parseInt(value.substring(rIndex + 1));
							if (start1 == -1 || end1 = -1) {
								System.err.println("Invalid format for coor2 " + value + ".");
								printUsage();
								return;
							}
						} else {
							System.err.println("Unkonwn option " + arg + ".");
							printUsage();
							return;
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
							return;
						}
					}
					
					continue;
				}
				
				
			}
		}
	}
}
