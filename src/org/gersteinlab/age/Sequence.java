package org.gersteinlab.age;

import java.io.*;

public class Sequence {
	private String _seq;
	private String _name;
	private int _start;
	private boolean _reverse;
	
	/* Next and previous sequences */
	private Sequence _next;
	private Sequence _prev;
	
	/**
	 * 
	 * @param seq
	 */
	public Sequence(String seq)
	{
		_seq = seq;
		_name = "";
		_start = 1;
		_reverse = false;
		_next = null;
		_prev = null;
		
		String ret = "";
		for (int i = 0; i < _seq.length(); i++) {
			char c = _seq.charAt(i);
			if (c == ' ')
				continue;
			else if (isNuc(c))
				ret = ret + c;
			else if (!isGap(c)) {
				System.err.println("Unrecognized nucleotile " + c + " in sequence.");
				_seq = "";
				return;
			}
		}
		
		_seq = ret;
	}
	
	/**
	 * 
	 * @param file
	 * @param start
	 * @param end
	 * @throws IOException
	 */
	public Sequence(String file, int start, int end) 
		throws IOException
	{
		if (start > 0 && end > 0 && end < start) {
			System.err.println("Invalid coordinates (" + start + ", " + end + ")." + 
		                       "Start is larger than end.");
			System.err.println("No sequence read!");
			return;
		}
		
		if (start >= 0)
			_start = start;
		
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line;
		
		while ((line = in.readLine()) != null)
			if (line.length() > 0)
				break;
		
		if (line.charAt(0) != '>') {
			System.err.println("Cannot find fasta header in file " + file);
			return;
		}
		
		_name = line.substring(1);
		
		int numPassed = 0;
		char c;
		while ((line = in.readLine()) != null) {
			String checked = "";
			for (int i = 0; i < line.length(); i++) {
				c = line.charAt(i);
				if (c == ' ')
					continue;
				else if (isNuc(c))
					checked = checked + c;
				else if (!isGap(c)) {
					System.err.println("Unrecognized nucleotide " + c + " in file " + file);
					_seq = "";
					return;
				}
			}
			
			int numChecked = checked.length();
			if (numChecked == 0)
				break;
			
			int newPassed = numPassed + numChecked;
			if (newPassed < start) {
				// Haven't reached start yet
				// Don't do anything
			} else if (numPassed < start && newPassed >= end) {
				// Spans whole region
				_seq = _seq + checked.substring(start - numPassed - 1, end - start + 1);
			} else if (numPassed < start && newPassed >= start) {
				// Spans start only
				_seq = _seq + checked.substring(start - numPassed - 1);
			} else if (end < start) {
				// No end
				_seq = _seq + checked;
			} else if (newPassed < end) {
				// Not reached end yet
				_seq = _seq + checked;
			} else if (numPassed < end && newPassed >= end) {
				// Spans end only
				_seq = _seq + checked.substring(0, end - numPassed);
			}
			
			numPassed = newPassed;
			line = "";
		}
		
		in.close();
	}
	
	/**
	 * 
	 * @param seq
	 * @param name
	 * @param start
	 * @param reverse
	 */
	public Sequence(String seq, String name, int start, boolean reverse)
	{
		_seq = seq;
		_name = name;
		_start = start;
		_reverse = reverse;
		_next = null;
		_prev = null;
	}
	
	/**
	 * Clone method.
	 * @return Sequence
	 */
	public Sequence clone()
	{
		return new Sequence(_seq, _name, _start, _reverse);
	}
	
	/**
	 * Getter for name
	 * @return name
	 */
	public String name()
	{
		return _name;
	}
	
	/**
	 * Getter for sequence
	 * @return sequence
	 */
	public String sequence()
	{
		return _seq;
	}
	
	/**
	 * Getter for start
	 * @return start position for sequence
	 */
	public int start()
	{
		return _start;
	}
	
	/**
	 * Getter for reverse flag
	 * @return whether sequence has been reversed
	 */
	public boolean reverse()
	{
		return _reverse;
	}
	
	/**
	 * Getter for next sequence
	 * @return next Sequence
	 */
	public Sequence next()
	{
		return _next;
	}
	
	/**
	 * Getter for previous Sequence
	 * @return previous Sequence
	 */
	public Sequence prev()
	{
		return _prev;
	}
	
	/**
	 * Returns complement of given nucleotide
	 * @param c
	 * @return complement of c
	 */
	public static char complement(char c)
	{
		if      (c == 'a') return 't';
		else if (c == 'c') return 'g';
		else if (c == 't') return 'a';
		else if (c == 'g') return 'c';
		else if (c == 'A') return 'T';
		else if (c == 'C') return 'G';
		else if (c == 'T') return 'A';
		else if (c == 'G') return 'C';
		else if (c == 'u') return 'a';
		else if (c == 'U') return 'A';
		
		return c;
	}
	
	/**
	 * Checks if character is a nucleotide
	 * @param c
	 * @return true if nucleotide, false otherwise
	 */
	public static boolean isNuc(char c)
	{
		if (c == 'A' || c == 'a' ||
			c == 'T' || c == 't' ||
			c == 'C' || c == 'c' ||
			c == 'G' || c == 'g' ||
			c == 'U' || c == 'u' ||
			c == 'X' || c == 'x' ||
			c == 'N' || c == 'n' ||
			c == 'R' || c == 'r' ||
			c == 'Y' || c == 'y' ||
			c == 'M' || c == 'm' ||
			c == 'K' || c == 'k' ||
			c == 'S' || c == 's' ||
			c == 'W' || c == 'w' ||
			c == 'B' || c == 'b' ||
			c == 'D' || c == 'd' ||
			c == 'H' || c == 'h' ||
			c == 'V' || c == 'v')
			return true;
		
		return false;
	}
	
	/**
	 * Character for gap
	 * @return
	 */
	public static char gap()
	{
		return '-';
	}
	
	/**
	 * 
	 * @param c
	 * @return
	 */
	public static boolean isGap(char c)
	{
		if (c == '-')
			return true;
		if (c == ' ')
			return true;
		
		return false;
	}
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean sameNuc(char a, char b)
	{
		char aa = Character.toUpperCase(a);
		char bb = Character.toUpperCase(b);
		
		if ((aa == 'A' || aa == 'C' || aa == 'T' || aa == 'G' || aa == 'U') &&
			(bb == 'A' || bb == 'C' || bb == 'T' || bb == 'G' || bb == 'U') && aa == bb)
			return true;
		
		return false;
	}
	
	/**
	 * Generates the reverse complement of the sequence, flips the _reverse
	 * flag, and stores the reverse complement as the object's sequence
	 */
	public void revcom()
	{
		if (_reverse)
			_start -= _seq.length() - 1;
		else
			_start += _seq.length() - 1;
		
		_reverse = !_reverse;
		String ret = "";
		for (int i = _seq.length() - 1; i >= 0; i--)
			ret  = ret + complement(_seq.charAt(i));
		
		_seq = ret;
	}
	
	/**
	 * Returns a sub-Sequence with the given start and end coordinates. Note 
	 * that we are using 1-based coordinates
	 * @param start
	 * @param end
	 * @return
	 */
	public Sequence substr(int start, int end)
	{
		if (start > 0 && end > 0 && start > end)
			return null;
		if (start <= 0)
			start = 1;
		if (end <= 0 || end > _seq.length())
			end = _seq.length();
		
		int s = start;
		int e = end;
		
		if (_reverse) {
			s = _seq.length() - e + 1;
			e = _seq.length() - s + 1;
		}
		
		String tmp = "";
		if (_seq.length() > 0)
			tmp = _seq.substring(s - 1, e - s + 1);
		
		return new Sequence(tmp, _name, _start + start - 1, _reverse);
	}
	
	/**
	 * Add a sequence before the current sequence
	 * @param newSequence
	 * @return
	 */
	public boolean addBefore(Sequence newSequence)
	{
		// Check if input is good
		if (newSequence == null)
			return false;
		
		// Check if the input is the same object
		if (newSequence == this)
			return false;
		
		// Check if the call comes from addAfter
		if (newSequence.next() == this && _prev == null) {
			_prev = newSequence;
			return true;
		}
		
		// Check if they are already paired
		if (newSequence.prev() == this && _next == newSequence)
			return true;
		
		// Check if it can be added
		if (_prev != null || newSequence.next() != null || newSequence.prev() != null)
			return false;
		
		// Set previous
		_prev = newSequence;
		
		// Update next for newSequence
		if (newSequence.addAfter(this) == false) {
			_prev = null;
			return false;
		}
		
		return true;
	}
	
	/**
	 * Adding atom after
	 * XXX Examine the relationship between this and addBefore
	 * @param sequence
	 * @return
	 */
	public boolean addAfter(Sequence newSequence)
	{
		// Check if not null
		if (newSequence == null)
			return false;
		
		// Check if same object
		if (newSequence == this)
			return true; // XXX Shouldn't this be false???
		
		// Check if the call comes from addAfter
		if (newSequence.prev() == this && _next == null) {
			_next = newSequence;
			return true;
		}
		
		// Check if they are already paired
		if (newSequence.prev() == this && _next == newSequence)
			return true;
		
		// Check if it can be added
		if (_next != null || newSequence.next() != null || newSequence.prev() != null)
			return false;
		
		_next = newSequence;
		
		if (newSequence.addBefore(this) == false) {
			_next = null;
			return false;
		}
		
		return true;
		
	}
	
	public static Sequence parseSequences(String file) 
		throws IOException
	{
		Sequence first = null;
		Sequence last  = null;
		BufferedReader in = null;
		
		try {
			in = new BufferedReader(new FileReader(file));
		} catch (IOException e) {
			System.err.println("Cannot open input file " + file);
			return null;
		}
		
		String line;
		String name = "";
		String seq = "";
		
		while ((line = in.readLine()) != null) {
			int n = line.length();
			if (n == 0)
				continue;
			
			if (line.charAt(0) == '>') {
				if (name.length() > 0) {
					Sequence s = new Sequence (seq, name, 1, false);
					if (first != null) {
						last.addAfter(s);
						last = s;
					} else {
						first = last = s;
					}
				}
				
				line = line.substring(1);
				seq = "";
				continue;
			}
			
			String checked = "";
			for (int i = 0; i < n - 1; i++) {
				char c = line.charAt(i);
				if (isNuc(c))
					checked = checked + c;
				else if (!isGap(c)) {
					System.err.println("Unrecognized nucleotide " + c + " in file " + file);
					seq = "";
					name = "";
				}
			}
			
			seq = seq + checked;
		}
		in.close();
		
		if (name.length() > 0) {
			Sequence s = new Sequence(seq, name, 1, false);
			if (first != null) {
				last.addAfter(s);
				last = s;
			} else {
				first = last = s;
			}
		}
		
		return first;
	}
	
	// Destructor method??
}
