package org.gersteinlab.age;

public class Util {
	
	/**
	 * Pads an input string with leading spaces to fit a set field width.
	 * 
	 * XXX Change to use StringBuilder. It would be more efficient that way. 
	 * Even better, look into java.text.MessageFormat
	 * 
	 * @param str - string
	 * @param n - field width
	 * @return - input string padded with spaces
	 */
	public static String addw(int n, String str)
	{
		String padding = "";
		for (int i = str.length(); i < n; i++) {
			padding += " ";
		}
		
		return padding + str;
	}
	
	/**
	 * Overloaded addw used for ints
	 * 
	 * @param n
	 * @param num
	 * @return String containing integer padded with appropriate spaces
	 */
	public static String addw(int n, int num)
	{
		return Util.addw(n, Integer.toString(num));
	}
}
