package org.gersteinlab.age;

public class AliFragment {
	
	private AliFragment _next;
	private AliFragment _prev;
	private String _ali1;
	private String _ali2;
	private int _start1;
	private int _start2;
	private int _end1;
	private int _end2;
	
	public static final int WIDTH = 60;
	public static final int MARGIN = 9;
	
	public AliFragment(String ali1, String ali2, 
	                   int start1, int start2, 
	                   int end1, int end2)
	{
		_next = null;
		_prev = null;
		_ali1 = ali1;
		_ali2 = ali2;
		_start1 = start1;
		_start2 = start2;
		_end1 = end1;
		_end2 = end2;
	}
	
	public AliFragment next(AliFragment b)
	{
		if (_next != null)
			return null;
		
		_next = b;
		b.prev(this);
		return _next;
	}
	
	public AliFragment prev(AliFragment b)
	{
		if (_prev != null)
			return null;
		
		_prev = b;
		b.next(this);
		return _prev;
	}
	
	public AliFragment next()
	{
		return _next;
	}
	
	public AliFragment prev()
	{
		return _prev;
	}
	
	public int start1()
	{
		return _start1;
	}
	
	public int start2()
	{
		return _start2;
	}
	
	public int end1()
	{
		return _end1;
	}
	
	public int end2()
	{
		return _end2;
	}
	
	public AliCount countAligned()
	{
		AliCount count = new AliCount(0, 0, 0);
		int len = _ali1.length();
		
		if (_ali2.length() < len)
			len = _ali2.length();
		
		for (int i = 0; i < len; i++) {
			count.nAli++;
			if (Sequence.sameNuc(_ali1.charAt(i), _ali2.charAt(i)) != false)
				count.nId++;
			if (Sequence.isGap(_ali1.charAt(i)) || Sequence.isGap(_ali2.charAt(i)))
				count.nGap++;
		}
		
		return count;
	}
	
	public void printAlignment()
	{
		int inc1 = 1;
		int inc2 = 1;
		
		if (_end1 < _start1)
			inc1 = -1;
		if (_end2 < _start2)
			inc2 = -1;
		
		StringBuffer margin = new StringBuffer();
		for (int i = 0; i <= MARGIN; i++)
			margin.append(" ");
		
		int n = _ali1.length();
		if (_ali2.length() < n)
			n = _ali2.length();
		int ind1 = _start1;
		int ind2 = _start2;
		for (int i = 0; i < n; i += WIDTH) {
			int st1 = ind1;
			int st2 = ind2;
			String a1 = _ali1.substring(i, WIDTH);
			String a2 = _ali2.substring(i, WIDTH);
			int nuc1 = 0;
			int nuc2 = 0;
			StringBuffer match = new StringBuffer();
			for (int j = 0; j < a1.length(); j++) {
				if (Sequence.sameNuc(a1.charAt(j), a2.charAt(j)))
					match.append("|");
				else if (!Sequence.isGap(a1.charAt(j)) && !Sequence.isGap(a2.charAt(j)))
					match.append(".");
				else
					match.append(" ");
				
				if (!Sequence.isGap(a1.charAt(j))) {
					nuc1++;
					ind1 += inc1;
				}
				if (!Sequence.isGap(a2.charAt(j))) {
					nuc2++;
					ind2 += inc2;
				}
			}
			
			System.out.println();
			
			if (nuc1 > 0) {
				System.out.print(Util.addw(MARGIN, st1) + " " + a1);
				System.out.println(" " + (ind1 - inc1));
			} else {
				System.out.println(margin.toString() + a1);
			}
			
			System.out.println(margin.toString() + match.toString());
			
			if (nuc2 > 0) {
				System.out.print(Util.addw(MARGIN, st2) + " " + a2);
				System.out.println(" " + (ind2 - inc2));
			} else {
				System.out.println(margin.toString() + a2);
			}
		}
	}
}
