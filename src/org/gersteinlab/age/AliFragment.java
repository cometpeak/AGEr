package org.gersteinlab.age;

public class AliFragment {
	
	public class AliCount {
		
	}
	
	private AliFragment _next;
	private AliFragment _prev;
	private String _ali1;
	private String _ali2;
	private int _start1;
	private int _start2;
	private int _end1;
	private int _end2;
	
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
	
	public void countAligned()
	{
		
	}
}
