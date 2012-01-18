package org.gersteinlab.age;

public class Block {
	private Block _next;
	private Block _prev;
	private int _start1;
	private int _start2;
	private int _length;
	private String _desc;
	
	public Block(int st1, int st2, int l, String d)
	{
		_next = null;
		_prev = null;
		_start1 = st1;
		_start2 = st2;
		_length = l;
		_desc = d;
	}
	
	public Block(int st1, int st2, int l)
	{
		_next = null;
		_prev = null;
		_start1 = st1;
		_start2 = st2;
		_length = l;
		_desc = "";
	}
	
	public String toString()
	{
		StringBuffer ret = new StringBuffer(_desc);
		if (ret.length() > 0)
			ret.append(" ");
		
		ret.append("block\t");
		ret.append(_start1).append("\t");
		ret.append(_start2).append("\t");
		ret.append(_length);
		
		return ret.toString();
	}
	
	public Block next(Block b)
	{
		if (_next != null)
			return null;
		
		_next = b;
		b.prev(this);
		return _next;
	}
	
	public Block prev(Block b)
	{
		if (_prev != null)
			return null;
		
		_prev = b;
		b.next(this);
		return _prev;
	}
	
	public Block next()
	{
		return _next;
	}
	
	public Block prev()
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
	
	public int length()
	{
		return length();
	}
}
