package org.gersteinlab.age;

public class Scorer {
	private short[][] _scores;
	private short _match;
	private short _mismatch;
	private short _gapOpen;
	private short _gapExtend;
	
	public Scorer(short match, short mismatch, short gapOpen, short gapExtend)
	{
		_match = match;
		_mismatch = mismatch;
		_gapOpen = gapOpen;
		_gapExtend = gapExtend;
		
		_scores = new short[256][256];
		for (int i1 = 0; i1 < 256; i1++) {
			char c1 = (char) i1;
			
			if (Sequence.complement(c1) == c1)
				continue;
			for (int i2 = 0; i2 < 256; i2++) {
				char c2 = (char) i2;
				if (Sequence.complement(c2) == c2)
					continue;
				if (Sequence.sameNuc(c1, c2))
					_scores[c1][c2] = match;
				else
					_scores[c1][c2] = mismatch;
			}
		}
	}
	
	public short getScore(char c1, char c2)
	{
		return _scores[c1][c2];
	}
	
	public short[] getScores(char c)
	{
		return _scores[c];
	}
	
	public short getMismatch()
	{
		return  _mismatch;
	}
	
	public short getMatch()
	{
		return _match;
	}
	
	public short getGapOpen()
	{
		return _gapOpen;
	}
	
	public short getGapExtend()
	{
		return _gapExtend;
	}
}
