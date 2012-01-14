package org.gersteinlab.age;

public class AGEAligner {
	public final static short DIAGONAL   = 1;
	public final static short HORIZONTAL = 2;
	public final static short VERTICAL   = 3;
	public final static short MASK       = 3;
	
	public final static short FS_DIAGONAL   = DIAGONAL;
	public final static short FS_HORIZONTAL = HORIZONTAL;
	public final static short FS_VERTICAL   = VERTICAL;
	public final static short FS_MASK       = MASK;
	
	public final static short RS_DIAGONAL   = DIAGONAL << 2;
	public final static short RS_HORIZONTAL = HORIZONTAL << 2;
	public final static short RS_VERTICAL   = VERTICAL << 2;
	public final static short RS_MASK       = MASK << 2;
	
	public final static short FM_DIAGONAL   = DIAGONAL << 4;
	public final static short FM_HORIZONTAL = HORIZONTAL << 4;
	public final static short FM_VERTICAL   = VERTICAL << 4;
	public final static short FM_MASK       = MASK << 4;
	
	public final static short RM_DIAGONAL   = DIAGONAL << 6;
	public final static short RM_HORIZONTAL = HORIZONTAL << 6;
	public final static short RM_VERTICAL   = VERTICAL << 6;
	public final static short RM_MASK       = MASK << 6;
	
	public final static int INDEL_FLAG        = 0x01;
	public final static int TDUPLICATION_FLAG = 0x02;
	public final static int INVERSION_FLAG    = 0x04;
	public final static int INVR_FLAG         = 0x08;
	public final static int INVL_FLAG         = 0x10;
	
	public final static int ALL_FLAGS = INDEL_FLAG | TDUPLICATION_FLAG |
			INVERSION_FLAG | INVR_FLAG | INVL_FLAG;
	
	private Sequence _s1;
	private Sequence _s1_rc;
	private Sequence _s2;
	private String _seq1;
	private String _seq1_rc;
	private String _seq2;
	private int _len1;
	private int _len2;
	private AliFragment _frags;
	private AliFragment _frags_alt;
	
	private int[] _f_score;
	private int[] _r_score;
	private short[] _trace;
	
	private int _score_n1;
	private int _score_n2;
	private int _score_size;
	private int _flag;
	private AGEAligner _aux_aligner;
	
	// Coordinates of the excised region: left1, left2, right1, right2
	private final static int MAX_BPOINTS = 100;
	
	private int[][] _bpoints;
	private int _n_bpoints;
	private int _max;
	
	private short _match;
	private short _mismatch;
	private short _gap_open;
	private short _gap_extend;
	
	/**
	 * 
	 * @author dzc
	 *
	 */
	protected class ExcisedRange {
		public int start;
		public int end;
	}
	
	/**
	 * 
	 * @author dzc
	 *
	 */
	protected class BPIdentity {
		public int left;
		public int right;
		public int range;
	}
	
	/**
	 * 
	 * @author dzc
	 *
	 */
	protected class TrackingMaxima {
		public int lg1;
		public int lg2;
		public int rg1;
		public int rg2;
	}
	
	/**
	 * 
	 * @param s1
	 * @param s2
	 */
	public AGEAligner(Sequence s1, Sequence s2)
	{
		_s1 = s1;
		_s2 = s2;
		_seq1 = _s1.sequence();
		_len1 = _seq1.length();
		_seq2 = _s2.sequence();
		_len2 = _seq2.length();
		_s1_rc = s1.clone();
		_seq1_rc = _s1_rc.sequence();
		_frags = null;
		_frags_alt = null;
		_max = 0;
		_flag = 0;
		_aux_aligner = null;
		_match = 0;
		_mismatch = 0;
		_gap_open = 0;
		_gap_extend = 0;
		_bpoints = new int[MAX_BPOINTS][4];
		
		// Constructing of reverse complement (RC) of first sequence
		// RC could be needed for alignment in regions of inversion
		_s1_rc.revcom();
		
		_score_n1 = _len1 + 2;
		_score_n2 = _len2 + 2;
		_score_size = _score_n1 * _score_n2;
		
		// Scoring and trace matrices
		_f_score = new int[(int) _score_size];
		_r_score = new int[(int) _score_size];
		_trace   = new short[(int) _score_size];
		for (int i = 0; i < (int) _score_size; i++) {
			_f_score[i] = 0;
			_r_score[i] = 0;
			_trace[i]   = 0;
		}
	}
	
	/**
	 * 
	 * @param scr
	 * @param flag
	 * @return
	 */
	public boolean align(Scorer scr, int flag)
	{
		// Erase previous break points
		_n_bpoints = 0;
		
		if (_len1 <= 0 || _len2 <= 0)
			return false;
		
		if ((flag & ALL_FLAGS) == 0)
			return true;
		
		// Deciding on whether we need to use auxiliary aligner
		int auxFlag = 0;
		_flag = flag;
		if ((flag & INVERSION_FLAG) != 0) {
			_flag   = INVL_FLAG;
			auxFlag = INVR_FLAG;
		}
		
		// XXX AGE_TIME
		
		_match      = scr.getMatch();
		_mismatch   = scr.getMismatch();
		_gap_open   = scr.getGapOpen();
		_gap_extend = scr.getGapExtend();
		
		_calcScores(scr);
		// _printMatrix(f_score);
		// _printMatrix(r_score);
		
		// XXX AGE_TIME
		
		_calcMaxima();
		// _printMatrix(f_score);
		// _printMatrix(r_score);
		
		// XXX AGE_TIME
		
		_max = _findBPs(true, MAX_BPOINTS / 2);
		_max = _findBPs(false, MAX_BPOINTS / 2);
		// System.out.println("# bpoints = " + _n_bpoints);
		
		// XXX AGE_TIME
		
		_findAlignment();
		
		// XXX AGE_TIME
		
		// Doing auxiliary alignment
		_aux_aligner = null;
		if (auxFlag != 0) {
			_aux_aligner = new AGEAligner(_s1, _s2);
			if (_aux_aligner.align(scr, auxFlag) == false) {
				_aux_aligner = null;
			}
		}
		
		return true;
	}
	
	/**
	 * 
	 * @return
	 */
	public int score()
	{
		if (_aux_aligner != null && _aux_aligner.score() > _max)
			return _aux_aligner.score();
		
		return _max;
	}
	
	/**
	 * 
	 */
	public void printAlignment()
	{
		if (_aux_aligner != null && _aux_aligner.score() > _max) {
			_aux_aligner.printAlignment();
			return;
		}
		
		final String EXCISED_MESSAGE  = "EXCISED REGION";
		final String EXCISED_MESSAGEs = "EXCISED REGION(S)";
		
		System.out.println();
		System.out.print("MATCH = " + _match + ", ");
		System.out.print("MISMATCH = " + _mismatch + ", ");
		System.out.print("GAP OPEN  = " + _gap_open + ", ");
		System.out.print("GAP EXTEND = " + _gap_extend + ", ");
		
		if ((_flag & INDEL_FLAG) != 0)
			System.out.print(", INDEL");
		else if ((_flag & INVL_FLAG) != 0)
			System.out.print(", INVERSION");
		else if ((_flag & INVR_FLAG) != 0)
			System.out.print(", INVERSION");
		else if ((_flag & TDUPLICATION_FLAG) != 0)
			System.out.print(", TDUPLICATION");
		
		int inc1 = 1;
		int inc2 = 1;
		if (_s1.reverse())
			inc1 = -1;
		if (_s2.reverse())
			inc2 = -1;
		
		int s1 = _s1.start();
		int s2 = _s2.start();
		
		int e1 = s1 + inc1 * (_seq1.length() - 1);
		int e2 = s2 + inc2 * (_seq2.length() - 1);
		
		System.out.print("First  seq [");
		System.out.print(Util.addw(_calcWidth(s1, s2), s1) + ",");
		System.out.print(Util.addw(_calcWidth(e1, e2), e1) + "] => ");
		System.out.print(Util.addw(9, Integer.toString(_seq1.length())) + " nucs '");
		System.out.println(_s1.name() + "'");
		
		System.out.print("Second seq [");
		System.out.print(Util.addw(_calcWidth(s1, s2), s2) + ",");
		System.out.print(Util.addw(_calcWidth(e1, e2), e2) + "] => ");
		System.out.print(Util.addw(9, _seq2.length()) + " nucs '");
		System.out.println(_s2.name() + "'");
		
		int nFrag = 0;
		for (AliFragment f = _frags; f != null; f = f.next())
			nFrag++;
		
		int[] nAli = new int[nFrag + 1];
		int[] nId  = new int[nFrag + 1];
		int[] nGap = new int[nFrag + 1];
		nAli[0] = nId[0] = nGap[0] = 0;
		int index = 1;
		for (AliFragment f = _frags; f != null; f = f.next()) {
			AliCount count = f.countAligned();
			nAli[index] = count.nAli;
			nId[index]  = count.nId;
			nGap[index] = count.nGap;
			nAli[0] += nAli[index];
			nId[0]  += nId[index];
			nGap[0] += nGap[index];
			index++;
		}
		
		int identic = 0;
		int gap = 0;
		if (nAli[0] > 0) {
			identic = (int) (100.0 * nId[0] / nAli[0] + 0.5);
			gap     = (int) (100.0 * nGap[0] / nAli[0] + 0.5);
		}
		
		System.out.println("Score: " + Util.addw(9, _max));
		System.out.println("Aligned: " + Util.addw(9, nAli[0]) + "        nucs");
		System.out.print("Identic: " + Util.addw(9, nId[0]));
		System.out.println(" (" + Util.addw(3, identic) + "%) nucs");
		if (nFrag > 1) {
			System.out.print(" =>");
			for (int fr = 1; fr < nFrag + 1; fr++) {
				int e = (int) (100.0 * nId[fr] / nAli[fr] + 0.5);
				System.out.print(" " + Util.addw(9, nId[fr]));
				System.out.print(" (" + Util.addw(3, e) + "%)");
			}
		}
		System.out.println();
		System.out.print("Gaps:    " + Util.addw(9, nGap[0]));
		System.out.println(" (" + Util.addw(3, gap) + "%) nucs\n");
		
		// Printing aligned region coordinates
		if (nAli[0] > 0) {
			System.out.println("Alignment:");
			System.out.println(" first  seq =>  ");
			for (AliFragment f = _frags; f != null; f = f.next()) {
				if (f.prev() != null)
					System.out.print(" " + EXCISED_MESSAGE + " ");
				int ws = _calcWidth(f.start1(), f.start2());
				int we = _calcWidth(f.end1(), f.end2());
				System.out.print("[" + Util.addw(ws, f.start1()));
				System.out.print("," + Util.addw(we, f.end1()) + "]");
			}
			System.out.println();
			System.out.print(" second seq =>  ");
			for (AliFragment f = _frags; f != null; f = f.next()) {
				if (f.prev() != null)
					System.out.print(" " + EXCISED_MESSAGE + " ");
				int ws = _calcWidth(f.start1(), f.start2());
				int we = _calcWidth(f.end1(), f.end2());
				System.out.print("[" + Util.addw(ws, f.start2()));
				System.out.print("," + Util.addw(we, f.end2()) + "]");
			}
			System.out.println();
		}
		
		if (nFrag > 1) {
			System.out.println("\n" + EXCISED_MESSAGEs + ":");
			_printExcised(_frags, inc1, inc2);
			
			if (_n_bpoints > 1)
				System.out.println("ALTERNATIVE REGION(S): " + (_n_bpoints - 1));
			for (int i = _n_bpoints - 1; i >= 1; i--) {
				if (_findAlignment(i, true) == false) // Finding alternative alignment
					continue;
				if (_frags_alt == null)
					continue;
				
				_printExcised(_frags_alt, inc1, inc2);
				break;
			}
			
			// Printing sequence identity around breakpoints
			System.out.println("Identity at breakpoints: ");
			_printIdentityAtBPs(inc1, inc2);
			
			System.out.println("Identity outside breakpoints: ");
			_printIdentityOutsideBPs(inc1, inc2);
			
			System.out.println("Identity inside breakpoints: ");
			_printIdentityInsideBPs(inc1, inc2);
		}
		
		// Printing actual alignment
		for (AliFragment f = _frags; f != null; f = f.next()) {
			if (f != _frags)
				System.out.println("\n" + EXCISED_MESSAGE);
			f.printAlignment();
		}
	}
	
	/**
	 * 
	 * @param frags
	 * @param inc1
	 * @param inc2
	 */
	private void _printExcised(AliFragment frags, int inc1, int inc2)
	{
		int s, e, len;
		
		for (AliFragment f = frags; f != null && f.next() != null; f = f.next()) {
			ExcisedRange range = _getExcisedRange(f);
			if (range == null)
				continue;

			s = range.start;
			e = range.end;
			System.out.print(" first  seq => ");
			len = Math.abs(s - e - inc1);
			System.out.print(Util.addw(9, len) + " nucs");
			if (len > 0)
				System.out.print(" [" + s + "," + e + "]");
			System.out.println();
			System.out.print(" second seq => ");
			s = f.end2() + inc2;
			e = f.next().start2() - inc2;
			len = Math.abs(s - e - inc2);
			System.out.print(Util.addw(9, len) + " nucs");
			if (len > 0)
				System.out.print(" [" + s + "," + e + "]");
			System.out.println();
		}
	}
	
	/**
	 * Prints identity at breakpoints. Called by printAlignment().
	 * 
	 * @param inc1
	 * @param inc2
	 */
	private void _printIdentityAtBPs(int inc1, int inc2)
	{
		int s, e, left, right;
		
		for (AliFragment f = _frags; f != null && f.next() != null; f = f.next()) {
			ExcisedRange range = _getExcisedRange(f, 0, 1);
			if (range == null)
				continue;
			
			s = range.start;
			e = range.end;
			int bp1 = inc1 * (s - _s1.start());
			int bp2 = inc1 * (e - _s1.start());
			BPIdentity identity = _calcIdentity(_seq1, bp1, bp2);
			int nHom;
			left  = identity.left;
			right = identity.right;
			nHom  = identity.range;
			System.out.print(" first  seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + (s + inc1 * left) + "," + (s + inc1 * (right - 1)) + "] to");
				System.out.print(" [" + (e + inc1 * left) + "," + (e + inc1 * (right - 1)) + "]");
			}
			System.out.println();
			
			s = f.end2() + inc2;
			e = f.next().start2();
			bp1 = inc2 * (s - _s2.start());
			bp2 = inc2 * (e - _s2.start());
			identity  = _calcIdentity(_seq2, bp1, bp2);
			left  = identity.left;
			right = identity.right;
			nHom  = identity.range;
			System.out.print(" second seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + (s + inc2 * left) + "," + (s + inc2 * (right - 1)) + "] to");
				System.out.print(" [" + (e + inc2 * left) + "," + (e + inc2 * (right - 1)) + "]");
			}
			System.out.println();
		}
	}
	
	/**
	 * Prints identity outside breakpoints. Called by printAlignment().
	 * 
	 * @param inc1
	 * @param inc2
	 */
	private void _printIdentityOutsideBPs(int inc1, int inc2)
	{
		int s, e;
		
		for (AliFragment f = _frags; f != null && f.next() != null; f = f.next()) {
			ExcisedRange range = _getExcisedRange(f, -1, -1);
			if (range == null)
				continue;
			
			s = range.start;
			e = range.end;
			int bp1 = inc1 * (s - _s1.start());
			int bp2 = inc1 * (e - _s1.start());
			int nHom = _calcOutsideIdentity(_seq1, bp1, bp2);
			System.out.print(" first  seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + (s + inc1 * (nHom - 1)) + "," + s + "] to");
				System.out.print(" [" + e + "," + (e + inc1 * (nHom - 1)) + "]");
			}
			System.out.println();
			
			s = f.end2();
			e = f.next().start2();
			bp1 = inc2 * (s - _s2.start());
			bp2 = inc2 * (e - _s2.start());
			nHom = _calcOutsideIdentity(_seq2, bp1, bp2);
			System.out.print(" second seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + (s - inc2 * (nHom - 1)) + "," + s + "] to");
				System.out.print(" [" + e + "," + (e + inc2 * (nHom - 1)) + "]");
			}
			System.out.println();
		}
	}
	
	/**
	 * Prints identity inside breakpoints. Called by printAlignment().
	 * 
	 * @param inc1
	 * @param inc2
	 */
	private void _printIdentityInsideBPs(int inc1, int inc2)
	{
		int s, e;
		
		for (AliFragment f = _frags; f != null && f.next() != null; f = f.next()) {
			ExcisedRange range = _getExcisedRange(f, -1, -1);
			if (range == null)
				continue;
			
			s = range.start;
			e = range.end;
			int bp1 = inc1 * (s - _s1.start());
			int bp2 = inc1 * (e - _s1.start());
			int nHom = _calcInsideIdentity(_seq1, bp1, bp2);
			System.out.print(" first  seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + s + "," + (s + inc1 * (nHom - 1)) + "] to");
				System.out.print(" [" + (e - inc1 * (nHom - 1)) + "," + e + "]");
			}
			System.out.println();
			
			s = f.end2() + inc2;
			e = f.next().start2() - inc2;
			bp1 = inc2 * (s - _s2.start());
			bp2 = inc2 * (e - _s2.start());
			nHom = _calcInsideIdentity(_seq2, bp1, bp2);
			System.out.print(" second seq => " + Util.addw(9, nHom) + " nucs");
			if (nHom > 0) {
				System.out.print(" [" + s + "," + (s + inc2 * (nHom - 1)) + "] to");
				System.out.print(" [" + (e - inc2 * (nHom - 1)) + "," + e + "]");
			}
			System.out.println();
		}
	}
	
	/**
	 * Note: bp1 and bp2 are zero-based
	 * @param seq
	 * @param bp1
	 * @param bp2
	 */
	private BPIdentity _calcIdentity(String seq, int bp1, int bp2)
	{
		if (bp1 >= bp2)
			return null;
		
		BPIdentity identity = new BPIdentity();
		int n = seq.length();
		int delta = bp2 - bp1;
		
		if (delta == 0)
			return null;
		
		int start = bp1 - 1;
		identity.left = 0;
		identity.right = 0;
		
		while (start >= 0 && Sequence.sameNuc(seq.charAt(start), seq.charAt(start + delta))) {
			start--;
			identity.left--;
		}
		
		int end = bp1;
		while (end < n && Sequence.sameNuc(seq.charAt(end), seq.charAt(end + delta))) {
			end++;
			identity.right++;
		}
		
		identity.range = identity.right - identity.left;
		return identity;
	}
	
	/**
	 * 
	 * @param seq
	 * @param bp1
	 * @param bp2
	 * @return
	 */
	private int _calcOutsideIdentity(String seq, int bp1, int bp2)
	{
		if (bp1 >= bp2)
			return 0;
		
		int n = seq.length();
		int d1 = bp1 + 1;
		int d2 = n - bp2;
		int nCheck = d1;
		if (d2 < d1)
			nCheck = d2;
		
		int start1 = bp1 - nCheck + 1;
		for (int i = 0; i < nCheck; i++) {
			int nSame = 0;
			int start2 = bp2 - i;
			for (int j = i; j < nCheck; j++) {
				if (Sequence.sameNuc(seq.charAt(start1 + j), seq.charAt(start2 + j)))
					nSame++;
				else
					break;
			}
			if (nSame == nCheck - i)
				return nSame;
		}
		
		return 0;
	}
	
	/**
	 * 
	 * @param seq
	 * @param bp1
	 * @param bp2
	 * @return
	 */
	private int _calcInsideIdentity(String seq, int bp1, int bp2)
	{
		if (bp1 >= bp2)
			return 0;
		
		int nCheck = (bp2 - bp1 + 1) / 2;
		int start2 = bp2 - nCheck + 1;
		for (int i = 0; i < nCheck; i++) {
			int nSame = 0;
			int start1 = bp1 - i;
			for (int j = 1; j < nCheck; j++) {
				if (Sequence.sameNuc(seq.charAt(start1 + j), seq.charAt(start2 + j)))
					nSame++;
				else
					break;
			}
			if (nSame == nCheck - i)
				return nSame;
		}
		
		return 0;
	}
	
	/**
	 * 
	 * @param f
	 * @param addStart
	 * @param addEnd
	 * @return
	 */
	private ExcisedRange _getExcisedRange(AliFragment f, int addStart, int addEnd)
	{
		ExcisedRange range = _getExcisedRange(f);
		if (range == null)
			return null;
		
		int inc1 = 1;
		if (_s1.reverse() != false)
			inc1 = -1;
		
		range.start += inc1 * addStart;
		range.end   += inc1 * addEnd;
		
		return range;
	}
	
	/**
	 * 
	 * @param f
	 * @return
	 */
	private ExcisedRange _getExcisedRange(AliFragment f)
	{
		if (f.next() == null)
			return null;
		
		ExcisedRange range = new ExcisedRange();
		int inc1 = 1;
		if (_s1.reverse() != false)
			inc1 = -1;
		if ((_flag & INDEL_FLAG) != 0) {
			range.start = f.end1() + inc1;
			range.end   = f.next().start1() - inc1;
		} else if ((_flag & TDUPLICATION_FLAG) != 0) {
			range.start = f.next().start1();
			range.end   = f.end1();
		} else if ((_flag & INVL_FLAG) != 0) {
			range.start = f.end1() + inc1;
			range.end   = f.next().start1();
		} else if ((_flag & INVR_FLAG) != 0) {
			range.start = f.end1();
			range.end   = f.next().start1() - inc1;
		} else {
			return null;
		}
		
		return range;
	}
	
	/**
	 * 
	 * @param v1
	 * @param v2
	 * @return
	 */
	private int _calcWidth(int v1, int v2)
	{
		int width1 = 0;
		if (v1 < 0)
			v1++;
		while (v1 != 0) {
			width1++;
			v1 /= 10;
		}
		
		int width2 = 0;
		if (v2 < 0)
			v2++;
		while (v2 != 0) {
			width2++;
			v2 /= 10;
		}
		
		if (width1 > width2)
			return width1;
		return width2;
	}
	
	/**
	 * 
	 * @param matr
	 */
	private void _printMatrix(int[] matr)
	{
		int pos = 0;
		for (int i2 = 0; i2 < _score_n2; i2++) {
			for (int i1 = 0; i1 < _score_n1; i1++) {
				System.out.print(Util.addw(4, matr[pos]) + " ");
				pos++;
			}
			System.out.println();
		}
		System.out.println();
	}
	
	/**
	 * 
	 * @param forward
	 * @param maxBPs
	 * @return
	 */
	private int _findBPs(boolean forward, int maxBPs)
	{
		int max = 0;
		if ((_flag & INDEL_FLAG) != 0)
			max = _findIndelBPs(forward, maxBPs);
		else if ((_flag & (INVL_FLAG | INVR_FLAG | TDUPLICATION_FLAG)) != 0)
			max = _findInversionTDuplicationBPs(forward, maxBPs);
		else
			return max;
		
		return max;
	}
	
	/**
	 * 
	 * @param forward
	 * @param maxBPs
	 * @return
	 */
	private int _findInversionTDuplicationBPs(boolean forward, int maxBPs)
	{
		int max = 0;
		int nAdd = 0;
		int[] fMaxima = _f_score;
		int[] rMaxima = _r_score;
		if (forward == true) {
			int posf = _score_n1 - 2;
			int posr = _score_n1 + 1;
			for (int i2 = 0; i2 < _score_n2 - 1; i2++) {
				int val = fMaxima[posf] + rMaxima[posr];
				if (val > max)
					max = val;
				posf += _score_n1;
				posr += _score_n1;
			}
			
			posf = _score_n1 - 2;
			posr = _score_n1 + 1;
			for (int i2 = 0; i2 < _score_n2 - 1; i2++) {
				if (fMaxima[posf] + rMaxima[posr] == max) {
					for (int j1 = 0; j1 < _score_n1 - 1; j1++) {
						int ind1 = posr - _score_n1 + j1 - 1;
						if ((_trace[ind1] & FM_MASK) != 0)
							continue;
						for (int j2 = 0; j2 < _score_n1 - 1; j2++) {
							int ind2 = posr + j2;
							if (fMaxima[ind1] + rMaxima[ind2] == max) {
								int lg1 = j1;
								int lg2 = i2;
								int rg1 = j2 + 1;
								int rg2 = i2 + 1;
								TrackingMaxima maxima = _traceBackMaxima(lg1, lg2, rg1, rg2); 
								lg1 = maxima.lg1;
								lg2 = maxima.lg2;
								rg1 = maxima.rg1;
								rg2 = maxima.rg2;
								
								int res = _addBpoints(lg1, lg2, rg1, rg2);
								if (res > 0)
									nAdd++;
								if (nAdd >= maxBPs || res < 0)
									return max;
							}
						}
					}
				}
				posf += _score_n1;
				posr += _score_n1;
			}
		} else { // Reverse
			int posf = _score_size - _score_n1 - 2;
			int posr = _score_size - _score_n1 + 1;
			for (int i2 = _score_n2 - 1; i2 > 0; i2--) {
				if (fMaxima[posf] + rMaxima[posr] == max) {
					for (int j2 = _score_n1 - 1; j2 > 0; j2--) {
						int ind2 = posr + j2 - 1;
						if ((_trace[ind2] & RM_MASK) != 0)
							continue;
						for (int j1 = _score_n1 - 1; j1 > 0; j1--) {
							int ind1 = posr - _score_n1 + j1 - 2;
							if (fMaxima[ind1] + rMaxima[ind2] == max) {
								int lg1 = j1 - 1;
								int lg2 = i2 - 1;
								int rg1 = j2;
								int rg2 = i2;
								TrackingMaxima maxima = _traceBackMaxima(lg1, lg2, rg1, rg2);
								lg1 = maxima.lg1;
								lg2 = maxima.lg2;
								rg1 = maxima.rg1;
								rg2 = maxima.rg2;
								
								int res = _addBpoints(lg1, lg2, rg1, rg2);
								if (res > 0)
									nAdd++;
								if (nAdd >= maxBPs || res < 0)
									return max;
							}
						}
					}
				}
				posf -= _score_n1;
				posr -= _score_n1;
			}
		}
		
		return max;
	}
	
	/**
	 * 
	 * @param forward
	 * @param maxBPs
	 * @return
	 */
	private int _findIndelBPs(boolean forward, int maxBPs)
	{
		int max = 0;
		int nAdd = 0;
		int[] fMaxima = _f_score;
		int[] rMaxima = _r_score;
		
		if (forward) {
			int[] maxf = _f_score;
			int[] maxr = _r_score;
			int k = 0;
			
			for (int i2 = 0; i2 < _score_n2 - 1; i2++) {
				for (int i1 = 0; i1 < _score_n1 - 1; i1++) {
					int val = maxf[k] + maxr[k + _score_n1 + 1];
					if (val > max)
						max = val;
					k++;
				}
				k++;
			}

			k = 0;
			short[] tr = _trace;
			for (int i2 = 0; i2 < _score_n2 - 1; i2++) {
				for (int i1 = 0; i1 < _score_n1 - 1; i1++) {
					if ((tr[k] & FM_MASK) == 0 && maxf[k] + maxr[k + _score_n1 + 1] == max) {
						int lg1 = i1;
						int lg2 = i2;
						int rg1 = i1 + 1;
						int rg2 = i2 + 1;
						TrackingMaxima maxima = _traceBackMaxima(lg1, lg2, rg1, rg2);
						lg1 = maxima.lg1;
						lg2 = maxima.lg2;
						rg1 = maxima.rg1;
						rg2 = maxima.rg2;
						int res = _addBpoints(lg1, lg2, rg1, rg2);
						if (res > 0)
							nAdd++;
						if (nAdd >= maxBPs || res < 0)
							return max;
					}
					k++;
				}
				k++;
			}
		} else { // Reverse
			int[] maxf = _f_score;
			int[] maxr = _r_score;
			int k = _score_size;
			for (int i2 = _score_n2 - 1; i2 > 0; i2--) {
				for (int i1 = _score_n1 - 1; i1 > 0; i1--) {
					int val = maxf[k - _score_n1 - 2] + maxr[k - 1];
					if (val > max)
						max = val;
					k--;
				}
				k--;
			}
			
			k = _score_size;
			short[] tr = _trace;
			for (int i2 = _score_n2 - 1; i2 > 0; i2--) {
				for (int i1 = _score_n1 - 1; i1 > 0; i1--) {
					if ((tr[k] & RM_MASK) == 0 && maxf[k - _score_n1 - 2] + maxr[k - 1] == max) {
						int lg1 = i1 - 1;
						int lg2 = i2 - 1;
						int rg1 = i1; 
						int rg2 = i2;
						TrackingMaxima maxima = _traceBackMaxima(lg1, lg2, rg1, rg2);
						lg1 = maxima.lg1;
						lg2 = maxima.lg2;
						rg1 = maxima.rg1;
						rg2 = maxima.rg2;
						int res = _addBpoints(lg1, lg2, rg1, rg2);
						if (res > 0)
							nAdd++;
						if (nAdd >= maxBPs || res < 0)
							return max;
					}
					k--;
				}
				k--;
			}
		}
		
		return max;
	}
	
	/**
	 * 
	 * @param lg1
	 * @param lg2
	 * @param rg1
	 * @param rg2
	 * @return
	 */
	private TrackingMaxima _traceBackMaxima(int lg1, int lg2, int rg1, int rg2)
	{
		// XXX pos is declared long in C++ AGE
		int pos = lg2 * _score_n1 + lg1;
		int tr = _trace[pos] & FM_MASK;
		TrackingMaxima maxima = new TrackingMaxima();
		
		while (tr != 0) {
			if (tr == FM_VERTICAL) {
				pos -= _score_n1;
				lg2--;
			} else if (tr == FM_HORIZONTAL) {
				pos--;
				lg1--;
			} else {
				System.err.println("Internal error (1) in tracking.");
			}
			tr = _trace[pos] & FM_MASK;
		}
		
		pos = rg2 * _score_n1 + rg1;
		tr = _trace[pos] & RM_MASK;
		while (tr != 0) {
			if (tr == RM_VERTICAL) {
				pos += _score_n1;
				rg2++;
			} else if (tr == RM_HORIZONTAL) {
				pos++;
				rg1++;
			} else {
				System.err.println("Internal error (2) in tracking.");
			}
			tr = _trace[pos] & RM_MASK;
		}
		
		maxima.lg1 = lg1;
		maxima.lg2 = lg2;
		maxima.rg1 = rg1;
		maxima.rg2 = rg2;
		
		return maxima;
	}
	
	/**
	 * 
	 * @param lg1
	 * @param lg2
	 * @param rg1
	 * @param rg2
	 * @return
	 */
	private int _addBpoints(int lg1, int lg2, int rg1, int rg2)
	{
		// System.out.println(lg1 + " " + lg2 + " " + rg1 + " " + rg2);
		if (lg1 == 0 || lg2 == 0)
			lg1 = lg2 = 0;
		if (rg1 == _score_n1 - 1 || lg2 == _score_n2 - 1) {
			rg1 = _score_n1 - 1;
			rg2 = _score_n2 - 1;
		}
		
		for (int i = 0; i < _n_bpoints; i++) 
			if (_bpoints[i][0] - lg1 == _bpoints[i][2] - rg1 &&
				_bpoints[i][1] - lg2 == _bpoints[i][3] - rg2)
				return 0;
		
		if (_n_bpoints < MAX_BPOINTS) {
			_bpoints[_n_bpoints][0] = lg1;
			_bpoints[_n_bpoints][1] = lg2;
			_bpoints[_n_bpoints][2] = rg1;
			_bpoints[_n_bpoints][3] = rg2;
			_n_bpoints++;
			return 1;
		}
		
		// System.err.println("Buffer for breakpoints exceeded. Skipping ...);
		return -1;
	}
	
	/**
	 * 
	 * @return
	 */
	private boolean _findAlignment()
	{
		return _findAlignment(0, false);
	}
	
	/**
	 * 
	 * @param bpIndex
	 * @param forAlt
	 * @return
	 */
	private boolean _findAlignment(int bpIndex, boolean forAlt)
	{
		if (forAlt) {
			_frags_alt = null;
		} else {
			_frags = null;
		}
		
		if (bpIndex > _n_bpoints)
			return false;
		
		final int MAX_N_FRAGS = 2;
		int lg1 = _bpoints[bpIndex][0];
		int lg2 = _bpoints[bpIndex][1];
		int rg1 = _bpoints[bpIndex][2];
		int rg2 = _bpoints[bpIndex][3];
		Block[] bs = {_findLeftBlocks(lg1, lg2), _findRightBlocks(rg1, rg2)};
		
		int start1 = 0;
		int start2 = 0;
		int end1   = 0;
		int end2   = 0;
		
		AliFragment frags = null;
		String ali1 = "";
		String ali2 = "";
		
		for (int b = 0; b < MAX_N_FRAGS; b++) {
			Block bls = bs[b];
			if (bls == null)
				continue;
			
			boolean useRc1 = ((_flag & INVL_FLAG) == 1 && b == 1) ||
					         ((_flag & INVR_FLAG) == 1 && b == 0);
			Sequence s1 = _s1;
			Sequence s2 = _s2;
			String seq1 = _seq1;
			String seq2 = _seq2;
			
			if (useRc1) {
				seq1 = _seq1_rc;
				s1 = _s1_rc;
			}
			
			int inc1 = 1;
			int inc2 = 1;
			if (s1.reverse())
				inc1 = -1;
			if (s2.reverse())
				inc2 = -1;
			
			start1 = s1.start() + inc1 * (bls.start1() - 1);
			start2 = s2.start() + inc2 * (bls.start2() - 1);
			int nucInd = 2 * b;
			int ind1 = bls.start1() - 1;
			int ind2 = bls.start2() - 1;
			while (bls != null) {
				end1 = s1.start() + inc1 * (bls.start1() + bls.length() - 2);
				end2 = s2.start() + inc2 * (bls.start2() + bls.length() - 2);
				int st1 = bls.start1() - 1;
				int st2 = bls.start2() - 1;
				int len = bls.length();
				while (ind1 < st1) {
					ali1 += seq1.charAt(ind1);
					ind1++;
					ali2 += Sequence.gap();
				}
				while (ind2 < st2) {
					ali1 += Sequence.gap();
					ali2 += seq2.charAt(ind2);
					ind2++;
				}
				for (int i = 0; i < len; i++) {
					ali1 += seq1.charAt(ind1);
					ali2 += seq2.charAt(ind2);
					ind1++;
					ind2++;
				}
				bls = bls.next();
			}
			if (ali1.length() > 0 && ali2.length() > 0 && ali1.length() == ali2.length()) {
				AliFragment f = new AliFragment(ali1, ali2, start1, start2, end1, end2);
				if (frags == null)
					frags = f;
				else
					frags = frags.next(f);
			}
			ali1 = "";
			ali2 = "";
		}
		
		while (frags != null && frags.prev() != null)
			frags = frags.prev();
		
		// Set everything in bs to null so garbage collector can do its thing.
		for (int i = 0; i < MAX_N_FRAGS; i++) {
			bs[i] = null;
		}
		bs = null;
		
		if (forAlt)
			_frags_alt = frags;
		else
			_frags = frags;
		
		return true;
	}
	
	private Block _findLeftBlocks(int left1, int left2)
	{
		Block bls = null;
		
		// Producing alignment blocks for left alignment
		int nSave = _len1;
		if (_len2 > nSave)
			nSave = _len2;
		
		int[] p1 = new int[nSave];
		int[] p2 = new int[nSave];
		nSave = 0;
		
		int pos = left2 * _score_n1 + left1;
		int tr = _trace[pos] & FS_MASK;
		while (tr != 0) {
			if (tr == FS_DIAGONAL) {
				p1[nSave] = left1;
				p2[nSave] = left2;
				left1--;
				left2--;
				nSave++;
				pos -= (_score_n1 + 1);
			} else if (tr == FS_HORIZONTAL) {
				left1--;
				pos--;
			} else if (tr == FS_VERTICAL) {
				left2--;
				pos -= _score_n1;
			} else {
				System.err.println("Internal error (3) in tracking.");
			}
			tr = _trace[pos] & FS_MASK;
		}
		
		int start1;
		int start2;
		if (nSave > 0) {
			int index = nSave - 1;
			int n = 1;
			start1 = p1[index];
			start2 = p2[index];
			index--;
			while (index >= 0) {
				if (p1[index] - p1[index + 1] == 1 &&
				    p2[index] - p2[index + 1] == 1)
					n++;
				else {
					Block b = new Block(start1, start2, n, "left");
					if (bls == null)
						bls = b;
					else
						bls = bls.next(b);
					start1 = p1[index];
					start2 = p2[index];
					n = 1;
				}
				index--;
			}
			Block b = new Block(start1, start2, n, "left");
			if (bls == null)
				bls = b;
			else
				bls = bls.next(b);
		}
		
		if (bls != null) {
			while (bls.prev() != null)
				bls = bls.prev();
		}
		return bls;
	}
	
	private Block _findRightBlocks(int right1, int right2)
	{
		Block bls = null;
		
		// Producing alignment blocks for right alignment
		int pos = right2 * _score_n1 + right1;
		int start1 = 0;
		int start2 = 0;
		int n = 0;
		int tr = _trace[pos] & RS_MASK;
		while (tr != 0) {
			if (tr == RS_DIAGONAL) {
				if (right1 <= _len1 && right2 <= _len2) {
					if (n == 0) {
						start1 = right1;
						start2 = right2;
						n = 1;
					} else {
						n++;
					}
				}
				right1++;
				right2++;
				pos += (_score_n1 + 1);
			} else if (tr == RS_HORIZONTAL) {
				right1++;
				pos++;
				if (n > 0) {
					Block b = new Block(start1, start2, n, "right");
					if (bls == null)
						bls = b;
					else
						bls = bls.next(b);
					n = 0;
				}
			} else if (tr == RS_VERTICAL) {
				right2++;
				pos += _score_n1;
				if (n > 0) {
					Block b = new Block(start1, start2, n, "right");
					if (bls == null)
						bls = b;
					else
						bls = bls.next(b);
					n = 0;
				}
			} else {
				System.err.println("Internal error (4) in tracking.");
			}
			tr = _trace[pos] & RS_MASK;
		}
		
		if (n > 0) {
			Block b = new Block(start1, start2, n, "right");
			if (bls == null)
				bls = b;
			else
				bls = bls.next(b);
		}
		
		return bls;
	}
	
	/**
	 * Finds the maxima for the score in the matrix in order to get the 
	 * optimal alignment
	 */
	private void _calcMaxima()
	{
		int thread1_i2 = -100;
		int thread2_i2 = -100;
		
		// In the C++ version, this is the first OMP section. We're not going to
		// do threads right now. I'm putting this section in the same block for
		// now though we're doing this completely serially for now.
		{
			int trIndex = 1;
			short[] tr = _trace;
			for (int i = 1; i < _score_n1; i++) {
				tr[trIndex] |= FM_HORIZONTAL;
				trIndex++;
			}
			for (int i = 1; i < _score_n2 - 1; i++) {
				tr[trIndex] |= FM_VERTICAL;
				trIndex += _score_n1;
			}
			
			// Assigning maxima for leading matrices (forward direction
			trIndex = _score_n1 + 1;
			int[] f  = _f_score;
			int[] f1 = _f_score;
			int[] f2 = _f_score;
			int fIndex  = _score_n1 + 1;
			int f1Index = _score_n1;
			int f2Index = 1;
			
			for (int i2 = 0; i2 < _len2; i2++) {
				thread1_i2 = i2; // XXX only used for OMP to mark row as in use
				
				// XXX We might not need this because it's only used in the C++
				// implementation to wait for the other OMP thread
				while (i2 == thread2_i2)
					;
				
				for (int i1 = 0; i1 < _len2; i1++) {
					short currTrace = 0;
					if (f2[f2Index] > f[fIndex] && i2 > 0) {
						f[fIndex] = f2[f2Index];
						currTrace = FM_VERTICAL;
					}
					if (f1[f1Index] > f[fIndex] && i1 > 0) {
						f[fIndex] = f1[f1Index];
						currTrace = FM_HORIZONTAL;
					}
					tr[trIndex] |= currTrace;
					fIndex++;
					f1Index++;
					f2Index++;
					trIndex++;
				}
				fIndex += 2;
				f1Index += 2;
				f2Index += 2;
				trIndex += 2;
			}
			thread1_i2 = -100;
		}
		
		// the other OMP section in the C++ implementation. Again, put in its
		// own scope for now.
		{
			int trIndex = _score_size - 2;
			short[] tr = _trace;
			for (int i = 1; i < _score_n1; i++) {
				tr[trIndex] |= FM_HORIZONTAL;
				trIndex--;
			}
			for (int i = 1; i < _score_n2; i++) {
				tr[trIndex] |= FM_VERTICAL;
				trIndex -= _score_n1;
			}
			
			// Assigning maxima for trailing matrices (reverse direction)
			trIndex = _score_size - _score_n1 - 2;
			int[] f  = _r_score;
			int[] f1 = _r_score;
			int[] f2 = _r_score;
			int fIndex  = _score_size - _score_n1 - 2;
			int f1Index = _score_size - _score_n1 - 1;
			int f2Index = _score_size - 2;
			
			int tempLen1 = _len1 - 1;
			int tempLen2 = _len2 - 1;
			for (int i2 = tempLen2; i2 >= 0; i2--) {
				thread2_i2 = i2; // XXX only used for OMP to mark row as in use
				int tmp = 0;
				
				// XXX We might not need this because it's only used in the C++
				// implementation to wait for the other OMP thread
				while (i2 == thread1_i2 && tmp < _score_size)
					tmp++;
				
				for (int i1 = tempLen1; i1 >= 0; i1--) {
					short currTrace = 0;
					if (f2[f2Index] > f[fIndex]) {
						f[fIndex] = f2[f2Index];
						currTrace = RM_VERTICAL;
					}
					if (f1[f1Index] > f[fIndex] && i1 < tempLen1) {
						f[fIndex] = f1[f1Index];
						currTrace = RM_HORIZONTAL;
					}
					tr[trIndex] |= currTrace;
					fIndex--;
					f1Index--;
					f2Index--;
					trIndex--;
				}
				fIndex -= 2;
				f1Index -= 2;
				f2Index -= 2;
				trIndex -= 2;
			}
			thread2_i2 = -100;
		}
	}
	
	/**
	 * 
	 * @param scr
	 */
	private void _calcScores(Scorer scorer)
	{
		int thread1_i2 = -100;
		int thread2_i2 = -100;
		
		// XXX As with _calcMaxima, we have OMP sections. This is the first
		// section. Just keep in it's own scope for now.
		{
			int[] fs1 = _f_score;
			int[] fs2 = _f_score;
			short[] tr1 = _trace;
			short[] tr2 = _trace;
			int fs1Index = _score_n1;
			int fs2Index = 0;
			int tr1Index = _score_n1;
			int tr2Index = 1;
			String c2 = _seq2;
			int c2Index = 0;
			for (int i2 = 0; i2 < _len2; i2++) {
				thread1_i2 = i2; // XXX only used for OMP to mark row as in use
				
				// XXX Might not be necessary anymore since we're not using OMP 
				while (i2 == thread2_i2)
					;
				
				String c1 = _seq1;
				if ((_flag & INVR_FLAG) == 1) 
					c1 = _seq1_rc;
				int c1Index = 0;
				
				short[] scores = scorer.getScores(c2.charAt(0));
				for (int i1 = 0; i1 < _len1; i1++) {
					int score = fs2[fs2Index] + scores[c1.charAt(c1Index)];
					fs2Index++;
					c1Index++;
					int sg2 = fs2[fs2Index];
					int tr = tr2[tr2Index] & FS_MASK;
					tr2Index++;
					if (tr == FS_HORIZONTAL || tr == FS_VERTICAL)
						sg2 += scorer.getGapExtend();
					else
						sg2 += scorer.getGapOpen();
					
					int sg1 = fs1[fs1Index];
					fs1Index++;
					tr = tr1[tr1Index] & FS_MASK;
					tr1Index++;
					if (tr == FS_HORIZONTAL || tr == FS_VERTICAL)
						sg1 += scorer.getGapExtend();
					else
						sg1 += scorer.getGapOpen();
					
					short currTrace = FS_DIAGONAL;
					if (sg1 >= score) {
						score = sg1;
						currTrace = FS_HORIZONTAL;
					}
					if (sg2 >= score) {
						score = sg2;
						currTrace = FS_VERTICAL;
					}
					if (score < 0) {
						score = 0;
						currTrace = 0;
					}
					tr1[tr1Index] |= currTrace;
					fs1[fs1Index] = score;
				}
				fs1Index += 2;
				fs2Index += 2;
				tr1Index += 2;
				tr2Index += 2;
				c2Index++;
			}
			thread1_i2 = -100;
		}
		
		// XXX The other OMP section
		{
			// Filling in the reverse direction
			int[] rs1 = _r_score;
			int[] rs2 = _r_score;
			short[] tr1 = _trace;
			short[] tr2 = _trace;
			int rs1Index = _score_size - 1 - _score_n1;
			int rs2Index = _score_size - 1;
			int tr1Index = _score_size - 1 - _score_n1;
			int tr2Index = _score_size - 2;
			
			String c2 = _seq2;
			int c2Index = _len2 - 1;
			
			for (int i2 = _len2 - 1; i2 >= 0; i2--) {
				thread2_i2 = i2; // XXX Indicate the row is in use
				
				int tmp = 0;
				// XXX Only used for waiting for other thread
				while (i2 == thread1_i2 && tmp < _score_size)
					tmp++;
				
				String c1 = _seq1;
				int c1Index = _len1 - 1;
				if ((_flag & INVL_FLAG) == 1)
					c1 = _seq1_rc;
				
				short[] scores = scorer.getScores(c2.charAt(c2Index));
				
				for (int i1 = _len1 - 1; i1 >= 0; i1--) {
					int score = rs2[rs2Index] + scores[c1.charAt(c1Index)];
					rs2Index--;
					c1Index--;
					
					int sg2 = rs2[rs2Index];
					int tr = tr2[tr2Index] & RS_MASK;
					tr2Index--;
					if (tr == RS_HORIZONTAL || tr == RS_VERTICAL)
						sg2 += scorer.getGapExtend();
					else
						sg2 += scorer.getGapOpen();
					
					int sg1 = rs1[rs1Index];
					rs1Index--;
					
					tr = tr1[tr1Index] & RS_MASK;
					tr1Index--;
					if (tr == RS_HORIZONTAL || tr == RS_VERTICAL)
						sg1 += scorer.getGapExtend();
					else
						sg1 += scorer.getGapOpen();
					
					short currTrace = RS_DIAGONAL;
					if (sg1 >= score) {
						score = sg1;
						currTrace = RS_HORIZONTAL;
					}
					if (sg2 >= score) {
						score = sg2;
						currTrace = RS_VERTICAL;
					}
					if (score < 0) {
						score = 0;
						currTrace = 0;
					}
					tr1[tr1Index] |= currTrace;
					rs1[rs1Index] = score;
				}
				rs1Index -= 2;
				rs2Index -= 2;
				tr1Index -= 2;
				tr2Index -= 2;
				c2Index--;
			}
			thread2_i2 = -100;
		}
	}
}
