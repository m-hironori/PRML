package prml.chap3;

/**
 * 線形回帰の既定関数用インターフェース
 */
public interface BasisFunction {
	public double phi(int i, double x);
}
