package prml.chap3;

/**
 * ‘½€®Šî’êŠÖ”
 */
public class PolynomialBasisFunction implements BasisFunction {

	@Override
	public double phi(int j, double x) {
		if (j == 0) {
			return 1;
		}
		return Math.pow(x, j);
	}

}
