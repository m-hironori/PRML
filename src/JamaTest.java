import Jama.Matrix;

public class JamaTest {
	public static void main(String[] args) {
		double[][] varA = { { 2., 2., 3. }, { 3., 5., 2. }, { 5., 3., 3. } };
		Matrix A = new Matrix(varA);
		// Matrix b = Matrix.random(3,1);
		double[][] varb = { { 15. }, { 19. }, { 20. } };
		Matrix b = new Matrix(varb);
		Matrix x = A.solve(b);
		Matrix Residual = A.times(x).minus(b);
		// double rnorm = Residual.normInf();
		printMatrix(A);
		printMatrix(b);
		printMatrix(x);
		printMatrix(Residual);
	}

	private static void printMatrix(Matrix x) {
		for (int j = 0; j < x.getColumnDimension(); j++) {
			System.out.print("\t");
			System.out.print("[" + j + "]");
		}
		System.out.println();
		for (int i = 0; i < x.getRowDimension(); i++) {
			System.out.print("[" + i + "]");
			for (int j = 0; j < x.getColumnDimension(); j++) {
				System.out.print("\t");
				System.out.print(x.get(i, j));
			}
			System.out.println();
		}
		return;
	}
}
