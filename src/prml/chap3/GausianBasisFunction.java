package prml.chap3;

/**
 * ガウス基底関数
 */
public class GausianBasisFunction implements BasisFunction {

	double[] u;
	double s;

	// パラメータを指定する
	// ガウスの平均と分散を指定。
	// ただし、u[0]は無視する。（u[0]はφ_0(x)に相当し、このときは必ずφ_0(x)=1であるため）
	public GausianBasisFunction(double[] u, double s) {
		this.u = u;
		this.s = s;
	}

	// 定義域とモデル数を指定。パラメータは自動決定
	public GausianBasisFunction(int m, double begin, double end) {
		this.u = makeAutoParamU(m, begin, end);
		this.s = makeAutoParamS(m, begin, end);;
	}

	// 基底関数の結果
	// ただし、j=0の場合、必ず１を返す。
	@Override
	public double phi(int j, double x) {
		if (j == 0) {
			return 1;
		}
		return Math.exp(-1 * Math.pow(x - u[j], 2) / (2 * s * s));
	}

	// パラメータuを自動的に調整する
	// 定義域を、M個で等分に分割する位置にする。
	// 戻り値： ret[1...m] = param u[1...m]
	// ex:M=10, 定義域：0.0 - 1.0
	// u = NaN, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
	private static double[] makeAutoParamU(int m, double begin, double end) {
		double[] u = new double[m];
		
		double s = (end - begin) / m;
		u[0] = Double.NaN; //ret[0]は使わない
		for (int i = 1; i < m; i++) {
			u[i] = begin + s * i;
		}
		return u;
	}

	// パラメータsを自動的に調整する（定義域をモデル数で等分割）
	// 戻り値：param s
	// ex:m=10, 定義域：0.0 - 1.0
	// s = 0.1
	private static double makeAutoParamS(int m, double begin, double end) {
		double s = (end - begin) / m;
		return s;
	}
}
