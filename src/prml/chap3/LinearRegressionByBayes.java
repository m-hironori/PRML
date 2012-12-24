package prml.chap3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import prml.util.XYGraph;

import Jama.Matrix;

/**
 * 線形回帰モデル：ベイズ（基底関数を選択）
 */
public class LinearRegressionByBayes {
	int numM;
	BasisFunction func;
	double alpha;
	double beta;
	List<double[]> training;
	Matrix mN;
	Matrix SNInv;
	Matrix PHI;
	
	private LinearRegressionByBayes() {}

	public LinearRegressionByBayes(int numM, BasisFunction func, double alpha,
			double beta, List<double[]> training) {
		this.numM = numM;
		this.func = func;
		this.alpha = alpha;
		this.beta = beta;
		this.training = training;

		this.mN = null;
		this.SNInv = null;
		this.PHI = null;
	}

	// wの事後分布のパラメータmを取得
	public double[] getParamM() {
		return this.mN.getColumnPackedCopy();
	}

	// wの事後分布のパラメータSを取得
	public double[][] getParamS() {
		return this.SNInv.inverse().getArrayCopy();
	}

	/**
	 * ベイズによる線形回帰。重みwの事後分布を推定。
	 */
	public void learn() {
		// 学習データ:(x,t)={(x_1,t_1),...,(x_n,t_n)
		// 曲線フィット：y(x,w) = w_0 * φ_0(x) + w_1 * φ_1(x) + w_2 * φ_2(x) ...
		//                  = ∑_j={0～M-1} w_j * φ_j(x)

		// w事前分布:
		// p(w)=ガウス(w|m_0, S_0)
		// m_0 = 0, S_0 = α * I
		// w事後分布：
		// p(w|t)=ガウス(w|m_N, S_N)
		// m_N = S_N ( S_0^-1 * m_0 + β * Φ^T * t )
		//     = β * S_N * Φ^T * t // m_0が全部0の場合
		// S_N = S_0 + β * Φ^T * Φ
		// (m:M次元(wのガウス分布の平均), S:M*M次元(wのガウス分布の分散))
		// (Φ:N*M次元, t:N*1次元, w:M*重み)
		// Φ_{i,j} = φ_j(x_i)

		// 計画行列Φを計算
		if(this.PHI == null){ //学習データが変更されたときだけ再計算
			this.PHI = new Matrix(training.size(), this.numM);
			for (int i = 0; i < training.size(); i++) {
				for (int j = 0; j < this.numM; j++) {
					double val_Pij = this.func.phi(j, training.get(i)[0]);
					PHI.set(i, j, val_Pij);
				}
			}
		}

		// 目標変数の訓練データベクトルt
		Matrix t = new Matrix(training.size(), 1);
		for (int i = 0; i < training.size(); i++) {
			double val_ti = training.get(i)[1];
			t.set(i, 0, val_ti);
		}

		// wの事前分布
		// wの事前分布の平均m_0
		// Matrix m0 = new Matrix(numM,1); //全部0
		// wの事前分布の分散の逆行列 S_0^-1
		Matrix I = new Matrix(this.numM, this.numM);
		for (int i = 0; i < this.numM; i++) {
			for (int j = 0; j < this.numM; j++) {
				if (i == j)
					I.set(i, j, 1.0);
			}
		}
		Matrix S0Inv = I.times(this.alpha);

		// 事後分布を計算
		this.SNInv = S0Inv.plus(PHI.transpose().times(PHI).times(this.beta));
		this.mN = this.SNInv.inverse().times(PHI.transpose()).times(t).times(this.beta);
	}
	
	//y(x)の予測分布の平均
	public double getPredictiveDistMean(double x){
		// Mean = m_N^T * φ(x)

		//phi(x)
		Matrix phi = new Matrix(this.numM, 1);
		for (int j = 0; j < this.numM; j++) {
			double val_Pj = this.func.phi(j, x);
			phi.set(j, 0, val_Pj);
		}
		
		Matrix y = this.mN.transpose().times(phi);
		
		return y.get(0, 0);
	}

	//y(x)の予測分布の分散の二乗
	public double getPredictiveDistVar(double x){
		// Var^2 = 1/beta + φ(x)^T * S_N * φ(x)

		//phi(x) 
		Matrix phi = new Matrix(this.numM, 1);
		for (int j = 0; j < this.numM; j++) {
			double val_Pj = this.func.phi(j, x);
			phi.set(j, 0, val_Pj);
		}
		
		Matrix y_Var = phi.transpose().times(this.SNInv.inverse()).times(phi);
		
		return 1.0 / this.beta + y_Var.get(0, 0);
	}

	
	public static void main(String args[]) throws IOException {
		// 学習データファイル名
		String filename = args[0];
		// モデルのパラメータ数
		int m = Integer.parseInt(args[1]);
		// 基底関数
		String BasisFuncName = args[2];
		// alpha, beta
		double alpha = Double.parseDouble(args[3]);
		double beta = Double.parseDouble(args[4]);

		// データをロード
		BufferedReader br = new BufferedReader(new FileReader(
				new File(filename)));
		String line;
		List<double[]> trainingData = new ArrayList<double[]>();
		try {
			while ((line = br.readLine()) != null) {
				String recStr[] = line.split(" ", 2);
				double[] rec = { Double.parseDouble(recStr[0]),
						Double.parseDouble(recStr[1]) };
				trainingData.add(rec);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// 基底関数選択
		BasisFunction func = null;
		if (BasisFuncName.equals("GAUSIAN")) {
			func = new GausianBasisFunction(m, 0.0, 1.0);
		} else if (BasisFuncName.equals("POLY")) {
			func = new PolynomialBasisFunction();
		}

		// 線形回帰ベイズ
		LinearRegressionByBayes lrBayes = new LinearRegressionByBayes(m, func,
				alpha, beta, trainingData);

		// wの事後分布計算(ベイズ)
		lrBayes.learn();

		// wの事後分布の平均m
		double[] w_Avg = lrBayes.getParamM();
		System.out.println("----");
		for (int i = 0; i < w_Avg.length; i++) {
			System.out.println(String.format("m[%d]\t%f", i, w_Avg[i]));
		}

		// wの事後分布の分散S
		double[][] w_Var = lrBayes.getParamS();
		System.out.println("----");
		for (int i = 0; i < w_Var.length; i++) {
			System.out.print(String.format("s[%d]", i));
			for (int j = 0; j < w_Var[0].length; j++) {
				System.out.print(String.format("\t%f", w_Var[i][j]));
			}
			System.out.println();
		}

		// グラフ描画用
		// 正解データ
		double[][] sineValues = makeSineValues();
		// 訓練データ
		double[][] trainValues = makeTrainValues(trainingData);
		// 学習結果
		double[][] resultValues = makeResultValues(lrBayes);
		// 分散
		double[][] resultPlusSigma = makeSigmaValues(resultValues, lrBayes, 1);
		// 分散
		double[][] resultMinusSigma = makeSigmaValues(resultValues, lrBayes, -1);
		
		// グラフ表示
		XYGraph xyGraph = new XYGraph(
				String.format("Fit Sine By Bayes(m=%d, %s, alpha=%.03f, beta=%.03f)", m, BasisFuncName, alpha, beta), "X", "Y");
		xyGraph.addDataValues("Sin(x)", sineValues, true);
		xyGraph.addDataValues("Training", trainValues, false);
		xyGraph.addDataValues("BayesFit", resultValues, true);
		xyGraph.addDataValues("BayesFit+sigma", resultPlusSigma, true);
		xyGraph.addDataValues("BayesFit-sigma", resultMinusSigma, true);
		xyGraph.rangeX(0.0, 1.0);
		xyGraph.rangeY(-1.0, 1.0);
		xyGraph.saveGraphAsPNG("sin.png", 500, 300); // viewメソッドの後に呼び出すと、動作がおかしいので注意
		xyGraph.view(700, 700);
	}

	private static double[][] makeSineValues() {
		double[][] ret = new double[2][101];
		// 0-1を100個のデータで埋める
		for (int i = 0; i <= 100; i++) {
			ret[0][i] = i / 100.0; // X
			ret[1][i] = Math.sin(2.0 * Math.PI * ret[0][i]); // Y
		}
		return ret;
	}

	private static double[][] makeTrainValues(List<double[]> trainingDataValues) {
		double[][] ret = new double[2][trainingDataValues.size()];
		int i = 0;
		for (double[] rec : trainingDataValues) {
			ret[0][i] = rec[0]; // X
			ret[1][i] = rec[1]; // Y
			i++;
		}
		return ret;
	}

	private static double[][] makeResultValues(LinearRegressionByBayes lrBayes) {
		double[][] ret = new double[2][101];
		// 0-1を100個のデータで埋める
		for (int i = 0; i <= 100; i++) {
			ret[0][i] = i / 100.0; // X
			ret[1][i] = lrBayes.getPredictiveDistMean(ret[0][i]); //Y
		}
		return ret;
	}

	private static double[][] makeSigmaValues(double[][] result, LinearRegressionByBayes lrBayes, double numS) {
		double[][] ret = new double[2][result[0].length];

		for (int i = 0; i < result[0].length; i++) {
			ret[0][i] = result[0][i]; // X
			ret[1][i] = result[1][i] + numS * Math.sqrt(lrBayes.getPredictiveDistVar(ret[0][i])); //result + numS * sigma
		}
		return ret;
	}


}
