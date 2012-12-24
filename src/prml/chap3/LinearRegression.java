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
 * 線形回帰モデル（基底関数版）
 */
public class LinearRegression {
	int numM;
	BasisFunction func;
	
	public LinearRegression(int numM, BasisFunction func) {
		this.numM = numM;
		this.func = func;
	}

/**
 * 二乗誤差による線形回帰
 *
 * 学習データ:(x,t)={(x_1,t_1),...,(x_n,t_n)
 * 曲線フィット：y(x,w) = w_0 * φ_0(x) + w_1 * φ_1(x) + w_2 * φ_2(x) ... 
 *               = ∑_j={0～M-1} w_j * φ_j(x)
 * Φ^T * Φ * w = Φ^T　* t (Φ:N*M次元, t:N*1次元, w:M*1次元)
 * Φ_{i,j} = φ_j(x_i)
 */
	public double[] getWeightsByMinSqrtErr(List<double[]> training) {
		//計画行列Φを計算
		Matrix PHI = new Matrix(training.size(), numM);
		for (int i = 0; i < training.size(); i++) {
			for (int j = 0; j < numM; j++) {
				double val_Pij = func.phi(j, training.get(i)[0]);
				PHI.set(i, j, val_Pij);
			}
		}
		
		//目的変数の訓練データベクトルt
		Matrix t = new Matrix(training.size(), 1);
		for (int i = 0; i < training.size(); i++) {
			double val_Ti = training.get(i)[1];
			t.set(i, 0, val_Ti);
		}

		//重みベクトルw
		Matrix w = PHI.transpose().times(PHI).inverse().times(PHI.transpose()).times(t);
				
		return w.getColumnPackedCopy();

	}

/**
 * 二乗誤差と正則化による線形回帰
　* 学習データ:(x,t)={(x_1,t_1),...,(x_n,t_n)
 * 曲線フィット：y(x,w) = w_0 * φ_0(x) + w_1 * φ_1(x) + w_2 * φ_2(x) ... 
 *                 = ∑_j={0～M-1} w_j * φ_j(x)
 *
 * ( λ * I - Φ^T * Φ ) * w = Φ^T　* t (Φ:N*M次元, t:N*1次元, w:M*1次元, λ：正則化係数)
 * Φ_{i,j} = φ_j(x_i)
 *
 * @param training
 * @param r 正則化係数
 * @return 重み
 */
	public double[] getWeightsByMinSqrtErrReg(List<double[]> training, double r) {
		//計画行列Φを計算
		Matrix PHI = new Matrix(training.size(), numM);
		for (int i = 0; i < training.size(); i++) {
			for (int j = 0; j < numM; j++) {
				double val_Pij = func.phi(j, training.get(i)[0]);
				PHI.set(i, j, val_Pij);
			}
		}
		
		//目的変数の訓練データベクトルt
		Matrix t = new Matrix(training.size(), 1);
		for (int i = 0; i < training.size(); i++) {
			double val_Ti = training.get(i)[1];
			t.set(i, 0, val_Ti);
		}

		// 行列λI
		Matrix rI = new Matrix(numM, numM);
		for (int i = 0; i < numM; i++) {
			for (int j = 0; j < numM; j++) {
				if (i == j)
					rI.set(i, j, r);
				else
					rI.set(i, j, 0.0);
			}
		}

		//重みベクトルw
		Matrix w = (rI.plus(PHI.transpose().times(PHI))).inverse().times(PHI.transpose()).times(t);
				
		return w.getColumnPackedCopy();
	}

	public static void main(String args[]) throws IOException {
		// 学習データファイル名
		String filename = args[0];
		// モデルのパラメータ数
		int m = Integer.parseInt(args[1]);
		//　基底関数
		String BasisFuncName = args[2];
		// ln λ
		int ln_r = Integer.parseInt(args[3]);

		//基底関数選択
		BasisFunction func = null;
		if(BasisFuncName.equals("GAUSIAN")){
			func = new GausianBasisFunction(m, 0.0, 1.0);
		} else if(BasisFuncName.equals("POLY")){
			func = new PolynomialBasisFunction();
		}

		// 線形回帰オブジェクト
		LinearRegression pcFitEM = new LinearRegression(m, func);

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

		// 重み計算(二乗和誤差最小)
		double[] w_MinSqrtErr = pcFitEM.getWeightsByMinSqrtErr(trainingData);
		System.out.println("----");
		for (int i = 0; i < w_MinSqrtErr.length; i++) {
			System.out.println(w_MinSqrtErr[i]);
		}

		// 重み計算(二乗和誤差+正則化最小)
		double r = Math.pow(Math.E, ln_r);
		double[] w_MinSqrtErrReg = pcFitEM.getWeightsByMinSqrtErrReg(trainingData, r);
		System.out.println("----");
		for (int i = 0; i < w_MinSqrtErrReg.length; i++) {
			System.out.println(w_MinSqrtErrReg[i]);
		}
		
		//グラフ描画用
		//正解データ
		double[][] sineValues = makeSineValues();
		//訓練データ
		double[][] trainValues = makeTrainValues(trainingData);
		//学習結果(二乗和)
		double[][] resultValues = makeResultValues(w_MinSqrtErr, m, func);
		//学習結果(二乗和+正則化)
		double[][] resultValuesReg = makeResultValues(w_MinSqrtErrReg, m, func);

		//グラフ表示
		XYGraph xyGraph = new XYGraph("Fit Sine(m="+ String.valueOf(m)+", Func="+ BasisFuncName+")", "X", "Y");
		xyGraph.addDataValues("Sin(x)", sineValues, true);
		xyGraph.addDataValues("Training", trainValues, false);
		xyGraph.addDataValues("MinSqrtErr", resultValues, true);
		xyGraph.addDataValues("MinSqrtErrReg(ln rabmda = "+String.valueOf(ln_r)+")", resultValuesReg, true);
		xyGraph.rangeX(0.0, 1.0);
		xyGraph.rangeY(-1.0, 1.0);
		xyGraph.saveGraphAsPNG("sin.png", 500, 300); //viewメソッドの後に呼び出すと、動作がおかしいので注意
		xyGraph.view(700, 700);
	}

	private static double[][] makeSineValues() {
		double[][] ret = new double[2][101];
		//0-1を100個のデータで埋める
		for(int i=0; i<=100; i++){
			ret[0][i] = i/100.0; //X
			ret[1][i] = Math.sin(2.0 * Math.PI * ret[0][i]); //Y
		}
		return ret;
	}

	private static double[][] makeTrainValues(List<double[]> trainingDataValues) {
		double[][] ret = new double[2][trainingDataValues.size()];
		int i=0;
		for(double[] rec: trainingDataValues){
			ret[0][i] = rec[0]; //X
			ret[1][i] = rec[1]; //Y
			i++;
		}
		return ret;
	}

	private static double[][] makeResultValues(double[] w, int m, BasisFunction func) {
		double[][] ret = new double[2][101];
		//0-1を100個のデータで埋める
		for(int i=0; i<=100; i++){
			ret[0][i] = i/100.0; //X
			for(int j=0; j<m; j++){ //Y
				ret[1][i] += w[j] * func.phi(j, ret[0][i]); 
			}
		}
		return ret;
	}
}
