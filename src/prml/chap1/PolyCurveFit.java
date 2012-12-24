package prml.chap1;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import prml.util.XYGraph;

import Jama.Matrix;

public class PolyCurveFit {

	ArrayList<double[]> training;

	public PolyCurveFit() {
		this.training = new ArrayList<double[]>();
	}

	public void addTrainingData(double[] data) {
		this.training.add(data);
	}

	public List<double[]> getTrainingDataValues() {
		return this.training;
	}

	public double[] getWeightsByMinSqrtErr(int m) {
		// 学習データ:(x,t)={(x_1,t_1),...,(x_n,t_n)
		// 曲線フィット：y(x,w) = w_0 + w_1 * x^1 + w_2 * x^2 ... = ∑_j=0～M w_j * x^j

		// Aw = T (A:m*m次元, T:m*1次元, w:m*1次元の重み)
		// A_ij = ∑_n=1～N (x_n)^{i+j}
		// T_i = ∑_n=1～N (x_n)^i * t_n
		// w重み

		//行列Aを直接計算
		Matrix A = new Matrix(m + 1, m + 1);
		for (int i = 0; i < m + 1; i++) {
			for (int j = 0; j < m + 1; j++) {
				double val_Aij = 0;
				for (int n = 0; n < this.training.size(); n++) {
					val_Aij += Math.pow(this.training.get(n)[0], i + j);
				}
				A.set(i, j, val_Aij);
			}
		}
		
		// 行列T
		Matrix T = new Matrix(m + 1, 1);
		for (int i = 0; i < m + 1; i++) {
			double val_Ti = 0;
			for (int n = 0; n < this.training.size(); n++) {
				val_Ti += this.training.get(n)[1]
						* Math.pow(this.training.get(n)[0], i);
			}
			T.set(i, 0, val_Ti);
		}

		Matrix w = A.solve(T);

		// for debug
//		printMatrix(A);
//		printMatrix(T);
//		printMatrix(w);
//		Matrix Residual = A.times(w).minus(T);
//		printMatrix(Residual);

		return w.getColumnPackedCopy();

	}


	public double[] getWeightsByMinSqrtErrReg(int m, double r) {
		// 学習データ:(x,t)={(x_1,t_1),...,(x_n,t_n)
		// 曲線フィット：y(x,w) = w_0 + w_1 * x^1 + w_2 * x^2 ... = ∑_j=0～M w_j * x^j

		// (A+λI)w = T (A:m*m次元, λ:正則化用係数, I:m*m次元単位行列, T:m*1次元, w:m*1次元の重み)
		// A_ij = ∑_n=1～N (x_n)^{i+j}
		// T_i = ∑_n=1～N (x_n)^i * t_n
		// w重み

		// 行列A
		Matrix A = new Matrix(m + 1, m + 1);
		for (int i = 0; i < m + 1; i++) {
			for (int j = 0; j < m + 1; j++) {
				double val_Aij = 0;
				for (int n = 0; n < this.training.size(); n++) {
					val_Aij += Math.pow(this.training.get(n)[0], i + j);
				}
				A.set(i, j, val_Aij);
			}
		}

		// 行列λI
		Matrix rI = new Matrix(m + 1, m + 1);
		for (int i = 0; i < m + 1; i++) {
			for (int j = 0; j < m + 1; j++) {
				if (i == j)
					rI.set(i, j, r);
				else
					rI.set(i, j, 0.0);
			}
		}

		// 行列T
		Matrix T = new Matrix(m + 1, 1);
		for (int i = 0; i < m + 1; i++) {
			double val_Ti = 0;
			for (int n = 0; n < this.training.size(); n++) {
				val_Ti += this.training.get(n)[1]
						* Math.pow(this.training.get(n)[0], i);
			}
			T.set(i, 0, val_Ti);
		}

		Matrix A_rI = A.plus(rI);
		Matrix w = A_rI.solve(T);

		// for debug
		// printMatrix(A);
		// printMatrix(T);
		// printMatrix(w);
		// Matrix Residual = A.times(w).minus(T);
		// printMatrix(Residual);

		return w.getColumnPackedCopy();

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

	public static void main(String args[]) throws IOException {
		// 学習データファイル名
		String filename = args[0];
		// 次数
		int m = Integer.parseInt(args[1]);
		// ln λ
		int ln_r = Integer.parseInt(args[2]);

		// フィッティングオブジェクト
		PolyCurveFit pcFitEM = new PolyCurveFit();

		// データをロードして学習データとして追加
		BufferedReader br = new BufferedReader(new FileReader(
				new File(filename)));
		String line;
		try {
			while ((line = br.readLine()) != null) {
				// System.out.println(line);
				String recStr[] = line.split(" ", 2);
				double[] rec = { Double.parseDouble(recStr[0]),
						Double.parseDouble(recStr[1]) };
				pcFitEM.addTrainingData(rec);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// 重み計算(二乗和誤差最小)
		double[] w_MinSqrtErr = pcFitEM.getWeightsByMinSqrtErr(m);
		System.out.println("----");
		for (int i = 0; i < w_MinSqrtErr.length; i++) {
			System.out.println(w_MinSqrtErr[i]);
		}

		// 重み計算(二乗和誤差+正則化最小)
		double r = Math.pow(Math.E, ln_r);
		double[] w_MinSqrtErrReg = pcFitEM.getWeightsByMinSqrtErrReg(m, r);
		System.out.println("----");
		for (int i = 0; i < w_MinSqrtErrReg.length; i++) {
			System.out.println(w_MinSqrtErrReg[i]);
		}
		
		//グラフ描画用
		//正解データ
		double[][] sineValues = makeSineValues();
		//訓練データ
		double[][] trainValues = makeTrainValues(pcFitEM.getTrainingDataValues());
		//学習結果(二乗和)
		double[][] resultValues = makeResultValues(w_MinSqrtErr, m);
		//学習結果(二乗和)
		double[][] resultValuesReg = makeResultValues(w_MinSqrtErrReg, m);

		//グラフ表示
		XYGraph xyGraph = new XYGraph("Fit Sine(m="+ String.valueOf(m)+")", "X", "Y");
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

	private static double[][] makeResultValues(double[] w, int m) {
		double[][] ret = new double[2][101];
		//0-1を100個のデータで埋める
		for(int i=0; i<=100; i++){
			ret[0][i] = i/100.0; //X
			for(int j=0; j<=m; j++){ //Y
				ret[1][i] += w[j] * Math.pow(ret[0][i],j); 
			}
		}
		return ret;
	}


}
