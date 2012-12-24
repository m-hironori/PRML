package prml.chap3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import prml.util.XYGraph;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * 線形回帰モデル：ベイズ（基底関数を選択）
 */
public class LinearRegressionEvidence extends LinearRegressionByBayes{
	double[] phiEV;
	
	public LinearRegressionEvidence(int numM, BasisFunction func, double alpha,
			double beta, List<double[]> training) {
		super(numM, func, alpha, beta, training);
		this.phiEV = null;
	}

	// alphaを取得
	public double getAlpha(){
		return this.alpha;
	}
	
	// betaを取得
	public double getBeta(){
		return this.beta;
	}

	/**
	 * alphaとbetaを自動調節。
	 */
	public void arrangeHyperParam() {
		// α = γ /(m_N^T * m_N)
		// β^-1 = 1/(N-γ) ∑_{n=1..N}{ t_n - m_N^T φ(x_n) }^2
		// γ = ∑_{i=0～M-1} λ_i / ( α + λ_i )
		// λ_i : β * Φ^T * Φ　の行列の固有値
		
		// λ_i : β * Φ^T * Φ　の行列の固有値(先に一回だけ計算すればＯＫ)
		// ここでは、β はかけずに Φ^T * Φ の固有値を計算する
		if(this.phiEV == null){
			this.phiEV = new double[this.numM];
			Matrix R = this.PHI.transpose().times(this.PHI).eig().getD();
			for(int j=0; j<this.numM; j++){
				this.phiEV[j] = R.get(j, j);
			}
		}

		// λ_i : β * Φ^T * Φ　の行列の固有値
		// ここでは、βをかける
		double[] lambdas = new double[this.numM];
		for(int j=0; j<this.numM; j++){
			lambdas[j] = this.beta * this.phiEV[j];
		}

		// γ = ∑_{i=0～M-1} λ_i / ( α + λ_i )
		double gamma = 0;
		for(int j=0; j<this.numM; j++){
			gamma += lambdas[j] / (this.alpha + lambdas[j]);
		}

		// α = γ /(m_N^T * m_N)
		this.alpha = gamma / this.mN.transpose().times(this.mN).get(0, 0);

		// β^-1 = 1/(N-γ) ∑_{n=1..N}{ t_n - m_N^T * φ(x_n) }^2 = 1/(N-γ) * || t - Φ * m_N ||^2
		double betaInv = 0;
		Matrix t = new Matrix(training.size(), 1);
		for (int i = 0; i < training.size(); i++) {
			double val_ti = training.get(i)[1];
			t.set(i, 0, val_ti);
		}
		double e = t.minus(this.PHI.times(this.mN)).normF();
		betaInv = 1.0/((double)this.training.size() - gamma) * Math.pow(e, 2);
		this.beta = 1.0 / betaInv;
	}
	
	//エビデンスを計算
	public double calcEvidence(){
		// ln p(t|α,β) = M/2 * ln α + N/2 * ln β - E(m_N) - 1/2 ln |A| - N/2 * ln(2π)
		// E(m_N) = β/2 * || t - Φ * m_N ||^2 + α/2 * m_N^T * m_N
		// A = α * I + β * Φ^T * Φ 
		
		double evidence = 0.0;

		//ヘッセ行列A
		Matrix I = new Matrix(this.numM, this.numM);
		for (int i = 0; i < this.numM; i++) {
			for (int j = 0; j < this.numM; j++) {
				if (i == j)
					I.set(i, j, 1.0);
			}
		}
		Matrix A = I.times(this.alpha).plus( this.PHI.transpose().times(this.PHI).times(this.beta));

		//E(m_N)
		Matrix t = new Matrix(training.size(), 1);
		for (int i = 0; i < training.size(); i++) {
			double val_ti = training.get(i)[1];
			t.set(i, 0, val_ti);
		}
		double e = t.minus(this.PHI.times(this.mN)).normF();
		double mNex = this.beta / 2.0 * Math.pow(e, 2) + this.alpha / 2.0 * this.mN.transpose().times(this.mN).get(0, 0);

		//エビデンス
		//ln p(t|α,β) = M/2 * ln α + N/2 * ln β - E(m_N) - 1/2 ln |A| - N/2 * ln(2π)
		evidence = ((double)this.numM) / 2.0 * Math.log(this.alpha)
					+ ((double)this.training.size()) / 2.0 * Math.log(this.beta)
					- mNex
					- 1.0 / 2.0 * Math.log(A.det())
					- ((double)this.training.size()) / 2.0 * Math.log(Math.PI);

//http://sage.math.canterbury.ac.nz/home/pub/95/
// A.detをA.normにすると教科書に大体あっているらしい　→　確かにそうなる
//		evidence = ((double)this.numM) / 2.0 * Math.log(this.alpha)
//				+ ((double)this.training.size()) / 2.0 * Math.log(this.beta)
//				- mNex
//				- 1.0 / 2.0 * Math.log(A.norm1())
//				- ((double)this.training.size()) / 2.0 * Math.log(Math.PI);

		return evidence;
	}
	

	
	public static void main(String args[]) throws IOException {
		// 学習データファイル名
		String filename = args[0];
		// モデルのパラメータ数
		int mMax = Integer.parseInt(args[1]);
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

		//モデル数を変えながらエビデンスを計算
		double[] evidences = new double[mMax+1];
		for(int m = 1; m <= mMax; m++){
			// 基底関数選択
			BasisFunction func = null;
			if (BasisFuncName.equals("GAUSIAN")) {
				func = new GausianBasisFunction(m, 0.0, 1.0);
			} else if (BasisFuncName.equals("POLY")) {
				func = new PolynomialBasisFunction();
			}
		
			// 線形回帰ベイズ
			LinearRegressionEvidence lrBayes = new LinearRegressionEvidence(m, func,
					alpha, beta, trainingData);
			// wの事後分布計算(ベイズ)
			lrBayes.learn();
			// エビデンス
			evidences[m] = lrBayes.calcEvidence();
		}
		
		// グラフ描画
		// エビデンスデータ
		double[][] evidenceValues = makeEvidenceValues(evidences);
		
		// グラフ表示
		XYGraph xyGraph = new XYGraph(
				String.format("Evidence(%s, alpha=%.03f, beta=%.03f)", BasisFuncName, alpha, beta), "M", "Evidence");
		xyGraph.addDataValues("Evidence", evidenceValues, true);
		xyGraph.rangeX(0.0, mMax);
		xyGraph.saveGraphAsPNG("evidence.png", 500, 300); // viewメソッドの後に呼び出すと、動作がおかしいので注意
		xyGraph.view(700, 700);


		//---------------------------------------------------
		//モデル数を固定し、aplhaとbetaを調節しながら、エビデンスを計算
		int cntMax = 10;
		double[] evidences2 = new double[cntMax+1];
		// 基底関数選択
		BasisFunction func = null;
		if (BasisFuncName.equals("GAUSIAN")) {
			func = new GausianBasisFunction(mMax, 0.0, 1.0);
		} else if (BasisFuncName.equals("POLY")) {
			func = new PolynomialBasisFunction();
		}
		// 線形回帰ベイズ
		LinearRegressionEvidence lrBayes = new LinearRegressionEvidence(mMax, func,
				alpha, beta, trainingData);
		System.out.println("------------------------");
		System.out.println("alpha\tbeta");
		for(int cnt = 1; cnt <= cntMax; cnt++){
			System.out.println(lrBayes.getAlpha() + "\t" + lrBayes.getBeta());
			// wの事後分布計算(ベイズ)
			lrBayes.learn();
			// エビデンス
			evidences2[cnt] = lrBayes.calcEvidence();
			//alphaとbetaを調節
			lrBayes.arrangeHyperParam();
		}
		// グラフ表示
		double[][] evidenceValues2 = makeEvidenceValues(evidences2);
		XYGraph xyGraph2 = new XYGraph(
				String.format("Evidence(%s, M=%d)", BasisFuncName, mMax), "cnt", "Evidence");
		xyGraph2.addDataValues("Evidence", evidenceValues2, true);
		xyGraph2.rangeX(0.0, mMax);
		xyGraph2.saveGraphAsPNG("evidence2.png", 500, 300); // viewメソッドの後に呼び出すと、動作がおかしいので注意
		xyGraph2.view(700, 700);
		
	}

	private static double[][] makeEvidenceValues(double[] evidences) {
		double[][] ret = new double[2][evidences.length - 1];
		for (int i = 1; i < evidences.length; i++) {
			//0番目は省く
			ret[0][i-1] = i;
			ret[1][i-1] = evidences[i];
		}
		return ret;
	}
}
