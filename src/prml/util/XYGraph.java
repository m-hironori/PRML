package prml.util;

import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;

public class XYGraph {

	private JFreeChart chart;
	
	public XYGraph(String graphName, String labelX, String labelY){
		this.chart = ChartFactory.createXYLineChart(
				graphName, // The chart title
				labelX, // x axis label
				labelY, // y axis label
				null, // The dataset for the chart
				PlotOrientation.VERTICAL, true, // Is a legend required?
				false, // Use tooltips
				false // Configure chart to generate URLs?
				);
	}

	//データを追加（double[0]=X座標リスト、double[1]=Y座標リスト）。線も表示するか否か
	public void addDataValues(String name, double[][] dataValuses, boolean line){
		DefaultXYDataset dataset = new DefaultXYDataset();
		dataset.addSeries(name, dataValuses);
		
		XYPlot plot = this.chart.getXYPlot();

		int datasetCnt = plot.getDatasetCount();
		plot.setDataset(datasetCnt, dataset);

		XYLineAndShapeRenderer renderer;
		if(line){
			renderer = new XYLineAndShapeRenderer(
				true, // 線を表示しない
				false // 点を表示する
				);
		}else{
			renderer = new XYLineAndShapeRenderer(
				false, // 線を表示しない
				true // 点を表示する
				);
		}
		
		plot.setRenderer(datasetCnt, renderer);
		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
	}

	//X座標の表示区間
	public void rangeX(double min, double max){
		ValueAxis xAxis = this.chart.getXYPlot().getDomainAxis();
		xAxis.setAutoRange(false);
		xAxis.setRange(min, max);
	}

	//Y座標の表示区間
	public void rangeY(double min, double max){
		ValueAxis yAxis = this.chart.getXYPlot().getRangeAxis();
		yAxis.setAutoRange(false);
		yAxis.setRange(min, max);
	}
	
	//表示（画面サイズを指定）
	public void view(int sizeX, int sizeY){
		String graphName = this.chart.getTitle().getText();
		JFrame frame = new JFrame(graphName);
		frame.getContentPane().add(new ChartPanel(this.chart));
		frame.setSize(sizeX, sizeY);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}

	//保存（ファイル名と画面サイズを指定）
	public void saveGraphAsPNG(String filename, int sizeX, int sizeY) throws IOException {
		File imageFile = new File(filename);
		ChartUtilities.saveChartAsPNG(imageFile, this.chart, sizeX, sizeY);
	}

	//テスト用メイン
	public static void main(String[] args) throws IOException{
		//サンプルデータ(sinとcosとy=x-1)
		final int NUM_VALUES = 60;
		double[][] sineValues = new double[2][NUM_VALUES];
		double[][] cosineValues = new double[2][NUM_VALUES];
		// X values
		for (int i = 0; i < NUM_VALUES; i++) {
			sineValues[0][i] = i / 10.0;
			cosineValues[0][i] = i / 10.0;
		}
		// Y values
		for (int i = 0; i < NUM_VALUES; i++) {
			sineValues[1][i] = Math.sin(sineValues[0][i]);
			cosineValues[1][i] = Math.cos(cosineValues[0][i]);
		}
		//サンプルデータ(y=a*x-1)
		double[][] scatterplotValues = new double[2][10];
		for (int i = 0; i < 10; i++) {
			scatterplotValues[0][i] = (double) i; // x
			scatterplotValues[1][i] = (double) i * (0.25 / 2.0) - 1.0; // y
		}

		//表示テスト
		XYGraph xyGraph = new XYGraph("Sine / Cosine Curves", "X", "Y");
		xyGraph.addDataValues("Sin(x)", sineValues, true);
		xyGraph.addDataValues("Cos(x)", cosineValues, true);
		xyGraph.addDataValues("y=a*x-1", scatterplotValues, false);
		xyGraph.saveGraphAsPNG("sin-cos.png", 500, 300); //先に保存しないと、グラフがおかしくなる
		xyGraph.view(700, 700);
	}

}
