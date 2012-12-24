//https://www.coderanch.com/how-to/java/JFreeChartDemo
import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;

public class XYChartExample extends JFrame {

	// how many values are there in the chart
	public static final int NUM_VALUES = 60;

	public static void main(String[] args) throws IOException {
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

		DefaultXYDataset dataset = new DefaultXYDataset();
		dataset.addSeries("Sine", sineValues);
		dataset.addSeries("Cosine", cosineValues);

		// ScatterPlotのデータ
		double[][] scatterplotValues = new double[2][10];
		for (int i = 0; i < 10; i++) {
			scatterplotValues[0][i] = (double) i; // x
			scatterplotValues[1][i] = (double) i * (0.25 / 2.0) - 1.0; // y
		}
		DefaultXYDataset scatdataset = new DefaultXYDataset();
		scatdataset.addSeries("sample", scatterplotValues);

		// Create the chart
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Sine / Cosine Curves", // The chart title
				"X", // x axis label
				"Y", // y axis label
				dataset, // The dataset for the chart
				PlotOrientation.VERTICAL, true, // Is a legend required?
				false, // Use tooltips
				false // Configure chart to generate URLs?
				);
		
		// scatterplot用データセットを追加
		XYPlot plot = chart.getXYPlot();
		int cnt = plot.getDatasetCount();
		plot.setDataset(1, scatdataset);
		XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer(
				false, // 線を表示しない
				true // 点を表示する
				);
		plot.setRenderer(1, renderer2);
		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

		// save it to a file
//		File imageFile = new File("sine-cosine.png");
//		ChartUtilities.saveChartAsPNG(imageFile, chart, 500, 300);

		// display it on the screen
		JFrame frame = new JFrame();
//		JFrame frame = new XYChartExample();
		frame.getContentPane().add(new ChartPanel(chart));
		frame.setSize(700, 500);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
}
