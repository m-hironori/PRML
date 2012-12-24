//http://thinkit.co.jp/cert/tech/4/6/3.htm


import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;

/**
 * 複数軸グラフのサンプル
 */
public class DualAxisChartSample {
	public static void main(String[] args) {

	    // Seriesの作成
		String barSeries1 = "First";
		String barSeries2 = "Second";
		String lineSeries1 = "Third";

		// カテゴリーの作成
		String category1 = "Category 1";
		String category2 = "Category 2";
		String category3 = "Category 3";
		String category4 = "Category 4";

		// 棒グラフのデータセットの生成
		DefaultCategoryDataset barDataset = new DefaultCategoryDataset();

		barDataset.addValue(1.0, barSeries1, category1);
		barDataset.addValue(4.0, barSeries1, category2);
		barDataset.addValue(3.0, barSeries1, category3);
		barDataset.addValue(5.0, barSeries1, category4);

		barDataset.addValue(5.0, barSeries2, category1);
		barDataset.addValue(7.0, barSeries2, category2);
		barDataset.addValue(6.0, barSeries2, category3);
		barDataset.addValue(8.0, barSeries2, category4);


		// 折れ線グラフのデータセットの生成
		DefaultCategoryDataset lineDataset = new DefaultCategoryDataset();

		lineDataset.addValue(15.0, lineSeries1, category1);
		lineDataset.addValue(24.0, lineSeries1, category2);
		lineDataset.addValue(31.0, lineSeries1, category3);
		lineDataset.addValue(25.0, lineSeries1, category4);

		// ベースとなるグラフを生成
		JFreeChart chart =
			ChartFactory.createBarChart(
				"Dual Axis Sample",
				"Category",
				"Value",
				barDataset,
				PlotOrientation.VERTICAL,
				true,
				true,
				false);

		CategoryPlot plot = chart.getCategoryPlot();
		// 追加する折れ線グラフのデータセットを追加
		plot.setDataset(1, lineDataset);
		plot.mapDatasetToRangeAxis(1, 1);

		// 追加する折れ線グラフの軸を設定
		CategoryAxis domainAxis = plot.getDomainAxis();
		final ValueAxis axis2 = new NumberAxis("Secondary");
		plot.setRangeAxis(1, axis2);

		// 追加する折れ線グラフの表示設定
		LineAndShapeRenderer renderer2 = new LineAndShapeRenderer();
		renderer2.setToolTipGenerator(new StandardCategoryToolTipGenerator());
		plot.setRenderer(1, renderer2);
		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

        // PNGに出力
        File outputFile = new File("./output/DualAxisChartSample.png");
        try {
            ChartUtilities.saveChartAsPNG(outputFile, chart, 500, 500);
        } catch (IOException ioEx) {
            ioEx.printStackTrace();
        }
	}
}
