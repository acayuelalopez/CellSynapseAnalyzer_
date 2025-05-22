import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import ij.IJ;

public class PointsByDistanceTest {
	public static List<Double> minDistance(List<Double> xRed, List<Double> yRed, Point2D myPoint) {
		List<Point2D> points = new ArrayList<Point2D>();
		List<Double> finalDistances = new ArrayList<Double>();
		for (int i = 0; i < xRed.size(); i++)
			points.add(new Point2D.Double(xRed.get(i), yRed.get(i)));

		//Collections.sort(points, createComparator(myPoint));

		for (Point2D p : points) {
			IJ.log(myPoint.distanceSq(p)+"----distance");
			finalDistances.add(myPoint.distanceSq(p));
		}
		// List<Point2D> result = points.subList(0, index);
		// IJ.log("The closest points with distance <=" + maxDistance + " are " +
		// result);
		return finalDistances;

	}

	private static Comparator<Point2D> createComparator(Point2D p) {
		final Point2D finalP = new Point2D.Double(p.getX(), p.getY());
		return new Comparator<Point2D>() {
			@Override
			public int compare(Point2D p0, Point2D p1) {
				double ds0 = p0.distanceSq(finalP);
				double ds1 = p1.distanceSq(finalP);
				return Double.compare(ds0, ds1);
			}

		};
	}

	public static Double findMin(List<Double> list) {

		// check list is empty or not
		if (list == null || list.size() == 0) {
			return Double.MAX_VALUE;
		}

		// create a new list to avoid modification
		// in the original list
		List<Double> sortedlist = new ArrayList<>(list);

		// sort list in natural order
		Collections.sort(sortedlist);

		// first element in the sorted list
		// would be minimum
		return sortedlist.get(0);
	}

}