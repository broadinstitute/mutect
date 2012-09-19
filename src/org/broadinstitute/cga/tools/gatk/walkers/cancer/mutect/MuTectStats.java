package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import java.util.Arrays;
import java.util.List;

/**
 * Collection of Statistical methods and tests used by MuTect
 */
public class MuTectStats {

    public static double calculateMAD(double[] dd, double median) {
        double[] dev = new double[dd.length];
        for(int i=0; i<dd.length; i++) {
            dev[i] = Math.abs(dd[i] - median);
        }
        return getMedian(dev);

    }

    public static double getMedian(double[] data) {
        Arrays.sort(data);
        Double result;

        if (data.length % 2 == 1) {
            // If the number of entries in the list is not even.

            // Get the middle value.
            // You must floor the result of the division to drop the
            // remainder.
            result = data[(int) Math.floor(data.length/2) ];

        } else {
            // If the number of entries in the list are even.

            // Get the middle two values and average them.
            Double lowerMiddle = data[data.length/2 ];
            Double upperMiddle = data[data.length/2 - 1 ];
            result = (lowerMiddle + upperMiddle) / 2;
        }

        return result;
    }

    public static double[] convertIntegersToDoubles(List<Integer> integers)
    {
        double[] ret = new double[integers.size()];
        for (int i=0; i < ret.length; i++)
        {
            ret[i] = integers.get(i);
        }
        return ret;
    }
}
