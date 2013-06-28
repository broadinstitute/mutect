package org.broadinstitute.cga.tools.gatk.filters;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.filters.ReadFilter;

/**
 * Only use reads from the specified library
 *
 * @author kcibul
 * @since Aug 15, 2012
 *
 */

public class LibraryReadFilter extends ReadFilter {
    @Argument(fullName = "library", shortName = "library", doc="The name of the library to keep, filtering out all others", required=true)
    private String LIBRARY_TO_KEEP = null;

    public boolean filterOut( final SAMRecord read ) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        return ( readGroup == null || !readGroup.getLibrary().equals( LIBRARY_TO_KEEP ) );
    }
}
