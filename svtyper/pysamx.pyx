cimport cython
from libc.stdlib cimport malloc, free

from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport AlignmentFile as AF
from pysam.libcalignedsegment cimport AlignedSegment, makeAlignedSegment
import pysam

# not meant to be used with multiprocessing or parallel code!
cdef class AlignmentFile(AF):
    def count(self,
              contig=None,
              start=None,
              stop=None,
              region=None,
              until_eof=False,
              read_callback="all",
              reference=None,
              end=None):
        cdef AlignedSegment read
        cdef long counter = 0
        if not self.is_open:
            raise ValueError("I/O operation on closed file")

        cdef int filter_method = 0
        if read_callback == "all":
            filter_method = 1
        elif read_callback == "nofilter":
            filter_method = 2

        for read in self.fetch(contig=contig,
                               start=start,
                               stop=stop,
                               reference=reference,
                               end=end,
                               region=region,
                               until_eof=until_eof):
            # apply filter
            if filter_method == 1:
                # filter = "all"
                if (read.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800)):
                    continue
            elif filter_method == 2:
                # filter = "nofilter"
                pass
            else:
                if not read_callback(read):
                    continue
            counter += 1

        return counter

    @cython.boundscheck(False)
    cdef int grab_filtered_reads(self, list data, long size, int tid, int start, int stop):
        cdef int retval = 0
        cdef bam1_t *b
        cdef hts_itr_t *iter = NULL
        cdef long i = 0
        cdef AlignedSegment read

        b = bam_init1()
        iter = sam_itr_queryi(self.index, tid, start, stop)
        retval = hts_itr_next(hts_get_bgzfp(self.htsfile), iter, b, self.htsfile)
        while retval >= 0:
            if (b.core.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800)):
                pass
            else:
                read = makeAlignedSegment(b, self)
                data[i] = read
                i += 1
            retval = hts_itr_next(hts_get_bgzfp(self.htsfile), iter, b, self.htsfile)

        hts_itr_destroy(iter)
        bam_destroy1(b)

        if retval == -2:
            return -2

        if i != size:
            return -3

        return 0

    @cython.boundscheck(False)
    def bulkfetch(self, contig=None, start=None, stop=None, region=None, tid=None, num_reads=None, reference=None, end=None):
        cdef int rtid, rstart, rstop, has_coord
        cdef long N = num_reads
        cdef long i
        cdef list reads = [None for i in range(N)]

        if contig is None:
            raise ValueError('Need to pass in a valid contig value!')
        if start is None:
            raise ValueError('Need to pass in a valid start value!')
        if stop is None:
            raise ValueError('Need to pass in a valid stop value!')
        if num_reads is None:
            raise ValueError('Need to pass in a valid num_reads value!')

        if not self.has_index():
            raise ValueError("didn't find the index file corresponding to {}".format(self.filename))

        if not self.is_open:
            raise IOError("Trying to do an I/O operation on a closed file")

        has_coord, rtid, rstart, rstop = self.parse_region(contig, start, stop, region, tid,
                                                           end=end, reference=reference)

        cdef int retval
        retval = self.grab_filtered_reads(reads, num_reads, rtid, start, stop)

        if retval == -2:
            raise IOError('truncated file')

        if retval == -3:
            raise IOError('Expectation ({}) and sam_itr_query fetch ({}) count mismatch!'.format(i, N))

        return reads

# python setup_pysamx.py build_ext --inplace
# >>> import pysam
# >>> import pysamx
# >>> foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r')
# >>> foo.count(contig='2', start=2798760, stop=2798780, read_callback='all')
# 18
# >>> foo.bulkfetch(contig='2', start=2798760, stop=2798780, num_reads=18)
# from timeit import timeit
# old
# >>> timeit("foo = pysam.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); list(foo.fetch(contig='2', start=2798760, stop=2798780))", 'import pysam; import pysamx', number=1000)
# new
# >>> timeit("foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); foo.bulkfetch(contig='2', start=2798760, stop=2798780, num_reads=18)", 'import pysam; import pysamx', number=1000)

# 2nd test (475 reads fetch)

# old
#timeit("foo = pysam.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); list(foo.fetch(contig='2', start=2798760, stop=7009000))", 'import pysam; import pysamx', number=1000)
#0.6631929874420166

# new
#timeit("foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); foo.bulkfetch(contig='2', start=2798760, stop=7009000, num_reads=475)", 'import pysam; import pysamx', number=1000)
#0.6480650901794434
