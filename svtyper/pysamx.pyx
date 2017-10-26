cimport cython
from libc.stdlib cimport malloc, free

from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport AlignmentFile as AF, IteratorRow, IteratorRowAll, IteratorRowAllRefs
from pysam.libcalignedsegment cimport AlignedSegment, makeAlignedSegment
import pysam

cdef class IteratorRowRegionFiltered(IteratorRow):
    """*(AlignmentFile samfile, int tid, int beg, int stop,
    int multiple_iterators=False)*
    iterate over mapped reads in a region.
    .. note::
        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.
    """
    cdef hts_itr_t * iter

    def __init__(self, AlignmentFile samfile,
                 int tid, int beg, int stop,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        if not samfile.has_index():
            raise ValueError("no index available for iteration")

        with nogil:
            self.iter = sam_itr_queryi(
                self.samfile.index,
                tid,
                beg,
                stop)

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        with nogil:
            self.retval = hts_itr_next(hts_get_bgzfp(self.htsfile),
                                       self.iter,
                                       self.b,
                                       self.htsfile)

    def __next__(self):
        self.cnext()
        while self.retval >= 0:
            if (self.b.core.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800)):
                # skip read if:
                # https://github.com/samtools/htslib/blob/master/htslib/sam.h
                # BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY
                pass
            else:
                break
            self.cnext()

        if self.retval >= 0:
            return makeAlignedSegment(self.b, self.samfile)
        elif self.retval == -2:
            # Note: it is currently not the case that hts_iter_next
            # returns -2 for a truncated file.
            # See https://github.com/pysam-developers/pysam/pull/50#issuecomment-64928625
            raise IOError('truncated file')
        else:
            raise StopIteration

    def __dealloc__(self):
        hts_itr_destroy(self.iter)

# not meant to be used with multiprocessing or parallel code!
cdef class AlignmentFile(AF):
    cdef long ccount(self, int tid, int start, int stop):
        cdef int retval = 0
        cdef bam1_t *b
        cdef hts_itr_t *iter = NULL
        cdef long i = 0

        b = bam_init1()
        iter = sam_itr_queryi(self.index, tid, start, stop)
        retval = hts_itr_next(hts_get_bgzfp(self.htsfile), iter, b, self.htsfile)
        while retval >= 0:
            if (b.core.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800)):
                # skip read if:
                # https://github.com/samtools/htslib/blob/master/htslib/sam.h
                # BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY
                pass
            else:
                i += 1
            retval = hts_itr_next(hts_get_bgzfp(self.htsfile), iter, b, self.htsfile)

        hts_itr_destroy(iter)
        bam_destroy1(b)

        if retval == -2:
            return -2

        return i

    def count(self,
              contig=None,
              start=None,
              stop=None,
              region=None,
              reference=None,
              tid=None,
              end=None):
        cdef AlignedSegment read
        cdef long counter = 0
        if not self.is_open:
            raise ValueError("I/O operation on closed file")

        cdef int rtid, rstart, rstop, has_coord

        if contig is None:
            raise ValueError('Need to pass in a valid contig value!')
        if start is None:
            raise ValueError('Need to pass in a valid start value!')
        if stop is None:
            raise ValueError('Need to pass in a valid stop value!')

        if not self.has_index():
            raise ValueError("didn't find the index file corresponding to {}".format(self.filename))

        if not self.is_open:
            raise IOError("Trying to do an I/O operation on a closed file")

        has_coord, rtid, rstart, rstop = self.parse_region(contig, start, stop, region, tid,
                                                           end=end, reference=reference)

        counter = self.ccount(rtid, start, stop)

        if counter == -2:
            raise IOError('truncated file')

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
                # skip read if:
                # https://github.com/samtools/htslib/blob/master/htslib/sam.h
                # BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY
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

    def fetch2(self,
              contig=None,
              start=None,
              stop=None,
              region=None,
              tid=None,
              until_eof=False,
              multiple_iterators=False,
              reference=None,
              end=None):
        """fetch reads aligned in a :term:`region`.
        See :meth:`AlignmentFile.parse_region` for more information
        on genomic regions.  :term:`reference` and `end` are also accepted for
        backward compatiblity as synonyms for :term:`contig` and `stop`,
        respectively.
        Without a `contig` or `region` all mapped reads in the file
        will be fetched. The reads will be returned ordered by reference
        sequence, which will not necessarily be the order within the
        file. This mode of iteration still requires an index. If there is
        no index, use `until_eof=True`.
        If only `reference` is set, all reads aligned to `reference`
        will be fetched.
        A :term:`SAM` file does not allow random access. If `region`
        or `contig` are given, an exception is raised.
        :class:`~pysam.FastaFile`
        :class:`~pysam.IteratorRow`
        :class:`~pysam.IteratorRow`
        :class:`~IteratorRow`
        :class:`IteratorRow`
        Parameters
        ----------
        until_eof : bool
           If `until_eof` is True, all reads from the current file
           position will be returned in order as they are within the
           file. Using this option will also fetch unmapped reads.
        multiple_iterators : bool
           If `multiple_iterators` is True, multiple
           iterators on the same file can be used at the same time. The
           iterator returned will receive its own copy of a filehandle to
           the file effectively re-opening the file. Re-opening a file
           creates some overhead, so beware.
        Returns
        -------
        An iterator over a collection of reads.
        Raises
        ------
        ValueError
            if the genomic coordinates are out of range or invalid or the
            file does not permit random access to genomic coordinates.
        """
        cdef int rtid, rstart, rstop, has_coord

        if not self.is_open:
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rstop = self.parse_region(contig, start, stop, region, tid,
                                                          end=end, reference=reference)

        # Turn of re-opening if htsfile is a stream
        if self.is_stream:
            multiple_iterators = False

        if self.is_bam or self.is_cram:
            if not until_eof and not self.is_remote:
                if not self.has_index():
                    raise ValueError(
                        "fetch called on bamfile without index")

            if has_coord:
                return IteratorRowRegionFiltered(
                    self, rtid, rstart, rstop,
                    multiple_iterators=multiple_iterators)
            else:
                if until_eof:
                    return IteratorRowAll(
                        self,
                        multiple_iterators=multiple_iterators)
                else:
                    # AH: check - reason why no multiple_iterators for
                    # AllRefs?
                    return IteratorRowAllRefs(
                        self,
                        multiple_iterators=multiple_iterators)
        else:
            if has_coord:
                raise ValueError(
                    "fetching by region is not available for SAM files")

            if multiple_iterators == True:
                raise ValueError(
                    "multiple iterators not implemented for SAM files")

            return IteratorRowAll(self,
                                  multiple_iterators=multiple_iterators)

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
#timeit("foo = pysam.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); foo.count(contig='2', start=2798760, stop=7009000, read_callback='all')", 'import pysam; import pysamx', number=1000)
#timeit("foo = pysam.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); list(foo.fetch(contig='2', start=2798760, stop=7009000))", 'import pysam; import pysamx', number=1000)
#0.6631929874420166

# new
#timeit("foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); foo.count(contig='2', start=2798760, stop=7009000)", 'import pysam; import pysamx', number=1000)
#timeit("foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); foo.bulkfetch(contig='2', start=2798760, stop=7009000, num_reads=475)", 'import pysam; import pysamx', number=1000)
#0.6480650901794434

# fetch2 addon
# timeit("foo = pysamx.AlignmentFile('../tests/data/NA12878.target_loci.sorted.bam', 'r'); list(foo.fetch2(contig='2', start=2798760, stop=2798780))", 'import pysam; import pysamx', number=1000)
