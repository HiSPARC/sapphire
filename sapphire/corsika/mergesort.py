# coding: utf-8
from __future__ import division
import os
import tables
import tempfile
from sapphire.utils import pbar
from heapq import merge


class TableMergeSort(object):

    """ Sort a PyTables HDF5 table either in memory or on-disk """

    _iterators = []
    _BUFSIZE = 100000
    hdf5_temp = None

    def __init__(self, key, inputfile, outputfile=None, tempfile=None,
                 tablename='groundparticles', destination=None,
                 overwrite=False, progress=True):

        """ Initialize the class

        :param key: the name of the column which is to be sorted.
        :param inputfile: PyTables HDF5 input file.
        :param outputfile: optional PyTables HDF5 output file. If None the
            inputfile will be used for output.
        :param tempfile: optional PyTables HDF5 tempfile. If not specified
             a temp file will be created and removed when finished.
        :param tablename: the name of the table to sort.
        :param destination: optional name of the sorted table.
        :param overwrite: if True, overwrite destination table.
        :param progress: if True, show verbose output and progress.

        """
        self.key = key
        self.hdf5_in = inputfile
        self.table = self.hdf5_in.get_node('/%s' % tablename)
        self.description = self.table._v_dtype
        self.nrows = len(self.table)
        self.overwrite = overwrite
        self.progress = progress
        self.tempfile = tempfile

        if outputfile is None:
            if destination is not None:
                self.hdf5_out = inputfile
                self.destination = destination
            else:
                raise RuntimeError("Must specify either an outputfile or a "
                                   "destination table")
        else:
            self.hdf5_out = outputfile
            if destination is not None:
                self.destination = destination
            else:
                self.destination = tablename  # PATH OR NAME ?

        try:
            self.hdf5_out.get_node('/%s' % destination)
            if self.overwrite:
                self.hd5_out.remove_nove('/', self.destination, recursive=True)
            else:
                raise RuntimeError("Destination table exists and overwrite "
                                   "is False")
        except tables.NoSuchNodeError:
            self.outtable = self.hdf5_out.create_table('/', self.destination,
                                                       self.description,
                                                       expectedrows=self.nrows)

        self._calc_nrows_in_chunk()

        if self.nrows > self.nrows_in_chunk:
            if self.tempfile is None:
                self.tempfile_path = self._create_tempfile_path()
                self.hdf5_temp = tables.open_file(self.tempfile_path, 'w')
            else:
                self.hdf5_temp = tempfile
            if self.progress:
                parts = int(len(self.table) / self.nrows_in_chunk) + 1
                print "On disk mergesort in %d parts." % parts
        else:
            if self.progress:
                print "Table can be sorted in memory."

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            self.tempfile_path
        except AttributeError:
            return
        if self.hdf5_temp is not None:
            self.hdf5_temp.close()
            os.remove(self.tempfile_path)

    def sort(self):
        """Sort the table"""

        chunk = self.nrows_in_chunk
        nrows = self.nrows
        parts = int(nrows / chunk) + 1
        if parts == 1:
            if self.progress:
                print "Sorting table in memory and writing to disk."
            self._sort_chunk(self.outtable, 0, nrows)
        else:
            if self.progress:
                print "Sorting in %d chunks of %d rows:" % (parts, chunk)

            for idx, start in pbar(enumerate(range(0, nrows, chunk)),
                                   length=parts, show=self.progress):
                table_name = 'temp_table_%d' % idx
                table = self.hdf5_temp.create_table('/', table_name,
                                                    self.description,
                                                    expectedrows=chunk)
                iterator = self._sort_chunk(table, start, start + chunk)
                self._iterators.append(iterator)

            rowbuf = self.outtable._get_container(self._BUFSIZE)
            idx = 0

            if self.progress:
                print "Merging:"

            for keyedrow in pbar(merge(*self._iterators), length=nrows,
                                 show=self.progress):
                x, row = keyedrow

                if idx == self._BUFSIZE:
                    self.outtable.append(rowbuf)
                    self.outtable.flush()
                    idx = 0

                rowbuf[idx] = row.fetch_all_fields()
                idx += 1

            # store last lines in buffer
            self.outtable.append(rowbuf[0:idx])
            self.outtable.flush()

    def _iter_chunk(self, table):
        """Iterate over sorted rows of particles

        return key first

        """
        for row in table.iterrows():
            yield row[self.key], (row)

    def _sort_table(self, tablechunk):
        """Sort the chunk

        Sort using mergesort, this is much faster than quicksort and
        heapsort for large chunks.

        """
        return tablechunk.sort(order=self.key, kind='mergesort')

    def _sort_chunk(self, table, start, stop):
        """Read, sort and store a chunk of the input"""

        tablechunk = self.table.read(start=start, stop=stop)
        self._sort_table(tablechunk)
        table.append(tablechunk)
        table.flush()

        return self._iter_chunk(table)

    def _calc_nrows_in_chunk(self):
        """Determine maximum in memory table (sort) size

        This should be based on available memory because larger chunks are
        much faster. It is hard to determine the RAM that will be available,
        and the total RAM that will be used by Python.

        Currently not implemented.

        Each CORSIKA groundparticles row is about 36 bytes, so 1e7 rows are
        about 350 MB.

        """
        self.nrows_in_chunk = int(5e7)

    def _create_tempfile_path(self, temp_dir=None):
        """Create a temporary file, close it, and return the path"""

        f, path = tempfile.mkstemp(suffix='.h5', dir=temp_dir)
        os.close(f)
        return path
