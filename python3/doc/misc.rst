.. _misc:

*************
Miscellaneous
*************

.. currentmodule:: HTSeq

.. doctest:: 
   :hide:

   >>> import HTSeq


``FileOrSequence``
==================

.. class:: FileOrSequence( filename_or_sequence )

   This class is a a canvenience wrapper around a file.
   
   The construcutor takes one argument, which may either be a string,
   which is interpreted as a file name (possibly with path), or a
   connection, by which we mean a text file opened for reading, or 
   any other object that can provide an iterator over strings 
   (lines of the file).
   
   The advantage of passing a file name instead of an already opened file
   is that if an iterator is requested several times, the file will be
   re-opened each time. If the file is already open, its lines can be read
   only once, and then, the iterator stays exhausted.      
   
   Furthermore, if a file name is passed that end in ".gz" or ".gzip"
   (case insensitive), it is transparently gunzipped.
   
   
   .. attribute:: FileOrSequence.fos
      
      The argument passed to the constructor, i.e., a filename or a sequence
      
   .. attribute:: FileOrSequence.line_no
   
      The line number (1-based) of the most recently read line. Initially None.
      
   .. method:: FileOrSequence.get_line_number_string( )
   
      Returns a string describing the position in the file. Useful for error messages.   


Version
=======

.. attribute:: __version__

   a string containing the current version

