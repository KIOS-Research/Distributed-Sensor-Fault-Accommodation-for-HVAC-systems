Name:   read_dxf
Author: Steven Michael (smichael@ll.mit.edu)
Date:   3/10/2005

############################################################

The following code allows MATLAB to read in a text DXF 
file.  It does not load any color or texture information


############################################################

Compilation:

The base distribution includes binary MATLAB functions for Linux and
Windows.  The functions were compiled with Matlab R14.  I have not
tried them with other versions, but they should work.  The Windows
compiler used is MS Visual Studio.NET 2003 (v7.1).  The linux compiler
is gcc version 2.96 (RedHat Linux 7.3).

To compile in Linux, simply type "make" in the base directory.
Some variables may need to be changed in the Makefile depending
upon MATLAB version and C compiler.

I've included the MS Visual Studio project for compiling the
"read_dxf.dll" file as well.

############################################################


Example Usage:

dxf = read_dxf('dxf_filename');
surfdxf(dxf);
