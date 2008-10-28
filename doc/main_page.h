/***************************************************************************
 *            main_page.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! 

\file main_page.h
\brief Main page of Doxygen documentation



\mainpage

\section Introduction

%Ariadne is a C++ package for set-based analysis of dynamical and control systems, including reachability analysis and verification.

\section Documentation

Look at the <a href="tutorial.html">Tutorial</a> for a guide on how to use %Ariadne to perform reachability analysis of dynamic systems using the Python interface to %Ariadne. 
A more detailed description of %Ariadne's capabilities can be found in the documentation for the individual <a href="modules.html">Modules</a>.
Some information on the mathematical foundations of %Ariadne can be found in the <a href="pages.html">Related Pages</a>.


\section Download

The homepage of %Ariadne is <a href="http://ariadne.parades.rm.cnr.it/">http://ariadne.parades.rm.cnr.it/</a>.

You can check out the latest version on the Subversion repository by typing:

  \code svn checkout http://svn.parades.rm.cnr.it/ariadne/ariadne/trunk/ ariadne/ \endcode

To make the code documentation, change to the ariadne/trunk/ directory and type:

  \code doxygen \endcode

\section Requirements

To compile %Ariadne, you will need:
  - A C++ compiler (we recommend g++ version 4.0.2 or higher, which can be downloaded from <a href="http://gcc.gnu.org/">http://gcc.gnu.org</a>) 
  - A version of the 'make' utility (such as GNU make <a href="http://www.gnu.org/software/make/">http://www.gnu.org/software/make/</a>)

  - The Boost C++ Libraries (version 1.35.0 or higher) <a href="http://www.boost.org/">http://www.boost.org/</a>.

\subsection Optional Optional packages

To enable optional support for multiple-precision arithmetic, you will need:
  - The GNU Multiple-Precision Library (version 4.1.2 or higher)  <a href="http://www.swox.com/gmp/">http://www.swox.com/gmp/</a>.
  - MPFR Library (version 2.2.1 or higher) <a href="http://www.mpfr.org/">http://www.mpfr.org/</a>.

To use the plotting facilities, you will need:
  - cairo (version 1.7 or higher) <a href="http://www.cairographics.org/">http://www.cairographics.org/</a>.

To enable pop-up draphics windows, you will need:
  - the GTK+ toolkit (version 2.0 or higher) <a href="http://www.gtk.org/">http://www.gtk.org/</a>.

For the Python interface, you will also need:
  - Python (version 2.4 or higher) <a href="http://www.python.org/">http://www.python.org/</a>.

To make the source code documentation, you will need:
  - Doxygen (version 1.4.6 or higher is recommended) 
       <a href="http://www.stack.nl/~dimitri/doxygen/">http://www.stack.nl/~dimitri/doxygen/</a>.

To build the source code from the Subversion repository, you will also need
  - GNU autotools
      - autoconf <a href="http://www.gnu.org/software/autoconf/">http://www.gnu.org/software/autoconf/</a>.
      - automake <a href="http://www.gnu.org/software/automake/">http://www.gnu.org/software/automake/</a>.
      - libtool <a href="http://www.gnu.org/software/libtool/">http://www.gnu.org/software/libtool/</a>.
  - GNU m4 <a href="http://www.gnu.org/software/m4/">http://www.gnu.org/software/m4/</a>.
  - Perl <a href="http://www.perl.org/">http://www.perl.org/</a>.

All the above packages are available for Linux/Unix, Apple OSX and Microsoft Windows operating systems. 
They are standard packages for most Linux distributions.

\section installation Installation from a source tarball

Unpack the source tarball and change to the main directory
\code 
  tar -xzf ariadne-x.y.z.tar.gz
  cd ariadne-x.y.z
\endcode
From main source directory, type
\code 
  ./configure
  make
  make install
\endcode
The default installation directory for the library is /usr/local/lib/ (or /usr/local/lib64 on 64-bit machines), and for the Python interface is /usr/lib/python/2.x/site-packages/.
These defaults can be changed by using the --prefix flag to ./configure, as in
\code 
  ./configure --prefix=$HOME
\endcode

By default, the C++ library is compiled with support for double-precision floating-point numbers (64 bit). The flags \c --enable-mp can be used to enable support for multiple-precision floating-point numbers.
For example, to configure with multiple-precision floating-point support enabled, use
\code 
  ./configure --enable-mp
\endcode

\section svninstallation Installation from the Subversion repository

If installing from the Subversion repository, change to the ariadne/trunk/ directory and type 
\code 
  ./bootstrap
\endcode
And follow the directions above.

*/
