
.. _package-sec:

=====================================
Appendix: Creating Binary Packages
=====================================

This page gives instructions for how to use cmake to create a binary .dmg 
installer file for Mac OS X and a .deb or .rpm binary package for different 
distributions of linux.  This information is not relevant to most users, 
and is provided only as reference for core developers. 

The first step in creating a package, for any operating system, is to follow 
the instruction given on the page about :ref:`install-compile-cmake-sec` for 
installing dependencies and obtaining the source code. The remaining 
instructions given here assume that all dependencies are already installed, 
and that a copy of the source code has been installed within the directory 
structure described in the instructions for compiling from source.

**Mac OS X**

On Mac OS X, after installing all dependencies and installing a copy of the
source code in a repository named git/, one must:

    * Change directory to the pscf/cmake directory.

    * From the directory, enter::

          > cmake -DBUILD_DMG=1 -DCMAKE_INSTALL_PREFIX=. ../git

    * Then enter::

          > make -j4
          > make package
          > make package

The instruction to run "make package" twice is not a typo: This appears to be 
necessary to get the packaging utility to install all of the shared libraries
correctly.  The absence of an "make install" command is also intentional - 
the install target is not supported when BUILD_DMG is defined.  This procedure 
creates a standard Mac .dmg installer file with a name of the form 
pscf<version>-Darwin.dmg in the pscf/cmake directory.

**Linux (Fedora or Ubuntu)**

On a linux system, one must:

    * Change directory to the pscf/cmake directory.

    * From the directory, enter::

          > cmake -DCMAKE_INSTALL_PREFIX=.. ../git

    * Then enter::

          > make -j4
          > make install
          > make package

The "make install" command will install the software in the users pscf/
directory, in subdirectories named bin/, lib/ and share/. The "make package"
command should then create either a file named pscf<version>-linux.rpm, 
if the procedure is performed on a system such as Fedora that uses redhat 
package manager (rpm) files, or a file named pscf<version>-linux.deb if 
built on a system such as Ubuntu that uses .deb files. RPM packages can
only be built on a system that use .rpm packages and .deb packages can 
only be built on a system that use .deb files.

On a system that uses .rpm files, to check the RPM for detailed 
information (Metadata, Dependencies, and File Contents), enter::

   > rpm --info -qpR -qlvp pscf-1.0.0-Linux.rpm 

On a system that uses .deb package files, to check the .deb file for 
semi-detailed information, enter::

    # This extracts multiple files
    ar -vx pscf-1.0.0-Linux.deb
    # See the files that would be installed
    tar tvfz data.tar.gz 

