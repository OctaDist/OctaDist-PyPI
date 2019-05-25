====================
Installing and Usage
====================

you can use OctaDist on any platform as long as Python 3 and
require packages (dependencies) are installed.

PyPI
----

OctaDist CLI version is also available on Python package index library, 
which can be found at https://pypi.org/project/octadist.

The end-user can use `pip`, a Python package-management system, 
to find and install OctaDist and other dependencies on OS at the same time.

The following commands might be useful:

- Install program::

   pip install octadist

- Upgrade to the latest version::

   pip install --upgrade octadist

- Upgrade/downgrade to a certain version, for example, version 2.5.3::

   pip install --upgrade octadist==2.5.3

More details on installing Python package can be found its official website: 
https://packaging.python.org/tutorials/installing-packages.

Anaconda 
--------

OctaDist is also available on Anaconda cloud server.
The channel of OctaDist is at https://anaconda.org/rangsiman/octadist.

It can be installed on system using command::

    conda install -c rangsiman octadist 


The platforms that OctaDist-Conda supported: [![Anaconda-Server Badge][Conda-platform-badge]][Conda-platform-link]


Checking
--------

You can check if the OctaDist is installed correctly, for example:

- Using ``pip``::

    pip search octadist

- Try importing ``octadist`` package::

    import octadist
    print(octadist.__version__)