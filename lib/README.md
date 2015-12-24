# matplotlib2tikz

[![Code Health](https://landscape.io/github/nschloe/matplotlib2tikz/master/landscape.png)](https://landscape.io/github/nschloe/matplotlib2tikz/master)
[![Documentation Status](https://readthedocs.org/projects/matplotlib2tikz/badge/?version=latest)](https://readthedocs.org/projects/matplotlib2tikz/?badge=latest)


This is matplotlib2tikz, a Python script for converting matplotlib figures into
native Pgfplots (TikZ) figures.

To download matplotlibtikz, go to its page on GitHub
https://github.com/nschloe/matplotlib2tikz.

### Installation
1. Place the matplotlib2tikz script in a directory where Python can find
    it (e.g., `$PYTHONPATH`).
    You can install it systemwide with

        sudo python setup.py install

   or place the script `matplotlib2tikz.py` into the directory where you
   intend to use it.

2. Make sure that your LaTeX installation includes the packages

  * TikZ (aka PGF, >=2.00), and
  *  Pgfplots (>=1.3).

### Usage

1. Generate your matplotlib plot as usual.

2. Instead of pyplot.show(), invoke matplotlib2tikz by

        tikz_save( 'myfile.tikz' );

   to store the TikZ file as myfile.tikz. Load the libary with:

        from matplotlib2tikz import save as tikz_save

      _Optional:_
      The scripts accepts several options, for example `height`, `width`,
      `encoding`, and some others. Invoke by

        tikz_save( 'myfile.tikz', figureheight='4cm', figurewidth='6cm' )


     IMPORTANT:
     Height and width must be set large enough; setting it too low it may
     result in a LaTeX compilation failure such as
        - Dimension Too Large, or
        - Arithmetic Overflow
      (see information about these errors in the manual of Pgfplots)

      To specify the dimension of the plot from within the LaTeX document, try
        
        tikz_save('myfile.tikz',
                  figureheight = '\\figureheight',
                  figurewidth = '\\figurewidth'
                  )

      and in the LaTeX source
 
        \newlength\figureheight
        \newlength\figurewidth
        \setlength\figureheight{4cm}
        \setlength\figurewidth{6cm}
        \input{myfile.tikz}

3. Add the contents of `myfile.tikz` into your LaTeX source code; a convenient
   way of doing so is to use `\input{/path/to/myfile.tikz}`. Also make sure
   that at the header of your document the packages TikZ and Pgfplots are
   included:

        \usepackage{tikz}
        \usepackage{pgfplots}

   Optionally, to use features of the latest Pgfplots package (as of
   Pgfplots 1.3), insert

        \pgfplotsset{compat=newest}

### License

matplotlib2tikz is published under the GNU Lesser General Public License v3.0.

### Contributing

If you experience bugs, would like to contribute, have nice examples of what
matplotlib2tikz can do, or if you are just looking for more information, then
please visit the web page of matplotlib2tikz
https://github.com/nschloe/matplotlib2tikz.
