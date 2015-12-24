#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2010--2014 Nico Schlömer
#
# This file is part of matplotlib2tikz.
#
# matplotlib2tikz is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# matplotlib2tikz is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# matplotlib2tikz.  If not, see <http://www.gnu.org/licenses/>.
#
from os import path
from matplotlib import pyplot as pp
import matplotlib2tikz
import testfunctions

# XXX: There seems to be an issue with the legends that do not appear in the
# pdf versions. The old legends are carried over from the previous run and are
# transfered to the next run. I believe this is a problem with the testing
# method, because if one adds print text in _draw_legend pyplot returns the
# correct legends when this test is run


def _main():
    # get command line arguments
    test_list = _parse_options()

    tex_file_path = './tex/acid.tex'

    # directory where all the generated files will end up
    data_dir = './data'

    # how to get from the LaTeX file to the data
    tex_relative_path_to_data = '../data'

    figure_width = '7.5cm'

    # open file for writing
    file_handle = open(tex_file_path, 'w')

    write_document_header(file_handle, figure_width)

    # Get all function names from the testfunctions module.
    function_strings = dir(testfunctions)
    # Remove all functions that start with '_' and that are not
    # in a certain exclude list (uuuugly).
    exclude_list = {'mpl', 'np', 'pp'}
    tmp = []
    for s in function_strings:
        if s[0] != '_' and s not in exclude_list:
            tmp.append(s)
    function_strings = tmp

    if not test_list is None:  # actually treat a sublist of test_functions
        # remove duplicates and sort
        test_list = sorted(set(test_list))
    else:
        # all indices
        test_list = xrange(0, len(function_strings))

    for k in test_list:
        print 'Test function %d (%s)...' % (k, function_strings[k]),
        pp.cla()
        pp.clf()
        # plot the test example
        comment = getattr(testfunctions, function_strings[k])()

        # plot reference figure
        pdf_path = data_dir + '/test' + repr(k) + '.pdf'
        pp.savefig(pdf_path)

        pdf_path = path.join(tex_relative_path_to_data,
                             path.basename(pdf_path)
                             )
        # Open figure, insert PDF
        file_handle.write('% test plot ' + str(k) + '\n'
                          '\\begin{figure}%\n'
                          '\\centering%\n'
                          '\\begin{tabular}{cc}\n'
                          '\includegraphics[width=\\figwidth]'
                          '{' + str(pdf_path) + '}%\n'
                          '&\n'
                          )
        # convert to TikZ
        tikz_path = data_dir + '/test%r.tex' % k
        tikz_tex_path = path.join(tex_relative_path_to_data,
                                  path.basename(tikz_path)
                                  )
        try:
            matplotlib2tikz.save(
                tikz_path,
                figurewidth=figure_width,
                tex_relative_path_to_data=tex_relative_path_to_data,
                show_info=False
                )
            file_handle.write('\\input{%s}\n' % tikz_tex_path)
        except:
            file_handle.write('% fail\n')

        # Close the figure
        file_handle.write('\\end{tabular}\n'
                          '\\caption{' + str(comment) + ' (test ID '
                          + str(k) + ').}%\n'
                          '\\end{figure}\\clearpage\n\n'
                          )
        print 'done.'
    write_document_closure(file_handle)
    file_handle.close()
    return


def write_document_header(file_handle, figure_width):
    '''Write the LaTeX document header to the file.
    '''
    file_handle.write('\\documentclass[landscape]{scrartcl}\n'
                      '\\usepackage{graphicx}\n'
                      '\\usepackage{pgfplots}\n'
                      '\\usepgfplotslibrary{groupplots}\n'
                      '\\pgfplotsset{compat=newest}\n\n'
                      '\\newlength\\figwidth%\n'
                      '\\setlength\\figwidth{' + figure_width + '}\n\n'
                      '\\begin{document}\n\n'
                      )
    return


def write_document_closure(file_handle):
    '''Write the LaTeX document closure to the file.
    '''
    file_handle.write('\\end{document}')
    return


def write_file_comparison_entry(file_handle,
                                pdf_path,
                                tikz_path,
                                test_id,
                                comment
                                ):
    '''Write the Tikz vs. PDF comparison figures to the LaTeX file.
    '''
    file_handle.write('% test plot ' + str(test_id) + '\n'
                      '\\begin{figure}%\n'
                      '\\centering%\n'
                      '\\begin{tabular}{cc}\n'
                      '\includegraphics[width=\\figwidth]'
                      '{' + str(pdf_path) + '}%\n'
                      '&\n'
                      '\input{' + str(tikz_path) + '}%\n'
                      '\\end{tabular}\n'
                      '\\caption{' + str(comment) + ' (test ID '
                      + str(test_id) + ').}%\n'
                      '\\end{figure}\\clearpage\n\n'
                      )
    return


def _parse_options():
    '''Parse input options.'''
    import argparse
    parser = argparse.ArgumentParser(description=
                                     'Acid test for matplotlib2tikz.'
                                     )
    parser.add_argument('--tests', '-t',
                        metavar='TEST_INDICES',
                        nargs='+',
                        type=int,
                        help='tests to perform'
                        )
    args = parser.parse_args()
    return args.tests


if __name__ == '__main__':
    # execute the test
    _main()
