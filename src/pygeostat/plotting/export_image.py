#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Function to export a matplotlib figure object as multiple files simply while minimizing
whitespace surrounding the final image"""
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import os
import fileinput
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from PIL import PngImagePlugin


def export_image(output_file=None, fltype=None, pad=0.03, dpi=300, custom=None, pdfpages=None, delim=None,
              Metadata=True, **kwargs):
    """
    This function exports a figure with the specified file name and type(s) to the specified
    location. Multiple file types can be exported at once. Avoids the use of plt.tight_layout()
    which can behave odd and minimizes whitespace on the edge of figures.

    Note:
        This function is typically called within plotting functions but can be used on its own.

    Extensions are not required in the ``output_file`` argument. They will be added according to what
    ``fltype`` is set to. The default output file types are png and eps. However, if extensions
    are provided, they will be used, provided that the argument ``fltype`` is not passed. The
    ``custom`` argument  provides extra flexibility if the default settings of this function are
    not desired. If the ``custom`` functionality is the only desired output, ``fltype`` can be set
    to ``False`` to prevent additional exports.

    PS and EPS files need to have their font definitions fixed so that they will be called properly
    which is done automatically if they are used.

    Figures can also be appended to an existing ``mpl.backends.backend_pdf.PdfPages`` object passed
    with the ``pdfpages`` argument. These objects are used to created multi-page PDF documents.

    Parameters:
        output_file (str or list): Details the file location and file name in one parameter or a list of
            files to export with or without file extensions. If not file extensions are provided,
            the parameter ``fltype`` will need to be specified if its defaults are not desired
        fltype (str, list, bool): The file extension or list of extensions. See plt.savefig() docs
            for which file types are supported. Can set to ``False`` to prevent normal
            functionality when ``custom`` settings are the only desired output.
        pad (float): The amount of padding around the figure
        dpi (int): The output file resolution
        custom (dict):  Indicates a custom dpi and file extension if an odd assortment of files are
            needed
        pdfpages (mpl.object): A multi-page PDF file object created by
            ``mpl.backends.backend_pdf.PdfPages``. If a PdfPages object is passed, the figure is
            exported to it in addition to the other files if specified. Use this method to generate
            PDF files with multiple pages.
        delim (str): delimiter in the output_file str passes that indicates different types of files. set
            to `None` to ensure filenames with spaces can be used.
        kwargs: Any other permissible keyword arguments to send to plt.savefig() (e.g., )

    Examples:
        A simple call using the default settings exporting two images at 300 dpi:

        .. code-block:: python

            gs.export_image(output_file='../Figures/histplt')
            # 'histplt.png' and 'histplt.eps' are exported


        A call specifying only a png file and altering the dpi setting and setting the background
        to transparent (via \**kwargs):

        .. code-block:: python

            gs.export_image(output_file='../Figures/histplt.png', dpi=250, transparent=True)
            #'histplt.png' is exported in '../Figures/'

        A call using only the custom argument:
        
        .. code-block:: python

            gs.export_image(output_file='../Figures/histplt', fltype=False, custom={600:'png', 200:'png'})
            # 'histplt_600.png' and 'histplt_200.png' are exported in '../Figures/'

        A call using a combination of arguments:

        .. code-block:: python

            gs.export_image(output_file='../Figures/histplt', custom={600:'jpg'})
            'histplt.png' and 'histplt.eps' at 300 dip in addition to histplt_600.jpg' are exported in '../Figures/'

        A call using a more complicated combination of arguments:

        .. code-block:: python

            gs.export_image(output_file=['../Figures/png/histplt', '../Figures/eps/histplt'], custom={600:'png'})
            #'histplt.png' @ 300 dpi and 'histplt_600.png' @ 600 dpi are placed in '../Figures/png/' while 'histplt.eps' is placed in '../Figures/eps/'

        Create a PDFPages matplotlib object and save the figure to it:

        .. code-block:: python
            :linenos:

            from matplotlib.backends.backend_pdf import PdfPages
            pdfpages = PdfPages('output_file.pdf')
            gs.location_plot(data_file)
            gs.export_image(pdfpages=pdfpages)
            pdfpages.close()

    """
    def exporttompl(flname, dpi, Metadata, **kwargs):
        """
        Small function to pass information to matplotlib for exporting images for easier
        iteration
        """
        _, flext = os.path.splitext(flname)

        # Set up a pdfpage automatically so we can insert metadata
        if flext.lower() == '.pdf':
            pdffig = PdfPages(flname)
            plt.savefig(pdffig, format='pdf', bbox_inches='tight', pad_inches=pad,
                        dpi=dpi, **kwargs)
            if Metadata:
                meta_dictionary = pdffig.infodict()
                for meta in Metadata:
                    meta_dictionary[meta] = Metadata[meta]
            pdffig.close()
        else:
            plt.savefig(flname, bbox_inches='tight', pad_inches=pad, dpi=dpi, **kwargs)
            # Correct EPS or PS font definition if needed

            if flext in ['.eps', '.ps']:
                # Grab a font for now
                font = mpl.rcParams['font.family']
                for line in fileinput.FileInput(flname, inplace=1):
                    # Get the font code used in the EPS or PS file
                    if "/FontName" in line:
                        dump = copy.deepcopy(line)
                        dump = dump.split('/')
                        dump = dump[2].split(' ')
                        font = dump[0]
                    line = line.replace("/b'%s'" % font, "/%s" % font)
                    print(line, end='')
            if flext == '.png' and Metadata is not None:
                im = Image.open(flname)
                meta = PngImagePlugin.PngInfo()

                for x in Metadata:
                    meta.add_text(x, Metadata[x])
                im.save(flname, "png", pnginfo=meta)
            # if flext.lower() is in ['.jpeg', '.jpeg']:

    # --------------------------------------------------------------------------------------
    # START OF MAIN FUNCTION
    # --------------------------------------------------------------------------------------
    # Sanity check
    if not output_file and not pdfpages:
        raise ValueError("One or both of `output_file` or `pdfpages` needs to be specified.")

    # set up basic metadata info
    if Metadata is True:
        Metadata = {'Author': os.getcwd()}

    if output_file:
        # If input oufl is a string, convert it to a list
        if delim is not None:
            if isinstance(output_file, str):
                if delim in output_file:
                    return ValueError("Output files need to be delimited by a single space not a comma")
                output_file = output_file.split(delim)
        # If input fltype is a string, convert it to a list
        if isinstance(fltype, str):
            if ',' in fltype:
                return ValueError("File types need to be delimited by a single space not a comma")
            fltype = fltype.split(' ')
        # Get the extensions from output_file if they exist
        exts = []
        output_files = []
        if not isinstance(output_file, list):
            output_file = [output_file]
        for fl in output_file:
            outf, ext = os.path.splitext(fl)
            if ext != '':
                exts.append(ext.replace(".", ""))
            output_files.append(outf)
        if len(exts) > 0 and fltype is None:
            fltype = exts
        output_file = output_files
        # Set fltype to defaults if True or None is passed and they aren't contained in output_file
        if fltype is None or (isinstance(fltype, bool) and fltype):
            fltype = ['png', 'eps']
        # Extend the output_file if needed
        if len(output_file) == 1 and fltype:
            output_file = output_file * len(fltype)
        # Export figures required through regular usage
        if not isinstance(fltype, bool) or fltype:
            for i, ext in enumerate(fltype):
                flname = '%s.%s' % (output_file[i], ext)
                exporttompl(flname, dpi, Metadata, **kwargs)
        # Export the custom files
        if custom is not None:
            for cus_dpi in custom:
                flname = '%s_%s.%s' % (output_file[0], cus_dpi, custom[cus_dpi])
                exporttompl(flname, cus_dpi, Metadata, **kwargs)

    # Export the the PDFPages object if required
    if pdfpages:
        # Try and use the variable fig if it exists (works for gs plotting functions as fig is set)
        try:
            pdfpages.savefig(fig, bbox_inches='tight', pad_inches=pad, dpi=dpi, **kwargs)
        except:
            # Grab the mpl figure if it isn't already set
            fig = plt.gcf()
            pdfpages.savefig(fig, bbox_inches='tight', pad_inches=pad, dpi=dpi, **kwargs)
