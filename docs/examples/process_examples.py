import time
import glob
import datetime
import os

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import HTMLExporter
import codecs
from numpy import random


def run_notebook_once(book, outpath="./"):
    """ Execute `book` optionally in a different path, overwriting `book` with the output """
    try:
        # read the notebook relative to this dir
        with open(book) as f:
            nb = nbformat.read(f, as_version=4)
        # ep = ExecutePreprocessor(timeout=99999, kernel_name="python3")
        # # run the notebook in the outdir
        # ep.preprocess(nb, {"metadata": {"path": outpath}, 
        #                                 "kernelspec": { "display_name": "Python",
        #                                                 "language": "python",
        #                                                 "name": "python3"
        #                                                 }
        #                                 }
        #                 )
        # write a temporary notebook
        tempbook = str(random.randint(1000, size=1)[0]) + book
        with open(tempbook, 'wt') as f:
            nbformat.write(nb, f)
        # read in the temporary notebook before removing
        output = nbformat.read(tempbook, as_version=4)
        os.remove(tempbook)
        # export to an HTML
        exporter = HTMLExporter()
        output, _ = exporter.from_notebook_node(output)
        i = book.index('.ipynb')
        html = book[:i] + '.html'
        codecs.open(html, 'w', encoding='utf-8').write(output)
    except:
        print("`rundemos` failed on %s" % book)
        raise


def main():
    """ Execute the demos in this directory """

    notebooks = glob.glob("*.ipynb")
    for book in notebooks:
        run_notebook_once(book)


if __name__ == "__main__":
    # invoke this with `python rundemos.py (test5)`
    start_time = time.time()
    main()
    elapsed_time = str(datetime.timedelta(seconds=(time.time() - start_time)))
    print("Demo execution time {}".format(elapsed_time))
