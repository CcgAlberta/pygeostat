#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''programs.py: Contains methods for running GSLIB and CCG programs'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import subprocess
import os
import time
import warnings
from ..utility.logging import printerr
from .program_utils import dedent_parstr
from .. pygeostat_parameters import Parameters

# Seed the random number generator for spacing parallel calls
#   and seeding ACORNI random number generators
import random
random.seed()


class Program(object):
    '''Base class containing routines for running GSLIB programs'''

    def __init__(self, program=None, parstr=None, parfile='temp', getpar=None,
                 nogetarg=False, defaultdict={}):
        self.program = self.__check_program(program)
        self.parstr = parstr
        # nogetarg is a flag for FORTRAN programs that call READ
        # instead of using GETARG - this requires that a PIPE be used
        # rather than traditional command line arguments
        self.nogetarg = nogetarg
        self.parfile = parfile
        self.defaultdict = defaultdict
        if getpar is None:
            getpar = Parameters['config.getpar']
        if getpar:
            self.getparfile()

    def __check_program(self, program):

        try:
            popenobject_ = subprocess.Popen(program)
            popenobject_.terminate()
            return program
        except FileNotFoundError:
            if Parameters['config.use_shipped_executable']:
                print('Unable to find the provided program! Trying to use the executable pool shipped with pygeostat!')
                exec_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..', r'executable'))
                _, program = os.path.split(program)
                return os.path.join(exec_dir, program)
            else:
                raise FileNotFoundError('The executable file does not exist. \nConsider setting the parameter `config.use_shipped_executable` to True to use the shipped executable files.')


    def run(self, parstr=None, parfile=None, program=None, nogetarg=None,
            filehandle=None, logfile=None, testfilename=None, pardict=None,
            quiet=False, liveoutput=None, chdirpath=None):
        """
        Runs a GSLIB style program using the subprocess module and prints the
        output. On an error, the output is printed and an exception is raised

        The only required parameters are `program` and `parstr`, all other
        parameters are optional.

        Parameters:
            parstr (str): parameters, taken from self if None
            parfile (str or callable): name of parameter file to create, or a callable function
                that returns a unique path/parfile.par string to use for this parfile.
            program (str): name of GSLIB/CCG program to run, taken from self if None
            nogetarg: uses a pipe + communicate for the call instead of arguments
            filehandle: handle for file to write program output to
            logfile (str): filename for a log file to write program output to. If file already
                exists it will overwrite the file. If filehandle is passed then logfile will
                be ignored.
            testfilename (list): name(s) to check for availability prior to
                executing the program
            quiet: if quiet, don't print to let the user know it is calling
            liveoutput (bool): live update the ipython notebook or calling
                script with the output of the called program
            chdirpath (str): Some programs have files they use which prevent
                running multiple instances of the same program from the same
                directory. This can be used to chang the path the program is
                called from to prevent conflicts like this (for example, with kt3d_lva)

        Examples:
            An example for setting up a few key words and running histplt

            >>> histpltpar = '''                  Parameters for HISTPLT
            ...           **********************
            ...
            ... START OF PARAMETERS:
            ... {datafl}          -file with data
            ... {varcol}   0                        -   columns for variable and weight
            ... -1.0     1.0e21              -   trimming limits
            ... {outfl}                   -file for PostScript output
            ...  0.0      -20.0               -attribute minimum and maximum
            ... -1.0                         -frequency maximum (<0 for automatic)
            ... 20                           -number of classes
            ... 0                            -0=arithmetic, 1=log scaling
            ... 0                            -0=frequency,  1=cumulative histogram
            ... 0                            -   number of cum. quantiles (<0 for all)
            ... 3                            -number of decimal places (<0 for auto.)
            ... {varname}                                                    -title
            ... 1.5                          -positioning of stats (L to R: -1 to 1)
            ... -1.1e21                      -reference value for box plot
            ... '''
            >>>
            >>> histplt = gs.Program(program='histplt', parfile='histplt.par')
            >>>
            >>> histplt.run(parstr=histpltpar.format(datafl=datafl.flname,
            ...                                      varcol=datafl.gscol('Bitumen'),
            ...                                      varname='Bitumen',
            ...                                      outfl='histplt_bitumen.ps'))

        """
        if program is None:
            program = self.program
        else:
            program = self.__check_program(program)
        if parstr is None:
            parstr = self.parstr
        if nogetarg is None:
            nogetarg = self.nogetarg
        if parfile is None:
            parfile = self.parfile
        # Update the working dictionary
        workingdict = self.defaultdict.copy()
        if pardict is not None:
            workingdict.update(pardict)

        # Logic Checks and figureing out defaults
        if liveoutput is None and filehandle is None and logfile is None:
            liveoutput = True
        elif liveoutput is None and logfile:
            liveoutput = True
        elif liveoutput is None and filehandle:
            liveoutput = False
        if filehandle:
            if logfile:
                warnings.warn('filehandle was passed so logfile is being ignored')
                logfile = None

        if isinstance(logfile, str):
            consoleout = open(logfile, 'w')
        elif logfile:
            raise TypeError('logfile must be string')

        # parfile could be a callable that generates tempfiles just before the program is run
        # e.g. in parallel
        try:
            parfile = parfile()
            shouldclean = True
        except TypeError:
            shouldclean = False
            pass

        program_call = [program, parfile]
        if not quiet:
            print('Calling: ', program_call)

        # remove leading whitespace on each line of parstr so it can be nested correctly in a
        # function or code-folded in the notebook!
        parstr = dedent_parstr(parstr)

        # Format the parstr first to catch KerError for printing if need be
        try:
            parstr = parstr.format(**workingdict)
        except KeyError as err:
            errormsg = 'The following Dictionary Key was not found: ' + str(err)
            raise KeyError(errormsg)

        # Create the parameter file
        try:
            if chdirpath is None:
                try:
                    with open(parfile, 'w', newline='', encoding='utf8') as pf:
                        pf.write(parstr)
                except:
                    with open(parfile, 'w') as pf:
                        pf.write(parstr)
            else:
                try:
                    with open(chdirpath + '/' + parfile, 'w', newline='', encoding='utf8') as pf:
                        pf.write(parstr)
                except:
                    with open(chdirpath + '/' + parfile, 'w') as pf:
                        pf.write(parstr)
        except:
            errormsg = 'Could not create parameter file: ' + parfile
            raise Exception(errormsg)
        if testfilename is not None:
            time.sleep(random.random())
            # Test if this file is open and available up to a max of 5 times
            # before just running anyway
            nattempts = 0
            while nattempts < 5:
                if not isinstance(testfilename, str):
                    try:
                        with open(testfilename, 'r'):
                            testfilename.readline()
                        break
                    except:
                        nattempts += 1
                        time.sleep(1.0 + random.random() * nattempts)
                else:
                    try:
                        for testfile in testfilename:
                            with open(testfile, 'r'):
                                testfile.readline()
                        break
                    except:
                        nattempts += 1
                        time.sleep(1.0 + random.random() * nattempts)

        if nogetarg:
            process = subprocess.Popen([program], stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            # message = process.communicate(input=bytes(parfile+'\n', 'ascii'))
            try:  # py3.x
                message = process.communicate(input=bytes(parfile + ' \n', 'ascii'))
            except:
                try:  # py2.7
                    message = process.communicate(input=unicode(parfile + ' \n', 'ascii'))
                except:
                    raise Exception('Something went wrong with calling ' + str(program))
            # print(message[0].decode('ascii'))
            if filehandle is None:
                if logfile is None:
                    print(''.join([msg.decode('ascii') for msg in message]))
                else:
                    consoleout.write(''.join([msg.decode('ascii') for msg in message]))
                    if liveoutput:
                        print(''.join([msg.decode('ascii') for msg in message]))
            else:
                filehandle.write(''.join([msg.decode('ascii') for msg in message]))
                if liveoutput:
                    print(''.join([msg.decode('ascii') for msg in message]))

        else:
            # The standandard for single-threaded calls to executables, allows for live output
            # and catches common errors
            try:
                p = subprocess.Popen(program_call, cwd=chdirpath,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, bufsize=1,
                                     universal_newlines=True)
                fulloutput = ''
                for line in p.stdout:    # still not working correctly output not
                    if liveoutput:
                        print(line, end='')  # constantly fed to the ipython notebook..
                    if "forrtl: The pipe has been ended." in line:
                        print("**************************************************************\n"
                              "WARNING: it looks like the program %s does not use `getarg`\n"
                              "Consider the argument `nogetarg=True` for the call to %s.run()!\n"
                              "**************************************************************\n" %
                              (self.program, self.program))
                        if shouldclean:
                            os.remove(parfile)
                        return
                    fulloutput += line
                if filehandle:
                    filehandle.write(fulloutput)
                elif logfile:
                    consoleout.write(fulloutput)
            except subprocess.CalledProcessError:
                # If the program crashed, raise an exception
                errormsg = 'ERROR executing: ' + str(program_call)
                if logfile:
                    consoleout.write(errormsg)
                    consoleout.close()
                raise Exception(errormsg)
            except FileNotFoundError:
                # Notify that perhaps this executable doesn't exist in the path ...
                errormsg = 'ERROR calling: ' + str(program_call[0]) + '. Check path to exe'
                if logfile:
                    consoleout.write(errormsg)
                    consoleout.close()
                raise Exception(errormsg)
            except:
                # This subprocess.check_output is required to call exe's in parallel with
                # gs.runparallel()
                # filehandles does not work in parallel
                try:
                    # Call the subprocess and buffer the output
                    subprocess_output = subprocess.check_output(program_call, cwd=chdirpath)
                except subprocess.CalledProcessError:
                    # If the program crashed, raise an exception
                    errormsg = 'ERROR executing: ' + str(program_call)
                    if logfile:
                        consoleout.write(errormsg)
                        consoleout.close()
                    raise Exception(errormsg)
                except FileNotFoundError:
                    # Notify that perhaps this executable doesnt exist in the path ...
                    errormsg = 'ERROR calling: ' + str(program_call[0]) + '. Check path to exe'
                    if logfile:
                        consoleout.write(errormsg)
                        consoleout.close()
                    raise Exception(errormsg)
                # Print out the decoded output from the subprocess if successful
                try:
                    decoded_output = subprocess_output.decode("utf-8")
                except UnicodeDecodeError:
                    print('Warning: UnicodeDecodeError decoding subprocess output with utf-8.')
                    decoded_output = ''
                if filehandle:
                    filehandle.write(decoded_output)
                if logfile:
                    consoleout.write(decoded_output)
                else:
                    print(decoded_output)
                # Check for signal of GSLIB crash in output in case the program
                # doesn't return a non-zero error code
                if 'FAILED' in decoded_output:
                    print('FAILED found in gslib output - logging.')
                    errormsg = 'ERROR executing: ' + ' '.join(program_call)
                    if logfile:
                        consoleout.write(errormsg)
                        consoleout.close()
                    raise Exception(errormsg)
        # Close logfile if it is open
        if logfile and not consoleout.closed:
            consoleout.close()
        if shouldclean:
            os.remove(parfile)

    def writepar(self, parstr=None, parfile=None, pardict=None):
        '''
        Writes out the parameter file without running the program, which can be
        helpful for checking

        Parameters:
            parstr (str): This is the parameter file string initiated with the program.
            parfile (str): file name or path to save the parameter file to.
            pardict (dict): Dictionary for the variables in the parameter file.
        '''
        # Keep the default parameter string if no new string provided
        if parstr is None:
            parstr = self.parstr
        if parfile is None:
            parfile = self.parfile
        # Update the working dictionary
        workingdict = self.defaultdict.copy()
        if pardict is not None:
            workingdict.update(pardict)
        # Create the parameter file
        try:
            with open(parfile, 'w', newline='', encoding='utf8') as pf:
                pf.write(parstr.format(**workingdict))
        except:
            errormsg = 'Could not create parameter file: ' + parfile
            raise Exception(errormsg)

    def getparfile(self, quiet=True):
        """
        Get the parfile from this program by copying to the clipboard or
        printing the parfile. This replaces the need to pre-execute a program
        to get the parfile, but relies on CCG programs being properly configured
        to generate the correct parfile upon first execution.

        This function requires `pyperclip`. This is a dependency of pygeostat, but
        can be installed with::

            > pip install pyperclip

        Parameters:
            quiet (bool) : The function will copy the parameter file to the
                clipboard with the block quotes as default (quiet = True)

        """
        import tempfile
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as dirname:
            # Run program to automatically generate a parfile
            try:
                p = subprocess.Popen(self.program, stdout=subprocess.PIPE,
                                     stdin=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     cwd=dirname)
                try:  # py3.x
                    message = p.communicate(input=bytes(' \n', 'ascii'))
                except:
                    try:
                        message = p.communicate(input=u' \n')
                    except:
                        raise Exception('Something went wrong with calling ' + self.program)
            except FileNotFoundError:
                raise FileNotFoundError('The program ' + self.program + ' was not found, '
                                        'check definitions!')
            # assume that the parfile is called "self.program + '.par' "
            if os.path.isfile(dirname + '/' + self.program + '.par'):
                parfl = dirname + '/' + self.program + '.par'
            else:  # grab the newest .par file that was created in this directory
                import glob
                try:
                    parfl = max(glob.iglob(dirname + '/*.par'), key=os.path.getctime)
                except ValueError:
                    raise ValueError('%s did not create a parfile! You may require the parfile '
                                     'from the ccg knowledgebase. Try: '
                                     'https://ccgsrv.geostats.ualberta.ca/ccgkb/doku.php?'
                                     'do=search&id=%s' % (self.program, self.program))
            # read the par file
            with open(parfl, 'r') as fl:
                parstr = fl.read()
            os.remove(parfl)

        # the panda's copy stopped working... no idea why. pyperclip seems to do a
        # better job anyway!
        if quiet:
            try:
                import pyperclip
                progstr = self.program
                if '/' in progstr:
                    progstr = progstr[progstr.rfind('/') + 1:]
                if '.' in progstr:
                    progstr = progstr[:progstr.rfind('.')]
                pyperclip.copy('parstr = """' + parstr + '"""')
                print(parfl + ' has been copied to the clipboard')
            except (ModuleNotFoundError, ImportError):
                print('Looks like pyperclip is not installed, use:\n',
                      '>>> pip install pyperclip\n',
                      'to install. In the meantime, here is the par file!\n')
                print(parstr)
        else:
            print(parstr)
