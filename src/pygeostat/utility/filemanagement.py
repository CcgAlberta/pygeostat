#!/usr/bin/env python
# -*- coding: utf-8 -*

'''filemanagement.py: Contains basic utility functions for creating/removing files and directories'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------

__author__ = 'pygeostat development team'
__date__ = '2018'
__version__ = '1.0.0'


def mkdir(dirnames, verbose=False):
    '''
    Make directory(s).

    Parameters:
        dirnames(str or list): directory(s) to create
        verbose(bool): print to screen if an error is encountered,
            including if the dirnames already exists
    '''
    from os import makedirs
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        try:
            makedirs(dirname)
        except:
            if verbose:
                print('{} not created - likely exists already!'.format(dirname))


def rmdir(dirnames, verbose=False):
    '''
    Remove directory(s).

    Parameters:
        dirnames(str or list): directory(s) to remove
        verbose(bool): print to screen if an error is encountered,
            including if the dirnames has already been removed
    '''
    from shutil import rmtree
    if isinstance(dirnames, str):
        dirnames = [dirnames]
    for dirname in dirnames:
        try:
            rmtree(dirname)
        except:
            if verbose:
                print('{} not removed - likely does not exist!'.format(dirname))


def rmfile(filenames, verbose=False):
    '''
    Remove file(s).

    Parameters:
        filenames(str): file to remove
        verbose(bool): print to screen if an error is encountered,
            including if the file has already been removed
    '''
    from os import remove
    if isinstance(filenames, str):
        filenames = [filenames]
    for filename in filenames:
        try:
            remove(filename)
        except:
            if verbose:
                print('{} not removed - likely does not exist!'.format(filename))


def get_executable(source='gslib', access_token=None, clean=False):

    '''
    Gets a collection of executable files from a protected repository using an access token. Note that in order to use this function, git needs to be installed on the target computer.
    
    Parameters:
        source (str): gslib or CCG as the source of software.
        access_token (str): An access token to authorize access to the target private repository for CCG software. Access token is available for CCG members and can be found at CCG knowledge base. 
        clean (bool): Option to clean the executable directory prior to upload the files from the target private repository. Note that choosing this option will delete the existing executable files.

    **Examples**

    Installing GSLIB executable files

    .. code-block:: python
        
        import pygeostat as gs
        gs.get_executable(source='GSLIB')
    
    .. image:: ./figures/GslibGetFiles.png
    |

    Installing CCG software

    .. code-block:: python
        
        import pygeostat as gs
        gs.get_executable(source='CCG', clean=True)

    .. image:: ./figures/CcgGetFiles.png
    |

    Using the installed software

    .. code-block:: python

        # Load example data included in pygeostat
        data_file = gs.ExampleData('cluster')

        # Initialize the Program object
        nscore = gs.Program('nscore', getpar=True)

        # Run normal score transformation
        parstr = """ Parameters for NSCORE
                     *********************
        START OF PARAMETERS:
        {datafile}                -  file with data
        3  1  2  3                -  number of variables and columns
        5                         -  column for weight, 0 if none
        0                         -  column for category, 0 if none
        0                         -  number of records if known, 0 if unknown
        -1.0e21   1.0e21          -  trimming limits
        0                         -transform using a reference distribution, 1=yes
        ../histsmth/histsmth.out  -file with reference distribution.
        1   2   0                 -  columns for variable, weight, and category
        101                       -maximum number of quantiles, 0 for all
        {outfl}                   -file for output
        {trnfl}                   -file for output transformation table
        """

        pars = dict(datafile=data_file.flname,
                    outfl = 'nscore.out',
                    trnfl = 'nscore.trn')
        nscore.run(parstr=parstr.format(**pars), liveoutput=True)
    '''

    import subprocess
    import os
    import shutil
    import sys
    import getpass

    if access_token is None and source.upper() == 'CCG':
        access_token = getpass.getpass('Please provide a valid access token')

    # check if git is available
    try:
        popenobject_ = subprocess.Popen('git')
        popenobject_.terminate()
    except Exception as exception_content:
        print ('Unsuccessful! Could not access git. Please make sure git is installed.')
        raise Exception(str(exception_content)) from None

    # Target direcory for executable
    target_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..', r'executable'))

    # Clean the directory first before uploading the executable files
    if clean:
        try:
            for file_name in os.listdir(target_dir):
                rmfile(os.path.join(target_dir,file_name))
        except Exception:
            print('Unable to clean the previously installed files!')

    # Create a temporary directory to get
    temp_dir = 'temp1618'
    mkdir(temp_dir)

    if source.lower() == 'gslib':
        print('The software is available under gslib license agreement (http://www.gslib.com)')
        
        command = 'git clone https://github.com/CcgAlberta/GslibRepository.git {target}'.format(access_token=access_token,target = temp_dir)

        if 'win' in sys.platform:
            source_dir = os.path.join(temp_dir, 'Windows64')
        else:
            source_dir = os.path.join(temp_dir, 'Linux64')
    elif source.lower() == 'ccg':
        print('The software is available under CCG software term of use (http://www.ccgalberta.com/software-terms-of-use)')
        if 'win' in sys.platform:
            command = 'git clone https://CcgUser:{access_token}@github.com/CcgAlberta/CcgSoftware.git {target}'.format(access_token=access_token,target = temp_dir)

            print('The software is suited for Microsoft platform')
        else:
            command = 'git clone https://CcgUser:{access_token}@github.com/CcgAlberta/CcgSoftware_Linux.git {target}'.format(access_token=access_token,target = temp_dir)

            print('The software is suited for Linux platform')

        source_dir = os.path.join(temp_dir, 'CcgFortranExecutable')
    else:
        raise ValueError('Wrong source was provided.')
        
    try:
        main_process = subprocess.Popen(command, cwd='.',
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                shell=True,
                                bufsize=1,
                                universal_newlines=True)
        rc=main_process.wait()
    except Exception as exception_content:
        try:
            logs = main_process.stdout.read().replace('\r', '')
            main_process.stdout.close()
        except UnboundLocalError:
            logs = str(exception_content)

        logs = logs.replace(access_token, 'retracted')
        
        __remove_temp_dir(temp_dir)
        
        raise Exception(logs) from None

    try:
        for file_name in os.listdir(source_dir):
            
            srcfile = os.path.join(source_dir,file_name)
            
            destfile = os.path.join(target_dir,file_name)
            
            shutil.copy(srcfile,destfile)

    except Exception as exception_content:
        logs = 'Unable to retrive data from the repository. Please make sure correct access token has been provided!' + str(exception_content)
        raise Exception(logs) from None

    finally:
        __remove_temp_dir(temp_dir)


def __remove_temp_dir(temp_dir):
    import os
    import sys
    import stat
    from os import path
    import shutil

    # Remove the temporary directory
    for root, dirs, files in os.walk(temp_dir):  
        for dir in dirs:
            os.chmod(path.join(root, dir), stat.S_IRWXU)
        for file in files:
            os.chmod(path.join(root, file), stat.S_IRWXU)
    shutil.rmtree(temp_dir)


def list_executable():
  """
  PRovides a list installed executable files under pygeostat
  """  
  import pathlib
  
  exec_dir = pathlib.Path(__file__).parent.parent.joinpath("executable")
  
  exe_files = [f.name for f in exec_dir.glob("*.exe")]
  
  return exe_files