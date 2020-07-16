#!/usr/bin/env python
# -*- coding: utf-8 -*

'''filemanagement.py: Contains basic utility functions for creating/removing files and directories'''
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function, unicode_literals

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


def get_executable(source='gslib', access_token=None, clean=True):

    '''
    Gets a collection of executable files from a protected repository using an access token.
    
    Parameters:
        source (str): gslib or CCG as the source of software
        access_token (str): An access token to authorize access to the target private repository
        clean (bool): Option to clean the executable directory prior to upload the files from the target private repository

    **Examples**

    .. code-block:: python
        
        import pygeostat as gs
        gs.get_executable(source='CCG', clean=True)

    |

    '''

    import subprocess
    import os
    import shutil
    import sys
    import getpass

    if access_token is None:
        access_token = getpass.getpass('Please provide a valid access token')

    # check if git is available
    try:
        popenobject_ = subprocess.Popen('git')
        popenobject_.terminate()
    except Exception:
        print ('Unsuccessful! Could not access git. Please make sure git is installed.')

    # Target direcory for executable
    target_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), r'..', r'executable'))

    # Clean the directory first before uploading the executable files
    if clean:
        for file_name in os.listdir(target_dir):
            rmfile(os.path.join(target_dir,file_name))

    # Create a temporary directory to get
    temp_dir = 'temp1618'
    mkdir(temp_dir)

    if source.lower() == 'gslib':
        print('The software is available under gslib license agreement (http://www.gslib.com)')
        
        command = 'git clone https://CcgUser:{access_token}@github.com/CcgAlberta/GslibRepository.git {target}'.format(access_token=access_token,target = temp_dir)

        if 'win' in sys.platform:
            source_dir = os.path.join(temp_dir, 'Windows64')
        else:
            source_dir = os.path.join(temp_dir, 'Linux64')
    else:
        print('The software is available under CCG software term of use (http://www.ccgalberta.com/software-terms-of-use)')
        command = 'git clone https://CcgUser:{access_token}@github.com/CcgAlberta/CcgSoftware.git {target}'.format(access_token=access_token,target = temp_dir)

        print('The software is suited for Microsoft platform')
        source_dir = os.path.join(temp_dir, 'CcgFortranExecutable')
        
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