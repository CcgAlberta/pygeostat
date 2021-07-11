# Contributing

When contributing to this repository, please first discuss the change you wish to make via github issues with the owners of this repository before making a change. 

If someone wish to report any bugs, this is encouraged to log the issue on github along with a zip file containing a notebook example of what is not working with clear explanations, and possibly screen shots that can be used to reproduce the issue. The owners will get notifications when an issue is raised.

This document outlines the conventions followed in pygeostat development.

## Pull Request Process

- Create a personal fork of the repository on your Github account
- Clone the fork on your local machine
- If you created your fork a while ago be sure to pull upstream changes into your local repository.
- Implement/fix your feature …
- Follow the conventions of the project
- If the project has tests run them
- Write or adapt tests as needed
- Add or change the documentation as needed
- From your fork open a pull request
- Your commit message should describe what the commit does to the code. The commit should mention whether your pull request closes an specific issue/bug (mentioning the issue) or makes an improvement
- The Pull Request will be merged if two of the owners have reviewed it and approved the pull request
- Once the pull request is approved and merged you can pull the changes from upstream to your local repo

## Naming Convention

 There are certain parameters/concepts in pygeostat that are frequently used for example, trimming limits. The purpose of pygeostat convention is to keep a consistent naming convention to refer to those parameters.

 The Python naming convention is used as the main guideline for naming modules, classes and functions. 

### 1. General

- Avoid using names that are too general or too wordy. Strike a good balance between the two.
- Bad: data_structure, my_list, info_map, dictionary_for_the_purpose_of_storing_data_representing_word_definitions
- Good: user_profile, menu_options, word_definitions
- Please don’t name things/variables “O”, “l”, or “I”
- When using CamelCase names, capitalize all letters of an abbreviation (e.g. HTTPServer)

### 2. Packages

- Package names should be all lower case
- When multiple words are needed, an underscore should separate them
- It is usually preferable to stick to 1 word names

### 3. Modules

- Module names should be all lower case
- When multiple words are needed, an underscore should separate them
- It is usually preferable to stick to 1 word names

### 4. Classes

- Class names should follow the UpperCaseCamelCase convention
- Python’s built-in classes, however are typically lowercase words
- Exception classes should end in “Error”

### 5. Global (module-level) Variables

- Global variables should be all lowercase
- Words in a global variable name should be separated by an underscore

### 6. Instance Variables

- Instance variable names should be all lower case
- Words in an instance variable name should be separated by an underscore
- Non-public instance variables should begin with a single underscore
- If an instance name needs to be mangled, two underscores may begin its name

### 7. Methods

- Method names should be all lower case
- Words in an method name should be separated by an underscore
- Non-public method should begin with a single underscore
- If a method name needs to be mangled, two underscores may begin its name

### 8. Method Arguments

- Instance methods should have their first argument named ‘self’.
- Class methods should have their first argument named ‘cls’

### 9. Functions

- Function names should be all lower case
- Words in a function name should be separated by an underscore

### 10. Constants

- Constant names must be fully capitalized
- Words in a constant name should be separated by an underscore

## File Header

 A list of file headers is provided here
       
    1. The first line should be shebang/hashbang line (specific to linux-like OS) to identify the absolute path to the bash interpreter. This makes it possible to run the file as a script invoking the interpreter implicitly. When a text file has a shebang, it is interpreted as an executable file. The following makes sure that the python executable will be found across unix-like distributions
    
        #!/usr/bin/env python

    2. Specifying the encoding of a Python file comes from PEP 0263 - Defining Python Source Code Encodings.

        # -*- coding: utf-8 -*-

    3. Next should be the docstring with a description. If the description is long, the first line should be a short summary that makes sense on its own, separated from the rest by a newline.
    
    All code, including import statements, should follow the docstring. Otherwise, the docstring will not be recognized by the interpreter, and you will not have access to it in interactive sessions (i.e. through obj.__doc__) or when generating documentation with automated tools

    4. Using __future__ to manage incompatible changes and introduction of new keywords

        from __future__ import absolute_import, division, print_function