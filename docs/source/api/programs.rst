.. _programs:
.. public

Scripting CCG Executables
=========================

Pygeostat can be used for scripting parameter files and using CCG/GSLIB software (executable files) to create advanced geostatistical workflows. This is done using **Program** class.


Program Class
*************
.. autoclass:: pygeostat.programs.programs.Program

Get Parameter File
++++++++++++++++++
.. automethod:: pygeostat.programs.programs.Program.getparfile

Run Program
+++++++++++
.. automethod:: pygeostat.programs.programs.Program.run

Write Parameter File
++++++++++++++++++++
.. automethod:: pygeostat.programs.programs.Program.writepar

Run Programs in Parallel
************************
.. autofunction:: pygeostat.programs.program_utils.runparallel

Misc Program Utilities
**********************
.. autofunction:: pygeostat.programs.program_utils.parallel_function
.. autofunction:: pygeostat.programs.program_utils.rseed
.. autofunction:: pygeostat.programs.program_utils.rseed_list
.. autofunction:: pygeostat.programs.program_utils.parstr_kwargs
.. autofunction:: pygeostat.programs.program_utils.dedent_parstr


