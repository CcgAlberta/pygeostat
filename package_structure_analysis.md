Enter the path to the package directory (e.g., /path/to/pygeostat/): 
Found 64 Python files to analyze
Analyzing pygeostat_parameters.py...
Analyzing __init__.py...
Analyzing datautils/utils.py...
Analyzing datautils/fastcomps.py...
Analyzing datautils/desurvey.py...
Analyzing datautils/__init__.py...
Analyzing datautils/labels.py...
Analyzing programs/program_utils.py...
Analyzing programs/__init__.py...
Analyzing programs/programs.py...
Analyzing plotting/correlation_matrix_plot.py...
Analyzing plotting/contour_plot.py...
Analyzing plotting/utils.py...
Analyzing plotting/set_style.py...
Analyzing plotting/qq_plot.py...
Analyzing plotting/accuracy_plot.py...
Analyzing plotting/variogram_plot.py...
Analyzing plotting/variogram_plot_simulation.py...
Analyzing plotting/cmaps.py...
Analyzing plotting/variogram_utils.py...
Analyzing plotting/histogram_plot_simulation.py...
Analyzing plotting/histogram_plot.py...
Analyzing plotting/validation_plot.py...
Analyzing plotting/pit_plot.py...
Analyzing plotting/gaussian_mv.py...
Analyzing plotting/location_plot.py...
Analyzing plotting/scatter_plot.py...
Analyzing plotting/slice_plot.py...
Analyzing plotting/loadings_plot.py...
Analyzing plotting/__init__.py...
Analyzing plotting/drill_plot.py...
Analyzing plotting/subplots.py...
Analyzing plotting/image_grid.py...
Analyzing plotting/mds_plot.py...
Analyzing plotting/grid_slice_plot.py...
Analyzing plotting/probability_plot.py...
Analyzing plotting/export_image.py...
Analyzing fortran/wrappers.py...
Analyzing fortran/compile_linux.py...
Analyzing fortran/sourcedefs.py...
Analyzing fortran/compile.py...
Analyzing fortran/__init__.py...
Analyzing fortran/resource/compile_lapack.py...
Analyzing fortran/resource/compile_lapack_linux.py...
Analyzing multivariate/utils.py...
Analyzing multivariate/__init__.py...
Analyzing utility/filemanagement.py...
Analyzing utility/logging.py...
Analyzing utility/__init__.py...
Analyzing utility/deprecation.py...
Analyzing data/data.py...
Analyzing data/grid_definition.py...
Analyzing data/h5_pytables.py...
Analyzing data/__init__.py...
Analyzing data/h5_io.py...
Analyzing data/iotools.py...
Analyzing transformations/transformations.py...
Analyzing transformations/rotations.py...
Analyzing transformations/__init__.py...
Analyzing statistics/utils.py...
Analyzing statistics/cdf.py...
Analyzing statistics/postsim.py...
Analyzing statistics/kde.py...
Analyzing statistics/__init__.py...

=== Package Analysis Summary ===
Analyzing package directory: ./pygeostat
Total functions found: 431
Total classes found: 16

=== Classes ===

File: data/data.py
  Class: DataFile
    Methods (33):
      - __init__
      - __str__
      - __len__
      - __getitem__
      - __setitem__
      - _get_specialattr
      - _get_shape
      - head
      - tail
      - _get_columns
      - _get_locations
      - _get_info
      - rename
      - drop
      - _get_xyz
      - check_for_duplicate_cols
      - setcol
      - setvarcols
      - _get_nvar
      - __repr__
      - setcatdict
      - check_datafile
      - write_file
      - gscol
      - describe
      - copy
      - gendict
      - applydict
      - unique_cats
      - truncatenans
      - addcoord
      - infergriddef
      - spacing
  Class: DictFile
    Methods (3):
      - __init__
      - read_dict
      - write_dict

File: data/grid_definition.py
  Class: GridDef
    Methods (26):
      - __init__
      - __str__
      - __repr__
      - copy
      - _parse_grid_string
      - count
      - origin
      - _get_grid_array
      - _get_block_volume
      - convert_to_2d
      - extents
      - outline_points
      - index3d_to_index1d
      - index1d_to_index3d
      - get_index
      - get_index3d
      - get_coordinates
      - get_vertical_indices
      - get_slice_coordinate
      - get_slice_index
      - change_blocknum
      - change_blocksize
      - generate_grid_points
      - pad_grid
      - random_points
      - random_indices

File: data/h5_io.py
  Class: H5Store
    Methods (11):
      - __init__
      - __str__
      - __getitem__
      - __setitem__
      - __enter__
      - __exit__
      - _loadfile
      - _get_paths
      - close
      - datasets
      - iteritems

File: data/h5_pytables.py
  Class: DataFile
    Methods (1):
      - retrieve_hdf5

File: datautils/desurvey.py
  Class: Drillhole
    Methods (6):
      - __init__
      - __str__
      - delx
      - dely
      - delz
      - getxyz

File: plotting/gaussian_mv.py
  Class: GmmUtility
    Methods (17):
      - __init__
      - __get_mixtures_from_list
      - __get_mixtures_from_file
      - pdf_marginal
      - summary_plot
      - __bivariate_plot
      - __univariate_plot
      - bivariate_plot
      - get_moments
      - conditional_moments
      - univariate_conditional_plot
      - univariate_conditional_pdf
      - get_conditional_pdf
      - univariate_conditional_cdf
      - univariate_pdf_from_mixture_plot
      - univariate_pdf_from_mixture
      - get_modality_measure
  Class: UnivariateNormal
    Methods (3):
      - __init__
      - pdf
      - __pdf
  Class: MultivariateNormal
    Methods (3):
      - __init__
      - pdf
      - __pdf

File: programs/program_utils.py
  Class: GennyWithAMemry
    Methods (3):
      - __init__
      - _grnd
      - __call__

File: programs/programs.py
  Class: Program
    Methods (5):
      - __init__
      - __check_program
      - run
      - writepar
      - getparfile

File: pygeostat_parameters.py
  Class: Parameters
    Methods (15):
      - __init__
      - __setitem__
      - __getitem__
      - __iter__
      - __repr__
      - __str__
      - __delitem__
      - find_all
      - restore_defaults
      - describe
      - save
      - load
      - reset_systemdefault
      - set_systemdefault
      - get_systemdefault
  Class: PlotStyle
    Methods (10):
      - __init__
      - set_style
      - update
      - _check_style_parameter
      - __setitem__
      - _generate
      - update_mplrcParams
      - restore_defaults
      - restore_mplrcParams
      - get_systemdefault

File: transformations/rotations.py
  Class: Arrow3D
    Methods (2):
      - __init__
      - draw

File: transformations/transformations.py
  Class: Arcball
    Methods (9):
      - __init__
      - place
      - setaxes
      - constrain
      - constrain
      - down
      - drag
      - next
      - matrix

File: utility/deprecation.py
  Class: PygeostatDeprecationWarning
    Methods (0):

=== Functions by File ===

File: data/data.py
  - __init__
    [No documentation]
  - __str__
  - __len__
  - __getitem__
  - __setitem__
    [No documentation]
  - _get_specialattr
  - _get_shape
    [No documentation]
  - head
  - tail
  - _get_columns
  - _get_locations
    [No documentation]
  - _get_info
  - rename
  - drop
  - _get_xyz
    [No documentation]
  - check_for_duplicate_cols
  - setcol
  - setvarcols
  - _get_nvar
    [No documentation]
  - __repr__
    [No documentation]
  - setcatdict
  - check_datafile
  - write_file
  - gscol
  - describe
  - copy
  - gendict
  - applydict
  - dictapply
    [No documentation]
  - unique_cats
  - truncatenans
  - addcoord
  - infergriddef
  - spacing
  - __init__
    [No documentation]
  - read_dict
  - write_dict
  - ExampleData

File: data/grid_definition.py
  - __init__
    [No documentation]
  - __str__
  - __repr__
  - copy
  - _parse_grid_string
  - count
  - origin
  - _get_grid_array
  - _get_block_volume
  - convert_to_2d
  - extents
  - outline_points
  - index3d_to_index1d
  - isall
    [No documentation]
  - index1d_to_index3d
  - get_index
  - get_index3d
  - get_coordinates
  - get_vertical_indices
  - get_slice_coordinate
  - get_slice_index
  - change_blocknum
  - change_blocksize
  - generate_grid_points
  - pad_grid
  - random_points
  - random_indices

File: data/h5_io.py
  - _fixh5path
  - h5_combine_data
  - write_h5
  - read_h5
  - ish5dataset
  - __init__
    [No documentation]
  - __str__
  - __getitem__
  - __setitem__
  - __enter__
  - __exit__
  - _loadfile
  - _get_paths
  - close
  - datasets
  - iteritems

File: data/h5_pytables.py
  - retrieve_hdf5
  - write_h5_pytables
  - open_h5
  - read_h5_pytables
  - convert_rlzns_to_h5

File: data/iotools.py
  - read_file
  - _test_file_open
  - read_gslib
  - read_csv
  - compile_pygsb
  - isbinary
  - is_binary_string
    [No documentation]
  - read_gsb
  - _data_trim
  - write_gslib
  - write_csv
  - write_gsb
  - write_vtk
  - _data_fillnan
  - write_hvtk
  - file_nlines
  - writeout_gslib_gmm
  - readvarg

File: datautils/desurvey.py
  - __init__
    [No documentation]
  - __str__
    [No documentation]
  - delx
  - dely
  - delz
  - getxyz
  - set_desurvey
  - get_desurvey
  - desurvey
  - set_comps
  - get_comps
  - upscale_worker
  - length_in_interval
  - serial_upscale
  - length_in_interval
  - serial_get_comps
  - varupscale
  - fast_comps
  - length_in_interval
  - write_vtp

File: datautils/fastcomps.py
  - get_composites
  - row_length_in_interval
  - comp_value_and_length

File: datautils/labels.py
  - insert_real_idx
  - make_labels

File: datautils/utils.py
  - round_sigfig
  - fileheader
  - corrmatstr
  - slice_grid
  - slicescatter
  - fixpath
  - is_numeric
  - ensure_dir
  - ensure_path
  - nearest_eucdist

File: fortran/compile.py
  - quote
  - ensure_dir
  - ensure_path
  - load_ifortenv
  - buildintelf2py
  - relinkintel
  - buildgnuf2py
  - localbuild
    [No documentation]
  - relinkgnu
  - build_lapack
  - _buildcall
  - _checkerror
  - _updatestatus
    [No documentation]
  - extensioncleaner
  - assert_environment
  - build_custom_fortran
  - build_pygeostat_fortran

File: fortran/compile_linux.py
  - ensure_dir
  - ensure_path
  - buildgnuf2py
  - localbuild
    [No documentation]
  - _buildcall
  - _checkerror
  - _updatestatus
    [No documentation]
  - build_pygeostat_fortran_linux

File: fortran/resource/compile_lapack.py
  - build_lapack_gnu
    [No documentation]
  - build_lapack_intel
    [No documentation]

File: fortran/resource/compile_lapack_linux.py
  - build_lapack_gnu
    [No documentation]

File: fortran/wrappers.py
  - subsample

File: multivariate/utils.py
  - mds

File: plotting/accuracy_plot.py
  - accuracy_plot

File: plotting/contour_plot.py
  - contour_plot

File: plotting/correlation_matrix_plot.py
  - correlation_matrix_plot

File: plotting/drill_plot.py
  - drill_plot

File: plotting/export_image.py
  - export_image
  - exporttompl

File: plotting/gaussian_mv.py
  - _tickoff
  - setup_plot
  - __init__
    [No documentation]
  - __get_mixtures_from_list
  - __get_mixtures_from_file
  - pdf_marginal
  - summary_plot
  - __bivariate_plot
  - __univariate_plot
  - bivariate_plot
  - get_moments
  - conditional_moments
  - univariate_conditional_plot
  - univariate_conditional_pdf
  - get_conditional_pdf
  - univariate_conditional_cdf
    [No documentation]
  - univariate_pdf_from_mixture_plot
  - univariate_pdf_from_mixture
  - get_modality_measure
  - get_density
    [No documentation]
  - __init__
    [No documentation]
  - pdf
    [No documentation]
  - __pdf
    [No documentation]
  - __init__
    [No documentation]
  - pdf
    [No documentation]
  - __pdf
    [No documentation]

File: plotting/grid_slice_plot.py
  - grid_slice_plot

File: plotting/histogram_plot.py
  - histogram_plot
  - singlecdf

File: plotting/histogram_plot_simulation.py
  - histogram_plot_simulation
  - _itersimulated_data
    [No documentation]

File: plotting/image_grid.py
  - image_grid

File: plotting/loadings_plot.py
  - loadings_plot

File: plotting/location_plot.py
  - location_plot

File: plotting/mds_plot.py
  - mds_plot

File: plotting/pit_plot.py
  - pit_plot

File: plotting/probability_plot.py
  - probability_plot

File: plotting/qq_plot.py
  - qq_plot

File: plotting/scatter_plot.py
  - scatter_plot
  - scatter_plots
  - scatter_plots_lu
  - _handle_variables_wt
  - _tickoff

File: plotting/set_style.py
  - set_plot_style
  - wrapper_plot_function
    [No documentation]

File: plotting/slice_plot.py
  - slice_plot
  - _color_handling_gridded

File: plotting/subplots.py
  - subplots
  - subplots_clean

File: plotting/utils.py
  - titleoverlap
  - ticktitlecheck
  - supaxislabel
  - getminmax
  - get_contcbarargs
  - get_palette
  - get_cmap
  - catcmapfromcontinuous
  - get_label
  - clrmplmem
  - addxticks
  - get_supaxislocs
  - _get_cmap
  - setup_plot
  - _get_mpl_cbar_callback
    [No documentation]
  - cax_pos_callback
    [No documentation]
  - get_cbar_axis
  - format_subplot_axis
  - scalebar
  - smart_annotate
  - get_statblk
  - _set_stat_fontsize
  - applytickbins
  - format_plot
  - _format_axis_xy
  - _format_grid
  - _spatial_labels
  - _spatial_pointdata
  - _spatial_griddata
  - _spatial_slice
  - _spatial_orient2fig
  - _spatial_aspect
  - _format_tick_labels
  - _tickoverlap
  - tickcheck

File: plotting/validation_plot.py
  - validation_plot

File: plotting/variogram_plot.py
  - variogram_plot

File: plotting/variogram_plot_simulation.py
  - variogram_plot_simulation

File: plotting/variogram_utils.py
  - get_uniquevarids
  - trimylim

File: programs/program_utils.py
  - runparallel
  - pbar_update
    [No documentation]
  - __init__
    [No documentation]
  - _grnd
    [No documentation]
  - __call__
  - rseed_list
  - parallel_function
  - pbar_update
    [No documentation]
  - parstr_kwargs
  - dedent_parstr
  - temp_gslib_file

File: programs/programs.py
  - __init__
    [No documentation]
  - __check_program
    [No documentation]
  - run
  - writepar
  - getparfile

File: pygeostat_parameters.py
  - _validate_string
  - _validate_string_or_None
    [No documentation]
  - _validate_float
  - _validate_float_or_None
    [No documentation]
  - _validate_int
  - _validate_int_or_None
    [No documentation]
  - _validate_color
  - _validate_list
  - _validate_list_or_None
    [No documentation]
  - _validate_list_int
  - _validate_list_int_or_None
    [No documentation]
  - _validate_bool
  - _validate_dict_or_string
  - _validate_dict
  - _validate_dict_or_None
    [No documentation]
  - _validate_kde_or_color
  - _validate_griddef_or_None
  - _write_pardict
    [No documentation]
  - _parse_pardict
    [No documentation]
  - __init__
    [No documentation]
  - __setitem__
    [No documentation]
  - __getitem__
    [No documentation]
  - __iter__
  - __repr__
    [No documentation]
  - __str__
    [No documentation]
  - __delitem__
    [No documentation]
  - find_all
  - restore_defaults
  - describe
  - save
  - load
  - reset_systemdefault
  - set_systemdefault
  - get_systemdefault
  - __init__
    [No documentation]
  - set_style
  - update
    [No documentation]
  - _check_style_parameter
    [No documentation]
  - __setitem__
  - _generate
  - update_mplrcParams
  - restore_defaults
  - restore_mplrcParams
  - get_systemdefault
  - _parameters_initialize
    [No documentation]
  - _plot_style_initilizer
    [No documentation]

File: statistics/cdf.py
  - cdf
  - percentile_from_cdf
  - z_percentile
  - variance_from_cdf
  - stdev_from_cdf
  - build_indicator_cdf

File: statistics/kde.py
  - kde_scipy
  - kde_statsmodels_u
  - kde_statsmodels_m
  - kde_sklearn

File: statistics/postsim.py
  - postsim_multfiles

File: statistics/utils.py
  - weighted_mean
  - weighted_variance
  - weighted_skew
  - weighted_kurtosis
  - weighted_covariance
  - weighted_correlation
  - weighted_correlation_rank
  - near_positive_definite
  - accsim
  - accmik
  - _interval_responses

File: transformations/rotations.py
  - get_rotation_matrix
  - principalvectors
  - azmdip
  - principaldirs
  - __init__
    [No documentation]
  - draw
    [No documentation]
  - plotprincipalvectors
  - drawgsvectorwidget
  - view_gslib_axes
  - drawellipsoid
  - drawgsaniswidget
  - view_gslib_axes

File: transformations/transformations.py
  - identity_matrix
  - translation_matrix
  - translation_from_matrix
  - reflection_matrix
  - reflection_from_matrix
  - rotation_matrix
  - rotation_from_matrix
  - scale_matrix
  - scale_from_matrix
  - projection_matrix
  - projection_from_matrix
  - clip_matrix
  - shear_matrix
  - shear_from_matrix
  - decompose_matrix
  - compose_matrix
  - orthogonalization_matrix
  - affine_matrix_from_points
  - superimposition_matrix
  - euler_matrix
  - euler_from_matrix
  - euler_from_quaternion
  - quaternion_from_euler
  - quaternion_about_axis
  - quaternion_matrix
  - quaternion_from_matrix
  - quaternion_multiply
  - quaternion_conjugate
  - quaternion_inverse
  - quaternion_real
  - quaternion_imag
  - quaternion_slerp
  - random_quaternion
  - random_rotation_matrix
  - __init__
  - place
  - setaxes
  - constrain
  - constrain
  - down
  - drag
  - next
  - matrix
  - arcball_map_to_sphere
  - arcball_constrain_to_axis
  - arcball_nearest_axis
  - vector_norm
  - unit_vector
  - random_vector
  - vector_product
  - angle_between_vectors
  - inverse_matrix
  - concatenate_matrices
  - is_same_transform
  - _import_module

File: utility/deprecation.py
  - _generate_deprecation_message
    [No documentation]
  - warn_deprecated
  - deprecated
  - deprecate
    [No documentation]
  - finalize
    [No documentation]
  - finalize
    [No documentation]
  - finalize
    [No documentation]
  - wrapper
    [No documentation]

File: utility/filemanagement.py
  - mkdir
  - rmdir
  - rmfile
  - get_executable
  - __remove_temp_dir
    [No documentation]
  - list_executable

File: utility/logging.py
  - printerr
  - log_progress

=== External Imports Used ===
- data.data
- data.grid_definition
- datautils.labels.make_labels
- datautils.utils
- matplotlib._cm_listed
- matplotlib.backends.backend_pdf
- matplotlib.collections
- matplotlib.colors
- matplotlib.lines
- matplotlib.patches
- matplotlib.path
- matplotlib.pyplot
- matplotlib.ticker
- mpl_toolkits.axes_grid1
- mpl_toolkits.mplot3d
- multivariate.utils
- pandas.api.types
- pyevtk.hl
- pygeostat.fortran.subsample
- rpy2.robjects
- rpy2.robjects.packages
- scipy.stats
- scriptnotifier.utils
- sklearn.neighbors
- statistics.cdf
- statistics.kde
- statistics.utils
- statsmodels.nonparametric.kde
- statsmodels.nonparametric.kernel_density
- utility.logging
- xml.etree.cElementTree


=== Fortran Usage Analysis ===

Fortran-related imports:
- from fortran import pygsb
- from fortran import pygsb
- from pygeostat.fortran.subsample import subsample
- from pygeostat.fortran.subsample import subsample

Functions potentially using Fortran:

File: data/h5_io.py
  - write_h5 (line 96)
    Doc: Write data to an HDF5 file using the python package H5PY. The file is appended to and in

File: data/iotools.py
  - compile_pygsb (line 210)
    Doc: Compiles 'pygeostat/fortran/src/pygsb.f90' using 'pygeostat/fortran/compile.py'
  - read_gsb (line 270)
    Doc: Reads in a CCG GSB (GSLIB-Binary) file.

File: fortran/compile.py
  - buildintelf2py (line 210)
    Doc: Compiles a single PYD using the intel compiler
  - relinkintel (line 319)
    Doc: Relinks a single PYD assuming the required code is already compiled and only the link step needs
  - buildgnuf2py (line 376)
    Doc: Compiles a single PYD using the gfortran compiler
  - relinkgnu (line 476)
    Doc: Relinks a single PYD assuming the required code is already compiled and only the link step needs
  - assert_environment (line 658)
    Doc: Ensure that the specified compiler(s) can be found. If compiler == `'gnu'`, first check if a
  - build_custom_fortran (line 732)
    Doc: This function is intended to allow arbitrary f2py extensions to be constructed using the tooling
  - build_pygeostat_fortran (line 802)
    Doc: This function builds the f2py extensions with minimal interaction from the user. The goal is to

File: fortran/compile_linux.py
  - buildgnuf2py (line 62)
    Doc: Compiles a single PYD using the gfortran compiler
  - build_pygeostat_fortran_linux (line 197)
    Doc: This function builds the f2py extensions with minimal interaction from the user. The goal is to

File: fortran/resource/compile_lapack.py
  - build_lapack_gnu (line 6)
  - build_lapack_intel (line 18)

File: fortran/resource/compile_lapack_linux.py
  - build_lapack_gnu (line 11)

File: fortran/wrappers.py
  - subsample (line 4)
    Doc: When datasets become too large to efficiently use, a sub-sample maybe extracted and used in

File: plotting/histogram_plot_simulation.py
  - histogram_plot_simulation (line 23)
    Doc: histogram_plot_simulation emulates the pygeostat histogram_plot program as a means of checking histogram

File: programs/programs.py
  - __init__ (line 30)

File: utility/filemanagement.py
  - get_executable (line 75)
    Doc: Gets a collection of executable files from a protected repository using an access token. Note that in order to use this function, git needs to be installed on the target computer.
