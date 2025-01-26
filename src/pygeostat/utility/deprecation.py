#
# ATTRIBUTATION TO MATPLOTLIB:
# The template for the PygeostatDeprecationWarning class and related functions within this
# deprecation module is adapted from the MatplotLibDeprecationWarning and related module,
# as it provides a flexible and extensible approach to deprecation warnings.
#
# The Matplotlib License is is available at: https://matplotlib.org/users/license.html,
# which was accessed on March 26, 2018. Relevant to this application are the items:
#
# 1. This LICENSE AGREEMENT is between the Matplotlib Development Team (“MDT”), and the
# Individual or Organization (“Licensee”) accessing and otherwise using matplotlib software
# in source or binary form and its associated documentation.
#
# 2. Subject to the terms and conditions of this License Agreement, MDT hereby grants Licensee
# a nonexclusive, royalty-free, world-wide license to reproduce, analyze, test, perform and/or
# display publicly, prepare derivative works, distribute, and otherwise use matplotlib 2.2.2 alone
# or in any derivative version, provided, however, that MDT’s License Agreement and MDT’s notice
# of copyright, i.e., “Copyright (c) 2012-2013 Matplotlib Development Team; All Rights Reserved”
# are retained in matplotlib 2.2.2 alone or in any derivative version prepared by Licensee.
#
# 3. In the event Licensee prepares a derivative work that is based on or incorporates
# matplotlib 2.2.2 or any part thereof, and wants to make the derivative work available
# to others as provided herein, then Licensee hereby agrees to include in any such work a
# brief summary of the changes made to matplotlib 2.2.2.
#
# Following this instruction, it should be noted that the name occurences of matplotlib
# were changed to pygeostat within this module. This was peformed to avoid confusion
# with the application of this deprecation warning within pygeostat, so it is not
# confused with deprecation relating to matplotlib. The underlying code, however,
# is entirely attributed to the Matplotlib Development team and is so accredited here.

import warnings
import functools


class PygeostatDeprecationWarning(UserWarning):
    """
    A class for issuing deprecation warnings for Pygeostat users.

    In light of the fact that Python builtin DeprecationWarnings are ignored
    by default as of Python 2.7 (see link below), this class was put in to
    allow for the signaling of deprecation, but via UserWarnings which are not
    ignored by default.

    https://docs.python.org/dev/whatsnew/2.7.html#the-future-for-python-2-x

    Attributation to Matplotlib:
    The template for the PygeostatDeprecationWarning class and related functions within this
    deprecation module is adapted from the MatplotLibDeprecationWarning and related module,
    as it provides a flexible and extensible approach to deprecation warnings. Refer to
    the header of deprecation.py for a full attributation and licence of Matplotlib.
    """
    pass


gsDeprecation = PygeostatDeprecationWarning


def _generate_deprecation_message(since, message='', name='',
                                  alternative='', pending=False,
                                  obj_type='attribute',
                                  addendum=''):

    if not message:

        if pending:
            message = (
                'The %(name)s %(obj_type)s will be deprecated in a '
                'future version.')
        else:
            message = (
                'The %(name)s %(obj_type)s was deprecated in version '
                '%(since)s.')

    altmessage = ''
    if alternative:
        altmessage = ' Use %s instead.' % alternative

    message = ((message % {
        'func': name,
        'name': name,
        'alternative': alternative,
        'obj_type': obj_type,
        'since': since}) +
        altmessage)

    if addendum:
        message += addendum

    return message


def warn_deprecated(
        since, message='', name='', alternative='', pending=False,
        obj_type='attribute', addendum=''):
    """
    Used to display deprecation warning in a standard way.

    Parameters
    ----------
    since : str
        The release at which this API became deprecated.

    message : str, optional
        Override the default deprecation message.  The format
        specifier `%(name)s` may be used for the name of the function,
        and `%(alternative)s` may be used in the deprecation message
        to insert the name of an alternative to the deprecated
        function.  `%(obj_type)s` may be used to insert a friendly name
        for the type of object being deprecated.

    name : str, optional
        The name of the deprecated object.

    alternative : str, optional
        An alternative function that the user may use in place of the
        deprecated function.  The deprecation warning will tell the user
        about this alternative if provided.

    pending : bool, optional
        If True, uses a PendingDeprecationWarning instead of a
        DeprecationWarning.

    obj_type : str, optional
        The object type being deprecated.

    addendum : str, optional
        Additional text appended directly to the final message.

    Examples
    --------

        Basic example::

            # To warn of the deprecation of "Pygeostat.name_of_module"
            warn_deprecated('1.4.0', name='Pygeostat.name_of_module',
                            obj_type='module')

    """
    message = _generate_deprecation_message(
        since, message, name, alternative, pending, obj_type)

    warnings.warn(message, gsDeprecation, stacklevel=1)


def deprecated(since, message='', name='', alternative='', pending=False,
               obj_type=None, addendum=''):
    """
    Decorator to mark a function or a class as deprecated.

    Parameters
    ----------
    since : str
        The release at which this API became deprecated.  This is
        required.

    message : str, optional
        Override the default deprecation message.  The format
        specifier `%(name)s` may be used for the name of the object,
        and `%(alternative)s` may be used in the deprecation message
        to insert the name of an alternative to the deprecated
        object.  `%(obj_type)s` may be used to insert a friendly name
        for the type of object being deprecated.

    name : str, optional
        The name of the deprecated object; if not provided the name
        is automatically determined from the passed in object,
        though this is useful in the case of renamed functions, where
        the new function is just assigned to the name of the
        deprecated function.  For example::

            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative object that the user may use in place of the
        deprecated object.  The deprecation warning will tell the user
        about this alternative if provided.

    pending : bool, optional
        If True, uses a PendingDeprecationWarning instead of a
        DeprecationWarning.

    addendum : str, optional
        Additional text appended directly to the final message.

    Examples
    --------

        Basic example::

            @deprecated('1.4.0')
            def the_function_to_deprecate():
                pass

    """

    def deprecate(obj, message=message, name=name, alternative=alternative,
                  pending=pending, addendum=addendum):
        import textwrap

        if not name:
            name = obj.__name__

        if isinstance(obj, type):
            obj_type = "class"
            old_doc = obj.__doc__
            func = obj.__init__

            def finalize(wrapper, new_doc):
                try:
                    obj.__doc__ = new_doc
                except (AttributeError, TypeError):
                    # cls.__doc__ is not writeable on Py2.
                    # TypeError occurs on PyPy
                    pass
                obj.__init__ = wrapper
                return obj
        else:
            obj_type = "function"
            if isinstance(obj, classmethod):
                func = obj.__func__
                old_doc = func.__doc__

                def finalize(wrapper, new_doc):
                    wrapper = functools.wraps(func)(wrapper)
                    wrapper.__doc__ = new_doc
                    return classmethod(wrapper)
            else:
                func = obj
                old_doc = func.__doc__

                def finalize(wrapper, new_doc):
                    wrapper = functools.wraps(func)(wrapper)
                    wrapper.__doc__ = new_doc
                    return wrapper

        message = _generate_deprecation_message(
            since, message, name, alternative, pending,
            obj_type, addendum)

        def wrapper(*args, **kwargs):
            warnings.warn(message, gsDeprecation, stacklevel=2)
            return func(*args, **kwargs)

        old_doc = textwrap.dedent(old_doc or '').strip('\n')
        message = message.strip()
        new_doc = (('\n.. deprecated:: %(since)s'
                    '\n    %(message)s\n\n' %
                    {'since': since, 'message': message}) + old_doc)
        if not old_doc:
            # This is to prevent a spurious 'unexected unindent' warning from
            # docutils when the original docstring was blank.
            new_doc += r'\ '

        return finalize(wrapper, new_doc)

    return deprecate
