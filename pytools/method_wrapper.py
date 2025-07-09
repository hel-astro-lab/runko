import inspect


class MethodWrapper:
    """
    Utility to have objects which method calls can be customized.

    Method wrapper is constructed from action function which takes
    name of the called method as first argument.
    Methods of this class return function objects which when invoked
    forward the invokation to the underlying action function.

    In addition to action function pre and post functions can be specified.
    They work same as action, but are called before and after action function.

    Example:

    >>> def pre(x, y): print("pre: ", x, y)
    >>> def action(x, y): print("action: ", x, y)
    >>> def post(x, y): print("post: ", x, y)
    >>> f = MethodWrapper(action, pre=pre, post=post)
    >>> f.foo(10)
    pre: foo, 20
    action: foo, 20
    post: foo, 20
    >>> f.bar('42')
    pre: bar, 42
    action: bar, 42
    post: bar, 42
    """

    def __init__(self, action, pre=None, post=None):
        self._action = action
        self._pre = pre
        self._post = post


    def __getattr__(self, name):
        """
        This is called when accessed attribute is missing.
        Instead of raising AttributeError, we call our customization functions.

        Attributes starting with __ are treated as proper attributes
        which are not customizable, so AttributeError is raised for them.
        """

        if name.startswith("__"):
            msg = f"MethodWrapper invalid method name: {name}\n"
            msg += "In order to prevent unindented effects, attributes starting with `__`"
            msg += "can not be customized."

            raise AttributeError(msg)


        def wrapped_method_call(*vargs, **kwargs):
            if self._pre:
                try:
                    sig = inspect.signature(self._pre)
                    sig.bind(name, *vargs, **kwargs)
                except TypeError:
                    msg = f"In wrapped `{name}` method, "
                    msg += f"pre can not handle arguments: {name}, *{vargs}, **{kwargs}"
                    raise TypeError(msg)

            try:
                   sig = inspect.signature(self._action)
                   sig.bind(name, *vargs, **kwargs)
            except TypeError:
                   msg = f"In wrapped `{name}` method, "
                   msg += f"action can not handle arguments: {name}, *{vargs}, **{kwargs}"
                   raise TypeError(msg)

            if self._post:
                try:
                    sig = inspect.signature(self._post)
                    sig.bind(name, *vargs, **kwargs)
                except TypeError:
                    msg = f"In wrapped `{name}` method, "
                    msg += f"post can not handle arguments: {name}, *{vargs}, **{kwargs}"
                    raise TypeError(msg)

            self._pre(name, *vargs, **kwargs) if self._pre else None
            self._action(name, *vargs, **kwargs)
            self._post(name, *vargs, **kwargs) if self._post else None

        return wrapped_method_call
