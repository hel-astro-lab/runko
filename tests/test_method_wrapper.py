import unittest
import pytools


class method_wrapper(unittest.TestCase):

    def test_simple_logging(self):

        log = []

        def action(x: str):
            log.append(f"actionA({x})")
            log.append(f"actionB({x})")

        e = pytools.MethodWrapper(action)

        e.foo()
        self.assertEqual(log, ["actionA(foo)", "actionB(foo)"])
        e.foo()
        e.foo()
        self.assertEqual(log, 3 * ["actionA(foo)", "actionB(foo)"])


    def test_simple_logging_with_pre_and_post(self):

        log = []

        def pre(x: str):
            log.append(f"pre({x})")

        def action(x: str):
            log.append(f"actionA({x})")
            log.append(f"actionB({x})")

        def post(x: str):
            log.append(f"post({x})")

        e = pytools.MethodWrapper(action, pre=pre, post=post)

        e.foo()
        self.assertEqual(log, ["pre(foo)", "actionA(foo)", "actionB(foo)", "post(foo)"])


    def test_logging_with_vargs_only(self):

        log = []

        def pre(x: str, *vargs):
            msg = f"pre({x},"

            for arg in vargs:
                msg += (f"{arg},")

            log.append(msg + ")")

        def action(x: str, *vargs):
            msg = f"action({x},"

            for arg in vargs:
                msg += (f"{arg},")

            log.append(msg + ")")

        def post(x: str, *vargs):
            msg = f"post({x},"

            for arg in vargs:
                msg += (f"{arg},")

            log.append(msg + ")")

        e = pytools.MethodWrapper(action, pre=pre, post=post)

        e.foo(42, "bar")
        self.assertEqual(log, ["pre(foo,42,bar,)", "action(foo,42,bar,)", "post(foo,42,bar,)"])


    def test_logging_with_kwargs_only(self):

        log = []

        def pre(x: str, **kwargs):
            log.append(f"pre({x},A={kwargs['A']})")

        def action(x: str, **kwargs):
            log.append(f"action({x},B={kwargs['B']})")

        def post(x: str, **kwargs):
            log.append(f"post({x},C={kwargs['C']})")

        e = pytools.MethodWrapper(action, pre=pre, post=post)

        e.foo(C=1, A=42, B="bar")
        self.assertEqual(log, ["pre(foo,A=42)", "action(foo,B=bar)", "post(foo,C=1)"])


    def test_logging_with_vargs_and_kwargs(self):

        log = []

        def pre(x: str, *vargs, **kwargs):
            log.append(f"pre({x},{vargs[0]},A={kwargs['A']})")

        def action(x: str, *vargs, **kwargs):
            log.append(f"action({x},{vargs[1]},B={kwargs['B']})")

        def post(x: str, *vargs, **kwargs):
            log.append(f"post({x},{vargs[2]},C={kwargs['C']})")

        e = pytools.MethodWrapper(action, pre=pre, post=post)

        e.foo('a' , 'b', 'c', C=1, A=42, B="bar")
        self.assertEqual(log, ["pre(foo,a,A=42)", "action(foo,b,B=bar)", "post(foo,c,C=1)"])


    def test_unhandeled_method_name_pre(self):

        def pre():
            pass

        def action(x):
            pass

        e = pytools.MethodWrapper(action, pre=pre)

        with self.assertRaisesRegex(TypeError, r"foo.*pre"):
            e.foo()


    def test_unhandeled_method_name_action(self):

        def action():
            pass

        e = pytools.MethodWrapper(action)

        with self.assertRaisesRegex(TypeError, r"foo.*action"):
            e.foo()


    def test_unhandeled_method_name_post(self):

        def action(x):
            pass

        def post():
            pass

        e = pytools.MethodWrapper(action, post=post)

        with self.assertRaisesRegex(TypeError, r"foo.*post"):
            e.foo()


    def test_unhandeled_pre_vargs(self):

        def pre(x: str):
            pass

        def action(x: str, y: int):
            pass

        e = pytools.MethodWrapper(action, pre=pre)

        with self.assertRaisesRegex(TypeError, r"foo.*pre.*42"):
            e.foo(42)


    def test_unhandeled_action_vargs(self):

        def action(x: str):
            pass

        e = pytools.MethodWrapper(action)

        with self.assertRaisesRegex(TypeError, r"foo.*action.*42"):
            e.foo(42)


    def test_unhandeled_post_vargs(self):

        def action(x: str, y: int):
            pass

        def post(x: str):
            pass

        e = pytools.MethodWrapper(action, post=post)

        with self.assertRaisesRegex(TypeError, r"foo.*post.*42"):
            e.foo(42)


    def test_unhandeled_pre_kwargs(self):

        def pre(x: str, y: int):
            pass

        def action(x: str, y: int, **kwargs):
            pass

        e = pytools.MethodWrapper(action, pre=pre)

        with self.assertRaisesRegex(TypeError, r"foo.*pre.*bar.*100"):
            e.foo(42, bar=100)


    def test_unhandeled_action_kwargs(self):

        def action(x: str, y: int):
            pass

        e = pytools.MethodWrapper(action)

        with self.assertRaisesRegex(TypeError, r"foo.*action.*bar.*100"):
            e.foo(42, bar=100)


    def test_unhandeled_post_vargs(self):

        def action(x: str, y: int, **kwargs):
            pass

        def post(x: str, y: int):
            pass

        e = pytools.MethodWrapper(action, post=post)

        with self.assertRaisesRegex(TypeError, r"foo.*post.*bar.*100"):
            e.foo(42, bar=100)


if __name__ == "__main__":
    unittest.main()
