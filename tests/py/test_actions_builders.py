import unittest
import runko


quote = runko.actions.quote
sequence = runko.actions.sequence
parallel = runko.actions.parallel
sequence_builder = runko.sequence_builder
parallel_builder = runko.parallel_builder


class actions_sequence_builder(unittest.TestCase):
    def test_empty(self):
        builder = sequence_builder()
        correct = (sequence,)
        self.assertEqual(builder.get(), correct)


    def test_one_long(self):
        builder = sequence_builder()
        builder.then(4, 2)
        correct = (sequence, (quote, (4, 2)))
        self.assertEqual(builder.get(), correct)


    def test_two_long(self):
        builder = sequence_builder()
        builder.then(4, 2)
        builder.then(0, 1, 2)
        correct = (sequence, (quote, (4, 2)), (quote, (0, 1, 2)))
        self.assertEqual(builder.get(), correct)


    def test_two_long_chain(self):
        builder = sequence_builder()
        builder.then(4, 2).then(0, 1, 2)
        correct = (sequence, (quote, (4, 2)), (quote, (0, 1, 2)))
        self.assertEqual(builder.get(), correct)


    def test_empty_is_benign(self):
        builder = sequence_builder()

        for _ in range(10):
            builder.then()

        correct = (sequence,)
        self.assertEqual(builder.get(), correct)




class actions_parallel_builder(unittest.TestCase):
    def test_foo_empty(self):
        builder = parallel_builder()
        correct = (parallel,)
        self.assertEqual(builder.get(), correct)


    def test_one_wide(self):
        builder = parallel_builder()
        builder.add(4, 2)
        correct = (parallel, (quote, (4, 2)))
        self.assertEqual(builder.get(), correct)


    def test_two_wide(self):
        builder = parallel_builder()
        builder.add(4, 2)
        builder.add(0, 1, 2)
        correct = (parallel, (quote, (4, 2)), (quote, (0, 1, 2)))
        self.assertEqual(builder.get(), correct)


    def test_two_wide_chain(self):
        builder = parallel_builder()
        builder.add(4, 2).add(0, 1, 2)
        correct = (parallel, (quote, (4, 2)), (quote, (0, 1, 2)))
        self.assertEqual(builder.get(), correct)


    def test_empty_is_benign(self):
        builder = parallel_builder()

        for _ in range(10):
            builder.add()

        correct = (parallel,)
        self.assertEqual(builder.get(), correct)




class actions_nested_builders(unittest.TestCase):
    def test_sequence_in_sequence(self):
        inner = sequence_builder().then(43)
        outer = sequence_builder().then(42).then(inner)

        first = (quote, (42,))
        second = (quote, (sequence, (quote, (43,))))
        correct = (sequence, first, second)
        print(outer.get())
        print(correct)
        self.assertEqual(outer.get(), correct)


    def test_lazy_sequence_in_sequence(self):
        a, b = sequence_builder(), sequence_builder()
        builder = sequence_builder()
        builder.then(a).then(b)
        a.then(1)
        b.then(2)

        first = (quote, (sequence, (quote, (1,))))
        second = (quote, (sequence, (quote, (2,))))
        correct = (sequence, first, second)
        self.assertEqual(builder.get(), correct)


    def test_sequence_in_parallel(self):
        builder = parallel_builder()
        builder.add(sequence_builder().then(1)).add(1)

        first = (quote, (sequence, (quote, (1,))))
        second = (quote, (1,))
        correct = (parallel, first, second)
        self.assertEqual(builder.get(), correct)


    def test_parallel_in_parallel(self):
        builder = parallel_builder()
        builder.add(parallel_builder().add(1)).add(1)

        first = (quote, (parallel, (quote, (1,))))
        second = (quote, (1,))
        correct = (parallel, first, second)
        self.assertEqual(builder.get(), correct)




if __name__ == "__main__":
    unittest.main()
