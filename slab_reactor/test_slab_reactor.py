import unittest

from slab_reactor import run

class TestSlabReactor(unittest.TestCase):

    def test_keffs(self):

        ref_keffs = [1.258249047, 1.007066854, 0.805372646,
                     1.329576914, 1.298737436, 0.681362819, 0.191909997]

        cases, keffs = run()

        tmpl = "Case {}, ref={} != {}"
        for i in range(len(cases)):
            self.assertAlmostEqual(ref_keffs[i], keffs[i],
                 msg=tmpl.format(cases[i], ref_keffs[i], keffs[i]))


if __name__ == '__main__':
    unittest.main()
