"""
Direct code testing of the reference customization script.
"""

import customize


class TestMain:
    def test_basic(self):
        genotype_result_path = "test/expected_output/test.genotype.json"
        custom_ref_out_path = "test/output/custom_reference"

        main_args = [
            "--genotype",
            genotype_result_path,
            "-o",
            custom_ref_out_path,
        ]

        # Should run without error.
        customize.main(main_args)
