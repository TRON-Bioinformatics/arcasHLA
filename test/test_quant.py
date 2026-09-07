"""
Test read quantification for expected outputs and behavior.
"""

import json
import os.path
import quant


def _load_quant_res_jsons(output_dir):
    """
    Convenience function to load both gene and allele quant results.
    """
    return {
        "genes": json.load(open(os.path.join(output_dir, "test.quant.genes.json"))),
        "alleles": json.load(open(os.path.join(output_dir, "test.quant.alleles.json"))),
    }


class TestMain:
    def test_basic(
        self, extract_reads, customize_reference, tmp_path, expected_output_dir
    ):
        # Should run without error.
        quant.do_quantification(
            file=extract_reads,
            sample="test",
            ref=customize_reference,
            outdir=str(tmp_path),
        )

        actual_results = _load_quant_res_jsons(tmp_path)
        expected_results = _load_quant_res_jsons(expected_output_dir)

        assert actual_results["genes"] == expected_results["genes"]
        assert actual_results["alleles"] == expected_results["alleles"]
