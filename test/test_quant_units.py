"""
Basic unit tests for the lower level quant.py functions, including
division-by-zero edge cases.
"""

import json

import pandas as pd

import quant


def _kallisto_df():
    return pd.DataFrame(
        {
            "target_id": ["A*01:01", "A*02:01", "B*07:02"],
            "length": [100.0, 120.0, 90.0],
            "eff_length": [90.0, 110.0, 80.0],
            "est_counts": [10.0, 0.0, 5.0],
            "tpm": [50.0, 0.0, 25.0],
        }
    )


class TestParseKallistoResults:
    def test_basic_sum(self):
        results = _kallisto_df()
        gene_indices = {"A1": {0}, "A2": {1}, "B1": {2}}

        lengths, counts, tpm = quant.parse_kallisto_results(results, gene_indices)

        assert counts["A1"] == 10.0
        assert counts["A2"] == 0.0
        assert counts["B1"] == 5.0
        assert lengths["A1"] == 100.0
        assert tpm["B1"] == 25.0

    def test_missing_gene_defaults_to_zero(self):
        """defaultdict should return 0.0 (not raise) for unseen genes."""
        results = _kallisto_df()
        lengths, counts, tpm = quant.parse_kallisto_results(results, {})

        assert counts["not_a_gene"] == 0.0
        assert lengths["not_a_gene"] == 0.0
        assert tpm["not_a_gene"] == 0.0

    def test_empty_indices_returns_zero(self):
        results = _kallisto_df()
        gene_indices = {"A1": set()}

        _, counts, tpm = quant.parse_kallisto_results(results, gene_indices)

        assert counts["A1"] == 0.0
        assert tpm["A1"] == 0.0


class TestComputeQuantMetrics:
    def test_basic_heterozygous(self):
        genes = {"A": ["A1", "A2"]}
        genotype = {"A1": "A*01:01", "A2": "A*02:01"}
        counts = {"A1": 10.0, "A2": 30.0}
        tpm = {"A1": 5.0, "A2": 15.0}

        gene_results, allele_results = quant.compute_quant_metrics(
            genes, genotype, counts, tpm
        )

        assert gene_results["A"]["count"] == 40
        assert allele_results["A"]["allele1_count"] == 10
        assert allele_results["A"]["allele2_count"] == 30
        # baf = min(10/40, 30/40) = 0.25
        assert allele_results["A"]["baf"] == 0.25

    def test_zero_counts_do_not_raise_division_by_zero(self):
        """When both alleles have zero counts, baf computation must not raise
        ZeroDivisionError and should default to a safe value."""
        genes = {"A": ["A1", "A2"]}
        genotype = {"A1": "A*01:01", "A2": "A*02:01"}
        counts = {"A1": 0.0, "A2": 0.0}
        tpm = {"A1": 0.0, "A2": 0.0}

        gene_results, allele_results = quant.compute_quant_metrics(
            genes, genotype, counts, tpm
        )

        assert allele_results["A"]["baf"] != allele_results["A"]["baf"]  # NaN
        assert gene_results["A"]["abundance"] == "0%"

    def test_total_hla_count_zero_across_genes(self):
        """Overall total_hla_count of zero (all alleles unexpressed) should not
        raise a ZeroDivisionError when computing gene abundance."""
        genes = {"A": ["A1", "A2"], "B": ["B1", "B2"]}
        genotype = {
            "A1": "A*01:01",
            "A2": "A*02:01",
            "B1": "B*07:02",
            "B2": "B*08:01",
        }
        counts = {"A1": 0.0, "A2": 0.0, "B1": 0.0, "B2": 0.0}
        tpm = {"A1": 0.0, "A2": 0.0, "B1": 0.0, "B2": 0.0}

        gene_results, allele_results = quant.compute_quant_metrics(
            genes, genotype, counts, tpm
        )

        assert gene_results["A"]["abundance"] == "0%"
        assert gene_results["B"]["abundance"] == "0%"
        assert allele_results["A"]["baf"] != allele_results["A"]["baf"]  # NaN
        assert allele_results["B"]["baf"] != allele_results["B"]["baf"]  # NaN


class TestComputeLohCorrectionDf:
    def _allele_results(self, a1=10.0, a2=30.0):
        return {
            "A": {
                "allele1": "A*01:01",
                "allele2": "A*02:01",
                "allele1_count": a1,
                "allele2_count": a2,
            }
        }

    def test_basic_no_loss(self):
        allele_results = self._allele_results(a1=20.0, a2=20.0)

        corrections_df = quant.compute_loh_correction_df(
            allele_results, ["A"], purity=1.0, ploidy=2.0
        )

        assert isinstance(corrections_df, pd.DataFrame)
        assert {"A_CN_1", "A_LOSS"}.issubset(corrections_df.columns)

    def test_zero_allele_counts_do_not_raise(self):
        """Both allele counts zero should not raise ZeroDivisionError when
        computing baf1/baf2."""
        allele_results = self._allele_results(a1=0.0, a2=0.0)

        # Should not raise.
        corrections_df = quant.compute_loh_correction_df(
            allele_results, ["A"], purity=1.0, ploidy=2.0
        )

        assert corrections_df is not None

    def test_zero_purity_does_not_raise(self):
        """A purity of zero previously caused a ZeroDivisionError when
        computing the copy-number correction factors."""
        allele_results = self._allele_results(a1=10.0, a2=30.0)

        # Should not raise.
        corrections_df = quant.compute_loh_correction_df(
            allele_results, ["A"], purity=0.0, ploidy=2.0
        )

        assert corrections_df is not None


class TestSaveResults(object):
    def test_save_allele_results_heterozygous(self, tmp_path):
        allele_results = {
            "A": {
                "allele1": "A*01:01",
                "allele2": "A*02:01",
                "allele1_count": 10,
                "allele2_count": 30,
                "allele1_tpm": 5,
                "allele2_tpm": 15,
                "baf": 0.25,
            }
        }
        tsv_path = tmp_path / "out.tsv"
        json_path = tmp_path / "out.json"

        quant.save_allele_results(allele_results, str(tsv_path), str(json_path))

        assert tsv_path.exists()
        with open(json_path) as fh:
            assert json.load(fh) == allele_results

    def test_save_gene_results(self, tmp_path):
        gene_results = {"A": {"count": 40, "tpm": 20, "abundance": "100.0%"}}
        tsv_path = tmp_path / "out.tsv"
        json_path = tmp_path / "out.json"

        quant.save_gene_results(gene_results, str(tsv_path), str(json_path))

        assert tsv_path.exists()
        with open(json_path) as fh:
            assert json.load(fh) == gene_results
