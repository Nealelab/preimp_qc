preimp_qc
=============

A Python package for performing GWAS QC.

Installation
------------

For now you can install preimp_qc using the command below. In the near future, it will be uploaded to pypi

.. code:: bash

   pip install -r requirements.txt
   python setup.py sdist
   pip install dist/preimp_qc-0.1.0.tar.gz

Usage
-----

.. code:: bash

   $ preimp_qc --dirname TestData/ --basename Basename --inputType plink --reference GRCh37
   in the example above, inside the directory TestData, there must be three PLINK file Basename.*{bed,bim,fam}

Arguments and Options
---------------------

+------------------------+--------------------------------------------+
| **Argument/ Option**   | **Description**                            |
+========================+============================================+
| ``--dirname``          | Path to where the data is                  |
+------------------------+--------------------------------------------+
| ``--basename``         | Data basename                              |
+------------------------+--------------------------------------------+
| ``--inputType``        | Input type, plink or vcf                   |
+------------------------+--------------------------------------------+
| ``--annotations``      | Annotations file to be used for annotating |
|                        | the VCF file (ONLY for VCF input)          |
+------------------------+--------------------------------------------+
| ``--reference``        | Reference genome build e.g. GRCh37, GRCh38 |
+------------------------+--------------------------------------------+
| ``--qc_round``         | The number of times QC has been performed  |
|                        | on the data. Use 1 if it's the first round |
+------------------------+--------------------------------------------+
| ``--pre_geno``         | include only SNPs with missing-rate < NUM  |
|                        | (before ID filter), important for post     |
|                        | merge of multiple platforms                |
+------------------------+--------------------------------------------+
| ``--mind``             | include only IDs with missing-rate < NUM   |
+------------------------+--------------------------------------------+
| ``--fhet_y``           | include only female IDs with fhet < NUM    |
+------------------------+--------------------------------------------+
| ``--fhet_x``           | include only male IDs with fhet > NUM      |
+------------------------+--------------------------------------------+
| ``--geno``             | include only SNPs with missing-rate < NUM  |
+------------------------+--------------------------------------------+
| ``--midi``             | include only SNPs with missing-rate        |
|                        | -difference ("case/control) < NUM          |
+------------------------+--------------------------------------------+
| ``--withpna``          | include monomorphic (invariant) SNPs       |
+------------------------+--------------------------------------------+
| ``--maf``              | include only SNPs with MAF >= NUM          |
+------------------------+--------------------------------------------+
| ``--hwe_th_con``       | HWE_controls < NUM                         |
+------------------------+--------------------------------------------+
| ``--hwe_th_cas``       | HWE_cases < NUM                            |
+------------------------+--------------------------------------------+