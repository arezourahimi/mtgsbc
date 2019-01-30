This repository contains R implementation of the algorithms proposed in "A multitask multiple kernel learning formulation for discriminating early- and late-stage cancers" (manuscript under review).

* run_mtmkl_step1_coclustering.R : shows how to produce similarity matrices based on the TCGA cohorts and Hallmark pathways
* run_mtmkl_step2_collect_coclustering_data.R : produces the similarity matrices by combining the results obtained in step 1
* run_mtmkl_step3_using_aggregated_similarity_matrices.R : shows how to replicate multitask experiments on the TCGA cohorts using the similarity matrices produced in step 2.
* run_mtmkl_step4_collect_classification_data.R	: collects the classification results

MTGSBC methods
------------
* classification_helper.R => helper functions
* solve_classification_models_cplex.R => support vector machine classification solver, and the cutting plane model using CPLEX optimization software
* solve_classification_models_mosek.R => support vector machine classification solver, and the cutting plane model using Mosek optimization software
* group_lasso_multitask_multiple_kernel_classification_train.R => training procedure for multitask group Lasso MKL
* group_lasso_multitask_multiple_kernel_classification_test.R => test procedure for multitask group Lasso MKL
* mtmkl_coclustering_algorithm.R => the heuristic algorithm used for the diversificaion phase in step 1

If you use any of the algorithms implemented in this repository, please cite the following paper (under review):

Arezou Rahimi and Mehmet Gonen. "A multitask multiple kernel learning formulation for discriminating early- and late-stage cancers"
