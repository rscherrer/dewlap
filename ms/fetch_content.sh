# Run this script to copy selected files from the analysis folders to the manuscript folder

# Fetch figures
cp ../analyses/04-machine-learning/plots/classif_svm_pca.png figures
cp ../analyses/07-ANOVA/figure_anova.png figures
cp ../maps/map.pdf figures

# Fetch supplementary figures
cp ../analyses/04-machine-learning/plots/classif_lda_pca.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_lda_pca_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_lda_refl.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_lda_refl_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_svm_pca_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_svm_refl.png suppfigures
cp ../analyses/04-machine-learning/plots/classif_svm_refl_pooled.png suppfigures
cp ../maps/detailed_map.pdf suppfigures
cp ../analyses/03-PCA/figure_brightness.png suppfigures
cp ../analyses/03-PCA/figure_brightness_pooled.png suppfigures
cp ../analyses/09-parallelism/figure_contrasts.png suppfigures
cp ../analyses/10-distances/figure_distances2.png suppfigures
cp ../analyses/10-distances/figure_distances.png suppfigures
cp ../analyses/06-pooled/figure_pooled.png suppfigures
cp ../analyses/02-reflectance/figure_reflectance.png suppfigures
cp ../analyses/02-reflectance/figure_reflectance_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_lda_pca.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_lda_pca_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_lda_refl.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_lda_refl_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_svm_pca.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_svm_pca_pooled.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_svm_refl.png suppfigures
cp ../analyses/04-machine-learning/plots/importance_svm_refl_pooled.png suppfigures

# Fetch tables
cp ../analyses/07-ANOVA/table_anova.tex tables

# Fetch supplementary tables
cp ../analyses/01-counts/table_counts.tex supptables
cp ../analyses/08-spatial-correlation/table_sites.tex supptables
cp ../analyses/03-PCA/table_expvar.tex supptables
cp ../analyses/03-PCA/table_brightness.tex supptables
cp ../analyses/05-assumptions/table_multinorm.tex supptables
cp ../analyses/05-assumptions/table_covariance.tex supptables
cp ../analyses/05-assumptions/table_normality.tex supptables
cp ../analyses/06-pooled/table_anova_pooled.tex supptables
cp ../analyses/04-machine-learning/table_classif_svm_pca.tex supptables
cp ../analyses/07-ANOVA/table_kw.tex supptables
cp ../analyses/08-spatial-correlation/table_spatial.tex supptables
