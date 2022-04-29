import pandas as pd
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
import sys
import statsmodels.formula.api as smf
import os.path

if __name__ == "__main__":

    phe_code_file = str(sys.argv[1])
    phe_code = phe_code_file.replace('Phecode_','').replace('.csv.gz','')
    
    pgs_index = int(sys.argv[2])-1

    # qc file from UKB
    ukb_qc_file = "/home/yx322/rds/rds-jmmh2-post_qc_data/uk_biobank/reference_files/genetic/reference_files/full_release/QC_documents/sampleQC_fromUKB_withHeaders.txt"
    # read ukb samples qc data to get array and PCs info
    df_ukb_qc = pd.read_csv(ukb_qc_file, skiprows=[i for i in range(0, 31)], sep=' ')
    
    
    #### select white british only ####
    df_ukb_qc = df_ukb_qc.loc[df_ukb_qc['in.white.British.ancestry.subset'] == 1]
    df_ukb_qc = df_ukb_qc[['#UKB_ID1', 'genotyping.array', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]

    
    pheno_file = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/PGSCatalog/PheWAS/_phenotyped/" + phe_code_file
    df_pheno = pd.read_csv(pheno_file,compression='gzip')
    df_pheno = df_pheno[['eid','genid','sex','PHECODE_AgeAsTimescale','PHECODE_AgeAsTimescale_Years']]
    df_pheno = df_pheno.rename(columns={'eid':'idno','PHECODE_AgeAsTimescale':'PHENOTYPE','PHECODE_AgeAsTimescale_Years':'CENSOR_AGE'})
    
    # Incomporate pcs and arrary info
    df_pheno_all = pd.merge(df_pheno, df_ukb_qc, left_on='genid', right_on='#UKB_ID1')
    
    omics_pgs_files = ['/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Metabolon_full/UKB_Metabolon.sscore.gz',\
        "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Olink_full/UKB_Olink.sscore.gz", \
        '/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Somalogic_full/UKB_Somalogic.sscore.gz', \
        '/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_Nightingale_5e8/UKB_Nightingale.sscore.gz']

    platforms = ['Metabolon', 'Olink', 'Somalogic', 'Nightingale']

    # add in gene expression PGS by chr
    for chr in range(1, 23):
        GE_pgs_file = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_GE_5e-8/chr" + str(chr) + "/UKB_GE.sscore.gz"
        omics_pgs_files.append(GE_pgs_file)
        platform_name = "GE_chr" + str(chr)
        platforms.append(platform_name)
    
    i = pgs_index
    platform = platforms[i]
    SOMA_PGS_file = omics_pgs_files[i]
    df_soma_pgs = pd.read_csv(SOMA_PGS_file, sep='\t', compression='gzip')
    score_cols = list(df_soma_pgs.columns[1:])
    
    write_file = '/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/UKB_phecode_association/raw_assocs/' + platform + '_PGS_UKB_' + phe_code + '_assoc_full_EU.txt'
    
    df_save = pd.DataFrame(columns=['Trait', 'HR', 'HR_low', 'HR_high', 'pvalue'])

    for col_name in score_cols:

        df_one_pgs = df_soma_pgs[['IID', col_name]]

        df_pheno_test = pd.merge(df_one_pgs, df_pheno_all, left_on='IID', right_on='idno')
        df_pheno_test = df_pheno_test[[col_name, 'sex', 'PHENOTYPE', 'CENSOR_AGE', 'genotyping.array', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]

        df_pheno_test = pd.get_dummies(df_pheno_test, drop_first=True)
        df_pheno_test[[col_name, 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9','PC10']] = StandardScaler().fit_transform(df_pheno_test[[col_name, 'PC1', 'PC2', 'PC3', 'PC4', 'PC5','PC6', 'PC7', 'PC8', 'PC9', 'PC10']])

        # adjust pgs for PCs
        #reg_trait = col_name + ' ~ PC1+ PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10'
        #adj_result = smf.ols(reg_trait, data=df_pheno_test, missing='drop').fit()
        #df_pheno_test[col_name] = adj_result.resid
        #df_pheno_test[[col_name]] = StandardScaler().fit_transform(df_pheno_test[[col_name]])

        df_pheno_test['temp_col'] = df_pheno_test[col_name]
        reg_trait =  'temp_col ~ PC1+ PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10'
        adj_result = smf.ols(reg_trait, data=df_pheno_test, missing='drop').fit()
        df_pheno_test[col_name] = adj_result.resid
        del df_pheno_test['temp_col']
        # remove nan rows
        index_nan = list(df_pheno_test[df_pheno_test.isnull().any(axis=1)].index)
        df_pheno_test.drop(index_nan,inplace=True)
        
        df_pheno_test[[col_name]] = StandardScaler().fit_transform(df_pheno_test[[col_name]])
        
        if 'sex_Male' in df_pheno_test.columns:
            try:
                cph1 = CoxPHFitter()
                cph1.fit(df_pheno_test, duration_col='CENSOR_AGE', event_col='PHENOTYPE', strata=['sex_Male'])
            except:
                cph1 = CoxPHFitter()
                cph1.fit(df_pheno_test, duration_col='CENSOR_AGE', event_col='PHENOTYPE', strata=['sex_Male'],step_size=0.1)
        else:
            try:
                cph1 = CoxPHFitter()
                cph1.fit(df_pheno_test, duration_col='CENSOR_AGE', event_col='PHENOTYPE')
            except:
                cph1 = CoxPHFitter()
                cph1.fit(df_pheno_test, duration_col='CENSOR_AGE', event_col='PHENOTYPE',step_size=0.1)

        #cph1 = CoxPHFitter()
        # cph1.fit(df_testing_analysis_dummies, duration_col='duration', event_col='chd_case',strata=['sex_Male'])
        #cph1.fit(df_pheno_test, duration_col='CENSOR_AGE', event_col='PHENOTYPE', strata=['sex_Male'])
        results = [col_name, cph1.hazard_ratios_[col_name], cph1.summary.loc[col_name, 'exp(coef) lower 95%'], cph1.summary.loc[col_name, 'exp(coef) upper 95%'], cph1.summary.loc[col_name, 'p']]

        #print(results)

        df_save = df_save.append({'Trait': results[0], 'HR': results[1], 'HR_low': results[2], 'HR_high': results[3], 'pvalue': results[4]}, ignore_index=True)

    df_save = df_save.sort_values('pvalue')
    df_save.to_csv(write_file, sep='\t', index=False)

