import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import PySimpleGUI as sg

sg.theme('DarkTeal9')  # Change the theme here

# Activate automatic conversion between R and Pandas dataframes
pandas2ri.activate()

# Load the edgeR package
robjects.r('library(edgeR)')

def create_layout():
    layout = [
        [sg.Text("Count Matrix File"), sg.Input(key="counts_file", size=(30, 1)), sg.FileBrowse()],
        [sg.Text("Clinical Data File"), sg.Input(key="clinical_file", size=(30, 1)), sg.FileBrowse()],
        [sg.Text("Min. Total Read Counts"), sg.Input(key="min_total_counts", default_text="10", size=(10, 1))],
        [sg.Text("Design Factors (comma separated)"), sg.Input(key="design_factors", size=(30, 1))],
        [sg.Text("Min. Log Fold Change"), sg.Input(key="min_lfc", default_text="1", size=(10, 1))],
        [sg.Text("Max. P-value"), sg.Input(key="max_pval", default_text="0.05", size=(10, 1))],
        [sg.Button("Run DEA"), sg.Button("Help"), sg.Button("Exit")],
        [sg.Output(size=(80, 20))],
        [sg.Text('DEA with edgeR', key='Application name', size=(None, 1), justification='left', font=("Alike", 11, "bold"))],
        [sg.Text('Thesis Project. Created by Mohit Panwar. Supervised by Julia Åkesson.', key='credits', size=(None, 1), justification='left', font=("Alike", 9))]
    ]
    return layout

def run_DEA(count_matrix, clinical_file, min_total_counts, design_factors, min_lfc, max_pval):
    """
    Run differential expression analysis using the edgeR package on a count matrix.

    Parameters:
    count_matrix (pandas.DataFrame): A DataFrame containing the count data, with genes as rows and samples as columns.
    clinical_file (str): The file path of the clinical data file.
    min_total_counts (int): The minimum total read counts.
    design_factors (list): A list of design factors.
    min_lfc (float): The minimum log fold change.
    max_pval (float): The maximum p-value.
    """
    try:
        # Convert the Pandas dataframe to an R matrix
        count_matrix_r = pandas2ri.py2rpy(count_matrix)

        # Load the clinical data
        clinical_data = pd.read_csv(clinical_file)
        clinical_data_r = pandas2ri.py2rpy(clinical_data)

        # Create the DGEList object
        dge = robjects.r['DGEList'](counts=count_matrix_r, genes=robjects.vectors.FactorVector(count_matrix.index))

        # Filter out low count genes
        keep = dge.sum(axis=1) >= min_total_counts
        dge_filtered = dge.rx(keep, True)

        # Define the experimental design
        design = clinical_data_r[design_factors]
        dge_design = robjects.r['estimateDisp'](dge_filtered, design)

        # Fit the model and perform the differential expression analysis
        fit = robjects.r['glmQLFit'](dge_design, design)
        qlf = robjects.r['glmQLFTest'](fit)
        res = robjects.r['summary'](qlf)

        # Convert the R dataframe to a Pandas dataframe
        DEGs = pandas2ri.rpy2py(res)
        DEGs = DEGs.loc[(DEGs['logFC'] > min_lfc) & (DEGs['PValue'] < max_pval), ]

        # Output results to the window
        print(DEGs)

        # Save the results to a file
        DEGs.to_csv(clinical_file.replace('.csv', '_DEA_results.csv'))

        sg.popup(f'Differential expression analysis complete. Results written to {clinical_file.replace(".csv", "_DEA_results.csv")}')

    except Exception as e:
        sg.popup(f'Error: {e}')

def show_help():
    help_text = """
    How to use:

    1. Click "Browse" next to "Count Matrix File" to provide a count matrix file (CSV format). The file should have rows representing genes and columns representing samples. The first row should contain the sample names.
    
    2. Click "Browse" next to "Clinical Data File" to provide a clinical data file (CSV format). The file should have rows representing samples and columns representing clinical variables. The first row should contain the variable names.
    
    3. Enter the minimum total read counts in the "Min. Total Read Counts" box to filter out genes with low expression. For example, enter "10" to keep genes with at least 10 total read counts across all samples.
    
    4. Enter design factors as comma-separated values in the "Design Factors" box. Design factors are the clinical variables you want to compare. For example, if you want to compare samples based on their condition and account for batch effects, input "condition,batch". You can use any clinical variables available in your clinical data file as design factors.
    
    5. Enter the minimum log fold change value in the "Min. Log Fold Change" box. Genes with a log fold change below this threshold will not be considered significant.
    
    6. Enter the maximum p-value in the "Max. P-value" box. Genes with a p-value above this threshold will not be considered significant.
    
    7. Click "Run PyDESeq2" to perform the differential expression analysis. The results will be displayed in a new window.

    File formats:

    Count matrix file (CSV):
        Rows represent genes, and columns represent samples. The first row should contain sample names.
        Example:
            Sample1,Sample2,Sample3
            Gene1,10,20,30
            Gene2,50,60,70

    Clinical data file (CSV):
        Rows represent samples, and columns represent clinical variables. The first row should contain variable names.
        Example:
            Sample,Condition,Batch
            Sample1,Control,1
            Sample2,Treated,1
            Sample3,Treated,2

    Example Scenarios:

    1. Suppose you are studying the effects of a drug on cancer cells. You have performed an RNA-seq experiment and obtained gene expression data for both treated and untreated cancer cells. Additionally, you have performed the experiment in two different labs, introducing a potential batch effect. Your count matrix file contains the gene expression data, while your clinical data file has information about the treatment conditions and the lab in which the experiment was performed. To perform differential expression analysis comparing treated and untreated samples while accounting for the batch effect, you would enter "condition,lab" in the "Design Factors" box. Enter the minimum total read counts, minimum log fold change, and maximum p-value as needed. Click "Run PyDESeq2" to view the results.

    2. Suppose you are investigating the impact of diet on gene expression in a mouse model. You have three groups of mice: one fed a high-fat diet, one fed a normal diet, and one fed a calorie-restricted diet. Furthermore, the mice are from two different genetic backgrounds. Your count matrix file contains gene expression data, and your clinical data file has information about the diet and genetic background of each mouse. To compare gene expression changes between the different diets while accounting for the genetic background, you would enter "diet,genetic background" in the "Design Factors" box. Enter the minimum total read counts, minimum log fold change, and maximum p-value as needed. Click "Run PyDESeq2" to view the results.

    3. In a study on the effects of aging on gene expression, you have samples from young and old individuals. The samples were collected at different time points, introducing a potential confounding factor. Your count matrix file contains gene expression data, while your clinical data file includes information about the age of the individuals and the time point of sample collection. To perform differential expression analysis comparing young and old individuals while accounting for the time point, you would enter "age,time_point" in the "Design Factors" box. Enter the minimum total read counts, minimum log fold change, and maximum p-value as needed. Click "Run PyDESeq2" to view the results.
    """
    sg.popup_scrolled("Help", help_text, size=(80, 25))

def main():
    # Create the PySimpleGUI window
    window = sg.Window('Differential Expression Analysis', create_layout())

    while True:
        event, values = window.read()

        if event == sg.WINDOW_CLOSED or event == 'Exit':
            break
        
        if event == 'Help':
            show_help()    

        if event == 'Run DEA':
            counts_file = values['counts_file']
            clinical_file = values['clinical_file']
            min_total_counts = int(values['min_total_counts'])
            design_factors = values['design_factors'].split(',')
            min_lfc = float(values['min_lfc'])
            max_pval = float(values['max_pval'])

            try:
                # Load the count data from the file
                count_data = pd.read_csv(counts_file, index_col=0)

                # Run differential expression analysis using the modified run_DEA function
                run_DEA(count_data, clinical_file, min_total_counts, design_factors, min_lfc, max_pval)

            except Exception as e:
                sg.popup(f'Error: {e}')

            window.close()

if __name__ == '__main__':
    main()
