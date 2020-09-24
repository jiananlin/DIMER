<!---Usage of python code for RNA methylation differential analysis.
% Edited by DNA - September 24, 2020--->
#  **Dimer**
Dimer (not to be confused with the [chemical molecule](https://en.wikipedia.org/wiki/Dimer_(chemistry))) stands for <b>Di</b>fferential <b>Me</b>thylation of <b>R</b>NA transcripts. Dimer is a Python-based package for testing differential methylation in MeRIP-Seq data.

## Dependencies
Dimer is written in Python and runs on Python (2.7 or lower). It requires the following bioinformatics and python packages installed to run successfully:
1. bedtools - http://bedtools.readthedocs.io/en/latest/
2. pybedtools - https://daler.github.io/pybedtools/
3. pandas - https://pandas.pydata.org/index.html

## Installation
The package can be downloaded from Github and extracted directly into any directory:
```
download code from github
```
After extracting the code, navigate to the `path_to_package/script/ folder and change permissions for the file` split_exon_ff to be executable:
```
chmod +x path_to_package/script/dimer
```
You may as well consider adding the `path_to_package/script/` to your `PATH` variable.
## Usage
Once the path to DIMER is added to your ```PATH``` variable, you can simply run DIMER using
```
dimer --options arguments
```
### Options
Dimer uses BAM files and a gene annotation file (GTF/GTF2) as input files. Testing differential methylation requires four sets of BAM files - IP and input files for treatment and control. To minimize the input variables, dimer takes all the BAM file paths through a single tab-delimited text file as argument to the ```--info``` variable.
* -h/--help&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;show this help message and exit.
*  -g/--gtf &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;input gtf file.
*  -i /--info &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;input PATH file. Please refer to the readme file for the format of this PATH file.
* -o /--outdir &ensp;&ensp;&ensp;&ensp;&ensp;output directory to save all ouput files. Default is the current folder.
* -w/--winsize&ensp;&ensp;&ensp;&ensp;&ensp;window size used to create bins. Default is 50 base pairs.
* -p/--pathin&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;The directory of the bam files. If you have the valid directory information already in the file of -i option, you don't need to
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;use this option. If you only specify the filename in the file of -i option, you need to provide this option for DIMER to find &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;the input files.

A example file list for three versus three sample comparison is provided below. A template file is provided in ```txt/file_name_list_template.txt```.
```
Group	Condition Replicate   File
Treatment	IP	1
Treatment	IP	2
Treatment	IP	3
Treatment	Input	1
Treatment	Input	2
Treatment	Input	3
Control	IP	1
Control	IP	2
Control	IP	3
Control	Input	1
Control	Input	2
Control	Input	3
```
The only user-specific parameter ```-win-size``` determines the window size to be used to assign read counts for all the genes specified in the GTF file provided. The testing procedure in the package does not support different window sizes for different BAM files. Finally, the user can specify the output file name using ```--outfile``` argument.

### Command
Dimer can be run using the following command
```
dimer --gtf path_to_gtf_file --info file_name_list --winsize window_size --outdir path_to_output_count_file
```


### Output format
Dimer mainly performs two tasks
1. Calculate window-based counts from the BAM files.
2. Test for differential methylation at the gene level and report the test statistic and the p-value.

Once both tasks are completed, dimer produces three sets of files
1. `win_winsize_count_gene_cqtest.txt` : A tab-separated text file consisting of three columns - gene name,  the test statistic and p-value.
2. `win_100_count.txt` : A tab-separated text file with the counts for each window for all the IP and input samples under both the conditions.
3. `win_X_meth_lv_Y_Z.bedgraph` : Bedgraph files constructed using the specified window size (**X**) for all conditions (**Y** = Treatment, Control) and all samples (**Z**) within each condition.


## References
Ayyala D. N., Lin, J., Ouyang, Z.;  **Differential RNA Methylation using Multivariate Statistical Methods**, *submitted*.
