# <div align="center">GSynoX</div>
## <div align="center">A lightweight gene ID translator that deals with synonyms and deprecated aliases</div>

GSynoX is a lightweight software tool designed to simplify a common task in bioinformatics: translating gene IDs between different formats such as HGNC gene symbols, Ensembl IDs, and Entrez IDs. To achieve this efficiently, GSynoX handles the various aliases and deprecated IDs associated with each gene symbol. Leveraging data from the HGNC database, GSynoX enables consistent use of gene symbols across datasets, even when different aliases are employed for the same genes. This is particularly useful for harmonizing data during integrative analyses. Unlike tools like BioMart that perform on-the-fly translations via internet connectivity, GSynoX operates entirely offline. It maintains a local copy of gene ID mappings, ensuring portability and optimal performance in environments such as high-performance computing (HPC) clusters, where internet access may be restricted for security reasons.

## Hello world translation

This tutorial introduces GSynoX and demonstrates its capabilities for various tasks, including handling gene synonyms, prioritizing preferred IDs, and integrating custom datasets. Follow along to learn how to use GSynoX effectively in your bioinformatics workflows. Let’s get started!

To begin, we need to import the necessary classes and initialize the GSynoX class, which serves as the primary interface for interacting with the tool.


```python
from gsynox import GSynoX
from gsynox.builder import Builder
import pandas as pd

g = GSynoX()
```

    Using internal database /home/jose/anaconda3/envs/gxref/lib/python3.10/site-packages/gsynox/resources/master.pkl


To evaluate the current content of GSynoX master database, we can use the function `get_info`.


```python
g.get_info()
```




    {'version': 0.1,
     'date': 'Wed Nov 13 17:18:22 2024',
     'dbs': ['symbol', 'ensembl_gene', 'entrez', 'hgnc', 'uniprot'],
     'log': ['Added uniprot on Wed Nov 13 17:19:55 2024']}



Once created, GSynoX instance can easily translate a list of provided ids by the user, in this case from gene SYMBOL to ENSEMBL_GENE formats.


```python
some_genes = ["FOXP2", "MYC", "CEBPA"]
g.symbol_to_ensembl_gene(some_genes)
```




    ['ENSG00000128573', 'ENSG00000136997', 'ENSG00000245848']



## Working with synonyms

Gene symbols often have synonyms or deprecated aliases, which can vary widely across datasets and introduce significant variability. GSynoX addresses this challenge by enabling you to retrieve and manage these synonyms, ensuring consistent and standardized naming throughout your analyses



```python
gene = "MYC"
g.synonyms(gene)
```




    ['MYCC', 'bHLHe39', 'c-Myc']



When GSynoX translates an external ID, it defaults to returning the official gene symbol. However, the package can also retrieve all registered synonyms, providing additional context and flexibility for downstream analyses. To enable this feature, simply set the all_synonyms parameter to True.


```python
ens = g.symbol_to_ensembl_gene(gene)
ens
```




    'ENSG00000136997'




```python
g.ensembl_gene_to_symbol(ens)
```




    'MYC'




```python
g.ensembl_gene_to_symbol(ens, all_synonyms=True)
```




    ['MYC', 'MYCC', 'bHLHe39', 'c-Myc']



This functionality also works when the user provides more than ID for translation, providing as a results a list of lists.


```python
some_ens = g.symbol_to_ensembl_gene(some_genes)
g.ensembl_gene_to_symbol(some_ens, all_synonyms=True)
```




    [['FOXP2', 'CAGH44', 'SPCH1', 'TNRC10'],
     ['MYC', 'MYCC', 'bHLHe39', 'c-Myc'],
     ['CEBPA', 'C/EBP-alpha', 'CEBP']]



## Using a preferred list of IDs

The recommended best practice for harmonizing gene SYMBOLs across two datasets is to convert both to their official ids. However, in cases where one dataset cannot be easily converted (e.g., due to a complex data model), you can align the SYMBOLs by translating one dataset's IDs to match the specific aliases used in the other. This can be achieved by leveraging the `preferred_ids` parameter to prioritize the desired aliases during the translation process.




```python
g.synonyms(some_genes)
```




    [['CAGH44', 'SPCH1', 'TNRC10'],
     ['MYCC', 'bHLHe39', 'c-Myc'],
     ['C/EBP-alpha', 'CEBP']]




```python
g.ensembl_gene_to_symbol(some_ens)
```




    ['FOXP2', 'MYC', 'CEBPA']




```python
g.ensembl_gene_to_symbol(some_ens, preferred_ids=["CAGH44", "MYCC", "C/EBP-alpha"])
```




    ['CAGH44', 'MYCC', 'C/EBP-alpha']



## Dealing with missing IDs

GSynoX can handle cases where some of the provided IDs are missing or invalid, making it robust for messy datasets, By default, missing IDs return None.



```python
g.ensembl_gene_to_symbol(["FAKE1", "FAKE2"])
```




    [None, None]



For more flexibility, you can specify a custom default value for missing IDs.


```python
g.default_null_id = "unknown"
g.ensembl_gene_to_symbol(["FAKE1", "FAKE2"])
```




    ['unknown', 'unknown']



## Shorcuts to most popular 

GSynoX supports direct translation between the most popular ID types, such as ENSEMBL_GENE, ENTREZ or HGNC IDs.



```python
g.symbol_to_ensembl_gene(some_genes)
```




    ['ENSG00000128573', 'ENSG00000136997', 'ENSG00000245848']




```python
g.symbol_to_entrez(some_genes)
```




    ['93986', '4609', '1050']




```python
g.symbol_to_hgnc(some_genes)
```




    ['HGNC:13875', 'HGNC:7553', 'HGNC:1833']



Also, at the protein level, such as UNIPROT


```python
g.symbol_to_uniprot(some_genes)
```




    [['A0A994J3Y0',
      'Q8N6B5',
      'A0A0U1RQY3',
      'F8WDL6',
      'Q75MZ5',
      'A0A0U1RQR8',
      'X5D2H2',
      'Q0PRL4',
      'O15409',
      'A0A994J6W1',
      'A8MUV4',
      'A0A0U1RQE0',
      'C9JQP8',
      'A0A0U1RQM2',
      'A8MTU2'],
     ['Q14899', 'A0A494C1T8', 'H0YBG3', 'A0A0B4J1R1', 'Q16591', 'P01106'],
     ['P49715']]



Although GSynoX is mainly oriented to translate from/to gene symbol, it allows users to translate between two external databases such as ENSEMBL_GENE and ENTREZ ids.


```python
some_ens
```




    ['ENSG00000128573', 'ENSG00000136997', 'ENSG00000245848']




```python
g.cross_id(some_ens, "ensembl_gene", "entrez", select_one=True)
```




    ['93986', '4609', '1050']



## Adding a new type of identifier to GSynoX

GSynoX can be easily extended with custom datasets, allowing you to incorporate additional ID mappings specific to your research. The Builder class lets the user to manage and add custom databases.


```python
# Initilizing db Builder
b = Builder()
b.master = g.master
```

As an example, we can integrate a toy database containing specific IDs for three existing gene symbols. The input file for this database should consist of two columns: the first column lists the specific gene symbols, while the second contains their corresponding IDs from the external database. If a gene symbol maps to multiple IDs, these IDs should be listed in the second column, separated by colons. This format ensures that all possible translations are properly captured and usable.


```python
toy_db_file = "toy_db.tsv"
toy_db = pd.read_csv(toy_db_file, sep="\t", dtype="str")
toy_db
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>symbol</th>
      <th>external</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>MYC</td>
      <td>a,b,c</td>
    </tr>
    <tr>
      <th>1</th>
      <td>FOXP2</td>
      <td>d,e</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CEBPA</td>
      <td>f,g,h,i</td>
    </tr>
  </tbody>
</table>
</div>



To add the database we just have to use the `add_db` function of the `Builder` class as follows:


```python
b.add_db(toy_db_file, "toy")
g.master = b.master
```

Then, to use the newly added database we will employ the generic functions `symbol_to_id` and `id_to_symbol`.


```python
g.symbol_to_id(some_genes, "toy", select_one=False)
```




    [['d', 'e'], ['a', 'b', 'c'], ['f', 'g', 'h', 'i']]




```python
g.id_to_symbol(["a", "d", "g"], "toy")
```




    ['MYC', 'FOXP2', 'CEBPA']



Finally, thanks to the internal structure of GSynoX, the newly added database can be also use for cross tranlation to other external database.


```python
ens = g.symbol_to_ensembl_gene("FOXP2")
g.cross_id(ens, "ensembl_gene", "toy")
```




    ['d', 'e']


