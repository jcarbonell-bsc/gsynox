#!/usr/bin/env python
# coding: utf-8

# # <div align="center">GSynoX</div>
# ## <div align="center">A lightweight gene ID translator that deals with synonyms and deprecated aliases</div>
# 
# GSynoX is a lightweight software tool designed to simplify a common task in bioinformatics: translating gene IDs between different formats such as HGNC gene symbols, Ensembl IDs, and Entrez IDs. To achieve this efficiently, GSynoX handles the various aliases and deprecated IDs associated with each gene symbol. Leveraging data from the HGNC database, GSynoX enables consistent use of gene symbols across datasets, even when different aliases are employed for the same genes. This is particularly useful for harmonizing data during integrative analyses. Unlike tools like BioMart that perform on-the-fly translations via internet connectivity, GSynoX operates entirely offline. It maintains a local copy of gene ID mappings, ensuring portability and optimal performance in environments such as high-performance computing (HPC) clusters, where internet access may be restricted for security reasons.

# ## Hello world translation
# 
# This tutorial introduces GSynoX and demonstrates its capabilities for various tasks, including handling gene synonyms, prioritizing preferred IDs, and integrating custom datasets. Follow along to learn how to use GSynoX effectively in your bioinformatics workflows. Let’s get started!
# 
# To begin, we need to import the necessary classes and initialize the GSynoX class, which serves as the primary interface for interacting with the tool.

# In[70]:


from gsynox import GSynoX
from gsynox.builder import Builder
import pandas as pd

g = GSynoX()


# To evaluate the current content of GSynoX master database, we can use the function `get_info`.

# In[71]:


g.get_info()


# Once created, GSynoX instance can easily translate a list of provided ids by the user, in this case from gene SYMBOL to ENSEMBL_GENE formats.

# In[72]:


some_genes = ["FOXP2", "MYC", "CEBPA"]
g.symbol_to_ensembl_gene(some_genes)


# ## Working with synonyms
# 
# Gene symbols often have synonyms or deprecated aliases, which can vary widely across datasets and introduce significant variability. GSynoX addresses this challenge by enabling you to retrieve and manage these synonyms, ensuring consistent and standardized naming throughout your analyses
# 

# In[73]:


gene = "MYC"
g.synonyms(gene)


# When GSynoX translates an external ID, it defaults to returning the official gene symbol. However, the package can also retrieve all registered synonyms, providing additional context and flexibility for downstream analyses. To enable this feature, simply set the all_synonyms parameter to True.

# In[74]:


ens = g.symbol_to_ensembl_gene(gene)
ens


# In[75]:


g.ensembl_gene_to_symbol(ens)


# In[76]:


g.ensembl_gene_to_symbol(ens, all_synonyms=True)


# This functionality also works when the user provides more than ID for translation, providing as a results a list of lists.

# In[77]:


some_ens = g.symbol_to_ensembl_gene(some_genes)
g.ensembl_gene_to_symbol(some_ens, all_synonyms=True)


# ## Using a preferred list of IDs
# 
# The recommended best practice for harmonizing gene SYMBOLs across two datasets is to convert both to their official ids. However, in cases where one dataset cannot be easily converted (e.g., due to a complex data model), you can align the SYMBOLs by translating one dataset's IDs to match the specific aliases used in the other. This can be achieved by leveraging the `preferred_ids` parameter to prioritize the desired aliases during the translation process.
# 
# 

# In[78]:


g.synonyms(some_genes)


# In[79]:


g.ensembl_gene_to_symbol(some_ens)


# In[80]:


g.ensembl_gene_to_symbol(some_ens, preferred_ids=["CAGH44", "MYCC", "C/EBP-alpha"])


# ## Dealing with missing IDs
# 
# GSynoX can handle cases where some of the provided IDs are missing or invalid, making it robust for messy datasets, By default, missing IDs return None.
# 

# In[81]:


g.ensembl_gene_to_symbol(["FAKE1", "FAKE2"])


# For more flexibility, you can specify a custom default value for missing IDs.

# In[82]:


g.default_null_id = "unknown"
g.ensembl_gene_to_symbol(["FAKE1", "FAKE2"])


# ## Shorcuts to most popular 
# 
# GSynoX supports direct translation between the most popular ID types, such as ENSEMBL_GENE, ENTREZ or HGNC IDs.
# 

# In[83]:


g.symbol_to_ensembl_gene(some_genes)


# In[84]:


g.symbol_to_entrez(some_genes)


# In[85]:


g.symbol_to_hgnc(some_genes)


# Also, at the protein level, such as UNIPROT

# In[86]:


g.symbol_to_uniprot(some_genes)


# Although GSynoX is mainly oriented to translate from/to gene symbol, it allows users to translate between two external databases such as ENSEMBL_GENE and ENTREZ ids.

# In[87]:


some_ens


# In[88]:


g.cross_id(some_ens, "ensembl_gene", "entrez", select_one=True)


# ## Adding a new type of identifier to GSynoX
# 
# GSynoX can be easily extended with custom datasets, allowing you to incorporate additional ID mappings specific to your research. The Builder class lets the user to manage and add custom databases.

# In[89]:


# Initilizing db Builder
b = Builder()
b.master = g.master


# As an example, we can integrate a toy database containing specific IDs for three existing gene symbols. The input file for this database should consist of two columns: the first column lists the specific gene symbols, while the second contains their corresponding IDs from the external database. If a gene symbol maps to multiple IDs, these IDs should be listed in the second column, separated by colons. This format ensures that all possible translations are properly captured and usable.

# In[90]:


toy_db_file = "toy_db.tsv"
toy_db = pd.read_csv(toy_db_file, sep="\t", dtype="str")
toy_db


# To add the database we just have to use the `add_db` function of the `Builder` class as follows:

# In[91]:


b.add_db(toy_db_file, "toy")
g.master = b.master


# Then, to use the newly added database we will employ the generic functions `symbol_to_id` and `id_to_symbol`.

# In[92]:


g.symbol_to_id(some_genes, "toy", select_one=False)


# In[93]:


g.id_to_symbol(["a", "d", "g"], "toy")


# Finally, thanks to the internal structure of GSynoX, the newly added database can be also use for cross tranlation to other external database.

# In[94]:


ens = g.symbol_to_ensembl_gene("FOXP2")
g.cross_id(ens, "ensembl_gene", "toy")


