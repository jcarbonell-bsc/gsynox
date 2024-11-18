#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import numpy as np
# import pkg_resources
import os

class GSynoX():
    """
    A class for managing and translating gene identifiers across various databases.
    
    Attributes:
    -----------
    master : dict
        A dictionary containing the database mappings.
    default_null_id : Any
        The default return value for missing translations.
    """

    master = {}
    default_null_id = None
    
    def __init__(self, master_file=None, init_empty=False):
        """
        Initialize the GSynox object.

        Parameters:
        -----------
        master_file : str, optional
            Path to the master database file (default is None).
        init_empty : bool, optional
            If True, initializes an empty master database (default is False).
        """
        
        if init_empty is False:
            if master_file is None:
                master_file = os.path.abspath(os.path.dirname(__file__) + '/resources/master.pkl')
                print("Using internal database " + master_file)
            
            with open(master_file, "rb") as g:
                self.master = pickle.load(g)

    
    def __translate(self, ids:list, db:str, domain:str, select_one=False, preferred_ids=[]) -> list:
        """
        Translate a list of IDs to/from a given database.
          
        Parameters:
        -----------
        ids : list
            List of IDs to be translated.
        db : str
            External database name involved in the translation (e.g., 'ensembl_gene').
        domain : str
            Defines the translation direction (e.g., 'symbol_to_id', 'id_to_symbol').
        select_one : bool, optional
            If True, selects the first translation result (default is False).
        preferred_ids : list, optional
            A list of IDs to prioritize during translation (default is empty).
          
        Returns:
        --------
        list
            A list of translated IDs.
        """
                
        # prepare ids
        if isinstance(ids, str):
            ids = [ids]
        
        # translate
        if len(preferred_ids) == 0:
            tx = list(map(lambda id: self.__translate_id(id, db, domain), ids))
        else:
            tx = list(map(lambda id: self.__translate_preferred_id(id, db, domain, preferred_ids), ids))
        
        # select
        if select_one is True:
            tx = list(map(lambda x: self.__select_first(x), tx))
        
        if len(tx)==1:
            tx = tx[0]
        
        return tx
        
        
    def __translate_id(self, id:str, db:str, domain:str) -> list:
        """
        Translate a single ID to/from a given database.

        Parameters:
        -----------
        id : str
            The ID to be translated.
        db : str
            External database name involved in the translation.
        domain : str
            Defines the translation direction (e.g., 'symbol_to_id', 'id_to_symbol').

        Returns:
        --------
        list
            A list of translated IDs or the default null ID if not found.
        """
        
        if id in self.master[db][domain]:
            return self.master[db][domain][id]
        else:
            return self.default_null_id
        
        
    def __select_first(self, x:list) -> str:
        """
        Select the first element of a list, if available.
        
        Parameters:
        -----------
        x : list
            A list of IDs or None.
        
        Returns:
        --------
        str
            The first element of the list or the default null ID.
        """

        if x is None:
            return None
        else:
            if isinstance(x, str):
                return x
            else: 
                if len(x)>0:
                    return x[0]
                else:
                    return self.default_null_id
            
        
    def __translate_preferred_id(self, id:str, db:str, domain:str, preferred_ids:list) -> list:
        """
        Translate a single ID with prioritization of preferred IDs.

        Parameters:
        -----------
        id : str
            The ID to be translated.
        db : str
            External database name involved in the translation.
        domain : str
            Defines the translation direction.
        preferred_ids : list
            List of IDs to prioritize during translation.

        Returns:
        --------
        list
            A list of translated IDs or the default null ID if not found.
        """
        
        if id in self.master[db][domain]:
            self.master[db][domain][id]
            tids = self.master[db][domain][id]
            cross = [tid in preferred_ids for tid in tids]
            matches = np.where(cross)[0]
            if len(matches) == 0:
                return tids[0]
            else:
                return [tids[x] for x in matches]
        else:
            return self.default_null_id

    
    # cross ids
    def cross_id(self, ids:str, db1:str, db2:str, preferred_ids=[], select_one=False) -> str | list:
        """
        Translate an ID between two databases (e.g., from ENSEMBL to ENTREZ).

        Parameters:
        -----------
        id : str
            The ID to be translated.
        db1 : str
            Source database.
        db2 : str
            Target database.
        preferred_ids : list
            List of IDs to prioritize during translation.
        select_one : bool, optional
            If True, selects the first translation result (default is True).

        Returns:
        --------
        str | list
            The translated ID(s).
        """
        
        ss = self.id_to_symbol(ids, db1)
        return self.symbol_to_id(ss, db=db2, preferred_ids=preferred_ids, select_one=select_one)
    

    # Shortcuts
       
    def id_to_symbol(self, ids:list, db:str, preferred_ids=[], all_synonyms=False) -> list:
        """
        Translate a list of IDs to their official gene symbols.

        Parameters:
        -----------
        ids : list
            List of database-specific IDs.
        db : str
            Source database name.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        all_synonyms : bool, optional
            If True, retrieves all synonyms (default is False).

        Returns:
        --------
        list
            A list of translated symbols.
        """
        
        if all_synonyms is True or len(preferred_ids)>0:
            return self.__translate(ids, db, "id_to_all_symbols", select_one=all_synonyms is False, preferred_ids=preferred_ids)
        else:
            return self.__translate(ids, db, "id_to_symbol", select_one=True, preferred_ids=preferred_ids)
            
    
    
    def symbol_to_id(self, ids:list, db:str, preferred_ids=[], select_one=True) -> list:
        """
        Translate a list of gene symbols to a specific database ID.

        Parameters:
        -----------
        ids : list
            List of gene symbols.
        db : str
            Target database name.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        select_one : bool, optional
            If True, selects the first translation result (default is True).

        Returns:
        --------
        list
            A list of translated IDs.
        """
        
        return self.__translate(ids, db, "symbol_to_id", select_one=select_one, preferred_ids=preferred_ids)
    
    
    # ENSEMBL GENE
    
    def ensembl_gene_to_symbol(self, ids:list, preferred_ids=[], all_synonyms=False) -> list:
        """
        Translate a list of ENSEMBL_GENE IDs to their official gene symbols.

        Parameters:
        -----------
        ids : list
            List of database-specific IDs.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        all_synonyms : bool, optional
            If True, retrieves all synonyms (default is False).

        Returns:
        --------
        list
            A list of translated symbols.
        """
        
        return self.id_to_symbol(ids, "ensembl_gene", preferred_ids=preferred_ids, all_synonyms=all_synonyms)
    
    
    def symbol_to_ensembl_gene(self, ids:list) -> list:
        """
        Translate a list of gene symbols to their corresponding ENSEMBL_GENE ids

        Parameters:
        -----------
        ids : list
            List of gene symbols.

        Returns:
        --------
        list
            A list of translated IDs.
        """
        
        return self.symbol_to_id(ids, "ensembl_gene")
    

    # ENTREZ
    
    def entrez_to_symbol(self, ids:list, preferred_ids=[], all_synonyms=False) -> list:
        """
        Translate a list of ENTREZ IDs to their official gene symbols.

        Parameters:
        -----------
        ids : list
            List of database-specific IDs.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        all_synonyms : bool, optional
            If True, retrieves all synonyms (default is False).

        Returns:
        --------
        list
            A list of translated symbols.
        """
        
        return self.id_to_symbol(ids, "entrez", preferred_ids=preferred_ids, all_synonyms=all_synonyms)
        
    
    def symbol_to_entrez(self, ids:list) -> list:
        """
        Translate a list of gene symbols to their corresponding ENTREZ ids

        Parameters:
        -----------
        ids : list
            List of gene symbols.

        Returns:
        --------
        list
            A list of translated IDs.
        """
        
        return self.symbol_to_id(ids, "entrez")
    
    
    # UNIPROT
    
    def uniprot_to_symbol(self, ids:list, preferred_ids=[], all_synonyms=False) -> list:
        """
        Translate a list of UNIPROT IDs to their official gene symbols.

        Parameters:
        -----------
        ids : list
            List of database-specific IDs.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        all_synonyms : bool, optional
            If True, retrieves all synonyms (default is False).

        Returns:
        --------
        list
            A list of translated symbols.
        """
        
        return self.id_to_symbol(ids, "uniprot", preferred_ids=preferred_ids, all_synonyms=all_synonyms)
        
    
    def symbol_to_uniprot(self, ids:list) -> list:
        """
        Translate a list of gene symbols to their corresponding UNIPROT ids

        Parameters:
        -----------
        ids : list
            List of gene symbols.

        Returns:
        --------
        list
            A list of translated IDs.
        """
        
        return self.symbol_to_id(ids, "uniprot", select_one=False)
    
    # HGNC
    
    def hgnc_to_symbol(self, ids:list, preferred_ids=[], all_synonyms=False) -> list:
        """
        Translate a list of HGNC IDs to their official gene symbols.

        Parameters:
        -----------
        ids : list
            List of database-specific IDs.
        preferred_ids : list, optional
            Preferred IDs for translation priority (default is empty).
        all_synonyms : bool, optional
            If True, retrieves all synonyms (default is False).

        Returns:
        --------
        list
            A list of translated symbols.
        """
        
        return self.id_to_symbol(ids, "hgnc", preferred_ids=preferred_ids, all_synonyms=all_synonyms)
        
    
    def symbol_to_hgnc(self, ids:list) -> list:
        """
        Translate a list of gene symbols to their corresponding HGNC ids

        Parameters:
        -----------
        ids : list
            List of gene symbols.

        Returns:
        --------
        list
            A list of translated IDs.
        """
        
        return self.symbol_to_id(ids, "hgnc")
    
    # OFFICAL SYMBOL
    
    def official_symbol(self, ids:list) -> list:
        """
        Retrieve the official gene symbols for the given IDs.
    
        Parameters:
        -----------
        ids : list
            A list of gene identifiers.
    
        Returns:
        --------
        list
            A list of official gene symbols.
        """
        
        return self.__translate(ids, "symbol", "official", select_one=True)
    
    
    def synonyms(self, ids:list) -> list:
        """
        Retrieve synonyms for the given gene IDs.
    
        Parameters:
        -----------
        ids : list
            A list of gene identifiers.
    
        Returns:
        --------
        list
            A list of lists containing synonyms for each ID.
        """
        return self.__translate(ids, "symbol", "synonyms", select_one=False)
    

    
    # Utilities
    
    def get_info(self) -> dict:
        """
        Get metadata about the current master database.

        Returns:
        --------
        dict
            A dictionary of metadata fields.
        """
        
        return self.master["info"]
    
    
    def random_symbols(self, n:int) -> list:
        """
        Select a random set of N gene symbols.

        Parameters:
        -----------
        n : int
            The number of random genes to select.

        Returns:
        --------
        list
            A list of random gene symbols.
        """
        
        official_symbols = list(self.master["symbol"]["official"].values())
        indexes = np.random.randint(0, len(official_symbols), n)
        return [official_symbols[i] for i in indexes]
    
    

    