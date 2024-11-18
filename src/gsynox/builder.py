#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import pickle
import datetime
import numpy as np
import re


class Builder:
    """
    A class to build and manage a symbol mapping master dictionary from HGNC data and other external databases.
    
    Attributes:
    -----------
    master : dict
        A dictionary containing the database mappings.
    """
    
    master = {}
    
    def __safe(self, x: str) -> str:
        """
        Removes '_id' and trailing characters from a string and trims whitespace to create a safe id.

        Parameters:
        -----------
            x : str
                The input string to clean.

        Returns:
        --------
            str: 
                The cleaned string.
        """
        return re.sub('_id(.)*', '', x).strip()

    def build(self, hgnc_file: str, key_symbol="symbol", alias_field="alias_symbol", prev_field="prev_symbol", fields=None) -> None:
        """
        Constructs the `master` dictionary from an HGNC data file.

        Parameters:
        -----------
            hgnc_file : str
                Path to the HGNC TSV file.
            key_symbol : str 
                Column name for the primary symbol. Default is 'symbol'.
            alias_field : str
                Column name for alias symbols. Default is 'alias_symbol'.
            prev_field : str
                Column name for previous symbols. Default is 'prev_symbol'.
            fields: list
                Additional fields to include. Default is `["ensembl_gene_id", "entrez_id", "hgnc_id"]`.

        Returns:
        --------
            None
        """
        fields = fields or ["ensembl_gene_id", "entrez_id", "hgnc_id"]

        # Load raw file
        hgnc = pd.read_csv(hgnc_file, sep="\t", dtype="str")
        hgnc = hgnc.fillna("")
        
        # Initialize master dictionary
        self.master = {}
        x = datetime.datetime.now()
        self.master["info"] = {
            "version": 0.1,
            "date": x.strftime("%c"),
            "dbs": [key_symbol] + list(map(self.__safe, fields)),
            "log": [],
        }
        self.master["symbol"] = {"official": {}, "synonyms": {}}
        for f in fields:
            self.master[self.safe(f)] = {
                "id_to_symbol": {},
                "id_to_all_symbols": {},
                "symbol_to_id": {},
            }

        # Fill master
        n = hgnc.shape[0]
        official_symbols = hgnc[["symbol"]].values
        step = 1000
        for i in range(n):
            if (i % step) == 0:
                print(i, " of ", n)

            symbol = hgnc.loc[i, key_symbol]
            # Extract alias symbols
            raw_alias_symbol = hgnc.loc[i, alias_field]
            alias_symbols = raw_alias_symbol.split("|") if raw_alias_symbol else []

            # Extract previous symbols
            raw_prev_symbol = hgnc.loc[i, prev_field]
            prev_symbols = raw_prev_symbol.split("|") if raw_prev_symbol else []

            synonyms = np.unique(alias_symbols + prev_symbols).tolist()

            for f in fields:
                raw = hgnc.loc[i, f]
                if raw:
                    ids = raw.split("|")
                    for u in ids:
                        self.master[self.safe(f)]["id_to_symbol"][u] = symbol
                        self.master[self.safe(f)]["id_to_all_symbols"][u] = [symbol] + synonyms
                    self.master[self.safe(f)]["symbol_to_id"][symbol] = ids
                    for s in synonyms:
                        if s and s not in official_symbols:
                            self.master[self.safe(f)]["symbol_to_id"][s] = ids

            self.master["symbol"]["official"][symbol] = symbol
            for s in synonyms:
                if s and s not in official_symbols:
                    self.master["symbol"]["official"][s] = symbol
            self.master["symbol"]["synonyms"][symbol] = synonyms

    def save(self, output_file: str) -> None:
        """
        Saves the `master` dictionary to a pickle file.

        Parameters:
        -----------
            output_file : str
                Path to the output file.

        Returns:
        --------
            None
        """
        with open(output_file, "wb") as outfile:
            pickle.dump(self.master, outfile)

    def add_db(self, db_file: str, name: str) -> None:
        """
        Adds a new external database to the `master` dictionary.

        Parameters:
        -----------
            db_file : str
                Path to the external database file.
            name : str
                Name of the new database to add.

        Returns:
        --------
            None
        """
        self.master[name] = {"id_to_symbol": {}, "id_to_all_symbols": {}, "symbol_to_id": {}}
        x = datetime.datetime.now()
        self.master["info"]["log"].append(f"Added {name} on {x.strftime('%c')}")
        self.master["info"]["dbs"].append(name)

        # Load data
        db = pd.read_csv(db_file, sep="\t", dtype="str")
        db = db.fillna("")

        for i in range(db.shape[0]):
            symbol = db.iloc[i, 0]
            osymbol = self.master["symbol"]["official"][symbol]
            synonyms = self.master["symbol"]["synonyms"][osymbol]

            raw_ids = db.iloc[i, 1]
            if raw_ids:
                ids = raw_ids.split(",")
                self.master[name]["symbol_to_id"][osymbol] = ids
                for s in synonyms:
                    self.master[name]["symbol_to_id"][s] = ids
                for id in ids:
                    self.master[name]["id_to_symbol"][id] = osymbol
                    self.master[name]["id_to_all_symbols"][id] = [osymbol] + synonyms
