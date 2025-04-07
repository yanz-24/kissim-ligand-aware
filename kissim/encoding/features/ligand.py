"""
kissim.encoding.feature.ligand

Defines the Ligand features.
"""

import numpy as np
import pandas as pd

from kissim.encoding.features.base import BaseFeature


class LigandFeature(BaseFeature):
    """
    Ligand-based features: distances from ligand center to pocket residues.

    Attributes
    ----------
    name : str or int
        Name for structure encoding by this feature.
    _residue_ids : list of int
        Residue IDs.
    _residue_ixs : list of int
        Residue indices.
    _distances_ctd : list of float
        Distance from pocket residues to ligand centroid.
    _distances_cst : list of float
        Distance from pocket residues to ligand closest heavy atom from the residue.
    _distances_fct : list of float
        Distance from pocket residues to ligand furthest heavy atom from the residue.
    _distances_ftf : list of float
        Distance from pocket residues to the ligand atom farthest from the furthest atom.
    """

    def __init__(self):
        self.name = None
        self._residue_ids = None
        self._residue_ixs = None
        self._distances_ctd = None
        self._distances_cst = None
        self._distances_fct = None
        self._distances_ftf = None

    @classmethod
    def from_pocket(cls, pocket):
        """
        Generate ligand-based features from pocket.

        Parameters
        ----------
        pocket : kissim.io.PocketBioPython
            Biopython-based pocket object.

        Returns
        -------
        kissim.encoding.features.LigandFeature
            Ligand feature object.
        """
        feature = cls()
        feature.name = pocket.name
        feature._residue_ids = pocket._residue_ids
        feature._residue_ixs = pocket._residue_ixs

        # Calculate distances from ligand to pocket residues
        distances = pocket.calculate_distance_ligand()
        # TODO: add distance calculations

        return feature

    @property
    def values(self):
        """
        Feature values.

        Returns
        -------
        list of float
            Concatenation of all ligand-residue distances.
        """
        return (
            self._distances_ctd
            + self._distances_cst
            + self._distances_fct
            + self._distances_ftf
        )

    @property
    def details(self):
        """
        Feature details.

        Returns
        -------
        pandas.DataFrame
            DataFrame of residue-wise distance metrics.
        """
        return pd.DataFrame({
            "residue_id": self._residue_ids,
            "ctd": self._distances_ctd,
            "cst": self._distances_cst,
            "fct": self._distances_fct,
            "ftf": self._distances_ftf,
        })


    def _calculate_distance_ligand(self):
        """
        Calculate distances between ligand's key geometric points and all 
        pocket residues (CA atoms).

        Parameters
        ----------
        pocket : kissim.io.PocketBioPython
            Biopython-based pocket object.

        Returns
        -------
        list of float
            Distances between ligand's key geometric points and all pocket residues.
        """

        # TODO: implement this function
        # fetch ligand coords
        # fetch pocket residues coords
        # calculate the centroid of ligand 
        # calculate distances from ligand centroid to all pocket residues
        # calculate distances from ligand closest heavy atom to all pocket residues
