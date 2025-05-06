"""
kissim.encoding.feature.ligand

Defines the Ligand features.
"""

import numpy as np
import pandas as pd

from kissim.encoding.features.base import BaseFeature
from kissim.io.dataframe import PocketDataFrame, LigandDataFrame


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
    _ligand : dict of (str: list of float)
        Distances from ligand to pocket residues. The keys are:
        - dist_ctd: distance to centroid
        - dist_cst: distance to closest atom
        - dist_fct: distance to farthest atom
        - dist_ftf: distance to farthest-to-farthest atom
    """

    def __init__(self):
        self.name = None
        self._residue_ids = None
        self._residue_ixs = None
        self._ligand = None

    @classmethod
    def from_pocket(cls, pocket: PocketDataFrame, ligand: LigandDataFrame):
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
        feature._ligand = feature._calculate_distance_ligand(pocket, ligand)

        return feature

    @property
    def values(self):
        """
        Ligand feature values.

        Returns
        -------
        list of float
            Concatenation of all ligand-residue distances.
        """
        return self._ligand

    @property
    def details(self):
        """
        Feature details for ligand features.

        Returns
        -------
        pandas.DataFrame
            DataFrame of residue-wise distance metrics.
        """
        return pd.DataFrame({
            "residue_id": self._residue_ids,
            "ctd": self._ligand['dist_ctd'],
            "cst": self._ligand['dist_cst'],
            "fct": self._ligand['dist_fct'],
            "ftf": self._ligand['dist_ftf'],
        })


    def _calculate_distance_ligand(self, pocket: PocketDataFrame, ligand: LigandDataFrame):
        """
        Calculate distances between ligand's key geometric points and all 
        pocket residues (CA atoms).

        Parameters
        ----------
        pocket : kissim.io.PocketDataFrame
            Pocket object.
        ligand_df : kissim.io.LigandDataFrame.df 
            pandas.DataFrame object.

        Returns
        -------
        dict of (str: list of float)
            Distances from ligand to pocket residues. 
        """
        pocket_coords = pocket.ca_atoms[["atom.x", "atom.y", "atom.z"]].to_numpy()
        ligand_coords = ligand.df[["atom.x", "atom.y", "atom.z"]].to_numpy()
        ligand_centroid = ligand_coords.mean(axis=0)

        centroid_distances = []
        closest_atom_distances = []
        farthest_atom_distances = []
        ftf_distances = []

        for res_coord in pocket_coords:
            if pd.isnull(res_coord).all():
                centroid_distances.append(np.nan) # ctd
                closest_atom_distances.append(np.nan) # cst
                farthest_atom_distances.append(np.nan)
                ftf_distances.append(np.nan)

            else:
                # Distance to centroid (ctd)
                d_centroid = np.linalg.norm(res_coord - ligand_centroid)

                # Distances from this residue to all ligand atoms
                dists = np.linalg.norm(ligand_coords - res_coord, axis=1)

                # Closest and farthest atom indices
                idx_closest = np.argmin(dists) # cst
                idx_farthest = np.argmax(dists) # fct
                fct_coord = ligand_coords[idx_farthest]

                # Now compute distances from fct to all other ligand atoms
                fct_to_others = np.linalg.norm(ligand_coords - fct_coord, axis=1)
                idx_ftf = np.argmax(fct_to_others)
                ftf_coord = ligand_coords[idx_ftf]

                # Distance from residue to ftf
                d_ftf = np.linalg.norm(res_coord - ftf_coord)

                # Append all distances
                centroid_distances.append(d_centroid) # ctd
                closest_atom_distances.append(dists[idx_closest]) # cst
                farthest_atom_distances.append(dists[idx_farthest])
                ftf_distances.append(d_ftf)

        # Add results to DataFrame
        ligand_distances = {}
        ligand_distances['dist_ctd'] = centroid_distances
        ligand_distances['dist_cst'] = closest_atom_distances
        ligand_distances['dist_fct'] = farthest_atom_distances
        ligand_distances['dist_ftf'] = ftf_distances

        return ligand_distances
