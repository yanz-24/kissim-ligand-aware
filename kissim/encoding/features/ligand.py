"""
kissim.encoding.feature.ligand

Defines the Ligand features.
"""

import numpy as np
import pandas as pd

from kissim.encoding.features.base import BaseFeature
from kissim.io.dataframe import PocketDataFrame


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
    def from_pocket(cls, pocket: PocketDataFrame, ligand): #TODO: add ligand type hint
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
        distances = feature._calculate_distance_ligand(pocket, ligand)
        feature._distances_ctd = distances["dist_ctd"].to_numpy()
        feature._distances_cst = distances["dist_cst"].to_numpy()
        feature._distances_fct = distances["dist_fct"].to_numpy()
        feature._distances_ftf = distances["dist_ftf"].to_numpy()

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


    def _calculate_distance_ligand(self, pocket: pd.DataFrame, ligand_df: pd.DataFrame):
        """
        Calculate distances between ligand's key geometric points and all 
        pocket residues (CA atoms).

        Parameters
        ----------
        pocket : kissim.io.PocketDataFrame
            Pocket object.
        TODO: add ligand

        Returns
        -------
        pandas.DataFrame
            DataFrame of distances from ligand to pocket residues. Each 
        """

        # TODO: implement this function
        # fetch ligand coords
        # fetch pocket residues coords
        # calculate the centroid of ligand 
        # calculate distances from ligand centroid to all pocket residues
        # calculate distances from ligand closest heavy atom to all pocket residues
        pocket_coords = pocket.ca_atoms[["atom.x", "atom.y", "atom.z"]].to_numpy()
        ligand_coords = ligand_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
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
        residues_df = pd.DataFrame()
        residues_df['dist_ctd'] = centroid_distances
        residues_df['dist_cst'] = closest_atom_distances
        residues_df['dist_fct'] = farthest_atom_distances
        residues_df['dist_ftf'] = ftf_distances

        return residues_df
