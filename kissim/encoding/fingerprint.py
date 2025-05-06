"""
kissim.encoding.fingerprint

Defines the kissim fingerprint.
"""

import logging
from typing import Union

from kissim.io import KlifsToKissimData, PocketBioPython, PocketDataFrame, LigandDataFrame
from kissim.encoding import FingerprintBase
from kissim.encoding.features import (
    SiteAlignFeature,
    SideChainOrientationFeature,
    SolventExposureFeature,
    SubpocketsFeature,
    LigandFeature,
)

logger = logging.getLogger(__name__)


class Fingerprint(FingerprintBase):
    @classmethod
    def from_structure_klifs_id(cls, structure_klifs_id, klifs_session=None):
        """
        Calculate fingerprint for a KLIFS structure (by structure KLIFS ID).

        Parameters
        ----------
        structure_klifs_id : int
            Structure KLIFS ID.
        klifs_session : opencadd.databases.klifs.session.Session or None
            Local or remote KLIFS session.
            If None (default), set up remote KLIFS session.

        Returns
        -------
        kissim.encoding.Fingerprint
            Fingerprint.
        """

        data = KlifsToKissimData.from_structure_klifs_id(structure_klifs_id, klifs_session)
        if data is None:
            logger.warning(f"{structure_klifs_id}: Empty fingerprint (data unaccessible).")
            fingerprint = None
        else:
            fingerprint = cls.from_text(
                data.text,
                data.extension,
                data.residue_ids,
                data.residue_ixs,
                data.structure_klifs_id,
                data.kinase_name,
                data.ligand_expo_id,
            )
        return fingerprint

    @classmethod
    def from_text(
        cls, 
        text: str, 
        extension: str, 
        residue_ids: Union[list, int], 
        residue_ixs: Union[list, int], 
        structure_name: str, 
        kinase_name: str,
        ligand_expo_id: str = None,
        ):
        """
        Calculate fingerprint for a KLIFS structure (by complex data as text and pocket residue
        IDs and indices).

        Parameters
        ----------
        text : str
            Structural complex data as string (file content).
        extension : str
            Structural complex data format (file extension).
        residue_ids : list of int
            Pocket residue IDs.
        residue_ixs : list of int
            Pocket residue indices.
        structure_name : str  # TODO or structure_klifs_id?
            Structure name.
        kinase_name : str
            Kinase name.

        Returns
        -------
        kissim.encoding.Fingerprint
            Fingerprint.
        """

        # BioPython-based and DataFrame-based pocket are both necessary for fingerprint features
        pocket_bp = PocketBioPython.from_text(
            text, extension, residue_ids, residue_ixs, structure_name
        )
        pocket_df = PocketDataFrame.from_text(
            text, extension, residue_ids, residue_ixs, structure_name
        )
        if ligand_expo_id:
            ligand_df = LigandDataFrame.from_text(
                text, extension, ligand_expo_id
            )
        if pocket_bp is None or pocket_df is None:
            logger.warning(f"{structure_name}: Empty fingerprint (pocket unaccessible).")
            fingerprint = None
        else:
            fingerprint = cls()
            fingerprint.structure_klifs_id = structure_name
            fingerprint.kinase_name = kinase_name
            fingerprint.residue_ids = pocket_bp._residue_ids
            fingerprint.residue_ixs = pocket_bp._residue_ixs
            values_dict = {}
            values_dict["physicochemical"] = fingerprint._get_physicochemical_features_dict(
                pocket_bp
            )
            values_dict["spatial"] = fingerprint._get_spatial_features_dict(pocket_df)
            values_dict["ligand"] = fingerprint._get_ligand_features_dict(pocket_df, ligand_df)
            fingerprint.values_dict = values_dict

        return fingerprint

    def _get_physicochemical_features_dict(self, pocket_bp):
        """
        Get physicochemical features.

        Parameters
        ----------
        pocket_bp : kissim.io.PocketBioPython
            Biopython-based pocket object.

        Returns
        -------
        dict of list of float
            Feature values (values) for physicochemical properties (keys).
        """

        # Set up physicochemical features
        features = {}
        # Add SiteAlign features
        for sitealign_feature_name in ["size", "hbd", "hba", "charge", "aromatic", "aliphatic"]:
            feature = SiteAlignFeature.from_pocket(pocket_bp, sitealign_feature_name)
            features[sitealign_feature_name] = feature.values
        # Add side chain orientation feature
        feature = SideChainOrientationFeature.from_pocket(pocket_bp)
        features["sco"] = feature.values
        features["sco.vertex_angle"] = feature._vertex_angles
        # Add solvent exposure feature
        feature = SolventExposureFeature.from_pocket(pocket_bp)
        features["exposure"] = feature.values
        features["exposure.ratio"] = feature._ratio

        return features

    def _get_spatial_features_dict(self, pocket_df):
        """
        Get spatial features (distances and moments).

        Parameters
        ----------
        pocket_df : kissim.io.PocketDataFrame
            DataFrame-based pocket object.

        Returns
        -------
        dict of list of float
            Per-subpocket feature values [and coordinates] (values) for distances and moments
            [and subpocket centers] (keys).
        """

        # Set up spatial features
        features = {}
        # Add subpockets features
        feature = SubpocketsFeature.from_pocket(pocket_df)
        features["distances"] = feature._distances
        features["moments"] = feature._moments
        features["subpocket_centers"] = feature._subpocket_centers

        return features

    def _get_ligand_features_dict(self, pocket_df, ligand_df):
        """
        Get ligand-based spatial features.

        Parameters
        ----------
        pocket_df : kissim.io.PocketDataFrame
            DataFrame-based pocket object.
        ligand_df : kissim.io.LigandDataFrame
            DataFrame-based ligand object.

        Returns
        -------
        dict of list of float
            Distance values (values) for ligand-related geometric points (keys).
        """

        # Set up physicochemical features
        features = {}
        # Add ligand features
        feature = LigandFeature.from_pocket(pocket_df, ligand_df)
        features["ctd"] = feature._ligand['dist_ctd']
        features["cst"] = feature._ligand['dist_cst']
        features["fct"] = feature._ligand['dist_fct']
        features["ftf"] = feature._ligand['dist_ftf']

        return features
